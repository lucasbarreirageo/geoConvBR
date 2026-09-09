#' Interactive map of a species' points, EOO and AOO (MapBiomas + fire)
#'
#' Builds a \pkg{leaflet} map with the occurrence points, EOO hull and occupied
#' AOO cells. Optionally overlays the MapBiomas land-cover raster and the
#' MapBiomas Fire accumulated layer (fire recurrence: number of years burned),
#' both clipped to the EOO, as toggleable layers with legends.
#'
#' @param assessment A \code{geoconv_assessment} from \code{\link{assess_species}}.
#' @param species Species name to plot (default: first assessed).
#' @param mapbiomas Logical; overlay the MapBiomas raster (default \code{TRUE}).
#' @param fire Logical; overlay the MapBiomas Fire accumulated layer
#'   (default \code{FALSE}).
#' @param src Optional MapBiomas LULC GeoTIFF path/URL.
#' @param fire_src Optional MapBiomas Fire GeoTIFF path/URL.
#' @param max_pixels Target maximum raster size (longest side, pixels) for the
#'   overlays (default \code{800}).
#' @param lang Legend language: \code{"pt"} (default) or \code{"en"}. Controls
#'   the language of the MapBiomas class labels and the fire-frequency legend.
#' @param clip Geometry the rasters are clipped to: \code{"eoo"} (default, the
#'   EOO hull), \code{"aoo"} (the union of occupied AOO cells) or \code{"all"}
#'   (the union of both, so the raster spans the whole range).
#' @param protected Logical; overlay protected areas. Uses the protected-area
#'   layer already stored by \code{assess_species(..., protected = TRUE)}, or
#'   fetches it on the fly from WDPA otherwise. Default \code{FALSE}.
#' @param pa_src Optional local protected-area vector file for the overlay
#'   (offline).
#' @param pa_occ_only Logical; when \code{TRUE} (default) the map draws only the
#'   UCs that actually contain occurrences of the species. Set \code{FALSE} to
#'   draw every UC the range overlaps. Does not affect the EOO/AOO overlap
#'   percentages, which are computed over the whole range.
#' @return A \code{leaflet} widget.
#' @export
map_species <- function(assessment, species = NULL, mapbiomas = TRUE,
                        fire = FALSE, src = NULL, fire_src = NULL,
                        max_pixels = 800, lang = c("pt", "en"),
                        clip = c("eoo", "aoo", "all"),
                        protected = FALSE, pa_src = NULL, pa_occ_only = TRUE) {
  lang <- match.arg(lang)
  clip <- match.arg(clip)
  lulc_col   <- if (lang == "en") "class_en" else "class_pt"
  lbl_mb     <- if (lang == "en") "Land cover %s" else "Cobertura %s"
  grp_lc     <- if (lang == "en") "Land cover" else "Cobertura"
  lbl_fire   <- if (lang == "en") "Fire frequency<br>(years, 1985-2024)"
                else "Frequencia de fogo<br>(anos, 1985-2024)"
  grp_fire   <- if (lang == "en") "Fire frequency" else "Frequencia de fogo"
  stopifnot(inherits(assessment, "geoconv_assessment"))
  if (!requireNamespace("leaflet", quietly = TRUE))
    stop("Package 'leaflet' is required.", call. = FALSE)
  d <- assessment$detail
  if (is.null(species)) species <- names(d)[1]
  obj <- d[[species]]
  if (is.null(obj)) stop("Species not found in assessment: ", species, call. = FALSE)
  st <- assessment$settings

  # Geometry the rasters are clipped to: the EOO hull, the union of AOO cells,
  # or (clip = "all") the union of both so the raster spans the whole range.
  clip_geom <- .clip_geometry(clip, obj$eoo$hull, obj$aoo$cells)

  m <- leaflet::leaflet()
  # Key-free basemaps: the CartoDB Positron basemap now requires an API key
  # (tiles render an "API KEY REQUIRED" watermark), so use OpenStreetMap for the
  # light layer and Esri World Imagery for satellite - both usable without a key.
  m <- leaflet::addProviderTiles(m, "OpenStreetMap.Mapnik", group = "Light")
  m <- leaflet::addProviderTiles(m, "Esri.WorldImagery", group = "Satellite")
  
  # --- optional MapBiomas LULC overlay (under everything) ---
  mb_on <- FALSE
  if (isTRUE(mapbiomas) && !is.null(clip_geom) &&
      requireNamespace("terra", quietly = TRUE)) {
    # Prefer the product actually used for this species (may differ from the
    # assessment default when the range fell back to global Sentinel-2/Esri).
    lc_ini  <- obj$initiative %||% st$initiative %||% "brazil"
    lc_year <- obj$year %||% st$year
    lc_coll <- obj$collection %||% st$collection
    rr <- tryCatch({
      r <- .mb_raster_display(clip_geom, lc_year, lc_coll, src,
                              max_pixels = max_pixels, crs = NULL,
                              initiative = lc_ini)
      leg <- mb_legend(lc_coll, lc_ini)
      vals <- terra::unique(r)[, 1]; vals <- vals[!is.na(vals)]
      list(r = r, keep = leg[leg$code %in% vals, , drop = FALSE], year = lc_year)
    }, error = function(e) {
      warning("MapBiomas layer skipped: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (!is.null(rr) && nrow(rr$keep) > 0) {
      pal <- leaflet::colorFactor(rr$keep$hex, levels = rr$keep$code,
                                  na.color = "transparent")
      m <- leaflet::addRasterImage(m, rr$r, colors = pal, opacity = 0.75,
                                   method = "ngb", project = TRUE,
                                   group = grp_lc)
      m <- leaflet::addLegend(m, position = "bottomright",
                              colors = rr$keep$hex, labels = rr$keep[[lulc_col]],
                              opacity = 0.75,
                              title = sprintf(lbl_mb, rr$year),
                              group = grp_lc)
      mb_on <- TRUE
    }
  }
  
  # --- optional MapBiomas Fire overlay (over LULC, under vectors) ---
  fire_on <- FALSE
  if (isTRUE(fire) && !is.null(clip_geom) &&
      requireNamespace("terra", quietly = TRUE)) {
    fr <- tryCatch(
      .fire_raster_display(clip_geom,
                           fire_collection = st$fire_collection %||% 4,
                           host_collection = st$fire_host_collection %||% 9,
                           src = fire_src, max_pixels = max_pixels),
      error = function(e) {
        warning("Fire layer skipped: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    fv <- if (!is.null(fr)) {
      v <- terra::values(fr); v[is.finite(v)]
    } else numeric(0)
    if (length(fv) > 0) {
      brks <- c(1, 2, 4, 7, 11, 21, 41)          # fire-frequency classes (years)
      pal  <- leaflet::colorBin(fire_palette(length(brks) - 1), domain = c(1, 40),
                                bins = brks, right = FALSE, na.color = "transparent")
      m <- leaflet::addRasterImage(m, fr, colors = pal, opacity = 0.8,
                                   method = "ngb", project = TRUE,
                                   group = grp_fire)
      m <- leaflet::addLegend(m, position = "bottomleft", pal = pal,
                              values = c(1, 40), opacity = 0.8,
                              title = lbl_fire,
                              group = grp_fire)
      fire_on <- TRUE
    }
  }
  
  # --- optional federal Protected areas (UCs) overlay ---
  grp_uc <- if (lang == "en") "Protected areas" else "Areas Protegidas"
  pa_on <- FALSE
  if (isTRUE(protected)) {
    pa_sf <- obj$pa$layer
    if (is.null(pa_sf) && !is.null(clip_geom))
      pa_sf <- tryCatch(protected_areas(clip_geom, src = pa_src),
                        error = function(e) NULL)
    if (!is.null(pa_sf) && nrow(pa_sf) > 0) {
      pa_sf <- sf::st_transform(pa_sf, 4326)
      if (isTRUE(pa_occ_only) && !is.null(obj$points)) {
        occ <- sf::st_transform(sf::st_geometry(obj$points), 4326)
        pa_sf <- pa_sf[lengths(suppressMessages(
          sf::st_intersects(pa_sf, occ))) > 0, , drop = FALSE]
      }
    }
    if (!is.null(pa_sf) && nrow(pa_sf) > 0) {
      lab <- sprintf("<b>%s</b><br>%s", pa_sf$pa_name,
                     ifelse(is.na(pa_sf$pa_category), "", pa_sf$pa_category))
      m <- leaflet::addPolygons(m, data = pa_sf, color = "#1f8d49", weight = 1.5,
                                fillColor = "#1f8d49", fillOpacity = 0.15,
                                popup = lab, label = pa_sf$pa_name,
                                group = grp_uc)
      pa_on <- TRUE
    }
  }

  if (!is.null(obj$eoo$hull))
    m <- leaflet::addPolygons(m, data = sf::st_transform(obj$eoo$hull, 4326),
                              color = "#1f4e79", weight = 2,
                              fillColor = "#1f4e79", fillOpacity = 0.08,
                              group = "EOO (hull)")
  if (!is.null(obj$aoo$cells))
    m <- leaflet::addPolygons(m, data = sf::st_transform(obj$aoo$cells, 4326),
                              color = "#9c0027", weight = 1,
                              fillColor = "#d4271e", fillOpacity = 0.25,
                              group = "AOO (2 km cells)")
  pts <- sf::st_transform(sf::st_geometry(obj$points), 4326)
  co <- sf::st_coordinates(pts)
  m <- leaflet::addCircleMarkers(m, lng = co[, 1], lat = co[, 2], radius = 4,
                                 color = "#111111", fillColor = "#f1c40f",
                                 fillOpacity = 0.9, weight = 1,
                                 group = "Occurrences")
  
  overlay <- c("Occurrences", "EOO (hull)", "AOO (2 km cells)")
  if (pa_on)   overlay <- c(grp_uc, overlay)
   if (fire_on) overlay <- c(grp_fire, overlay)
  if (mb_on)   overlay <- c(grp_lc, overlay)
  m <- leaflet::addLayersControl(
    m, baseGroups = c("Light", "Satellite"), overlayGroups = overlay,
    options = leaflet::layersControlOptions(collapsed = FALSE))
  # Species name box: placed top-left (below the zoom control) and allowed to
  # wrap, so a long italic name is never clipped by the right edge of the map as
  # it was when anchored top-right with white-space:nowrap.
  m <- leaflet::addControl(
    m,
    html = sprintf(
      "<div style='display:inline-block;width:max-content;max-width:280px;white-space:normal;word-break:keep-all;overflow-wrap:normal;background:rgba(255,255,255,.9);padding:2px 8px;border-radius:6px;font-size:14px;line-height:1.25'><b>%s</b></div>",
      .sp_html(species)),
    position = "topleft")
  m
}

#' Habitat composition chart (natural / altered / water / other) for EOO and AOO
#'
#' Draws, for one species, two horizontal stacked bars (EOO and AOO) showing the
#' MapBiomas composition as a percentage of the mapped (observed) area: natural,
#' altered (anthropic), water and other. Segment labels are placed inside each
#' bar and the legend sits below the plot, so nothing overlaps. The headline
#' \emph{terrestrial} converted percentage is shown as a subtitle.
#'
#' When \pkg{ggplot2} is available (a hard dependency) the function returns a
#' \code{ggplot} object, which \code{\link{mas_plotly}} can turn into an
#' interactive chart; otherwise it draws a base-graphics fallback and returns
#' the percentage matrix invisibly.
#'
#' @param assessment A \code{geoconv_assessment}.
#' @param species Species name (default: first).
#' @param lang Label language: \code{"en"} (default) or \code{"pt"}.
#' @return A \code{ggplot} object (ggplot2 available) or, in the base-graphics
#'   fallback, invisibly the plotted percentage matrix.
#' @seealso \code{\link{mas_plotly}} for the interactive version.
#' @export
plot_conversion <- function(assessment, species = NULL, lang = c("en", "pt")) {
  lang <- match.arg(lang)
  stopifnot(inherits(assessment, "geoconv_assessment"))
  d <- assessment$detail
  if (is.null(species)) species <- names(d)[1]
  obj <- d[[species]]
  if (is.null(obj)) stop("Species not found: ", species, call. = FALSE)

  gl <- .mb_group_labels(lang)
  grp <- unname(gl$labels[c("natural", "anthropic", "water", "other")])
  cols <- stats::setNames(unname(gl$hex[c("natural", "anthropic", "water", "other")]), grp)
  other_lab <- grp[4]
  if (lang == "en") {
    xlab <- "% of mapped (observed) area"
    main_suffix_fmt <- " - composition (land cover %s)"
    sub_fmt  <- "Converted (terrestrial): EOO %s  |  AOO %s"
  } else {
    xlab <- "% da area mapeada (observada)"
    main_suffix_fmt <- " - composicao (cobertura %s)"
    sub_fmt  <- "Convertido (terrestre): EOO %s  |  AOO %s"
  }

  comp <- function(cv) {
    if (is.null(cv)) return(stats::setNames(rep(NA_real_, 4), grp))
    a <- stats::setNames(
      c(cv$natural_km2, cv$anthropic_km2, cv$water_km2, cv$other_km2), grp)
    tot <- sum(a, na.rm = TRUE)
    if (!is.finite(tot) || tot <= 0) return(stats::setNames(rep(NA_real_, 4), grp))
    100 * a / tot
  }
  M <- rbind(AOO = comp(obj$aoo_conversion), EOO = comp(obj$eoo_conversion))

  fmt <- function(cv) { v <- if (is.null(cv)) NA_real_ else cv$converted_pct
    if (is.finite(v)) sprintf("%.0f%%", v) else "-" }
  main_suffix <- sprintf(main_suffix_fmt, obj$year %||% assessment$settings$year)
  title      <- .sp_title_expr(species, main_suffix)   # ggplot (plotmath)
  title_html <- .sp_title_html(species, main_suffix)   # plotly / HTML
  subtitle <- sprintf(sub_fmt, fmt(obj$eoo_conversion), fmt(obj$aoo_conversion))

  # --- ggplot2 (interactive-ready) ---
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(
      range = factor(rep(c("AOO", "EOO"), each = length(grp)),
                     levels = c("AOO", "EOO")),
      group = factor(rep(grp, 2), levels = grp),
      pct   = c(M["AOO", ], M["EOO", ]),
      stringsAsFactors = FALSE
    )
    df$lab <- ifelse(is.finite(df$pct) & df$pct >= 6,
                     sprintf("%.0f%%", df$pct), "")
    df$lab_col <- ifelse(as.character(df$group) == other_lab, "black", "white")
    p <- ggplot2::ggplot(
      df, ggplot2::aes(x = .data[["pct"]], y = .data[["range"]],
                       fill = .data[["group"]])) +
      ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[["lab"]], colour = .data[["lab_col"]]),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3.4, fontface = "bold", show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = cols[grp], drop = FALSE, name = NULL) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
      ggplot2::labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
      .mas_theme() +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    attr(p, "pct") <- M   # keep the percentage matrix accessible
    attr(p, "mas") <- list(kind = "bar", df = df, palette = as.list(cols[grp]),
                           levels = grp, xlab = xlab,
                           title = title_html, subtitle = subtitle)
    return(p)
  }

  # --- base-graphics fallback ---
  op <- graphics::par(mar = c(5, 5, 4, 9), xpd = NA); on.exit(graphics::par(op), add = TRUE)
  bp <- graphics::barplot(t(M), horiz = TRUE, col = cols[grp], border = "white",
    xlim = c(0, 100), las = 1, xlab = xlab, main = title)

  for (r in seq_len(nrow(M))) {
    x <- M[r, ]; cm <- cumsum(c(0, x))
    for (k in seq_along(x)) if (is.finite(x[k]) && x[k] >= 6)
      graphics::text(cm[k] + x[k] / 2, bp[r], sprintf("%.0f%%", x[k]),
        cex = 0.8, col = if (grp[k] == other_lab) "black" else "white")
  }
  graphics::mtext(subtitle, side = 3, line = 0.2, cex = 0.85)
  usr <- graphics::par("usr")
  graphics::legend(
    x = usr[2] + 0.02 * diff(usr[1:2]), y = mean(usr[3:4]), yjust = 0.5,
    legend = grp, fill = cols[grp], border = "white",
    bty = "n", cex = 0.9, xpd = NA)
  invisible(M)
}