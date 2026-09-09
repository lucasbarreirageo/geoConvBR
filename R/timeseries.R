#' MapBiomas land-cover composition over time for a polygon
#'
#' For a given area (e.g. an EOO hull or the union of AOO cells), computes the
#' area of each MapBiomas class (or conservation group) for several years,
#' returning a tidy long table suitable for a stacked-area chart of percentage
#' of area versus year (as in the MapBiomas coverage figures).
#'
#' Each year is read independently (local windowed read or Earth Engine), so a
#' long span of annual data can take a while. By default a 5-year step is used;
#' pass \code{years = mb_years()} for the full annual series.
#'
#' @param geom An \code{sf}/\code{sfc} polygon (any CRS).
#' @param years Integer vector of years. If \code{NULL} (default) a 5-year step
#'   across the initiative's span (plus the last year) is used.
#' @param collection Integer MapBiomas collection. \code{NULL} uses the
#'   initiative default.
#' @param initiative One of \code{"brazil"} (default), \code{"amazonia"} or
#'   \code{"colombia"} (see \code{\link{mb_initiatives}}).
#' @param backend \code{"local"} (default, windowed GeoTIFF read) or
#'   \code{"gee"} (Earth Engine; Brazil only).
#' @param src Optional GeoTIFF path/URL for the local backend (overrides the
#'   public per-year URL).
#' @param by \code{"class"} (default; one series per MapBiomas class, with the
#'   official colours) or \code{"group"} (natural / anthropic / water / other).
#' @param include_not_observed Logical; keep the "Not Observed" class
#'   (default \code{FALSE}).
#' @param verbose Logical; print progress (default \code{TRUE}).
#' @param max_pixels Pixel budget for each year's local area read (default
#'   \code{5e7}). Large or widely-spread ranges are read per polygon part and/or
#'   decimated with majority resampling so the read stays bounded; class
#'   proportions are preserved. Ignored by the GEE backend.
#' @return A \code{data.frame} (long) with columns \code{year}, \code{label},
#'   \code{hex}, \code{group}, \code{area_km2} and \code{pct} (percentage of the
#'   mapped area in that year). For \code{by = "class"} it also has \code{code}
#'   and \code{class_en}.
#' @seealso \code{\link{plot_timeseries}}, \code{\link{timeseries_for_species}}
#' @export
cover_timeseries <- function(geom, years = NULL, collection = NULL,
                             initiative = "brazil",
                             backend = c("local", "gee"), src = NULL,
                             by = c("class", "group"),
                             include_not_observed = FALSE, verbose = TRUE,
                             max_pixels = 5e7) {
  backend <- match.arg(backend)
  by <- match.arg(by)
  # Planar (GEOS) geometry, as in assess_species(): avoids the S2 "degenerate
  # edge" rejection of a wide-ranging densified hull during the windowed reads.
  if (requireNamespace("sf", quietly = TRUE)) {
    .old_s2 <- suppressMessages(sf::sf_use_s2())
    on.exit(suppressMessages(sf::sf_use_s2(.old_s2)), add = TRUE)
    suppressMessages(sf::sf_use_s2(FALSE))
  }
  ini <- .mb_resolve_initiative(initiative)
  initiative <- ini$key
  if (is.null(collection)) collection <- ini$collection
  if (backend == "gee" && initiative != "brazil") backend <- "local"
  g <- sf::st_geometry(geom)
  if (length(g) == 0 || all(sf::st_is_empty(g)))
    stop("`geom` is empty; a polygon is required (EOO needs >= 3 unique points).",
         call. = FALSE)

  if (is.null(years)) {
    yy <- mb_years(collection, initiative)
    years <- sort(unique(c(seq(min(yy), max(yy), by = 5), max(yy))))
  }
  years <- sort(unique(as.integer(years)))
  leg <- mb_legend(collection, initiative)

  out <- vector("list", length(years))
  for (i in seq_along(years)) {
    y <- years[i]
    # Release the previous year's raster memory before the next read so the peak
    # footprint stays at one year, not the whole series (avoids the OS killing
    # the R process - "Terminated" - on memory-limited machines).
    if (i > 1L) invisible(gc(FALSE))
    if (verbose) message(sprintf("  coverage %d (%d/%d)...", y, i, length(years)))
    ca <- tryCatch({
      if (backend == "gee") {
        mb_class_areas_gee(g, year = y, collection = collection)
      } else {
        .mb_class_areas_local(g, year = y, collection = collection,
                              initiative = initiative, src = src,
                              max_pixels = max_pixels)
      }
    }, error = function(e) {
      warning(sprintf("year %d failed: %s", y, conditionMessage(e)), call. = FALSE)
      NULL
    })
    if (is.null(ca) || !nrow(ca)) next

    m <- merge(ca, leg, by = "code", all.x = TRUE)
    m$group[is.na(m$group)] <- "other"
    unk <- is.na(m$class_pt); m$class_pt[unk] <- paste0("Class ", m$code[unk])
    m$class_en[is.na(m$class_en)] <- m$class_pt[is.na(m$class_en)]
    m$hex[is.na(m$hex)] <- "#bdbdbd"
    if (!include_not_observed) m <- m[m$group != "not_observed", , drop = FALSE]
    if (!nrow(m)) next
    m$year <- y
    out[[i]] <- m
  }
  invisible(gc(FALSE))

  df <- do.call(rbind, out)
  if (is.null(df) || !nrow(df))
    stop("No data could be computed for the requested years.", call. = FALSE)

  if (by == "group") {
    gl <- .mb_group_labels("pt")   # cover_timeseries stores PT labels; plot_timeseries translates
    grp_lab <- gl$labels
    grp_hex <- gl$hex
    agg <- stats::aggregate(area_km2 ~ year + group, data = df, FUN = sum)
    agg$label <- unname(grp_lab[agg$group])
    agg$hex   <- unname(grp_hex[agg$group])
    res <- agg[, c("year", "group", "label", "hex", "area_km2")]
  } else {
    agg <- stats::aggregate(
      area_km2 ~ year + code + class_pt + class_en + group + hex,
      data = df, FUN = sum)
    names(agg)[names(agg) == "class_pt"] <- "label"
    res <- agg[, c("year", "code", "label", "class_en", "group", "hex", "area_km2")]
  }
  res[] <- lapply(res, function(x) if (is.factor(x)) as.character(x) else x)

  tot <- stats::aggregate(area_km2 ~ year, data = res, FUN = sum)
  names(tot)[2] <- "year_total"
  res <- merge(res, tot, by = "year")
  res$pct <- ifelse(res$year_total > 0, 100 * res$area_km2 / res$year_total, NA_real_)
  res$year_total <- NULL

  # Make the series rectangular: every label present in every year (0 where
  # absent). Without this, a class missing in some years leaves gaps (and
  # single-point groups) that make ggplot2::geom_area() fail.
  res <- .complete_years(res)

  res$area_km2 <- round(res$area_km2, 4)
  res$pct <- round(res$pct, 3)
  res[order(res$year, res$label), , drop = FALSE]
}

#' @keywords internal
#' @noRd
.complete_years <- function(res) {
  yrs <- sort(unique(res$year))
  labs <- unique(res$label)
  if (!length(yrs) || !length(labs)) return(res)
  full <- expand.grid(year = yrs, label = labs,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  meta_cols <- setdiff(names(res), c("year", "area_km2", "pct"))
  meta <- res[!duplicated(res$label), meta_cols, drop = FALSE]
  out <- merge(full, meta, by = "label", all.x = TRUE)
  out <- merge(out, res[, c("year", "label", "area_km2", "pct"), drop = FALSE],
               by = c("year", "label"), all.x = TRUE)
  out$area_km2[is.na(out$area_km2)] <- 0
  out$pct[is.na(out$pct)] <- 0
  out
}

#' Land-cover time series for one assessed species (EOO or AOO)
#'
#' Convenience wrapper that pulls the stored EOO hull or AOO geometry from an
#' assessment and runs \code{\link{cover_timeseries}} on it, reusing the
#' assessment's collection and backend by default.
#'
#' @param assessment A \code{geoconv_assessment} from \code{\link{assess_species}}.
#' @param species Species name (default: first assessed species).
#' @param range \code{"eoo"} (default) or \code{"aoo"}.
#' @param years,by,src,verbose Passed to \code{\link{cover_timeseries}}.
#' @param max_pixels Pixel budget per year, passed to
#'   \code{\link{cover_timeseries}}. Defaults to \code{8e6} (~8 million) - much
#'   lower than the single-snapshot assessment, because a composition
#'   \emph{trend} is a ratio and needs far less resolution, and every year is a
#'   separate read: the coarser budget makes a multi-year series over a large
#'   range several times faster with no meaningful change to the percentages.
#' @return The long data frame from \code{\link{cover_timeseries}}, with
#'   attributes \code{species} and \code{range} set.
#' @export
timeseries_for_species <- function(assessment, species = NULL,
                                   range = c("eoo", "aoo"), years = NULL,
                                   by = c("class", "group"), src = NULL,
                                   verbose = TRUE, max_pixels = 8e6) {
  stopifnot(inherits(assessment, "geoconv_assessment"))
  range <- match.arg(range)
  by <- match.arg(by)
  d <- assessment$detail
  if (is.null(species)) species <- names(d)[1]
  obj <- d[[species]]
  if (is.null(obj)) stop("Species not found in assessment: ", species, call. = FALSE)

  geom <- if (range == "eoo") obj$eoo$hull else .st_union_quiet(obj$aoo$cells)
  if (is.null(geom) || length(sf::st_geometry(geom)) == 0)
    stop(sprintf("No %s geometry for '%s'.", toupper(range), species), call. = FALSE)

  s <- assessment$settings
  # Prefer the product actually used for this species (Sentinel-2 / Esri when the
  # range fell back outside MapBiomas coverage).
  lc_ini  <- obj$initiative %||% s$initiative %||% "brazil"
  lc_coll <- obj$collection %||% s$collection
  ts <- cover_timeseries(geom, years = years, collection = lc_coll,
                         initiative = lc_ini,
                         backend = s$backend, src = src, by = by,
                         verbose = verbose, max_pixels = max_pixels)
  attr(ts, "species") <- species
  attr(ts, "range") <- toupper(range)
  ts
}

#' Stacked-area chart of land-cover percentage over time
#'
#' Draws percentage of area versus year as stacked areas, one band per class (or
#' group), using the official MapBiomas colours, in the style of the MapBiomas
#' coverage figures. Uses \pkg{ggplot2} when available, otherwise base graphics.
#'
#' @param ts A long data frame from \code{\link{cover_timeseries}} or
#'   \code{\link{timeseries_for_species}}.
#' @param title Optional plot title. If \code{NULL}, a title is built from the
#'   \code{species}/\code{range} attributes when present.
#' @param legend Logical; draw the class legend (default \code{TRUE}).
#' @param lang Legend language: \code{"en"} (default) or \code{"pt"}. When
#'   \code{"en"}, the band labels use the MapBiomas \code{class_en} column (when
#'   present in \code{ts}) and the axis/legend/title text is in English.
#' @return A \pkg{ggplot} object (invisibly) when ggplot2 is available; otherwise
#'   \code{NULL} after drawing a base plot.
#' @export
plot_timeseries <- function(ts, title = NULL, legend = TRUE,
                            lang = c("en", "pt")) {
  lang <- match.arg(lang)
  stopifnot(is.data.frame(ts),
            all(c("year", "label", "hex", "pct") %in% names(ts)))

  # Swap the (Portuguese) `label` for the English MapBiomas class names when
  # requested and the `class_en` column is available (by = "class" carries it;
  # by = "group" does not, so group labels are translated further below).
  if (lang == "en" && "class_en" %in% names(ts)) ts$label <- ts$class_en

  # UI strings by language
  lab <- if (lang == "en") {
    list(x = "Year", y = "Percentage (%)",
         fill = "Land-cover classes", legend = "Classes",
         title = "land cover over time")
  } else {
    list(x = "Ano", y = "Porcentagem (%)",
         fill = "Classes de cobertura do solo", legend = "Classes",
         title = "cobertura ao longo do tempo")
  }

  # guarantee a rectangular series (every label in every year) so a class that
  # is absent in some years cannot break ggplot2::geom_area()
  yrs_u <- sort(unique(ts$year)); labs_u <- unique(ts$label)
  if (nrow(ts) < length(yrs_u) * length(labs_u)) {
    sp_a <- attr(ts, "species"); rg_a <- attr(ts, "range")
    vals <- ts[, c("year", "label", "pct")]
    keep <- setdiff(names(ts), c("year", "pct"))
    meta <- ts[!duplicated(ts$label), keep, drop = FALSE]
    full <- expand.grid(year = yrs_u, label = labs_u,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ts <- merge(full, meta, by = "label", all.x = TRUE)
    ts <- merge(ts, vals, by = c("year", "label"), all.x = TRUE)
    ts$pct[is.na(ts$pct)] <- 0
    attr(ts, "species") <- sp_a; attr(ts, "range") <- rg_a
  }

  title_html <- title
  if (is.null(title)) {
    sp <- attr(ts, "species"); rg <- attr(ts, "range")
    if (!is.null(sp)) {
      suffix <- paste0(if (!is.null(rg)) paste0(" (", rg, ")") else "",
                       " - ", lab$title)
      title      <- .sp_title_expr(sp, suffix)   # ggplot (plotmath)
      title_html <- .sp_title_html(sp, suffix)   # plotly / HTML
    }
  }

  # canonical ordering: by MapBiomas hierarchy when codes are present
  if ("code" %in% names(ts)) {
    leg <- mb_legend()
    key <- match(ts$code, leg$code)
    ord <- order(key)
    lvls <- unique(ts$label[ord])
  } else {
    # group series (no code): relabel from the internal `group` key via the
    # single-source helper, so PT and EN never drift apart.
    gl <- .mb_group_labels(lang)
    if ("group" %in% names(ts)) ts$label <- unname(gl$labels[ts$group])
    pref <- unname(gl$labels[c("natural", "anthropic", "water", "other", "not_observed")])
    lvls <- c(intersect(pref, unique(ts$label)), setdiff(unique(ts$label), pref))
  }
  cols <- vapply(lvls, function(L) ts$hex[ts$label == L][1], character(1))
  ny <- length(unique(ts$year))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ts$label <- factor(ts$label, levels = lvls)
    layer <- if (ny < 2) {
      ggplot2::geom_col(position = "stack", width = 0.6)
    } else {
      ggplot2::geom_area(position = "stack", colour = NA)
    }
    p <- ggplot2::ggplot(
      ts, ggplot2::aes(x = .data[["year"]], y = .data[["pct"]],
                       fill = .data[["label"]])) +
      layer +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE,
                                 name = lab$fill) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = lab$x, y = lab$y, title = title) +
      .mas_theme()
    if (ny >= 2) p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
    if (!legend) p <- p + ggplot2::theme(legend.position = "none")
    attr(p, "mas") <- list(
      kind = "area",
      df = data.frame(year = ts$year, label = as.character(ts$label),
                      pct = ts$pct, stringsAsFactors = FALSE),
      palette = as.list(cols), levels = lvls, ny = ny,
      xlab = lab$x, ylab = lab$y, title = title_html)
    return(p)
  }

  # ---- base-graphics fallback (stacked polygons / bar for a single year) ----
  yrs <- sort(unique(ts$year))
  M <- matrix(0, nrow = length(yrs), ncol = length(lvls),
              dimnames = list(as.character(yrs), lvls))
  for (i in seq_len(nrow(ts))) {
    M[as.character(ts$year[i]), ts$label[i]] <- ts$pct[i]
  }

  op <- graphics::par(mar = c(4, 4, if (is.null(title)) 1 else 3,
                              if (legend) 11 else 1), xpd = NA)
  on.exit(graphics::par(op), add = TRUE)

  if (ny < 2) {
    graphics::barplot(matrix(M[1, ], ncol = 1), col = cols, border = "white",
                      ylim = c(0, 100), names.arg = rownames(M)[1],
                      ylab = lab$y, xlab = lab$x, main = title)
  } else {
    cum <- if (ncol(M) == 1) M else t(apply(M, 1, cumsum))
    plot(range(yrs), c(0, 100), type = "n", xlab = lab$x,
         ylab = lab$y, main = title)
    for (j in seq_along(lvls)) {
      lower <- if (j == 1) rep(0, length(yrs)) else cum[, j - 1]
      upper <- cum[, j]
      graphics::polygon(c(yrs, rev(yrs)), c(lower, rev(upper)),
                        col = cols[j], border = NA)
    }
  }
  if (legend) {
    graphics::legend("topright", inset = c(-0.32, 0), legend = lvls, fill = cols,
                     border = NA, bty = "n", cex = 0.75, title = lab$legend,
                     xpd = NA)
  }
  invisible(NULL)
}

#' Regression trendline of one land-cover class through time
#'
#' Fits a regression of a single land-cover class (or conservation group)
#' percentage of area against year and plots the observed points with the
#' fitted curve, its equation, R-squared and p-value. When the
#' \pkg{ggtrendline} package is installed
#' (\url{https://CRAN.R-project.org/package=ggtrendline}) it is used to draw the
#' trendline and annotate the statistics; otherwise the function falls back to a
#' plain \pkg{ggplot2} linear fit carrying the same statistics, and finally to
#' base graphics.
#'
#' @param ts A long data frame from \code{\link{cover_timeseries}} or
#'   \code{\link{timeseries_for_species}}.
#' @param class_label The \code{label} (class or group name) to analyse. When
#'   \code{NULL} (default), the class with the largest mean percentage across
#'   the series is used.
#' @param model \pkg{ggtrendline} model, one of \code{"line2P"} (linear,
#'   default), \code{"line3P"} (quadratic), \code{"log2P"}, \code{"exp2P"},
#'   \code{"exp3P"}, \code{"power2P"} or \code{"power3P"}. Ignored (a linear fit
#'   is used) when \pkg{ggtrendline} is not installed.
#' @param title Optional plot title. When \code{NULL}, a title is built from the
#'   class name and the \code{species}/\code{range} attributes when present.
#' @param lang Axis/title language: \code{"en"} (default) or \code{"pt"}. When
#'   \code{"en"} the English MapBiomas class names (\code{class_en}) are used
#'   when available.
#' @return A \pkg{ggplot} object when \pkg{ggplot2} is available; otherwise
#'   \code{NULL} (invisibly) after drawing a base plot.
#' @seealso \code{\link{plot_timeseries}}, \code{\link{timeseries_for_species}}
#' @export
plot_class_trendline <- function(ts, class_label = NULL, model = "line2P",
                                 title = NULL, lang = c("en", "pt")) {
  lang <- match.arg(lang)
  stopifnot(is.data.frame(ts),
            all(c("year", "label", "pct") %in% names(ts)))

  # Display labels follow plot_timeseries: English MapBiomas class names when
  # requested and available, otherwise the native (Portuguese) labels.
  disp <- if (lang == "en" && "class_en" %in% names(ts)) ts$class_en else ts$label

  # Resolve the requested class against the display name, the native label or
  # the English name, so a caller may pass any of them; default to the class
  # with the largest mean share.
  resolve <- function(x) {
    if (is.null(x)) return(NULL)
    if (x %in% disp) return(x)
    if (x %in% ts$label) return(disp[match(x, ts$label)])
    if ("class_en" %in% names(ts) && x %in% ts$class_en)
      return(disp[match(x, ts$class_en)])
    NULL
  }
  class_label <- resolve(class_label)
  if (is.null(class_label)) {
    mp <- tapply(ts$pct, disp, mean, na.rm = TRUE)
    class_label <- names(mp)[which.max(mp)]
  }

  keep <- disp == class_label
  d <- ts[keep, c("year", "pct"), drop = FALSE]
  d <- d[is.finite(d$year) & is.finite(d$pct), , drop = FALSE]
  d <- d[order(d$year), , drop = FALSE]

  lab <- if (lang == "en")
    list(x = "Year", y = "Percentage (%)", t = "trend")
  else
    list(x = "Ano", y = "Porcentagem (%)", t = "tendencia")

  hexc <- {
    h <- if ("hex" %in% names(ts)) ts$hex[keep] else NULL
    if (length(h)) h[1] else "#1f8d49"
  }

  ttl <- title
  if (is.null(ttl)) {
    sp <- attr(ts, "species"); rg <- attr(ts, "range")
    ttl <- paste0(class_label,
                  if (!is.null(rg)) paste0(" (", rg, ")") else "",
                  " - ", lab$t,
                  if (!is.null(sp)) paste0(" - ", sp) else "")
  }

  if (!nrow(d)) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      return(ggplot2::ggplot() + ggplot2::labs(title = ttl) + .mas_theme())
    }
    return(invisible(NULL))
  }

  # Preferred path: ggtrendline (equation + R2 + p-value baked into the plot).
  if (nrow(d) >= 3 &&
      requireNamespace("ggtrendline", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    p <- tryCatch(
      ggtrendline::ggtrendline(d$year, d$pct, model = model,
                               linecolor = "#233d2c", pointcolor = hexc) +
        ggplot2::labs(x = lab$x, y = lab$y, title = ttl) +
        .mas_theme(),
      error = function(e) NULL)
    if (!is.null(p)) return(p)
  }

  # Fallback: ggplot2 linear fit with lm statistics annotated.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["year"]],
                                         y = .data[["pct"]])) +
      ggplot2::geom_point(colour = hexc, size = 2)
    if (nrow(d) >= 2) {
      fit <- stats::lm(pct ~ year, data = d)
      cf  <- stats::coef(fit)
      sm  <- suppressWarnings(summary(fit))
      pv  <- tryCatch(
        stats::pf(sm$fstatistic[1L], sm$fstatistic[2L], sm$fstatistic[3L],
                  lower.tail = FALSE), error = function(e) NA_real_)
      eqn <- sprintf("y = %.4g %s %.4g x\nR^2 = %.3f%s",
                     cf[[1L]], if (cf[[2L]] >= 0) "+" else "-", abs(cf[[2L]]),
                     sm$r.squared,
                     if (is.finite(pv)) sprintf(", p = %.3g", pv) else "")
      p <- p +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                             colour = "#233d2c", fill = "#233d2c",
                             alpha = 0.15) +
        ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.05,
                          vjust = 1.3, label = eqn, size = 3.4)
    }
    return(p + ggplot2::labs(x = lab$x, y = lab$y, title = ttl) + .mas_theme())
  }

  # Base-graphics fallback.
  plot(d$year, d$pct, xlab = lab$x, ylab = lab$y, main = ttl, pch = 19,
       col = hexc)
  if (nrow(d) >= 2) graphics::abline(stats::lm(pct ~ year, data = d))
  invisible(NULL)
}
