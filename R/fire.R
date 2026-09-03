## R/fire.R — MapBiomas Fire (Fogo) integration for mappingAS
## Mirrors the LULC functions: read a windowed crop, tabulate, summarise,
## time series. Reuses mb_raster_local() (GDAL /vsicurl/ + on-disk cache).

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' Build a MapBiomas Fire (Fogo) GeoTIFF URL
#'
#' \code{"frequency"} is the fire-frequency layer (pixel = number of years
#' burned over \code{period}, 1-40); \code{"accumulated"} is binary (burned at
#' least once); \code{"annual"} is the per-year burned area.
#'
#' @param year Integer year (ignored when \code{product = "accumulated"}).
#' @param product \code{"accumulated"} (default) or \code{"annual"}.
#' @param fire_collection Fire collection number (default \code{4}).
#' @param host_collection Initiative folder hosting the fire data
#'   (default \code{9}; Fire Col.4 lives under \code{collection_9/}).
#' @param period Two-integer accumulated span (default \code{c(1985, 2024)}).
#' @return A length-1 character URL.
#' @examples
#' mb_fire_url(product = "accumulated")
#' mb_fire_url(2024, "annual")
#' @export
mb_fire_url <- function(year = 2024,
                        product = c("accumulated", "annual", "frequency"),
                        fire_collection = 4, host_collection = 9,
                        period = c(1985, 2024)) {
  product <- match.arg(product)
  base <- sprintf(paste0("https://storage.googleapis.com/mapbiomas-public/",
                         "initiatives/brasil/collection_%d/fire-col%d"),
                  as.integer(host_collection), as.integer(fire_collection))
  fc <- as.integer(fire_collection)
  p1 <- as.integer(period[1]); p2 <- as.integer(period[2])
  switch(product,
    annual      = sprintf("%s/annual_burned/mapbiomas_fire_col%d_br_annual_burned_%d.tif",
                          base, fc, as.integer(year)),
    frequency   = sprintf("%s/fire_frequency_v1/mapbiomas_fire_col%d_br_fire_frequency_%d_%d.tif",
                          base, fc, p1, p2),
    accumulated = sprintf("%s/accumulated_burned_v1/mapbiomas_fire_col%d_br_accumulated_burned_%d_%d.tif",
                          base, fc, p1, p2))
}

#' Crop a MapBiomas Fire raster to an area of interest (local backend)
#'
#' Thin wrapper over \code{\link{mb_raster_local}}; reuses its windowed
#' \code{/vsicurl/} read and on-disk cache.
#'
#' @inheritParams mb_fire_url
#' @param aoi An \code{sf}/\code{sfc} polygon (e.g. an EOO hull or AOO union).
#' @param src Optional GeoTIFF path/URL overriding the built fire URL.
#' @param mask,cache,cache_dir Passed to \code{\link{mb_raster_local}}.
#' @return A \pkg{terra} \code{SpatRaster} of fire values clipped to the AOI.
#' @export
fire_raster_local <- function(aoi, year = 2024,
                              product = c("accumulated", "annual", "frequency"),
                              fire_collection = 4, host_collection = 9,
                              src = NULL, mask = TRUE, cache = TRUE,
                              cache_dir = NULL) {
  product <- match.arg(product)
  src <- src %||% mb_fire_url(year, product, fire_collection, host_collection)
  r <- mb_raster_local(aoi, src = src, mask = mask, cache = cache,
                       cache_dir = cache_dir)
  names(r) <- "fire"
  r
}

#' Area (km^2) per fire value from a cropped fire raster
#'
#' Fire analogue of \code{\link{mb_class_areas_raster}}. For the accumulated
#' product, \code{freq} is the fire frequency (years burned).
#'
#' @param r A \pkg{terra} \code{SpatRaster} from \code{\link{fire_raster_local}}.
#' @return A \code{data.frame} with columns \code{freq} and \code{area_km2}.
#' @export
fire_areas <- function(r) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)
  a <- terra::cellSize(r, unit = "km", mask = TRUE)
  z <- terra::zonal(a, r, fun = "sum", na.rm = TRUE)
  names(z) <- c("freq", "area_km2")
  z$freq <- as.integer(z$freq)
  z[is.finite(z$area_km2) & z$area_km2 > 0, , drop = FALSE]
}

#' Summarise fire pressure inside a range
#'
#' Fire analogue of \code{\link{summarise_conversion}}. Pass the known range
#' area (\code{range_km2}, e.g. \code{eoo$area_km2}) as denominator so the
#' percentage is exact and robust to the raster's NoData handling.
#'
#' @param fire_areas A \code{data.frame} from \code{\link{fire_areas}}.
#' @param range_km2 Total range area (km^2). If \code{NULL}, the sum of all
#'   tabulated cells is used (valid only when unburned pixels are coded 0).
#' @return A list with \code{burned_km2}, \code{total_km2}, \code{burned_pct},
#'   \code{freq_mean}, \code{freq_max} and \code{by_freq}.
#' @export
summarise_fire <- function(fire_areas, range_km2 = NULL) {
  stopifnot(all(c("freq", "area_km2") %in% names(fire_areas)))
  bf <- fire_areas[fire_areas$freq > 0, , drop = FALSE]
  burned <- sum(bf$area_km2)
  total  <- range_km2 %||% sum(fire_areas$area_km2)
  w <- if (nrow(bf)) bf$area_km2 / sum(bf$area_km2) else numeric(0)
  list(
    burned_km2 = burned,
    total_km2  = total,
    burned_pct = if (total > 0) 100 * burned / total else NA_real_,
    freq_mean  = if (nrow(bf)) sum(bf$freq * w) else NA_real_,
    freq_max   = if (nrow(bf)) max(bf$freq) else NA_real_,
    by_freq    = bf[order(bf$freq), , drop = FALSE]
  )
}

#' @keywords internal
#' @noRd
## Accumulated-fire summary for one geometry (used by assess_species()).
.fire_for <- function(geom, range_km2, fire_collection, host_collection,
                      src, verbose, label) {
  if (is.null(geom) || length(geom) == 0 || all(sf::st_is_empty(geom))) {
    if (verbose) message(sprintf("  %s: no polygon; fire = NA.", label))
    return(NULL)
  }
  tryCatch({
    r <- fire_raster_local(geom, product = "accumulated", src = src,
                           fire_collection = fire_collection,
                           host_collection = host_collection)
    summarise_fire(fire_areas(r), range_km2 = range_km2)
  }, error = function(e) {
    warning(sprintf("%s fire failed: %s", label, conditionMessage(e)),
            call. = FALSE)
    NULL
  })
}

#' @keywords internal
#' @noRd
## Downsampled accumulated-fire raster for display (burned pixels only).
.fire_raster_display <- function(aoi, fire_collection = 4, host_collection = 9,
                                 src = NULL, max_pixels = 800, crs = NULL,
                                 cache = TRUE, product = "frequency") {
  r <- fire_raster_local(aoi, product = product, src = src,
                         fire_collection = fire_collection,
                         host_collection = host_collection, cache = cache)
  r[r <= 0] <- NA                       # drop unburned (freq/accum coded 0)
  d <- max(dim(r)[1:2])
  if (d > max_pixels)
    r <- terra::aggregate(r, fact = ceiling(d / max_pixels),
                          fun = "max", na.rm = TRUE)
  if (!is.null(crs)) r <- terra::project(r, crs, method = "near")
  names(r) <- "fire"
  r
}

#' Burned-area time series for a polygon (% of range x year)
#'
#' Fire analogue of \code{\link{cover_timeseries}}: reads the annual burned
#' layer for each year and returns burned area (km^2) and percentage of the
#' range. Each year is one windowed read, so a full annual series can be slow.
#'
#' @param geom An \code{sf}/\code{sfc} polygon.
#' @param years Integer vector (default \code{\link{mb_years}()}).
#' @param range_km2 Range area for the percentage denominator (optional).
#' @inheritParams mb_fire_url
#' @param verbose Logical; print progress.
#' @return A \code{data.frame} with \code{year}, \code{burned_km2},
#'   \code{burned_pct}.
#' @export
fire_timeseries <- function(geom, years = NULL, fire_collection = 4,
                            host_collection = 9, range_km2 = NULL,
                            verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)
  # Planar (GEOS) geometry, as in assess_species(): the S2 engine rejects a
  # densified, wide-ranging hull ("Edge is degenerate") during the reads.
  if (requireNamespace("sf", quietly = TRUE)) {
    .old_s2 <- suppressMessages(sf::sf_use_s2())
    on.exit(suppressMessages(sf::sf_use_s2(.old_s2)), add = TRUE)
    suppressMessages(sf::sf_use_s2(FALSE))
  }
  years <- years %||% mb_years()
  rows <- lapply(years, function(y) {
    # Free each year's raster memory on exit so the peak stays at one year
    # (avoids the OS killing the R process on memory-limited machines).
    on.exit(invisible(gc(FALSE)), add = TRUE)
    if (verbose) message("  fire ", y)
    r <- tryCatch(
      fire_raster_local(geom, y, "annual", fire_collection, host_collection),
      error = function(e) NULL)
    if (is.null(r)) return(NULL)
    a <- terra::cellSize(r, unit = "km", mask = TRUE)
    burned <- as.numeric(terra::global(a * (r > 0), "sum", na.rm = TRUE)[1, 1])
    total  <- range_km2 %||% as.numeric(terra::global(a, "sum", na.rm = TRUE)[1, 1])
    data.frame(year = y, burned_km2 = round(burned, 4),
               burned_pct = if (total > 0) round(100 * burned / total, 3) else NA_real_)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out[order(out$year), , drop = FALSE]
}

#' Burned-area time series for one assessed species (EOO or AOO)
#'
#' Convenience wrapper around \code{\link{fire_timeseries}} that pulls the
#' stored EOO hull or AOO geometry (and its area) from an assessment.
#'
#' @param assessment A \code{geoconv_assessment} from \code{\link{assess_species}}.
#' @param species Species name (default: first assessed).
#' @param range \code{"eoo"} (default) or \code{"aoo"}.
#' @param years,verbose Passed to \code{\link{fire_timeseries}}.
#' @return The data frame from \code{\link{fire_timeseries}} with attributes
#'   \code{species} and \code{range}.
#' @export
fire_timeseries_for_species <- function(assessment, species = NULL,
                                        range = c("eoo", "aoo"), years = NULL,
                                        verbose = TRUE) {
  stopifnot(inherits(assessment, "geoconv_assessment"))
  range <- match.arg(range)
  d <- assessment$detail
  if (is.null(species)) species <- names(d)[1]
  obj <- d[[species]]
  if (is.null(obj)) stop("Species not found: ", species, call. = FALSE)
  geom <- if (range == "eoo") obj$eoo$hull else .st_union_quiet(obj$aoo$cells)
  rng_km2 <- if (range == "eoo") obj$eoo$area_km2 else obj$aoo$area_km2
  if (is.null(geom) || length(sf::st_geometry(geom)) == 0)
    stop(sprintf("No %s geometry for '%s'.", toupper(range), species), call. = FALSE)
  s <- assessment$settings
  ts <- fire_timeseries(geom, years = years,
                        fire_collection = s$fire_collection %||% 4,
                        range_km2 = rng_km2, verbose = verbose)
  attr(ts, "species") <- species
  attr(ts, "range") <- toupper(range)
  ts
}

#' Sequential palette for fire frequency (YlOrRd)
#'
#' @param n Number of breaks (default \code{6}).
#' @return A character vector of hex colours.
#' @export
fire_palette <- function(n = 6) {
  grDevices::colorRampPalette(
    c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"))(n)
}

#' @keywords internal
#' @noRd
## Discrete fire-frequency classes (number of years burned, 1985-2024).
.fire_freq_bins <- function() {
  data.frame(lower = c(1, 2, 4, 7, 11, 21),
             upper = c(1, 3, 6, 10, 20, 40),
             label = c("1", "2-3", "4-6", "7-10", "11-20", "21-40"),
             hex   = fire_palette(6),
             stringsAsFactors = FALSE)
}

#' Bar/line chart of burned area (%) over time
#'
#' @param ts A data frame from \code{\link{fire_timeseries}} /
#'   \code{\link{fire_timeseries_for_species}}.
#' @param title Optional title (built from attributes when \code{NULL}).
#' @param lang Label language: \code{"en"} (default) or \code{"pt"}.
#' @return A \pkg{ggplot} object when ggplot2 is available; otherwise \code{NULL}
#'   after a base plot.
#' @export
plot_fire_timeseries <- function(ts, title = NULL, lang = c("en", "pt")) {
  lang <- match.arg(lang)
  stopifnot(is.data.frame(ts), all(c("year", "burned_pct") %in% names(ts)))
  ylab      <- if (lang == "en") "% of area burned" else "% da area queimada"
  suffix_fmt <- if (lang == "en") "%s - area burned per year (MapBiomas Fire)"
                else "%s - area queimada por ano (MapBiomas Fogo)"
  title_html <- title
  if (is.null(title)) {
    sp <- attr(ts, "species"); rg <- attr(ts, "range")
    if (!is.null(sp)) {
      suffix <- sprintf(suffix_fmt,
                        if (!is.null(rg)) paste0(" (", rg, ")") else "")
      title      <- .sp_title_expr(sp, suffix)   # ggplot (plotmath)
      title_html <- .sp_title_html(sp, suffix)   # plotly / HTML
    }
  }
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(ts, ggplot2::aes(x = .data[["year"]],
                                     y = .data[["burned_pct"]])) +
      ggplot2::geom_col(fill = "#fd8d3c") +
      ggplot2::geom_line(colour = "#bd0026", linewidth = 0.5) +
      ggplot2::labs(x = NULL, y = ylab, title = title) +
      .mas_theme()
    attr(p, "mas") <- list(
      kind = "fire",
      df = data.frame(year = ts$year, burned_pct = ts$burned_pct),
      ylab = ylab, title = title_html)
    p
  } else {
    graphics::barplot(ts$burned_pct, names.arg = ts$year, col = "#fd8d3c",
                      ylab = ylab, main = title %||% "")
    invisible(NULL)
  }
}
