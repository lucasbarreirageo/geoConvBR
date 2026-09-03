#' Build the MapBiomas land-use GeoTIFF URL for a given year and initiative
#'
#' Returns the public Google Cloud Storage URL of the MapBiomas
#' land-use/land-cover annual mosaic for the chosen initiative. All supported
#' initiatives publish annual single-year Cloud-Optimized GeoTIFFs on the same
#' public bucket, so they can be streamed with GDAL's \code{/vsicurl/} driver
#' (no Google Earth Engine, no Google Drive). See \code{\link{mb_initiatives}}
#' for the full list of initiatives, their default collections and year spans.
#'
#' @param year Integer year.
#' @param collection Integer collection number. \code{NULL} (default) uses the
#'   initiative's default collection.
#' @param initiative One of the keys of \code{\link{mb_initiatives}}
#'   (\code{"brazil"} by default).
#' @return A length-1 character URL.
#' @examples
#' mb_source_url(2024)
#' mb_source_url(2023, initiative = "amazonia")
#' mb_source_url(2024, initiative = "peru")
#' @export
mb_source_url <- function(year = 2024, collection = NULL, initiative = "brazil") {
  ini <- .mb_resolve_initiative(initiative)
  if (identical(ini$provider, "esri")) {
    stop("mb_source_url() builds MapBiomas URLs; the global Sentinel-2/Esri ",
         "layer is tiled - use s2_source_url() instead.", call. = FALSE)
  }
  if (is.null(collection)) collection <- ini$collection
  collection <- as.integer(collection)
  year <- as.integer(year)
  sprintf("https://storage.googleapis.com/mapbiomas-public/initiatives/%s",
          ini$build(year, collection))
}

#' @keywords internal
#' @noRd
.mb_set_gdal <- function() {
  # Read tuning for /vsicurl/ streaming. These only affect I/O speed (bigger
  # block cache, larger chunked range requests, HTTP/2 multiplexing, threaded
  # decompression) - never the pixel values returned, so area statistics are
  # unchanged.
  keys <- c("GDAL_DISABLE_READDIR_ON_OPEN", "CPL_VSIL_CURL_USE_HEAD",
            "GDAL_HTTP_MAX_RETRY", "GDAL_HTTP_RETRY_DELAY", "VSI_CACHE",
            "GDAL_CACHEMAX", "GDAL_NUM_THREADS", "GDAL_HTTP_MULTIPLEX",
            "GDAL_HTTP_VERSION", "CPL_VSIL_CURL_CHUNK_SIZE",
            "CPL_VSIL_CURL_CACHE_SIZE")
  # GDAL_CACHEMAX and the /vsicurl/ cache are kept modest (128 MB / 64 MB) so a
  # windowed read does not reserve hundreds of MB of block cache - important on
  # memory-limited machines/containers where an over-large cache plus the raster
  # itself gets the process OOM-killed ("Terminated"). Values, not just speed.
  vals <- c("EMPTY_DIR", "NO", "3", "1", "TRUE",
            "128", "ALL_CPUS", "YES",
            "2", "1048576",
            "67108864")
  old <- stats::setNames(lapply(keys, terra::getGDALconfig), keys)
  for (i in seq_along(keys)) terra::setGDALconfig(keys[i], vals[i])
  old
}

#' Crop a MapBiomas LULC raster to an area of interest (local backend)
#'
#' Reads the MapBiomas national mosaic and crops/masks it to \code{aoi} using
#' \pkg{terra}. By default it streams a \emph{windowed} read of only the AOI's
#' bounding box from the public GeoTIFF via GDAL's \code{/vsicurl/} driver, so
#' there is no Google Earth Engine account and no full-country download. For
#' offline or repeated use, point \code{src} to a local GeoTIFF instead.
#'
#' For very large ranges (e.g. continental EOOs) the windowed read can still
#' transfer a lot of data; in that case prefer \code{\link{mb_class_areas_gee}}.
#'
#' @param aoi An \code{sf}/\code{sfc} polygon (any CRS) defining the area to
#'   extract, e.g. an EOO hull or the union of AOO cells.
#' @param year Integer year (default \code{2024}).
#' @param collection Integer collection number. \code{NULL} uses the
#'   \code{initiative} default (10/6/3 for Brazil/Amazonia/Colombia).
#' @param initiative One of \code{"brazil"} (default), \code{"amazonia"} or
#'   \code{"colombia"} (see \code{\link{mb_initiatives}}).
#' @param src Optional path or URL to a MapBiomas GeoTIFF. If \code{NULL}
#'   (default) the public URL for \code{year}/\code{collection}/\code{initiative}
#'   is used through \code{/vsicurl/}.
#' @param mask Logical; if \code{TRUE} (default) pixels outside the polygon are
#'   set to \code{NA}. If \code{FALSE} only a rectangular crop is returned.
#' @param cache Logical; if \code{TRUE} (default) the rectangular windowed crop
#'   is cached on disk (keyed by source + bounding box), so reading the same
#'   area again (e.g. the EOO during assessment and later when mapping, or a
#'   re-run) does not re-download. Masking is always applied in memory after the
#'   cached crop is loaded.
#' @param cache_dir Directory for the windowed-crop cache (default a
#'   \code{mappingAS_mb_cache} folder under \code{tempdir()}).
#' @return A \pkg{terra} \code{SpatRaster} of MapBiomas pixel codes restricted to
#'   the AOI, in the raster's native CRS.
#' @export
mb_raster_local <- function(aoi, year = 2024, collection = NULL,
                            initiative = "brazil",
                            src = NULL, mask = TRUE,
                            cache = TRUE, cache_dir = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.", call. = FALSE)
  }
  ini <- .mb_resolve_initiative(initiative)
  # Global Sentinel-2 / Esri fallback: read the tiled global mosaic instead.
  if (identical(ini$provider, "esri")) {
    return(s2_raster_local(aoi, year = year, src = src, mask = mask,
                           cache = cache, cache_dir = cache_dir))
  }
  if (is.null(collection)) collection <- ini$collection
  aoi <- sf::st_geometry(aoi)
  if (is.na(sf::st_crs(aoi))) sf::st_crs(aoi) <- 4326

  # robust remote reads; restore the user's GDAL settings on exit
  old <- .mb_set_gdal()
  on.exit(for (k in names(old)) terra::setGDALconfig(k, old[[k]]), add = TRUE)

  if (is.null(src)) src <- mb_source_url(year, collection, ini$key)
  path <- if (grepl("^https?://", src)) paste0("/vsicurl/", src) else src

  # suppressWarnings: a failed open also emits a GDAL "no such file" *warning*
  # alongside the error; we surface a single, clear error below, so the extra
  # warning is only noise.
  r <- tryCatch(
    suppressWarnings(terra::rast(path)),
    error = function(e) {
      stop("Could not open MapBiomas raster at:\n  ", src,
           "\nIf you are offline or the host is blocked, download the GeoTIFF ",
           "and pass it via `src=`.\nOriginal error: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  aoi_r <- sf::st_transform(aoi, terra::crs(r))
  vect <- terra::vect(aoi_r)
  rc <- .mb_window(r, terra::ext(vect), src, year, collection, cache, cache_dir)
  if (mask) rc <- terra::mask(rc, vect)
  names(rc) <- "mapbiomas_class"
  rc
}

#' @keywords internal
#' @noRd
.mb_window <- function(r, e, src, year, collection, cache, cache_dir) {
  rc <- terra::crop(r, e, snap = "out")          # lazy windowed read
  if (!isTRUE(cache)) return(rc)

  if (is.null(cache_dir)) cache_dir <- file.path(tempdir(), "mappingAS_mb_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  key <- paste(basename(src), year, collection,
               paste(round(as.vector(e), 5), collapse = "_"), sep = "__")
  f <- file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_]+", "-", key), ".tif"))
  if (file.exists(f)) return(terra::rast(f))

  terra::writeRaster(rc, f, datatype = "INT1U", overwrite = TRUE,
                     gdal = "COMPRESS=LZW")       # materialise the read once
  terra::rast(f)
}

# Downsampled, optionally reprojected raster for *display* (not area calc).
# Tries a fast decimated read first (GDAL -outsize uses the COG overviews, so a
# large window is not streamed at native 30 m); on any failure it falls back to
# the native windowed read + aggregate, so behaviour is never lost.
#' @keywords internal
#' @noRd
.mb_raster_display <- function(aoi, year, collection, src, max_pixels = 600,
                               crs = NULL, cache = TRUE, initiative = "brazil") {
  if (.s2_is_esri(initiative)) {
    return(.s2_raster_display(aoi, year = year, src = src,
                              max_pixels = max_pixels, crs = crs, cache = cache))
  }
  r <- tryCatch(
    .mb_display_read(aoi, year, collection, src, max_pixels, initiative),
    error = function(e) NULL)
  if (is.null(r)) {
    r <- mb_raster_local(aoi, year = year, collection = collection,
                         initiative = initiative, src = src,
                         mask = TRUE, cache = cache)
    d <- max(dim(r)[1:2])
    if (d > max_pixels)
      r <- terra::aggregate(r, fact = ceiling(d / max_pixels),
                            fun = "modal", na.rm = TRUE)
  }
  if (!is.null(crs)) r <- terra::project(r, crs, method = "near")
  names(r) <- "mapbiomas_class"
  r
}

# Fast, display-only read: GDAL decimates on read via `-outsize` (using the
# raster overviews when present), avoiding a native-resolution stream over a
# large window. Used only for the map overlay - never for area statistics.
#' @keywords internal
#' @noRd
.mb_display_read <- function(aoi, year, collection, src, max_pixels = 600,
                             initiative = "brazil") {
  if (!requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE))
    stop("Packages 'terra' and 'sf' are required.", call. = FALSE)
  ini <- .mb_resolve_initiative(initiative)
  if (is.null(collection)) collection <- ini$collection
  aoi <- sf::st_geometry(aoi)
  if (is.na(sf::st_crs(aoi))) sf::st_crs(aoi) <- 4326
  if (is.null(src)) src <- mb_source_url(year, collection, ini$key)
  path <- if (grepl("^https?://", src)) paste0("/vsicurl/", src) else src

  old <- .mb_set_gdal()
  on.exit(for (k in names(old)) terra::setGDALconfig(k, old[[k]]), add = TRUE)

  r0    <- terra::rast(path)                       # metadata only (no pixels)
  aoi_r <- sf::st_transform(aoi, terra::crs(r0))
  bb    <- as.numeric(sf::st_bbox(aoi_r))          # xmin ymin xmax ymax
  # cap the output width at the window's native width to avoid upsampling
  win_w <- max(1L, as.integer(ceiling((bb[3] - bb[1]) / terra::xres(r0))))
  out_w <- min(as.integer(max_pixels), win_w)

  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)
  sf::gdal_utils(
    "translate", source = path, destination = tmp,
    options = c("-projwin", sprintf("%.10f", bb[1]), sprintf("%.10f", bb[4]),
                sprintf("%.10f", bb[3]), sprintf("%.10f", bb[2]),
                "-outsize", as.character(out_w), "0", "-r", "nearest"))
  rc <- terra::rast(tmp)
  rc <- terra::mask(rc, terra::vect(aoi_r)) + 0L   # detach from the temp file
  names(rc) <- "mapbiomas_class"
  rc
}

#' Tabulate area per MapBiomas class from a cropped raster
#'
#' Computes the true (geodesic, latitude-corrected) area of each MapBiomas pixel
#' class inside a cropped/masked raster. Works for rasters in geographic CRS,
#' where pixel area varies with latitude, by using \code{terra::cellSize()}.
#'
#' @param r A \pkg{terra} \code{SpatRaster} of MapBiomas codes (e.g. from
#'   \code{\link{mb_raster_local}}).
#' @return A \code{data.frame} with columns \code{code} and \code{area_km2}.
#' @export
mb_class_areas_raster <- function(r) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.", call. = FALSE)
  }
  area <- terra::cellSize(r, unit = "km", mask = TRUE)
  z <- terra::zonal(area, r, fun = "sum", na.rm = TRUE)
  names(z) <- c("code", "area_km2")
  z$code <- as.integer(z$code)
  z[is.finite(z$area_km2) & z$area_km2 > 0, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Budget-aware per-class area reader (local backend)
# ---------------------------------------------------------------------------
# A native-resolution read of the bounding box of a very large or widely-spread
# range can reach billions of pixels (a continental EOO at 30 m), exhausting
# memory/time. These helpers keep the read bounded:
#   * small ranges  -> a single native read (identical results to before);
#   * sparse ranges -> one native read per polygon part, summed (e.g. the AOO
#     cells, which cover little true area inside a huge bounding box), so the
#     empty space between parts is never streamed and resolution is preserved;
#   * one large part -> a decimated read that samples via the GeoTIFF overviews
#     (nearest-neighbour, so only the decimated pixels are read - memory stays
#     small), keeping the total within the budget.
# Sampling millions of pixels preserves the class proportions - and therefore
# the % converted / natural - so the headline Criterion B screening metrics are
# unchanged in practice.

#' Default pixel budget for the area (statistics) read: ~50 million pixels.
#' @keywords internal
#' @noRd
.area_max_pixels <- function() 5e7

#' Native pixel width/height/count of an extent (xmin,xmax,ymin,ymax) at a
#' raster's resolution.
#' @keywords internal
#' @noRd
.bbox_native_px <- function(r, ext) {
  w <- max(1, ceiling((ext[2] - ext[1]) / terra::xres(r)))
  h <- max(1, ceiling((ext[4] - ext[3]) / terra::yres(r)))
  c(w = w, h = h, n = as.numeric(w) * as.numeric(h))
}

#' Split a (multi)polygon geometry into its individual POLYGON parts as a list
#' of single-feature \code{sf} objects (keeping the CRS). Returned as \code{sf}
#' - not bare \code{sfc} - so \code{terra::vect()} accepts them across versions.
#' @keywords internal
#' @noRd
.geom_parts <- function(geom) {
  g <- sf::st_geometry(geom)
  g <- tryCatch(suppressWarnings(sf::st_cast(g, "POLYGON")),
                error = function(e) g)
  lapply(seq_along(g), function(i) sf::st_sf(geometry = g[i]))
}

#' Per-class MapBiomas areas over an AOI, bounded by a pixel budget (local
#' backend). Mirrors \code{mb_class_areas_raster(mb_raster_local(...))} for
#' ranges that fit the budget, and stays bounded for ones that do not.
#' @keywords internal
#' @noRd
.mb_class_areas_local <- function(aoi, year, collection, initiative, src = NULL,
                                  max_pixels = .area_max_pixels(),
                                  cache = TRUE, cache_dir = NULL) {
  ini <- .mb_resolve_initiative(initiative)
  if (identical(ini$provider, "esri"))
    return(.s2_class_areas(aoi, year = year, src = src, max_pixels = max_pixels,
                           cache = cache, cache_dir = cache_dir))
  if (is.null(collection)) collection <- ini$collection
  aoi <- sf::st_geometry(aoi)
  if (is.na(sf::st_crs(aoi))) sf::st_crs(aoi) <- 4326

  old <- .mb_set_gdal()
  on.exit(for (k in names(old)) terra::setGDALconfig(k, old[[k]]), add = TRUE)

  if (is.null(src)) src <- mb_source_url(year, collection, ini$key)
  path <- if (grepl("^https?://", src)) paste0("/vsicurl/", src) else src
  r0 <- tryCatch(
    suppressWarnings(terra::rast(path)),
    error = function(e)
      stop("Could not open MapBiomas raster at:\n  ", src,
           "\nOriginal error: ", conditionMessage(e), call. = FALSE))

  aoi_r <- sf::st_sf(geometry = sf::st_transform(aoi, terra::crs(r0)))
  full_px <- .bbox_native_px(r0, as.vector(terra::ext(terra::vect(aoi_r))))

  # Fast path: the whole window fits the budget -> native single read.
  if (full_px[["n"]] <= max_pixels)
    return(mb_class_areas_raster(mb_raster_local(
      aoi, year = year, collection = collection, initiative = initiative,
      src = src, mask = TRUE, cache = cache, cache_dir = cache_dir)))

  parts <- .geom_parts(aoi_r)
  part_px <- vapply(parts, function(p)
    .bbox_native_px(r0, as.vector(terra::ext(terra::vect(p))))[["n"]], numeric(1))

  if (length(parts) > 1L && sum(part_px) <= max_pixels) {
    # Sparse geometry (e.g. AOO cells): native read per part, exact.
    tabs <- Map(function(p, np) {
      v  <- terra::vect(p)
      rc <- terra::mask(terra::crop(r0, terra::ext(v), snap = "out"), v)
      mb_class_areas_raster(rc)
    }, parts, part_px)
  } else {
    # One or a few large contiguous parts (e.g. EOO hull): decimate each.
    tabs <- lapply(parts, function(p) .mb_area_read_part(path, p, r0, max_pixels))
  }
  .sum_class_areas(tabs)
}

#' Read one polygon part's class areas, decimating with majority resampling when
#' its native window exceeds the budget.
#' @keywords internal
#' @noRd
.mb_area_read_part <- function(path, part, r0, max_pixels) {
  ext <- as.vector(terra::ext(terra::vect(part)))   # xmin xmax ymin ymax
  px  <- .bbox_native_px(r0, ext)
  if (px[["n"]] <= max_pixels) {
    v  <- terra::vect(part)
    rc <- terra::mask(terra::crop(r0, terra::ext(v), snap = "out"), v)
    return(mb_class_areas_raster(rc))
  }
  fact  <- ceiling(sqrt(px[["n"]] / max_pixels))
  out_w <- max(1L, as.integer(floor(px[["w"]] / fact)))
  out_h <- max(1L, as.integer(floor(px[["h"]] / fact)))
  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)
  # Decimate on read with nearest-neighbour sampling: GDAL uses the GeoTIFF
  # overviews and reads ONLY the decimated pixels, so memory stays tiny even for
  # a continental window. (Majority/"mode" resampling would instead read the
  # whole native window to compute each cell's dominant class - hundreds of
  # millions of pixels - which OOM-kills the process on small machines.) For a
  # class-composition percentage over a large area, sampling millions of pixels
  # gives essentially the same proportions.
  ok <- tryCatch({
    sf::gdal_utils(
      "translate", source = path, destination = tmp,
      options = c("-projwin",
                  sprintf("%.10f", ext[1]), sprintf("%.10f", ext[4]),
                  sprintf("%.10f", ext[2]), sprintf("%.10f", ext[3]),
                  "-outsize", as.character(out_w), as.character(out_h),
                  "-r", "nearest"))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    # Fallback: native crop + modal aggregate (still bounded on output).
    v  <- terra::vect(part)
    rc <- terra::crop(r0, terra::ext(v), snap = "out")
    d  <- max(dim(rc)[1:2]); f <- max(1L, ceiling(d / sqrt(max_pixels)))
    if (f > 1L) rc <- terra::aggregate(rc, fact = f, fun = "modal", na.rm = TRUE)
    return(mb_class_areas_raster(terra::mask(rc, v)))
  }
  rc <- terra::rast(tmp)
  rc <- terra::mask(rc, terra::vect(part))
  mb_class_areas_raster(rc)
}

#' Sum a list of code/area_km2 tables into one, dropping empty pieces.
#' @keywords internal
#' @noRd
.sum_class_areas <- function(tabs) {
  tabs <- Filter(function(d) !is.null(d) && nrow(d) > 0, tabs)
  if (!length(tabs)) return(data.frame(code = integer(0), area_km2 = numeric(0)))
  all <- do.call(rbind, tabs)
  agg <- stats::aggregate(area_km2 ~ code, data = all, FUN = sum)
  agg$code <- as.integer(agg$code)
  agg[is.finite(agg$area_km2) & agg$area_km2 > 0, , drop = FALSE]
}