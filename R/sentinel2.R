# ----------------------------------------------------------------------------
# Sentinel-2 / Esri global 10 m land cover (fallback outside MapBiomas)
# ----------------------------------------------------------------------------
# When a species' occurrences fall outside every MapBiomas product shipped with
# mappingAS (i.e. outside the supported South-American countries), the package
# can fall back to the Esri / Impact Observatory / Microsoft "10 m Annual Land
# Use Land Cover" (9-class) series derived from ESA Sentinel-2 - the global
# product behind the ArcGIS Living Atlas Land Cover Explorer
# (https://livingatlas.arcgis.com/landcoverexplorer/).
#
# It follows the same philosophy as the MapBiomas backend: public
# Cloud-Optimized GeoTIFFs streamed with GDAL's /vsicurl/ driver, no Google
# Earth Engine and no account. The global mosaic is tiled by MGRS grid-zone
# designator (e.g. "47P") on a public Azure bucket, so an area of interest is
# read from the one or few tiles it intersects; per-class areas are computed in
# each tile's native UTM projection and summed, which keeps the geodesic area
# statistics correct and avoids double counting across tiles.

#' @keywords internal
#' @noRd
.s2_years <- function() 2017:2023

#' @keywords internal
#' @noRd
.s2_default_base_url <- function() {
  # Public AWS Open Data bucket of the Esri / Impact Observatory 9-class
  # 10 m Annual LULC (COGs, no token). Objects are flat: "<GZD>_<YYYY>.tif"
  # (e.g. "15S_2023.tif"). Registry: https://registry.opendata.aws/io-lulc/
  "https://io-10m-annual-lulc.s3.us-west-2.amazonaws.com"
}

#' Internal descriptor for the Sentinel-2 / Esri initiative
#'
#' Mirrors the shape of a MapBiomas registry entry (\code{key}, \code{label},
#' \code{collection}, \code{years}) so it can flow through the same code paths,
#' but carries \code{provider = "esri"} to route reading and legend selection.
#' @keywords internal
#' @noRd
.s2_initiative <- function() {
  list(key = "sentinel2",
       label = "Sentinel-2 / Esri 10 m Land Cover (Impact Observatory)",
       collection = NA_integer_, years = .s2_years(),
       provider = "esri", build = NULL)
}

#' @keywords internal
#' @noRd
.s2_is_esri <- function(initiative) {
  prov <- tryCatch(.mb_resolve_initiative(initiative)$provider,
                   error = function(e) NA_character_)
  identical(prov, "esri")
}

#' Clamp a year to the Sentinel-2 / Esri available range
#' @keywords internal
#' @noRd
.s2_clamp_year <- function(year) {
  y <- suppressWarnings(as.integer(year))
  yy <- .s2_years()
  if (length(y) != 1L || is.na(y)) return(max(yy))
  min(max(y, min(yy)), max(yy))
}

#' Approximate national bounding boxes (lon/lat) of the MapBiomas products
#'
#' Used by \code{\link{.auto_initiative}} to pick the right MapBiomas country
#' product from occurrence coordinates. Each is \code{c(xmin, xmax, ymin, ymax)}.
#' @keywords internal
#' @noRd
.mb_country_bboxes <- function() list(
  brazil    = c(-74.0, -34.0, -34.0,   5.3),
  argentina = c(-74.0, -53.0, -55.5, -21.0),
  bolivia   = c(-70.0, -57.0, -23.0,  -9.5),
  chile     = c(-76.0, -66.0, -56.0, -17.0),
  colombia  = c(-79.5, -66.8,  -4.3,  13.5),
  ecuador   = c(-81.2, -75.0,  -5.1,   1.7),
  peru      = c(-81.4, -68.6, -18.4,   0.1),
  venezuela = c(-73.4, -59.8,   0.6,  12.3),
  paraguay  = c(-62.7, -54.2, -27.6, -19.3),
  uruguay   = c(-58.5, -53.0, -35.0, -30.0)
)

#' Auto-select the land-cover product from occurrence coordinates
#'
#' If the occurrences fall (mostly) within South America, returns the most
#' likely MapBiomas country product (so MapBiomas is prioritised there); each
#' point is assigned to the smallest-area national box that contains it, and the
#' modal country wins. Otherwise returns \code{"sentinel2"} so the global
#' Esri / Sentinel-2 layer is used directly. Points in South America but outside
#' every national box (e.g. the Guianas) default to \code{"brazil"}, which the
#' assessment's data-driven fallback still switches to Sentinel-2 if that
#' product has no data there.
#' @param pts An \code{sf}/\code{sfc} of occurrence points (any CRS).
#' @return A single initiative key (a MapBiomas country or \code{"sentinel2"}).
#' @keywords internal
#' @noRd
.auto_initiative <- function(pts) {
  if (!requireNamespace("sf", quietly = TRUE)) return("sentinel2")
  g <- sf::st_geometry(pts)
  if (is.na(sf::st_crs(g))) sf::st_crs(g) <- 4326
  g <- suppressWarnings(sf::st_transform(g, 4326))
  co <- sf::st_coordinates(g)
  if (!nrow(co)) return("sentinel2")
  lon <- co[, 1]; lat <- co[, 2]
  in_sa <- lon >= -82 & lon <= -34 & lat >= -56 & lat <= 13
  if (mean(in_sa, na.rm = TRUE) < 0.5) return("sentinel2")  # mostly outside SA
  boxes <- .mb_country_bboxes()
  areas <- vapply(boxes, function(b) (b[2] - b[1]) * (b[4] - b[3]), numeric(1))
  assign1 <- function(lo, la) {
    hit <- names(boxes)[vapply(boxes, function(b)
      lo >= b[1] && lo <= b[2] && la >= b[3] && la <= b[4], logical(1))]
    if (!length(hit)) return(NA_character_)
    hit[which.min(areas[hit])]
  }
  cc <- mapply(assign1, lon[in_sa], lat[in_sa])
  cc <- cc[!is.na(cc)]
  if (!length(cc)) return("brazil")
  names(sort(table(cc), decreasing = TRUE))[1]
}

#' Years available for the global Sentinel-2 / Esri land-cover series
#'
#' @return An integer vector of the years covered by the Esri / Impact
#'   Observatory 10 m Annual Land Use Land Cover product (currently
#'   2017-2023).
#' @examples
#' range(s2_years())
#' @seealso \code{\link{esri_legend}}, \code{\link{s2_source_url}},
#'   \code{\link{assess_species}}
#' @export
s2_years <- function() as.integer(.s2_years())

#' Sentinel-2 / Esri global 10 m land-cover legend
#'
#' Returns the 9-class legend of the Esri / Impact Observatory / Microsoft
#' "10 m Annual Land Use Land Cover" product (derived from ESA Sentinel-2, the
#' data behind the ArcGIS Living Atlas Land Cover Explorer), with each pixel
#' class mapped to the same conservation groups used by \pkg{mappingAS} for
#' MapBiomas (\code{"natural"}, \code{"anthropic"}, \code{"water"},
#' \code{"other"}, \code{"not_observed"}). This lets a range assessed with the
#' global fallback flow through the exact same conversion, plotting and report
#' pipeline as a MapBiomas one (see \code{\link{summarise_conversion}}).
#'
#' Group mapping: Trees, Rangeland and Flooded Vegetation are \code{natural};
#' Crops and Built Area are \code{anthropic} (converted habitat); Water is
#' \code{water}; Bare Ground and Snow/Ice are \code{other} (ambiguous, excluded
#' from the conversion denominator by default); Clouds are \code{not_observed}.
#'
#' @return A \code{data.frame} with columns \code{code}, \code{class_en},
#'   \code{class_pt}, \code{hex} (official Esri colour), \code{level1} and
#'   \code{group}.
#' @examples
#' esri_legend()
#' subset(esri_legend(), group == "anthropic")$class_en
#' @seealso \code{\link{mb_legend}}, \code{\link{summarise_conversion}}
#' @export
esri_legend <- function() {
  data.frame(
    code     = c(1L, 2L, 4L, 5L, 7L, 8L, 9L, 10L, 11L),
    class_en = c("Water", "Trees", "Flooded Vegetation", "Crops",
                 "Built Area", "Bare Ground", "Snow/Ice", "Clouds",
                 "Rangeland"),
    class_pt = c("Agua", "Arvores", "Vegetacao Alagavel", "Agricultura",
                 "Area Construida", "Solo Exposto", "Neve/Gelo", "Nuvens",
                 "Vegetacao Campestre"),
    hex      = c("#1A5BAB", "#358221", "#87D19E", "#FFDB5C",
                 "#ED022A", "#EDE9E4", "#F2FAFF", "#C8C8C8",
                 "#C6AD8D"),
    level1   = c("Water", "Forest", "Herbaceous/Shrubby", "Farming",
                 "NonVegetated", "NonVegetated", "Water", "NotObserved",
                 "Herbaceous/Shrubby"),
    group    = c("water", "natural", "natural", "anthropic",
                 "anthropic", "other", "other", "not_observed",
                 "natural"),
    stringsAsFactors = FALSE
  )
}

#' Build the Sentinel-2 / Esri land-cover GeoTIFF URL for a tile and year
#'
#' Returns the public Cloud-Optimized GeoTIFF URL of one MGRS grid-zone tile of
#' the Esri / Impact Observatory 10 m Annual Land Use Land Cover product. The
#' global mosaic is tiled by MGRS grid-zone designator (e.g. \code{"47P"}); an
#' area of interest usually intersects one or a few tiles, which the local
#' backend resolves automatically (see \code{\link{s2_raster_local}}).
#'
#' @param tile MGRS grid-zone designator (UTM zone number + latitude band
#'   letter), e.g. \code{"47P"}. The public bucket zero-pads the zone to two
#'   digits (e.g. \code{"01C"}).
#' @param year Integer year (2017-2023).
#' @param base_url Base URL of the public bucket holding the per-tile COGs.
#'   Defaults to the public AWS Open Data bucket \code{io-10m-annual-lulc};
#'   kept overridable so a mirror or an updated layout can be used without
#'   changing the package.
#' @return A length-1 character URL of the form \code{<base>/<tile>_<year>.tif}.
#' @examples
#' s2_source_url("47P", 2023)
#' @seealso \code{\link{s2_raster_local}}, \code{\link{esri_legend}}
#' @export
s2_source_url <- function(tile, year = 2023, base_url = NULL) {
  if (is.null(base_url)) base_url <- .s2_default_base_url()
  y <- as.integer(year)
  sprintf("%s/%s_%d.tif", sub("/+$", "", base_url), tile, y)
}

# ---- MGRS grid-zone helpers -------------------------------------------------

#' @keywords internal
#' @noRd
.s2_bands <- function() {
  c("C", "D", "E", "F", "G", "H", "J", "K", "L", "M",
    "N", "P", "Q", "R", "S", "T", "U", "V", "W", "X")
}

#' UTM zone number (1-60) for a longitude
#' @keywords internal
#' @noRd
.s2_utm_zone <- function(lon) {
  lon <- ((lon + 180) %% 360) - 180            # wrap to [-180, 180)
  z <- as.integer(floor((lon + 180) / 6) + 1L)
  ((z - 1L) %% 60L) + 1L
}

#' MGRS latitude band letter (C-X) for a latitude
#' @keywords internal
#' @noRd
.s2_lat_band <- function(lat) {
  bands <- .s2_bands()
  lat <- max(min(lat, 84 - 1e-9), -80)
  idx <- as.integer(floor((lat + 80) / 8) + 1L)
  idx <- max(1L, min(idx, length(bands)))
  bands[idx]
}

#' Geographic bounds (lon/lat) of one MGRS grid-zone cell
#' @keywords internal
#' @noRd
.s2_cell_bounds <- function(zone, band) {
  bands <- .s2_bands()
  i <- match(band, bands)
  south <- -80 + (i - 1L) * 8
  north <- if (identical(band, "X")) 84 else south + 8
  west <- -180 + (zone - 1L) * 6
  east <- west + 6
  c(xmin = west, ymin = south, xmax = east, ymax = north)
}

#' Enumerate the MGRS tiles intersecting an AOI (in lon/lat)
#' @keywords internal
#' @noRd
.s2_tiles_for_aoi <- function(aoi_ll) {
  bb <- as.numeric(sf::st_bbox(aoi_ll))          # xmin ymin xmax ymax
  step <- 3                                       # sample bbox to catch all zones
  xs <- unique(c(seq(bb[1], bb[3], by = step), bb[3]))
  ys <- unique(c(seq(bb[2], bb[4], by = step), bb[4]))
  grid <- expand.grid(lon = xs, lat = ys)
  keys <- unique(mapply(function(lo, la)
    paste(.s2_utm_zone(lo), .s2_lat_band(la), sep = "|"),
    grid$lon, grid$lat))
  lapply(strsplit(keys, "|", fixed = TRUE), function(zb) {
    zone <- as.integer(zb[1]); band <- zb[2]
    list(zone = zone, band = band, cell = .s2_cell_bounds(zone, band))
  })
}

#' Candidate tile URLs (padded and unpadded zone) or an explicit source
#'
#' \code{src} may be an explicit single GeoTIFF (ends in \code{.tif}), used
#' as-is, or a base-URL prefix that overrides the default bucket. The tile's
#' MGRS designator is tried both zero-padded (\code{"01C"}, Sentinel-2 style)
#' and unpadded (\code{"1C"}) so the reader is tolerant of either layout.
#' @keywords internal
#' @noRd
.s2_tile_urls <- function(tile, year, src, base_url) {
  if (!is.null(src) && grepl("\\.tif+$", src, ignore.case = TRUE)) return(src)
  bu <- if (!is.null(src)) src
        else if (!is.null(base_url)) base_url
        else .s2_default_base_url()
  padded <- sprintf("%02d%s", tile$zone, tile$band)
  unpad  <- sprintf("%d%s", tile$zone, tile$band)
  urls <- s2_source_url(padded, year, bu)
  if (!identical(unpad, padded)) urls <- c(urls, s2_source_url(unpad, year, bu))
  unique(urls)
}

#' Clip an AOI to one MGRS cell so tiles never double-count overlaps
#' @keywords internal
#' @noRd
.s2_clip_to_cell <- function(aoi_ll, cell) {
  bb <- sf::st_bbox(c(xmin = unname(cell["xmin"]), ymin = unname(cell["ymin"]),
                      xmax = unname(cell["xmax"]), ymax = unname(cell["ymax"])),
                    crs = 4326)
  cellpoly <- sf::st_as_sfc(bb)
  inter <- tryCatch(
    suppressWarnings(sf::st_intersection(sf::st_geometry(aoi_ll), cellpoly)),
    error = function(e) NULL)
  if (is.null(inter) || length(inter) == 0 || all(sf::st_is_empty(inter)))
    return(NULL)
  inter
}

#' Read (and crop/mask) one Sentinel-2/Esri tile to the AOI
#' @keywords internal
#' @noRd
.s2_read_tile <- function(tile, aoi_ll, year, src, base_url,
                          cache = TRUE, cache_dir = NULL) {
  clip <- .s2_clip_to_cell(aoi_ll, tile$cell)
  if (is.null(clip)) return(NULL)

  urls <- .s2_tile_urls(tile, year, src, base_url)
  r <- NULL; used <- NA_character_
  for (u in urls) {
    path <- if (grepl("^https?://", u)) paste0("/vsicurl/", u) else u
    r <- tryCatch(terra::rast(path), error = function(e) NULL)
    if (!is.null(r)) { used <- u; break }
  }
  if (is.null(r)) return(NULL)

  clip_r <- sf::st_transform(clip, terra::crs(r))
  vect <- terra::vect(clip_r)
  rc <- .mb_window(r, terra::ext(vect), used, year, 0L, cache, cache_dir)
  rc <- terra::mask(rc, vect)
  names(rc) <- "esri_class"
  rc
}

#' Crop a Sentinel-2 / Esri LULC mosaic to an area of interest
#'
#' Global analogue of \code{\link{mb_raster_local}} for the Esri / Impact
#' Observatory 10 m Annual Land Use Land Cover product. Determines the MGRS
#' grid-zone tile(s) intersecting \code{aoi}, streams a windowed read of each
#' via GDAL's \code{/vsicurl/} driver (no Google Earth Engine, no account),
#' masks to the polygon and, when the AOI spans more than one tile, merges the
#' pieces (reprojected to geographic coordinates for display). For accurate area
#' statistics the assessment path tabulates each tile in its native projection
#' and sums, rather than reprojecting; use \code{\link{summarise_conversion}} on
#' the class areas for that.
#'
#' @param aoi An \code{sf}/\code{sfc} polygon (any CRS), e.g. an EOO hull or the
#'   union of AOO cells.
#' @param year Integer year (2017-2023, default \code{2023}).
#' @param src Optional path/URL to a single Esri LULC GeoTIFF (ends in
#'   \code{.tif}, used as-is), or a base-URL prefix that replaces the default
#'   bucket. If \code{NULL} (default) the public bucket is used through
#'   \code{/vsicurl/}.
#' @param base_url Base URL of the public bucket (see \code{\link{s2_source_url}}).
#' @param mask Logical; if \code{TRUE} (default) pixels outside the polygon are
#'   set to \code{NA}.
#' @param cache Logical; cache the windowed crops on disk (default \code{TRUE}).
#' @param cache_dir Directory for the crop cache (default under \code{tempdir()}).
#' @return A \pkg{terra} \code{SpatRaster} of Esri land-cover codes for the AOI.
#' @seealso \code{\link{esri_legend}}, \code{\link{s2_source_url}},
#'   \code{\link{mb_raster_local}}
#' @export
s2_raster_local <- function(aoi, year = 2023, src = NULL, base_url = NULL,
                            mask = TRUE, cache = TRUE, cache_dir = NULL) {
  if (!requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    stop("Packages 'terra' and 'sf' are required.", call. = FALSE)
  }
  aoi <- sf::st_geometry(aoi)
  if (is.na(sf::st_crs(aoi))) sf::st_crs(aoi) <- 4326
  aoi_ll <- sf::st_transform(aoi, 4326)

  old <- .mb_set_gdal()
  on.exit(for (k in names(old)) terra::setGDALconfig(k, old[[k]]), add = TRUE)

  tiles <- .s2_tiles_for_aoi(aoi_ll)
  parts <- lapply(tiles, .s2_read_tile, aoi_ll = aoi_ll, year = year,
                  src = src, base_url = base_url, cache = cache,
                  cache_dir = cache_dir)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    stop("Could not read any Sentinel-2/Esri land-cover tile for this area ",
         "(year ", year, ").\nIf you are offline or the host is blocked, ",
         "download the GeoTIFF and pass it via `src=`.", call. = FALSE)
  }
  parts <- lapply(parts, function(r)
    if (terra::same.crs(r, "EPSG:4326")) r
    else terra::project(r, "EPSG:4326", method = "near"))
  if (length(parts) == 1L) {
    r <- parts[[1]]
  } else {
    # Independently projected tiles have different origins; project each onto a
    # single common geographic grid so they align, then merge (tiles are
    # disjoint MGRS cells, so there is nothing to double-count).
    exs  <- lapply(parts, function(r) as.vector(terra::ext(r)))
    tmpl <- terra::rast(
      terra::ext(min(vapply(exs, `[`, 0, 1)), max(vapply(exs, `[`, 0, 2)),
                 min(vapply(exs, `[`, 0, 3)), max(vapply(exs, `[`, 0, 4))),
      resolution = terra::res(parts[[1]]), crs = "EPSG:4326")
    parts <- lapply(parts, function(r) terra::project(r, tmpl, method = "near"))
    r <- do.call(terra::merge, parts)
  }
  names(r) <- "esri_class"
  if (isTRUE(mask)) r <- terra::mask(r, terra::vect(aoi_ll))
  r
}

#' Per-class Esri land-cover areas over an AOI (native-projection sum)
#'
#' Reads each intersecting MGRS tile in its native UTM projection, tabulates the
#' geodesic area per pixel class and sums across tiles. Because each AOI pixel is
#' clipped to exactly one MGRS cell, tiles never double-count. Returned in the
#' \code{code}/\code{area_km2} shape expected by \code{\link{summarise_conversion}}.
#' @keywords internal
#' @noRd
.s2_class_areas <- function(aoi, year = 2023, src = NULL, base_url = NULL,
                            cache = TRUE, cache_dir = NULL,
                            max_pixels = .area_max_pixels()) {
  if (!requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    stop("Packages 'terra' and 'sf' are required.", call. = FALSE)
  }
  aoi <- sf::st_geometry(aoi)
  if (is.na(sf::st_crs(aoi))) sf::st_crs(aoi) <- 4326
  aoi_ll <- sf::st_transform(aoi, 4326)

  old <- .mb_set_gdal()
  on.exit(for (k in names(old)) terra::setGDALconfig(k, old[[k]]), add = TRUE)

  tiles <- .s2_tiles_for_aoi(aoi_ll)
  tabs <- lapply(tiles, .s2_tile_class_areas, aoi_ll = aoi_ll, year = year,
                 src = src, base_url = base_url, max_pixels = max_pixels)
  .sum_class_areas(tabs)
}

#' Per-class Esri areas for one MGRS tile, bounded by a pixel budget.
#'
#' Reads the AOI-clipped part of the tile in its native UTM projection. As for
#' MapBiomas, sparse geometry (e.g. AOO cells) is read one part at a time at
#' native 10 m resolution, while a single large part (e.g. an EOO hull spanning
#' the tile) is decimated with majority resampling so the read stays bounded.
#' @keywords internal
#' @noRd
.s2_tile_class_areas <- function(tile, aoi_ll, year, src, base_url, max_pixels) {
  clip <- .s2_clip_to_cell(aoi_ll, tile$cell)
  if (is.null(clip)) return(NULL)

  urls <- .s2_tile_urls(tile, year, src, base_url)
  r <- NULL; used <- NA_character_
  for (u in urls) {
    path <- if (grepl("^https?://", u)) paste0("/vsicurl/", u) else u
    r <- tryCatch(terra::rast(path), error = function(e) NULL)
    if (!is.null(r)) { used <- path; break }
  }
  if (is.null(r)) return(NULL)

  clip_r <- sf::st_transform(clip, terra::crs(r))
  parts <- .geom_parts(clip_r)
  part_px <- vapply(parts, function(p)
    .bbox_native_px(r, as.vector(terra::ext(terra::vect(p))))[["n"]], numeric(1))

  if (sum(part_px) <= max_pixels) {
    tabs <- lapply(parts, function(p) {
      v  <- terra::vect(p)
      rc <- terra::mask(terra::crop(r, terra::ext(v), snap = "out"), v)
      mb_class_areas_raster(rc)
    })
  } else {
    tabs <- lapply(parts, function(p) .mb_area_read_part(used, p, r, max_pixels))
  }
  .sum_class_areas(tabs)
}

#' Downsampled Sentinel-2/Esri raster for display (not area calc)
#' @keywords internal
#' @noRd
.s2_raster_display <- function(aoi, year, src, max_pixels = 600,
                               crs = NULL, cache = TRUE, base_url = NULL) {
  r <- s2_raster_local(aoi, year = year, src = src, base_url = base_url,
                       mask = TRUE, cache = cache)
  d <- max(dim(r)[1:2])
  if (d > max_pixels)
    r <- terra::aggregate(r, fact = ceiling(d / max_pixels),
                          fun = "modal", na.rm = TRUE)
  if (!is.null(crs)) r <- terra::project(r, crs, method = "near")
  names(r) <- "esri_class"
  r
}
