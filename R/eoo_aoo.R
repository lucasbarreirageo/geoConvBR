#' Extent of Occurrence (EOO) via minimum convex polygon 
#'
#' Computes the EOO as the area of the convex hull (minimum convex polygon)
#' enclosing all occurrence points, following the IUCN Red List guidelines and
#' the approach used by the \pkg{ConR} package. The hull is built in
#' longitude/latitude, its edges are densified with \code{st_segmentize()} so
#' that they follow great circles, and the area is measured directly on the
#' WGS84 ellipsoid (a "spheroid" measurement).
#'
#' Degenerate cases are handled so the function never errors:
#' \itemize{
#'   \item With fewer than three unique points a polygon cannot be built and the
#'     EOO is \code{NA}.
#'   \item When the unique points are perfectly collinear (they would form a
#'     line of area zero) a small amount of noise is added to the coordinates
#'     (\code{jitter}), exactly as in ConR, so that a valid polygon can still be
#'     produced.
#'   \item If, for any other reason, a valid polygon still cannot be built, the
#'     EOO is returned as \code{NA} rather than throwing an error.
#' }
#'
#' @param points An \code{sf} of POINT geometries (one species). Use
#'   \code{\link{read_occurrences}} to produce it.
#' @param segmentize_km Numeric edge length, in kilometres, used to densify the
#'   hull edges along great circles (default \code{20}, matching ConR).
#' @return A list with \code{area_km2} (numeric), \code{n_records} (integer),
#'   \code{n_unique} (distinct coordinates), \code{hull} (the hull as an
#'   \code{sfc} in WGS84, or \code{NULL} when the EOO is undefined) and
#'   \code{crs_laea} (a data-centred equal-area proj string, kept for mapping).
#' @examples
#' f <- system.file("extdata", "example_occurrences.csv", package = "mappingAS")
#' occ <- read_occurrences(f)
#' sp1 <- occ[occ$species == occ$species[1], ]
#' calc_eoo(sp1)$area_km2
#' @export
calc_eoo <- function(points, segmentize_km = 20) {
  .assert_points(points, "points")
  pts <- sf::st_geometry(sf::st_transform(points, 4326))
  co <- unique(round(sf::st_coordinates(pts)[, c("X", "Y"), drop = FALSE], 8))
  n_rec <- length(pts)
  n_uni <- nrow(co)

  # A data-centred equal-area CRS is kept in the output for downstream mapping
  # (the EOO area itself is measured on the ellipsoid, ConR-style).
  crs_laea <- tryCatch(laea_crs(pts), error = function(e) NA_character_)

  out <- list(area_km2 = NA_real_, n_records = n_rec, n_unique = n_uni,
              hull = NULL, crs_laea = crs_laea)
  if (n_uni < 3) {
    message("EOO undefined for < 3 unique locations (returning NA).")
    return(out)
  }

  # Break a perfectly straight line by jittering, as ConR does, so that a
  # polygon (rather than a degenerate line of area 0) can be built.
  if (.is_collinear(co)) {
    message("EOO: collinear points; adding noise (jitter) before hull.")
    guard <- 0L
    repeat {
      co[, 1] <- jitter(co[, 1], factor = 0.1)
      co[, 2] <- jitter(co[, 2], factor = 0.1)
      guard <- guard + 1L
      if (!.is_collinear(co) || guard > 100L) break
    }
  }

  hull <- tryCatch(
    .convex_hull_spheroid(co, segmentize_km),
    error = function(e) NULL
  )

  # Final safety net: if no valid polygon could be built, return NA (no error).
  if (is.null(hull) ||
      !any(as.character(sf::st_geometry_type(hull)) %in%
           c("POLYGON", "MULTIPOLYGON"))) {
    message("EOO undefined: could not build a valid polygon; returning NA.")
    return(out)
  }

  out$area_km2 <- .spheroid_area_km2(hull)  # m^2 -> km^2, S2-robust
  out$hull <- hull
  out
}

#' @keywords internal
#' @noRd
.convex_hull_spheroid <- function(co, segmentize_km) {
  idx <- grDevices::chull(co[, 1], co[, 2])
  idx <- c(idx, idx[1])
  ring <- .close_ring(.dedupe_vertices(co[idx, , drop = FALSE]))

  poly <- sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326)
  # Densify edges along great circles, then measure on the ellipsoid (ConR).
  poly <- sf::st_segmentize(poly, units::set_units(segmentize_km, "km"))
  # Segmentizing very long geodesic edges can emit consecutive coincident
  # vertices, which the S2 engine later rejects ("Edge is degenerate (duplicate
  # vertex)"). Drop them before validating so a large, widely-spread range still
  # yields a usable hull.
  poly <- .dedupe_polygon_rings(poly)
  tryCatch(sf::st_make_valid(poly), error = function(e) poly)
}

#' Remove consecutive (near-)duplicate vertices from a coordinate matrix.
#' Tolerance is ~1e-7 degrees (~1 cm), matching the S2 snap scale that
#' otherwise flags a zero-length edge as a degenerate duplicate.
#' @keywords internal
#' @noRd
.dedupe_vertices <- function(m, tol = 1e-7) {
  if (is.null(m) || nrow(m) < 2L) return(m)
  d <- abs(diff(m[, 1:2, drop = FALSE]))
  keep <- c(TRUE, (d[, 1] > tol) | (d[, 2] > tol))
  m[keep, , drop = FALSE]
}

#' Ensure a coordinate ring is explicitly closed (first vertex repeated last).
#' @keywords internal
#' @noRd
.close_ring <- function(m) {
  if (!is.null(m) && nrow(m) >= 1L &&
      any(m[1, ] != m[nrow(m), ]))
    m <- rbind(m, m[1, , drop = FALSE])
  m
}

#' De-duplicate the vertices of every ring of a single-POLYGON sfc, dropping any
#' ring left with too few points to be valid.
#' @keywords internal
#' @noRd
.dedupe_polygon_rings <- function(poly) {
  crs <- sf::st_crs(poly)
  g <- tryCatch(sf::st_geometry(poly)[[1]], error = function(e) NULL)
  if (is.null(g) || !length(g)) return(poly)
  rings <- lapply(g, function(r) .close_ring(.dedupe_vertices(r)))
  rings <- Filter(function(r) !is.null(r) && nrow(r) >= 4L, rings)
  if (!length(rings)) return(poly)
  tryCatch(sf::st_sfc(sf::st_polygon(rings), crs = crs),
           error = function(e) poly)
}

#' Ellipsoidal area (km^2) of a geographic geometry, robust to the S2 engine
#' rejecting a degenerate edge on very large hulls. On such a failure it retries
#' with S2 disabled, where \code{st_area()} on a longlat geometry falls back to
#' \pkg{lwgeom}'s geodesic (ellipsoidal) area - the same ConR-style measurement.
#' @keywords internal
#' @noRd
.spheroid_area_km2 <- function(g) {
  a <- tryCatch(sum(as.numeric(sf::st_area(g))) / 1e6, error = function(e) NA_real_)
  if (length(a) == 1L && is.finite(a)) return(a)
  old <- suppressMessages(sf::sf_use_s2())
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  g2 <- tryCatch(sf::st_make_valid(g), error = function(e) g)
  tryCatch(sum(as.numeric(sf::st_area(g2))) / 1e6, error = function(e) NA_real_)
}

#' @keywords internal
#' @noRd
.is_collinear <- function(co) {
  if (nrow(co) < 3) return(TRUE)
  r <- suppressWarnings(stats::cor(co[, 1], co[, 2]))
  is.na(r) || round(abs(r), 6) == 1
}

#' Area of Occupancy (AOO) via occupied grid cells 
#'
#' Computes the AOO as the number of occupied grid cells multiplied by the cell
#' area, following the IUCN Red List guidelines (reference scale: 2 x 2 km cells
#' = 4 km^2 each) and the approach used by the \pkg{ConR} package. Points are
#' snapped to a regular grid defined on a data-centred equal-area projection.
#'
#' Because the occupied-cell count depends on where the grid is placed, this
#' function tests several randomly translated grids and returns the
#' \emph{smallest} count, as ConR does. The returned \code{cells} correspond to
#' the grid placement that produced the reported value.
#'
#' @param points An \code{sf} of POINT geometries (one species).
#' @param cell_km Grid cell side length in kilometres (default \code{2}).
#' @param n_rep Integer; number of randomly translated grids to test (default
#'   \code{30}, matching the typical ConR setting). The minimum occupied-cell
#'   count over the replicates is returned.
#' @return A list with \code{area_km2}, \code{n_cells}, \code{cell_km},
#'   \code{cells} (occupied cells as an \code{sfc} polygon set in WGS84) and
#'   \code{crs_laea}.
#' @examples
#' f <- system.file("extdata", "example_occurrences.csv", package = "mappingAS")
#' occ <- read_occurrences(f)
#' sp1 <- occ[occ$species == occ$species[1], ]
#' calc_aoo(sp1)$area_km2
#' @export
calc_aoo <- function(points, cell_km = 2, n_rep = 30) {
  .assert_points(points, "points")
  pts <- sf::st_geometry(sf::st_transform(points, 4326))
  crs_laea <- laea_crs(pts)
  pts_p <- sf::st_transform(pts, crs_laea)

  size <- cell_km * 1000
  co <- sf::st_coordinates(pts_p)[, c("X", "Y"), drop = FALSE]

  best <- .aoo_min_grid(co, size, max(1L, as.integer(n_rep)))

  cells <- unique(data.frame(ix = best$ix, iy = best$iy))
  n_cells <- nrow(cells)

  # build occupied-cell polygons (for mapping / habitat extraction)
  polys <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    x0 <- best$ox + cells$ix[i] * size
    y0 <- best$oy + cells$iy[i] * size
    polys[[i]] <- sf::st_polygon(list(rbind(
      c(x0, y0), c(x0 + size, y0), c(x0 + size, y0 + size),
      c(x0, y0 + size), c(x0, y0)
    )))
  }
  cells_p <- sf::st_sfc(polys, crs = crs_laea)

  list(
    area_km2 = n_cells * (cell_km^2),
    n_cells = n_cells,
    cell_km = cell_km,
    cells = sf::st_transform(cells_p, 4326),
    crs_laea = crs_laea
  )
}

#' @keywords internal
#' @noRd
.aoo_min_grid <- function(co, size, n_rep) {
  best_n <- Inf
  best <- NULL
  for (h in seq_len(n_rep)) {
    # First replicate uses the origin grid; the rest are random translations.
    if (h == 1L) {
      ox <- 0
      oy <- 0
    } else {
      ox <- stats::runif(1) * size
      oy <- stats::runif(1) * size
    }
    ix <- floor((co[, "X"] - ox) / size)
    iy <- floor((co[, "Y"] - oy) / size)
    n  <- nrow(unique(cbind(ix, iy)))
    if (n < best_n) {
      best_n <- n
      best <- list(ix = ix, iy = iy, ox = ox, oy = oy)
      if (best_n == 1L) break  # cannot do better than a single cell
    }
  }
  best
}

