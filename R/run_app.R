#' Launch the mappingAS Shiny application
#'
#' Starts the interactive app for uploading occurrences, computing EOO/AOO and
#' MapBiomas habitat conversion, and exploring/exporting the results.
#'
#' @param launch.browser Logical; open in the default browser (default
#'   \code{TRUE}).
#' @param max_upload_mb Maximum size, in megabytes, of a file that can be
#'   uploaded in the app (default \code{500}). Shiny's own default is only 5 MB,
#'   which rejects large occurrence tables before they can be read; raise this
#'   for very large datasets. Set to \code{NULL} to leave any pre-existing
#'   \code{shiny.maxRequestSize} option untouched.
#' @param ... Passed to \code{shiny::runApp()}.
#' @return Invisibly \code{NULL}; called for its side effect.
#' @examples
#' if (interactive()) {
#'   run_app()
#'   # allow uploads up to 2 GB
#'   run_app(max_upload_mb = 2048)
#' }
#' @export
run_app <- function(launch.browser = TRUE, max_upload_mb = 500, ...) {
  for (p in c("shiny", "bslib", "leaflet", "DT", "htmltools", "htmlwidgets",
              "ggplot2", "plotly")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required to run the app.", call. = FALSE)
    }
  }
  if (!is.null(max_upload_mb)) {
    mb <- suppressWarnings(as.numeric(max_upload_mb)[1])
    if (!is.finite(mb) || mb <= 0)
      stop("`max_upload_mb` must be a single positive number.", call. = FALSE)
    options(shiny.maxRequestSize = mb * 1024^2)
  }
  app_dir <- system.file("shiny", package = "mappingAS")
  if (app_dir == "" || !file.exists(file.path(app_dir, "app.R"))) {
    stop("Could not locate the Shiny app. Reinstall mappingAS.", call. = FALSE)
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
