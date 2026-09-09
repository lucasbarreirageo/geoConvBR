#' Standalone HTML species factsheet with user metadata and photos
#'
#' Builds a single, self-contained HTML file that presents one species as a
#' printable/shareable factsheet - the kind hosted on a supplementary website.
#' It combines two sources of information:
#' \enumerate{
#'   \item \strong{Everything the assessment already computes} - range metrics
#'     (EOO/AOO and provisional Criterion B category), habitat conversion, fire
#'     and protected-area overlap - rendered as a key-metrics grid, the standard
#'     charts (composition, protection, and any land-cover / fire time series
#'     supplied), a "top anthropic activities" chart derived from the per-class
#'     land-cover breakdown, and the interpretive narrative of
#'     \code{\link{assessment_report}}.
#'   \item \strong{Information the package cannot know}, supplied by the user:
#'     the taxonomy (Family, Genus, Authority), the supporting-information block
#'     (Countries, System, Habitat, Biome, Vegetation), free-text land use and
#'     conservation units, a list of examined vouchers, a taxonomic reference
#'     (a Reflora / POWO link, or - for a newly described species - the article
#'     citation) and up to four photographs.
#' }
#' The whole page is self-contained (fonts fall back to the system stack, images
#' and charts are embedded as base64 data URIs), so the returned file can be
#' opened offline or published as-is (e.g. on GitHub Pages). Each photograph
#' carries a watermark in its lower-right corner with the owner name given in
#' \code{photo_credit}.
#'
#' @param assessment A \code{geoconv_assessment} from \code{\link{assess_species}}.
#' @param species Species name (default: the first assessed).
#' @param lang Factsheet language: \code{"en"} (default) or \code{"pt"}.
#' @param file Optional output path. When supplied the HTML is written there and
#'   the path is returned invisibly; otherwise the HTML is returned as a
#'   length-one character string.
#' @param family,genus,authority Taxonomy fields. \code{genus} defaults to the
#'   first word of \code{species}; the others are user-supplied (optional).
#' @param countries,system,habitat,biome,vegetation The "Supporting information"
#'   block, all optional free text.
#' @param land_use,conservation_units Optional overrides. When left \code{NULL}
#'   they are filled in automatically from the assessment: \code{land_use} from
#'   the anthropic land-cover classes within the EOO (the Conversion module) and
#'   \code{conservation_units} from the overlapping protected areas (the
#'   Protected-areas module). Pass a string to override the derived value; line
#'   breaks are preserved.
#' @param vouchers Optional examined material: a character vector (one voucher
#'   per element) or a single string with one voucher per line.
#' @param reference Optional taxonomic reference: a Reflora / POWO URL (rendered
#'   as a link) or, for a newly described taxon, the article citation.
#' @param photos Optional character vector of image file paths (PNG/JPEG/...);
#'   at most the first four are used.
#' @param photo_credit Watermark text (photo owner). A single value is applied
#'   to every photo; a vector is matched per photo.
#' @param cover_series,fire_series Optional land-cover / burned-area time series
#'   for this species (as for \code{\link{assessment_report}}); when supplied
#'   they add the temporal-trend narrative and their charts (a land-cover
#'   composition-over-time chart, and a burned-area-per-year chart).
#' @param map Logical; embed the distribution map (occurrence points, EOO, AOO
#'   and land cover). Default \code{TRUE}; set \code{FALSE} to skip it (e.g.
#'   offline or for speed).
#' @param map_interactive Logical; when \code{TRUE} (default) embed the
#'   \emph{interactive} Leaflet map (\code{\link{map_species}}, the same widget
#'   the app downloads via "Download map (HTML)") as a self-contained document
#'   inside an iframe, so the whole factsheet remains one portable file. When
#'   \code{FALSE}, or when the widget cannot be built (no \pkg{leaflet} /
#'   \pkg{htmlwidgets} / \pkg{pandoc}), the publication-ready static map
#'   (\code{\link{map_static}}) is embedded instead.
#' @param top_n_threats Number of anthropic land-cover classes to show in the
#'   "top anthropic activities" chart (default \code{5}).
#' @return The HTML string, or (when \code{file} is given) the file path,
#'   invisibly.
#' @seealso \code{\link{assessment_report}} for the narrative-only report.
#' @examples
#' \donttest{
#' occ <- read_occurrences(system.file("extdata", "example_occurrences.csv",
#'                                     package = "mappingAS"))
#' res <- assess_species(occ, year = 2024, verbose = FALSE)  # reads MapBiomas
#' out <- file.path(tempdir(), "factsheet.html")
#' factsheet_html(res, file = out, map = FALSE,
#'                family = "Gentianaceae", countries = "Brazil",
#'                system = "Terrestrial", biome = "Atlantic Forest")
#' }
#' @export
factsheet_html <- function(assessment, species = NULL,
                           lang = c("en", "pt"),
                           file = NULL,
                           family = NULL, genus = NULL, authority = NULL,
                           countries = NULL, system = NULL, habitat = NULL,
                           biome = NULL, vegetation = NULL,
                           land_use = NULL, conservation_units = NULL,
                           vouchers = NULL, reference = NULL,
                           photos = NULL, photo_credit = NULL,
                           cover_series = NULL, fire_series = NULL,
                           map = TRUE, map_interactive = TRUE,
                           top_n_threats = 5L) {
  lang <- match.arg(lang)
  stopifnot(inherits(assessment, "geoconv_assessment"))
  if (is.null(assessment$summary) || nrow(assessment$summary) == 0)
    stop("The assessment has no results to report.", call. = FALSE)
  if (is.null(species)) species <- assessment$summary$species[1]

  html <- .factsheet_build(
    assessment, species, lang,
    family = family, genus = genus, authority = authority,
    countries = countries, system = system, habitat = habitat,
    biome = biome, vegetation = vegetation,
    land_use = land_use, conservation_units = conservation_units,
    vouchers = vouchers, reference = reference,
    photos = photos, photo_credit = photo_credit,
    cover_series = cover_series, fire_series = fire_series,
    map = map, map_interactive = map_interactive,
    top_n_threats = top_n_threats)

  if (!is.null(file) && nzchar(file)) {
    writeLines(html, file, useBytes = TRUE)
    return(invisible(file))
  }
  html
}

# ---------------------------------------------------------------------------
# Small self-contained helpers (no extra package dependency)
# ---------------------------------------------------------------------------

# Vectorised, base-R base64 encoder for a raw vector. Used to inline images and
# rendered charts as data URIs so the factsheet is a single portable file.
.base64_encode <- function(raw) {
  if (!length(raw)) return("")
  alpha <- c(LETTERS, letters, as.character(0:9), "+", "/")
  n <- length(raw)
  pad <- (3L - (n %% 3L)) %% 3L
  x <- c(as.integer(raw), integer(pad))            # zero-pad to a multiple of 3
  m <- matrix(x, nrow = 3L)
  v <- m[1L, ] * 65536L + m[2L, ] * 256L + m[3L, ] # 24-bit groups (< 2^31)
  i1 <- v %/% 262144L
  i2 <- (v %/% 4096L) %% 64L
  i3 <- (v %/% 64L) %% 64L
  i4 <- v %% 64L
  chars <- character(4L * ncol(m))
  chars[c(TRUE,  FALSE, FALSE, FALSE)] <- alpha[i1 + 1L]
  chars[c(FALSE, TRUE,  FALSE, FALSE)] <- alpha[i2 + 1L]
  chars[c(FALSE, FALSE, TRUE,  FALSE)] <- alpha[i3 + 1L]
  chars[c(FALSE, FALSE, FALSE, TRUE)]  <- alpha[i4 + 1L]
  out <- paste(chars, collapse = "")
  if (pad > 0L) out <- paste0(substr(out, 1L, nchar(out) - pad), strrep("=", pad))
  out
}

.guess_mime <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         png = "image/png",
         jpg = "image/jpeg", jpeg = "image/jpeg",
         gif = "image/gif", webp = "image/webp",
         svg = "image/svg+xml", bmp = "image/bmp",
         "application/octet-stream")
}

# Read a file and return a `data:<mime>;base64,...` URI, or NULL on failure.
.file_data_uri <- function(path, mime = NULL) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path))
    return(NULL)
  sz <- file.info(path)$size
  if (!is.finite(sz) || sz <= 0) return(NULL)
  raw <- tryCatch(readBin(path, "raw", n = sz), error = function(e) NULL)
  if (is.null(raw)) return(NULL)
  if (is.null(mime)) mime <- .guess_mime(path)
  sprintf("data:%s;base64,%s", mime, .base64_encode(raw))
}

# Render a ggplot to a PNG data URI (headless-safe; NULL if it cannot be drawn).
.gg_data_uri <- function(g, width = 6.6, height = 3.4, dpi = 150) {
  if (!inherits(g, "ggplot")) return(NULL)
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({
    grDevices::png(tf, width = width, height = height, units = "in",
                   res = dpi, bg = "white")
    print(g)
    grDevices::dev.off()
    TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
  if (!isTRUE(ok) || !file.exists(tf)) return(NULL)
  on.exit(unlink(tf), add = TRUE)
  .file_data_uri(tf, "image/png")
}

# Build the interactive Leaflet map (map_species()) as a self-contained HTML
# document and return it as a `data:text/html;base64,...` URI to embed in an
# iframe. Returns NULL when the widget cannot be built or saved self-contained
# (e.g. leaflet / htmlwidgets / pandoc unavailable, or the species has no
# geometry) - the caller then falls back to the static map.
.leaflet_map_uri <- function(assessment, species, lang, st) {
  if (!requireNamespace("leaflet", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) return(NULL)
  m <- tryCatch(
    map_species(assessment, species = species,
                lang = if (lang == "pt") "pt" else "en",
                protected = isTRUE(st$protected)),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  tf <- tempfile(fileext = ".html")
  ok <- tryCatch({
    htmlwidgets::saveWidget(m, tf, selfcontained = TRUE, title = species)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok) || !file.exists(tf)) return(NULL)
  on.exit(unlink(tf), add = TRUE)
  .file_data_uri(tf, "text/html")
}

# Minimal HTML escaping for user-supplied text.
.esc_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

# Escape, then turn newlines into <br> (for multi-line free-text fields).
.esc_multiline <- function(x) {
  if (is.null(x)) return("")
  x <- paste(as.character(x), collapse = "\n")
  gsub("\n", "<br>", .esc_html(x), fixed = TRUE)
}

.blank <- function(x) is.null(x) ||
  (length(x) == 1L && (is.na(x) || !nzchar(trimws(as.character(x))))) ||
  (length(x) == 0L)

# The anthropic land-cover classes for a range, from the per-class table.
# Returns a data frame with label, area_km2, pct (of terrestrial area) and hex,
# ordered by area (descending). `top_n = Inf` keeps them all.
.factsheet_threats <- function(det, lang = "en", range = "eoo", top_n = 5L) {
  conv <- if (identical(range, "aoo")) det$aoo_conversion else det$eoo_conversion
  if (is.null(conv) || is.null(conv$by_class) || !is.data.frame(conv$by_class))
    return(NULL)
  byc <- conv$by_class
  if (!all(c("group", "area_km2") %in% names(byc))) return(NULL)
  a <- byc[!is.na(byc$group) & byc$group == "anthropic" &
           is.finite(byc$area_km2) & byc$area_km2 > 0, , drop = FALSE]
  if (!nrow(a)) return(NULL)
  a <- a[order(-a$area_km2), , drop = FALSE]
  if (is.finite(top_n)) a <- utils::head(a, top_n)
  lab <- if (lang == "pt" && "class_pt" %in% names(a)) a$class_pt
         else if ("class_en" %in% names(a)) a$class_en
         else as.character(a$code)
  terr <- sum(byc$area_km2[byc$group %in% c("natural", "anthropic")],
              na.rm = TRUE)
  data.frame(label = lab, area_km2 = a$area_km2,
             pct = if (terr > 0) 100 * a$area_km2 / terr else NA_real_,
             hex = if ("hex" %in% names(a)) a$hex else "#b34540",
             stringsAsFactors = FALSE)
}

# Long table of the top anthropic classes (ranked by EOO area) with their share
# in BOTH the EOO and the AOO, for a grouped bar chart.
.factsheet_threats_both <- function(det, lang = "en", top_n = 5L) {
  e <- .factsheet_threats(det, lang, "eoo", top_n)
  if (is.null(e)) e <- .factsheet_threats(det, lang, "aoo", top_n)
  if (is.null(e)) return(NULL)
  labs <- e$label
  ea <- .factsheet_threats(det, lang, "eoo", Inf)
  aa <- .factsheet_threats(det, lang, "aoo", Inf)
  look <- function(tab, l) {
    if (is.null(tab)) return(0)
    v <- tab$pct[match(l, tab$label)]
    if (length(v) != 1L || is.na(v)) 0 else v
  }
  rows <- lapply(labs, function(l) data.frame(
    label = l, range = c("EOO", "AOO"),
    pct = c(look(ea, l), look(aa, l)), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

.factsheet_threats_plot <- function(df, lang = "en") {
  if (is.null(df) || !nrow(df) ||
      !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  df$label <- factor(df$label, levels = rev(unique(df$label)))
  df$range <- factor(df$range, levels = c("EOO", "AOO"))
  xlab <- if (lang == "pt") "% da area terrestre" else "% of terrestrial range"
  ttl  <- if (lang == "pt") "Principais atividades antropicas (EOO e AOO)"
          else "Top anthropic activities (EOO and AOO)"
  cols <- c(EOO = "#1f8d49", AOO = "#e08a1e")
  dodge <- ggplot2::position_dodge(width = 0.72)
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["pct"]], y = .data[["label"]],
                                   fill = .data[["range"]])) +
    ggplot2::geom_col(width = 0.66, position = dodge,
                      colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", .data[["pct"]])),
                       position = dodge, hjust = -0.12, size = 2.9,
                       colour = "#1b2b23") +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.22))) +
    ggplot2::labs(x = xlab, y = NULL, title = ttl) +
    .mas_theme()
}

# Clean habitat-composition chart (Converted vs Natural, terrestrial %) for the
# EOO and AOO, from the summary percentages - so both bars are always complete
# and the labels always sit on their own colour.
.factsheet_composition_plot <- function(r, lang = "en") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  v <- c(r$eoo_converted_pct, r$eoo_natural_pct,
         r$aoo_converted_pct, r$aoo_natural_pct)
  if (all(is.na(v))) return(NULL)
  conv <- if (lang == "pt") "Convertido" else "Converted"
  nat  <- "Natural"
  # Factor levels are reversed relative to the desired visual order because
  # position_stack draws the first level at the right end.
  df <- data.frame(
    range = factor(rep(c("EOO", "AOO"), each = 2), levels = c("AOO", "EOO")),
    group = factor(rep(c(conv, nat), 2), levels = c(nat, conv)),
    pct   = c(r$eoo_converted_pct, r$eoo_natural_pct,
              r$aoo_converted_pct, r$aoo_natural_pct),
    stringsAsFactors = FALSE)
  df$lab <- ifelse(is.finite(df$pct) & df$pct >= 4, sprintf("%.0f%%", df$pct), "")
  cols <- stats::setNames(c("#1f8d49", "#d62839"), c(nat, conv))
  ttl  <- if (lang == "pt") "Composicao do habitat (cobertura)"
          else "Habitat composition (land cover)"
  xlab <- if (lang == "pt") "% da area terrestre" else "% of terrestrial area"
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["pct"]], y = .data[["range"]],
                                   fill = .data[["group"]])) +
    ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = .data[["lab"]]),
                       position = ggplot2::position_stack(vjust = 0.5),
                       colour = "white", fontface = "bold", size = 3.6) +
    ggplot2::scale_fill_manual(values = cols, name = NULL,
                               breaks = c(nat, conv)) +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::labs(x = xlab, y = NULL, title = ttl) +
    .mas_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# Protection chart (Natural-in-UC / Altered-in-UC / Outside) for EOO and AOO,
# built from the stored protected-area percentages, with each label guaranteed
# to sit on its own segment/colour.
.factsheet_protection_plot <- function(pa, lang = "en") {
  if (is.null(pa) || !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  num <- function(x) if (is.null(x) || is.na(x)) 0 else x
  has_nat <- !is.null(pa$eoo_nat_uc_pct) && is.finite(pa$eoo_nat_uc_pct)
  if (has_nat) {
    natl <- "Natural"
    altl <- if (lang == "pt") "Alterado (em UC)" else "Altered (in UC)"
    outl <- if (lang == "pt") "Fora de UC" else "Outside UC"
    lv <- c(natl, altl, outl)
    seg <- function(n, a) c(num(n), num(a), max(0, 100 - num(n) - num(a)))
    eoo <- seg(pa$eoo_nat_uc_pct, pa$eoo_alt_uc_pct)
    aoo <- seg(pa$aoo_nat_uc_pct, pa$aoo_alt_uc_pct)
    cols <- stats::setNames(c("#1f8d49", "#f2c14e", "#d9d9d9"), lv)
    lcol <- stats::setNames(c("white", "black", "black"), lv)
  } else {
    inl <- if (lang == "pt") "Em UC" else "In UC"
    outl <- if (lang == "pt") "Fora de UC" else "Outside UC"
    lv <- c(inl, outl)
    eoo <- c(num(pa$eoo_pct), 100 - num(pa$eoo_pct))
    aoo <- c(num(pa$aoo_pct), 100 - num(pa$aoo_pct))
    cols <- stats::setNames(c("#1f8d49", "#d9d9d9"), lv)
    lcol <- stats::setNames(c("white", "black"), lv)
  }
  df <- data.frame(
    range = factor(rep(c("EOO", "AOO"), each = length(lv)),
                   levels = c("AOO", "EOO")),
    group = factor(rep(lv, 2), levels = rev(lv)),
    pct   = c(eoo, aoo), stringsAsFactors = FALSE)
  df$lab <- ifelse(is.finite(df$pct) & df$pct >= 5, sprintf("%.0f%%", df$pct), "")
  df$lcol <- lcol[as.character(df$group)]
  ttl  <- if (lang == "pt") "Protecao por areas protegidas"
          else "Protection by protected areas"
  xlab <- if (lang == "pt") "% da area terrestre" else "% of terrestrial range"
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["pct"]], y = .data[["range"]],
                                   fill = .data[["group"]])) +
    ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = .data[["lab"]],
                                    colour = .data[["lcol"]]),
                       position = ggplot2::position_stack(vjust = 0.5),
                       fontface = "bold", size = 3.4, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cols, name = NULL, breaks = lv) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::labs(x = xlab, y = NULL, title = ttl) +
    .mas_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# Free-text land use derived from the conversion module (the anthropic classes
# present in the EOO), and conservation units from the protected-area list.
.derive_land_use <- function(det, lang = "en") {
  t <- .factsheet_threats(det, lang, "eoo", Inf)
  if (is.null(t) || !nrow(t)) return(NULL)
  paste0(t$label, " (", sprintf("%.1f%%", t$pct), ")", collapse = "; ")
}
.derive_cons_units <- function(det) {
  pa <- det$pa
  lst <- if (!is.null(pa)) pa$list else NULL
  nm <- NULL
  if (!is.null(lst) && is.data.frame(lst) && "pa_name" %in% names(lst)) {
    if ("n_occ" %in% names(lst)) {
      occ <- lst[!is.na(lst$n_occ) & lst$n_occ > 0, , drop = FALSE]
      nm <- if (nrow(occ)) occ$pa_name else lst$pa_name
    } else nm <- lst$pa_name
  } else if (!is.null(pa$layer) && "pa_name" %in% names(pa$layer)) {
    nm <- pa$layer$pa_name
  }
  nm <- unique(nm[!is.na(nm) & nzchar(nm)])
  if (!length(nm)) return(NULL)
  paste(nm, collapse = "; ")
}

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

.factsheet_build <- function(assessment, species, lang,
                             family, genus, authority,
                             countries, system, habitat, biome, vegetation,
                             land_use, conservation_units,
                             vouchers, reference,
                             photos, photo_credit,
                             cover_series, fire_series, map, map_interactive,
                             top_n_threats) {
  L <- function(en, pt) if (lang == "en") en else pt
  s <- assessment$summary
  r <- s[s$species == species, , drop = FALSE]
  if (nrow(r) == 0)
    stop("Species not found in assessment: ", species, call. = FALSE)
  r <- r[1, , drop = FALSE]
  st  <- assessment$settings
  det <- assessment$detail[[species]]
  sp_parts <- .sp_parts(species)
  if (.blank(genus)) genus <- strsplit(trimws(species), "\\s+")[[1]][1]
  # Land use and conservation units are derived from the package (the Conversion
  # and Protected areas modules) unless the caller supplied an explicit override.
  if (.blank(land_use)) land_use <- .derive_land_use(det, lang)
  if (.blank(conservation_units))
    conservation_units <- .derive_cons_units(det)

  fmt_km  <- function(x) if (.blank(x) || is.na(x))
    "&mdash;" else paste0(formatC(x, format = "f", big.mark = ",", digits = 2),
                          " km<sup>2</sup>")
  fmt_pct <- function(x) if (.blank(x) || is.na(x)) "&mdash;" else
    sprintf("%.1f%%", x)
  fmt_int <- function(x) if (.blank(x) || is.na(x)) "&mdash;" else
    formatC(x, format = "d", big.mark = ",")

  badge <- .iucn_badge(r$provisional_cat)

  # --- Photos (up to four) with watermark ---------------------------------
  photo_html <- ""
  if (!is.null(photos)) {
    photos <- as.character(photos)
    photos <- photos[!is.na(photos) & nzchar(photos)]
    photos <- utils::head(photos, 4L)
    credits <- if (.blank(photo_credit)) rep("", length(photos)) else {
      pc <- as.character(photo_credit)
      if (length(pc) == 1L) rep(pc, length(photos)) else
        rep(pc, length.out = length(photos))
    }
    cards <- character(0)
    for (i in seq_along(photos)) {
      uri <- .file_data_uri(photos[i])
      if (is.null(uri)) next
      wm <- if (nzchar(credits[i]))
        sprintf("<span class='fs-wm'>&copy; %s</span>", .esc_html(credits[i]))
        else ""
      cards <- c(cards, sprintf(
        "<figure class='fs-photo'><img src='%s' alt='%s'>%s</figure>",
        uri, .esc_html(species), wm))
    }
    if (length(cards))
      photo_html <- sprintf("<div class='fs-gallery fs-n%d'>%s</div>",
                            length(cards), paste(cards, collapse = ""))
  }

  # --- Taxonomy + supporting-information rows -----------------------------
  row <- function(lbl, val, italic = FALSE) {
    if (.blank(val)) return("")
    v <- if (isTRUE(italic)) sprintf("<i>%s</i>", .esc_html(val))
         else .esc_multiline(val)
    sprintf("<div class='fs-row'><dt>%s</dt><dd>%s</dd></div>", lbl, v)
  }
  tax_html <- paste0(
    row(L("Family", "Familia"), family),
    row(L("Genus", "Genero"), genus, italic = TRUE),
    row(L("Authority", "Autoria"),
        if (.blank(authority)) sp_parts$rest else authority))

  supp_html <- paste0(
    row(L("Countries", "Paises"), countries),
    row(L("System", "Sistema"), system),
    row(L("Habitat", "Habitat"), habitat),
    row(L("Biome", "Bioma"), biome),
    row(L("Vegetation", "Vegetacao"), vegetation))

  extra_html <- paste0(
    row(L("Land use", "Uso do solo"), land_use),
    row(L("Conservation units", "Unidades de conservacao"), conservation_units))

  # --- Vouchers -----------------------------------------------------------
  vouch_html <- ""
  if (!.blank(vouchers)) {
    vv <- if (length(vouchers) == 1L)
      strsplit(as.character(vouchers), "\r?\n")[[1]] else as.character(vouchers)
    vv <- trimws(vv); vv <- vv[nzchar(vv)]
    if (length(vv))
      vouch_html <- sprintf("<ul class='fs-vouchers'>%s</ul>",
        paste(sprintf("<li>%s</li>", .esc_html(vv)), collapse = ""))
  }

  # --- Reference (Reflora/POWO link or article citation) ------------------
  ref_html <- ""
  if (!.blank(reference)) {
    reference <- trimws(as.character(reference)[1])
    ref_html <- if (grepl("^https?://", reference, ignore.case = TRUE))
      sprintf("<a href='%s' target='_blank' rel='noopener'>%s</a>",
              .esc_html(reference), .esc_html(reference))
      else .esc_multiline(reference)
  }

  # --- Key metrics (already computed by the package) ----------------------
  metric <- function(lbl, val)
    sprintf("<div class='fs-metric'><span class='k'>%s</span><span class='v'>%s</span></div>",
            lbl, val)
  metrics <- c(
    metric(L("Provisional category", "Categoria provisoria"),
           sprintf("<span class='fs-cat' style='background:%s;color:%s'>%s</span>",
                   badge$bg, badge$fg,
                   if (.blank(r$provisional_cat)) "&mdash;" else
                     .esc_html(r$provisional_cat))),
    metric("EOO", fmt_km(r$eoo_km2)),
    metric("AOO", sprintf("%s (%s %s)", fmt_km(r$aoo_km2),
                          fmt_int(r$aoo_cells), L("cells", "celulas"))),
    metric(L("Occurrence records", "Registros de ocorrencia"),
           sprintf("%s (%s %s)", fmt_int(r$n_records), fmt_int(r$n_unique),
                   L("unique", "unicos"))))
  if (!is.null(r$n_subpop) && !is.na(r$n_subpop))
    metrics <- c(metrics,
      metric(L("Subpopulations (est.)", "Subpopulacoes (est.)"),
             fmt_int(r$n_subpop)))
  if (!is.null(r$n_locations) && !is.na(r$n_locations))
    metrics <- c(metrics,
      metric(L("Locations (est.)", "Localidades (est.)"),
             fmt_int(r$n_locations)))
  if (isTRUE(st$mapbiomas))
    metrics <- c(metrics,
      metric(L("Natural habitat (EOO / AOO)", "Habitat natural (EOO / AOO)"),
             sprintf("%s / %s", fmt_pct(r$eoo_natural_pct),
                     fmt_pct(r$aoo_natural_pct))),
      metric(L("Converted (EOO / AOO)", "Convertido (EOO / AOO)"),
             sprintf("%s / %s", fmt_pct(r$eoo_converted_pct),
                     fmt_pct(r$aoo_converted_pct))))
  if (isTRUE(st$fire) && !is.null(r$eoo_burned_pct))
    metrics <- c(metrics,
      metric(L("Burned at least once (EOO / AOO)",
               "Queimado ao menos uma vez (EOO / AOO)"),
             sprintf("%s / %s", fmt_pct(r$eoo_burned_pct),
                     fmt_pct(r$aoo_burned_pct))))
  if (isTRUE(st$protected) && !is.null(r$occ_in_uc_pct))
    metrics <- c(metrics,
      metric(L("Protected-area overlap (EOO / AOO)",
               "Sobreposicao com AP (EOO / AOO)"),
             sprintf("%s / %s", fmt_pct(r$eoo_uc_pct), fmt_pct(r$aoo_uc_pct))))
  metrics_html <- sprintf("<div class='fs-metrics'>%s</div>",
                          paste(metrics, collapse = ""))

  # --- Distribution map ---------------------------------------------------
  # Prefer the interactive Leaflet map the package already produces (same as
  # the app's "Download map (HTML)"), embedded as a self-contained document in
  # an iframe so the whole factsheet stays a single portable file. Fall back to
  # the static publication map if the widget cannot be built (no leaflet /
  # htmlwidgets / pandoc, or no geometry).
  map_cap <- L("Species range: occurrence points, EOO, AOO and land cover.",
               "Distribuicao da especie: pontos de ocorrencia, EOO, AOO e cobertura.")
  map_html <- ""
  if (isTRUE(map)) {
    map_uri <- if (isTRUE(map_interactive))
      .leaflet_map_uri(assessment, species, lang, st) else NULL
    if (!is.null(map_uri)) {
      map_html <- sprintf(
        "<iframe class='fs-mapframe' src='%s' loading='lazy' title='%s'></iframe><p class='fs-mapcap'>%s</p>",
        map_uri, .esc_html(species), map_cap)
    } else {
      png_uri <- .gg_data_uri(tryCatch(
        map_static(assessment, species = species,
                   lang = if (lang == "pt") "pt" else "en",
                   protected = isTRUE(st$protected)),
        error = function(e) NULL), width = 7.4, height = 6.2, dpi = 150)
      if (!is.null(png_uri))
        map_html <- sprintf(
          "<figure class='fs-fig fs-map'><img src='%s' alt='%s'><figcaption>%s</figcaption></figure>",
          png_uri, .esc_html(species), map_cap)
    }
  }

  # --- Charts (complementary information) ---------------------------------
  cover_list <- .as_series_list(cover_series)
  fire_list  <- .as_series_list(fire_series)
  figs <- list()
  add_fig <- function(cap, uri, w = "") {
    if (is.null(uri)) return(invisible())
    figs[[length(figs) + 1L]] <<- list(cap = cap, uri = uri)
  }
  # Top anthropic activities within BOTH the EOO and the AOO.
  thr_both <- .factsheet_threats_both(det, lang = lang, top_n = top_n_threats)
  add_fig(L("Top anthropic activities within the EOO and AOO.",
            "Principais atividades antropicas na EOO e na AOO."),
          .gg_data_uri(.factsheet_threats_plot(thr_both, lang),
                       width = 6.8, height = 3.4))
  if (isTRUE(st$mapbiomas))
    add_fig(L("Habitat composition of the EOO and AOO (land cover).",
              "Composicao do habitat na EOO e na AOO (cobertura)."),
            .gg_data_uri(.factsheet_composition_plot(r, lang),
                         width = 6.8, height = 2.7))
  if (isTRUE(st$protected))
    add_fig(L("Range protection by protected areas (EOO and AOO).",
              "Protecao da distribuicao por areas protegidas (EOO e AOO)."),
            .gg_data_uri(.factsheet_protection_plot(if (!is.null(det)) det$pa else NULL, lang),
                         width = 6.8, height = 2.7))
  for (ts in cover_list)
    add_fig(sprintf(L("Land-cover composition over time - %s.",
                      "Composicao da cobertura ao longo do tempo - %s."),
                    toupper(attr(ts, "range") %||% "")),
            .gg_data_uri(tryCatch(plot_timeseries(ts, lang = lang),
                         error = function(e) NULL), width = 6.6, height = 3.4))
  for (ts in fire_list)
    add_fig(sprintf(L("Burned area per year - %s.",
                      "Area queimada por ano - %s."),
                    toupper(attr(ts, "range") %||% "")),
            .gg_data_uri(tryCatch(plot_fire_timeseries(ts, lang = lang),
                         error = function(e) NULL), width = 6.6, height = 3.1))
  figs_html <- if (length(figs))
    sprintf("<div class='fs-figs'>%s</div>", paste(vapply(figs, function(f)
      sprintf("<figure class='fs-fig'><img src='%s' alt=''><figcaption>%s</figcaption></figure>",
              f$uri, f$cap), character(1)), collapse = "")) else ""

  # --- Narrative (reuse the assessment report), rendered as a single block of
  # bullet points. The Habitat-conversion section additionally lists the main
  # anthropic activities that altered the EOO and the AOO. --------------------
  b <- tryCatch(.report_build(assessment, species, lang, cover_series,
                              fire_series, figures = FALSE),
                error = function(e) NULL)
  narrative_html <- ""
  refs_pkg <- character(0)
  if (!is.null(b)) {
    conv_heading <- L("Habitat conversion", "Conversao de habitat")
    fmt_top <- function(t) if (is.null(t) || !nrow(t)) NULL else
      paste0(t$label, " (", sprintf("%.1f%%", t$pct), ")", collapse = "; ")
    top_e <- fmt_top(.factsheet_threats(det, lang, "eoo", top_n_threats))
    top_a <- fmt_top(.factsheet_threats(det, lang, "aoo", top_n_threats))
    conv_bullets <- character(0)
    if (!is.null(top_e)) conv_bullets <- c(conv_bullets, sprintf(L(
      "Main anthropic activities altering the EOO: %s.",
      "Principais atividades antropicas que alteraram a EOO: %s."), top_e))
    if (!is.null(top_a)) conv_bullets <- c(conv_bullets, sprintf(L(
      "Main anthropic activities altering the AOO: %s.",
      "Principais atividades antropicas que alteraram a AOO: %s."), top_a))
    secs <- vapply(b$sections, function(sec) {
      paras <- sec$p
      if (identical(sec$h, conv_heading) && length(conv_bullets))
        paras <- c(paras, conv_bullets)
      sprintf("<h3>%s</h3><ul class='fs-notes'>%s</ul>", sec$h,
              paste(sprintf("<li>%s</li>", paras), collapse = ""))
    }, character(1))
    narrative_html <- paste(secs, collapse = "")
    refs_pkg <- b$references
  }

  # --- References ---------------------------------------------------------
  ref_items <- character(0)
  if (nzchar(ref_html))
    ref_items <- c(ref_items, ref_html)
  if (length(refs_pkg))
    ref_items <- c(ref_items, vapply(refs_pkg, .esc_html, character(1)))
  refs_block <- if (length(ref_items))
    sprintf("<ol class='fs-refs'>%s</ol>",
            paste(sprintf("<li>%s</li>", ref_items), collapse = "")) else ""

  # --- Section wrappers ---------------------------------------------------
  section <- function(title, body)
    if (nzchar(body)) sprintf(
      "<section class='fs-card'><h2>%s</h2>%s</section>", title, body) else ""

  title_html <- sprintf("<i>%s</i>%s", .esc_html(sp_parts$italic),
    if (!.blank(authority)) paste0(" ", .esc_html(authority)) else
      if (nzchar(sp_parts$rest)) paste0(" ", .esc_html(sp_parts$rest)) else "")

  subtitle <- sprintf(L(
    "Preliminary IUCN Red List Criterion B screening (provisional) - generated on %s with the mappingAS R package.",
    "Triagem preliminar pelo Criterio B da Lista Vermelha da IUCN (provisoria) - gerada em %s com o pacote R mappingAS."),
    format(Sys.Date()))

  body_html <- paste0(
    "<header class='fs-header'>",
    sprintf("<div class='fs-title'><h1>%s</h1><p class='fs-sub'>%s</p></div>",
            title_html, subtitle),
    sprintf("<span class='fs-cat fs-cat-lg' style='background:%s;color:%s'>%s</span>",
            badge$bg, badge$fg,
            if (.blank(r$provisional_cat)) "&mdash;" else .esc_html(badge$code)),
    "</header>",
    photo_html,
    "<div class='fs-grid'>",
    section(L("Taxonomy", "Taxonomia"),
            if (nzchar(tax_html)) sprintf("<dl class='fs-dl'>%s</dl>", tax_html) else ""),
    section(L("Supporting information", "Informacoes de apoio"),
            if (nzchar(supp_html)) sprintf("<dl class='fs-dl'>%s</dl>", supp_html) else ""),
    section(L("Land use and conservation units",
              "Uso do solo e unidades de conservacao"),
            if (nzchar(extra_html)) sprintf("<dl class='fs-dl'>%s</dl>", extra_html) else ""),
    section(L("Examined vouchers", "Vouchers examinados"), vouch_html),
    "</div>",
    section(L("Key metrics", "Metricas principais"), metrics_html),
    section(L("Distribution map", "Mapa de distribuicao"), map_html),
    section(L("Complementary information (charts)",
              "Informacoes complementares (graficos)"), figs_html),
    section(L("Assessment notes", "Notas da avaliacao"), narrative_html),
    section(b$refs_heading %||% L("References", "Referencias"), refs_block))

  paste0(
    "<!DOCTYPE html>\n<html lang='", lang, "'>\n<head>\n",
    "<meta charset='utf-8'>\n",
    "<meta name='viewport' content='width=device-width, initial-scale=1'>\n",
    sprintf("<title>%s - %s</title>\n", .esc_html(species),
            L("Factsheet", "Ficha")),
    "<style>\n", .factsheet_css(), "\n</style>\n</head>\n",
    "<body>\n<main class='fs-page'>\n", body_html,
    "\n<footer class='fs-footer'>", .esc_html(sprintf(L(
      "Screening aid only - not a formal IUCN Red List assessment. mappingAS %s.",
      "Apenas auxilio de triagem - nao e uma avaliacao formal da Lista Vermelha da IUCN. mappingAS %s."),
      tryCatch(as.character(utils::packageVersion("mappingAS")),
               error = function(e) ""))),
    "</footer>\n</main>\n</body>\n</html>\n")
}

# Self-contained stylesheet for the factsheet (MapBiomas-inspired green).
.factsheet_css <- function() {
  "
  :root{--green:#1f8d49;--green-deep:#0e5c3a;--ink:#1b2b23;--muted:#6b776f;
        --line:#e3e8e3;--bg:#eef1ee;--surface:#ffffff;--radius:14px;}
  *{box-sizing:border-box;}
  body{margin:0;background:var(--bg);color:var(--ink);
       font-family:'Inter',system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;
       line-height:1.55;}
  .fs-page{max-width:960px;margin:0 auto;padding:28px 20px 48px;}
  .fs-header{display:flex;align-items:flex-start;justify-content:space-between;
             gap:16px;border-bottom:3px solid var(--green);padding-bottom:14px;
             margin-bottom:20px;}
  .fs-title{min-width:0;}
  .fs-title h1{margin:0 0 .15rem;font-size:1.9rem;font-weight:800;
               color:var(--green-deep);overflow-wrap:break-word;}
  .fs-title h1 i{font-style:italic;}
  .fs-sub{margin:0;color:var(--muted);font-size:.9rem;}
  .fs-cat{display:inline-block;padding:.15rem .5rem;border-radius:8px;
          font-weight:700;font-size:.85rem;}
  .fs-cat-lg{font-size:1.5rem;padding:.35rem .8rem;border-radius:12px;
             white-space:nowrap;box-shadow:0 2px 6px rgba(0,0,0,.12);}
  .fs-gallery{display:grid;gap:10px;margin:0 0 22px;}
  .fs-gallery.fs-n1{grid-template-columns:1fr;}
  .fs-gallery.fs-n2{grid-template-columns:repeat(2,1fr);}
  .fs-gallery.fs-n3,.fs-gallery.fs-n4{grid-template-columns:repeat(2,1fr);}
  .fs-photo{position:relative;margin:0;border-radius:var(--radius);
            overflow:hidden;background:#000;box-shadow:0 1px 3px rgba(0,0,0,.15);}
  .fs-photo img{display:block;width:100%;height:100%;object-fit:cover;
                aspect-ratio:4/3;}
  .fs-wm{position:absolute;right:8px;bottom:8px;padding:2px 8px;
         font-size:.72rem;font-weight:600;color:#fff;
         background:rgba(0,0,0,.5);border-radius:6px;
         text-shadow:0 1px 2px rgba(0,0,0,.6);backdrop-filter:blur(2px);}
  .fs-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:16px;
           margin-bottom:16px;}
  .fs-card{background:var(--surface);border:1px solid var(--line);
           border-radius:var(--radius);padding:16px 18px;
           min-width:0;overflow-wrap:break-word;}
  .fs-card h2{margin:0 0 .6rem;font-size:1.05rem;color:var(--green);
              border-bottom:1px solid var(--line);padding-bottom:.35rem;}
  .fs-card h3{margin:1rem 0 .3rem;font-size:.98rem;color:var(--green-deep);}
  .fs-card p{margin:.35rem 0;font-size:.92rem;overflow-wrap:break-word;}
  .fs-dl{margin:0;}
  .fs-row{display:grid;grid-template-columns:minmax(110px,36%) 1fr;gap:8px;
          padding:.3rem 0;border-bottom:1px dashed var(--line);}
  .fs-row:last-child{border-bottom:none;}
  .fs-row dt{margin:0;font-weight:600;color:var(--muted);font-size:.86rem;
             min-width:0;overflow-wrap:break-word;}
  .fs-row dd{margin:0;font-size:.92rem;min-width:0;overflow-wrap:break-word;}
  .fs-vouchers{margin:.2rem 0;padding-left:1.1rem;font-size:.9rem;}
  .fs-vouchers li{margin:.2rem 0;}
  .fs-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
              gap:10px;}
  .fs-metric{background:var(--bg);border:1px solid var(--line);
             border-radius:10px;padding:.5rem .7rem;min-width:0;}
  .fs-metric .k{display:block;font-size:.78rem;color:var(--muted);
                font-weight:600;overflow-wrap:break-word;}
  .fs-metric .v{display:block;font-size:1.05rem;font-weight:700;margin-top:2px;
                overflow-wrap:break-word;}
  .fs-figs{display:grid;grid-template-columns:1fr;gap:18px;}
  .fs-fig{margin:0;}
  .fs-fig img{display:block;width:100%;height:auto;border:1px solid var(--line);
              border-radius:10px;background:#fff;}
  .fs-fig figcaption{margin-top:.35rem;font-size:.82rem;color:var(--muted);
                     font-style:italic;}
  .fs-mapframe{width:100%;height:540px;border:1px solid var(--line);
               border-radius:10px;background:#fff;display:block;}
  .fs-mapcap{margin:.35rem 0 0;font-size:.82rem;color:var(--muted);
             font-style:italic;}
  .fs-notes{margin:.2rem 0 .6rem;padding-left:1.15rem;}
  .fs-notes li{margin:.3rem 0;font-size:.92rem;line-height:1.5;}
  .fs-refs{margin:.2rem 0;padding-left:1.2rem;font-size:.84rem;color:var(--ink);}
  .fs-refs li{margin:.3rem 0;line-height:1.4;word-break:break-word;}
  .fs-footer{margin-top:26px;padding-top:12px;border-top:1px solid var(--line);
             font-size:.8rem;color:var(--muted);text-align:center;}
  sup{font-size:.7em;}
  @media (max-width:640px){.fs-grid{grid-template-columns:1fr;}
    .fs-gallery.fs-n2,.fs-gallery.fs-n3,.fs-gallery.fs-n4{grid-template-columns:1fr;}
    .fs-header{flex-direction:column;}}
  @media print{body{background:#fff;}.fs-card{break-inside:avoid;}
    .fs-fig{break-inside:avoid;}}
  "
}
