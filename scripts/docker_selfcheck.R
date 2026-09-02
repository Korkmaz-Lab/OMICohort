# Does this tree carry everything the app opens, and will it run here?
# (2026-09-01, step 124.)
#
# WHY THIS EXISTS. The image is built from a WHITELIST (.dockerignore excludes `*`, then
# re-includes ~1.76 GB of a ~19 GB tree). Docker's re-include rules are subtle, and every
# way of getting them wrong fails the same silent way: the image builds, starts, serves a
# page, and the missing piece only surfaces when a visitor queries the cohort whose .h5 was
# dropped. So the build RUNs this, and a missing file fails the BUILD, by name.
#
# It is not a Docker script. It takes the tree it is pointed at and asks whether the app
# could run there, so it is equally the check to run before uploading a data bundle. The
# required set is DERIVED -- the R closure by parsing source() calls out of app.R, the
# matrices and regulons out of config/cohorts.tsv -- because a hand-written list of 53
# cohorts' files is exactly the artifact that goes stale when cohort 54 arrives.
#
# EXIT 0 = this tree can serve. Anything else names what is missing and why it matters.

ROOT <- getwd()
ok <- TRUE
note <- function(...) cat(sprintf(...), "\n", sep = "")
fail <- function(...) { ok <<- FALSE; cat("  MISSING: ", sprintf(...), "\n", sep = "") }

# --- 1. .Rprofile was read -----------------------------------------------------------------
# Not a formality. .Rprofile sets shiny.autoload.r = FALSE, and R only reads it from the
# WORKING DIRECTORY at startup -- so an image whose WORKDIR is wrong, or a CMD that passes
# --vanilla, silently re-enables Shiny's autoload. Shiny would then source all of R/ on
# startup, including the CLI drivers whose commandArgs() guard only holds under Rscript.
# The app would die at boot with an error naming a pipeline script, not the real cause.
note("[1] .Rprofile")
if (!identical(getOption("shiny.autoload.r"), FALSE))
  fail(".Rprofile was not read (shiny.autoload.r is %s, expected FALSE). R must be STARTED in %s.",
       paste(deparse(getOption("shiny.autoload.r")), collapse = ""), ROOT)

# --- 2. the R files the app actually sources -----------------------------------------------
# Transitive closure by parsing, not a list: source("R/x.R") and source(file.path(ROOT,
# "R", "x.R")) are both used in this project and both have to be followed.
note("[2] R source closure")
closure <- local({
  seen <- character(0); stack <- "app.R"
  while (length(stack)) {
    f <- stack[[1]]; stack <- stack[-1]
    if (f %in% seen) next
    seen <- c(seen, f)
    p <- file.path(ROOT, f)
    if (!file.exists(p)) next
    txt <- readLines(p, warn = FALSE)
    m <- regmatches(txt, gregexpr('(?:sys\\.)?source\\([^)]*?"([^"]+\\.R)"', txt, perl = TRUE))
    for (hit in unlist(m)) {
      g <- sub('.*"([^"]+\\.R)".*', "\\1", hit)
      cand <- if (file.exists(file.path(ROOT, g))) g else file.path("R", basename(g))
      stack <- c(stack, cand)
    }
  }
  sort(seen)
})
for (f in closure) if (!file.exists(file.path(ROOT, f))) fail("%s (sourced by the app)", f)
note("    %d files, all present", length(closure))

# --- 3. registry, and the matrices it implies ----------------------------------------------
# Sourcing data_access.R resolves ROOT and loads the registry; it also fails loudly on its
# own if config/cohorts.tsv is absent, which is the check for that file.
note("[3] registry and feature matrices")
suppressWarnings(suppressMessages(source(file.path(ROOT, "R", "data_access.R"))))
reg <- COHORTS
note("    %d cohorts, %d cancer types", nrow(reg), length(unique(reg$cancer_type)))
if (!file.exists(file.path(PROC, "clinical.db"))) fail("data/processed/clinical.db")

# Every cohort must have expr_; viper_ only where a network backs it. Rather than encode
# that rule twice, require expr_ and REPORT the viper_ count, then assert it against the
# registry's own network_dir column below.
n_expr <- 0L
for (co in reg$cohort) {
  f <- file.path(PROC, sprintf("expr_%s.h5", tolower(co)))
  if (file.exists(f)) n_expr <- n_expr + 1L else fail("data/processed/expr_%s.h5", tolower(co))
}
n_viper <- sum(file.exists(file.path(PROC, sprintf("viper_%s.h5", tolower(reg$cohort)))))
note("    expr_*.h5 %d/%d, viper_*.h5 %d", n_expr, nrow(reg), n_viper)

# --- 4. ARACNe regulons --------------------------------------------------------------------
# guide.R stops the app if a declared network is missing, so a short count here is a boot
# failure there. Same source of truth (network_dir), reached without calling guide.R.
note("[4] ARACNe regulons")
dirs <- unique(reg$network_dir[nzchar(reg$network_dir)])
for (d in dirs) {
  f <- file.path(ROOT, "networks", d, sprintf("regulon_%s.rds", d))
  if (!file.exists(f)) fail("networks/%s/regulon_%s.rds (the Guide counts these)", d, d)
}
note("    %d declared, %d on disk", length(dirs),
     sum(file.exists(file.path(ROOT, "networks", dirs, sprintf("regulon_%s.rds", dirs)))))

# --- 5. scan tables ------------------------------------------------------------------------
# app.R globs results/<cancer_type>/*.csv. An empty tissue is not an error there -- the tab
# just shows nothing -- which is why it is checked here instead.
note("[5] scan tables")
for (ct in sort(unique(reg$cancer_type))) {
  n <- length(Sys.glob(file.path(ROOT, "results", ct, "*.csv")))
  if (n == 0L) fail("results/%s/ has no scan CSV (the Multiple-query tab would be empty)", ct)
}
note("    %d CSVs across %d tissues",
     length(Sys.glob(file.path(ROOT, "results", "*", "*.csv"))),
     length(unique(reg$cancer_type)))

# --- 6. the two raw inputs that are runtime reads -------------------------------------------
# Both are annotation layers, not survival input; both are read straight off data/raw.
note("[6] reference inputs")
for (f in c("data/raw/depmap/CRISPRInferredCommonEssentials.csv",
            "data/raw/depmap/depmap_screened_genes.txt",
            "data/raw/depmap/SOURCE.txt",
            "data/raw/cptac/cptac_proteomics.tsv.gz"))
  if (!file.exists(file.path(ROOT, f))) fail("%s", f)

# --- 7. branding assets --------------------------------------------------------------------
# Read out of R/about.R and R/brand.R rather than listed, so a new logo is covered.
note("[7] www assets")
# EVALUATED, not scraped. A first version of this grepped R/about.R for *.png|*.jpg and
# demanded every hit -- and its first run failed on `korkmaz_lab.jpg`, a filename that
# appears ONLY inside the comment at about.R:157 warning that this exact misspelling is
# not the file on disk. Reading prose as a requirement turns a comment about a bug into
# a bug. ABOUT_LOGOS is the structure the page renders from, so ask it.
suppressWarnings(suppressMessages(source(file.path(ROOT, "R", "about.R"))))
assets <- c(vapply(ABOUT_LOGOS, `[[`, character(1), "file"), "favicon.png")
for (a in assets) if (!file.exists(file.path(ROOT, "www", a)))
  fail("www/%s (rendered by the About page)", a)
note("    %d referenced, all present", length(assets))

# --- 8. packages ---------------------------------------------------------------------------
# The nine the running app loads. `viper` is deliberately NOT here: VIPER scores are
# precomputed into viper_*.h5, so the image never needs that dependency chain.
note("[8] packages")
for (p in c("shiny", "DT", "shinythemes", "rhdf5", "RSQLite", "DBI",
            "data.table", "metafor", "survival")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) fail("R package %s is not installed", p) else note("    %-12s %s", p, v)
}

# --- 9. graphics devices -------------------------------------------------------------------
# A headless Linux container is exactly where this breaks. png() with no cairo and no X11
# is the failure that renders every KM curve as an error box, and it does not stop the app
# from starting. Probed by OPENING the device, for the reason R/plots.R gives: on some
# machines capabilities("cairo") reports TRUE while the device dies on a missing library.
note("[9] graphics devices")
probe <- function(open) {
  f <- tempfile(); before <- grDevices::dev.cur()
  got <- tryCatch({ suppressWarnings(open(f)); TRUE }, error = function(e) FALSE)
  if (got && !identical(grDevices::dev.cur(), before)) grDevices::dev.off()
  unlink(f); got
}
if (!probe(function(f) grDevices::png(f, width = 400, height = 300)))
  fail("png() will not open -- every KM curve would render as an error (need cairo or X11)")
if (!probe(function(f) grDevices::cairo_pdf(f, width = 6, height = 4)))
  note("    note: cairo_pdf unavailable; PDF downloads fall back to pdf() (fonts by name, not embedded)")
note("    cairo capability: %s", capabilities("cairo"))

# -------------------------------------------------------------------------------------------
if (!ok) {
  cat("\nFAILED. The tree at ", ROOT, " cannot serve the app.\n", sep = "")
  quit(status = 1L)
}
cat("\nOK. Everything the app opens is present at ", ROOT, ".\n", sep = "")
