# Data-access layer for OMICohort.
# Pulls clinical/survival from clinical.db and feature vectors from expr_<cohort>.h5.
# All feature matrices are gene-wise z-scored WITHIN cohort, so a value of ~1 unit
# ~= 1 SD in that cohort -> Cox HRs are per-SD and comparable across cohorts.

suppressPackageStartupMessages({
  library(rhdf5); library(RSQLite); library(DBI)
})

# ROOT is the project directory, and it is RESOLVED rather than written down. From the
# beginning of the project until 2026-08-31 (step 119) it was a hard-coded ABSOLUTE PATH to
# one home directory, under the project's pre-rename name -- which meant this file, the one
# every test and the app load first, could only work on one machine, in one directory, under
# one name. Two consequences that were both invisible while the literal happened to be right:
#
#   1. A CLONE OF THIS REPOSITORY COULD NOT RUN. Nothing errored at source() time; the first
#      failure was a missing clinical.db, several frames deep, saying nothing about why.
#   2. FORTY-FOUR TESTS OPEN WITH `ROOT <- getwd()` AND WERE SILENTLY OVERRULED. Sourcing this
#      file replaced their derived root with the literal. They passed because the two agreed
#      on this machine; the derivation they appear to perform was doing nothing.
#
# Note that pkg/OMICohort had already solved this for the packaged fork -- omicohort_root()
# resolves option, then env var, then a fallback -- and the source tree it was extracted from
# never picked the fix up. The order below is the same one, so the two agree.
#
# getwd() is the honest default because the project already REQUIRES R to be started in the
# project directory (.Rprofile sets options(shiny.autoload.r = FALSE) and is only read at
# startup, in the working directory). The check that follows is the point: a wrong root must
# fail HERE, naming what it looked for, rather than four calls later as a missing file.
ROOT <- local({
  r <- getOption("omicohort.root",
                 default = Sys.getenv("OMICOHORT_ROOT", unset = getwd()))
  r <- tryCatch(normalizePath(r, winslash = "/", mustWork = TRUE), error = function(e) r)
  if (!file.exists(file.path(r, "config", "cohorts.tsv")))
    stop(sprintf(paste0(
      "OMICohort: cannot locate the project root.\n",
      "  Looked in: %s\n",
      "  Expected:  %s\n",
      "R must be STARTED in the project directory (setwd() does not read .Rprofile), or the ",
      "root must be given explicitly:\n",
      "  options(omicohort.root = \"/path/to/omicohort\")   # per session\n",
      "  OMICOHORT_ROOT=/path/to/omicohort Rscript ...        # per process"),
      r, file.path(r, "config", "cohorts.tsv")))
  r
})
PROC <- file.path(ROOT, "data", "processed")

# The cohort registry (config/cohorts.tsv) is the single source of truth for which
# cohorts exist and their tissue/stratifier. Sourced here — AFTER ROOT is defined
# (registry.R needs ROOT) and since registry.R does NOT source this file, there is no
# cycle. Default breast cohorts resolve from the registry (CPTAC excluded: not in it).
#
# local = TRUE is REQUIRED, not stylistic: the CLI drivers (discovery_validation.R,
# step_b_scan.R) load this file with sys.source(envir = lib) to keep its CLI block
# dormant. A plain source() would evaluate registry.R in globalenv, where ROOT does
# not exist -> "object 'ROOT' not found". local = TRUE evaluates it in whichever env
# this file is being loaded into, so both the plain-source and sys.source paths work.
source(file.path(ROOT, "R", "registry.R"), local = TRUE)
SURV_COHORTS <- cohorts_for("breast")

.db_con <- function() dbConnect(RSQLite::SQLite(), file.path(PROC, "clinical.db"))

# TCGA patients can have 2-3 samples (multiple vials/portions). Counting a patient
# more than once is pseudo-replication in a Cox/KM model, so for survival we keep ONE
# sample per patient: primary tumor (type 01) preferred, then lowest vial letter.
# Barcode = TCGA-XX-XXXX(1:12)-<type 14:15><vial 16>. TCGA rows are identified by the
# raw barcode pattern (^TCGA-), independent of the dataset label, so this works for
# TCGA_BRCA, TCGA_OV, etc. Non-TCGA cohorts (METABRIC MB-####, SCAN-B F####) are one
# row per patient, so untouched.
.dedup_tcga <- function(d) {
  is_t <- grepl("^TCGA-", d$sample_id)
  if (!any(is_t)) return(d)
  t <- d[is_t, ]; rest <- d[!is_t, ]
  patient <- substr(t$sample_id, 1, 12)
  stype   <- substr(t$sample_id, 14, 15)
  vial    <- substr(t$sample_id, 16, 16)
  prio    <- ifelse(stype == "01", 0L, 1L)          # primary tumor first
  ord     <- order(patient, prio, vial)
  t <- t[ord, ]
  t <- t[!duplicated(patient[ord]), ]
  rbind(t, rest)
}

# Clinical + chosen endpoint as generic time/event columns.
# OS = overall survival; DFS = disease-free/progression (METABRIC RFS, TCGA PFI);
# DSS = disease-specific survival (TCGA CDR only).
get_clinical <- function(cohorts = SURV_COHORTS, endpoint = c("OS", "DFS", "DSS")) {
  endpoint <- match.arg(endpoint)
  cohorts <- validate_cohorts(cohorts)
  con <- .db_con(); on.exit(dbDisconnect(con))
  q <- sprintf(
    "SELECT s.sample_id, s.dataset_id, s.pam50, s.stage, s.idh,
            v.os_time, v.os_event, v.dfs_time, v.dfs_event, v.dss_time, v.dss_event
     FROM samples s JOIN survival v USING(sample_id)
     WHERE s.dataset_id IN (%s)",
    paste(sprintf("'%s'", cohorts), collapse = ","))
  d <- dbGetQuery(con, q)
  d$pam50[!is.na(d$pam50) & d$pam50 == ""] <- NA   # TCGA has "" placeholders
  d$stage[!is.na(d$stage) & d$stage == ""] <- NA   # empty stage -> NA (generic stratifier)
  d <- .dedup_tcga(d)                              # one sample per TCGA patient
  d$time  <- switch(endpoint, OS = d$os_time,  DFS = d$dfs_time,  DSS = d$dss_time)
  d$event <- switch(endpoint, OS = d$os_event, DFS = d$dfs_event, DSS = d$dss_event)
  d
}

.h5_path <- function(cohort) file.path(PROC, sprintf("expr_%s.h5", tolower(cohort)))

# One feature (gene/protein/TF) as a named numeric vector across the given cohorts.
# Names are sample_id; cohort prefixes are NOT added (sample_ids are globally unique).
# kind "rppa" reads rppa_<cohort>.h5, whose features are ANTIBODIES (Akt_pS473), not
# gene symbols -- one gene can carry several. R/rppa.R owns the gene->antibody map.
get_feature <- function(feature, cohorts = SURV_COHORTS,
                       kind = c("expr", "viper", "immune", "cna", "rppa")) {
  kind <- match.arg(kind)
  cohorts <- validate_cohorts(cohorts)
  out <- list()
  n_with_h5 <- 0L   # cohorts that actually carry this kind of data (had the h5)
  for (coh in cohorts) {
    f <- file.path(PROC, sprintf("%s_%s.h5", kind, tolower(coh)))
    if (!file.exists(f)) next
    n_with_h5 <- n_with_h5 + 1L
    feats <- as.character(h5read(f, "/features"))
    idx <- match(feature, feats)
    if (is.na(idx)) { warning(sprintf("%s not in %s", feature, basename(f))); next }
    samps <- as.character(h5read(f, "/samples"))
    # rhdf5 reads the matrix as samples x features (transposed vs h5py's
    # features x samples), so a feature is a COLUMN: index = list(all rows, idx).
    col <- as.numeric(h5read(f, "/matrix", index = list(NULL, idx)))
    names(col) <- samps
    out[[coh]] <- col
  }
  # Empty result: name the aggregate cause ONCE, distinguishing "no requested cohort even
  # carries this kind of data" (a wrong kind/cohort pairing -- previously FULLY silent, no
  # h5 means no per-cohort warning fired either) from "the feature is in none of the cohorts
  # that do carry it" (a typo/unknown symbol). Still returns numeric(0) rather than stopping:
  # get_feature is a low-level accessor and an empty return is LEGITIMATE (the vocabulary
  # union -- the ovarian panel has no cna h5 -- and batch drivers like make_plots() rely on
  # it to skip a missing feature). The hard fail-loud lives at get_survival(), the analytical
  # boundary, which rejects an empty score outright so it can never reach a silent result.
  if (!length(out)) {
    if (n_with_h5 == 0)
      warning(sprintf("none of the requested cohorts carry %s data (no %s_*.h5): %s",
                      kind, kind, paste(cohorts, collapse = ", ")))
    else
      warning(sprintf("feature '%s' not found in any of the %d requested cohort(s) with %s data",
                      feature, n_with_h5, kind))
    return(numeric(0))
  }
  do.call(c, unname(out))   # concatenate, preserving sample_id names (no prefix)
}

# Convenience alias for gene expression.
get_expression <- function(gene, cohorts = SURV_COHORTS)
  get_feature(gene, cohorts, kind = "expr")
