# CPTAC-BRCA protein detectability: how many of the 122 mass-spec tumors actually
# carried a quantified value for this gene's product.
#
# This is a statement about the ASSAY, not about the patient. Shotgun MS is
# semi-stochastic -- whether a protein is quantified in a given run is driven by its
# abundance and how well its peptides ionize -- so a low count means "hard to measure
# here", never "absent in these tumors". It is context for a gene, in the same family
# as R/tumor_normal.R, and it is deliberately NOT a survival result: CPTAC-BRCA has 2
# OS deaths over a 19.7-month window and is absent from config/cohorts.tsv, so no
# survival path can reach it.
#
# Read from the RAW quant file, not expr_cptac.h5: that matrix is z-scored per protein
# (mean 0 by construction) and zero-imputes missing values, which erases the exact
# coverage signal wanted here.
#
# ---- WHY THIS IS A FUNCTION IN R/ AND NOT THREE LINES IN app.R (2026-08-18) ----------
# The quant file is keyed by gene SYMBOL but its rows are Ensembl protein isoforms:
# the upstream umich matrix has a (Name, Database_ID) MultiIndex and scripts/export_cptac.py
# flattens it to level 0, so one symbol can own several ENSP rows. Measured on the file
# as staged: 12882 rows, 12020 symbols, 707 of them multi-row.
#
# app.R built the lookup as setNames(rowSums(!is.na(vals)), genes) and indexed it by name.
# R's name lookup silently returns the FIRST match, so for a multi-row symbol the badge
# reported one arbitrary isoform's coverage as the gene's. Measured: 281 symbols where the
# first row disagrees with the truth, worst case MKI67 -- ENSP00000357642 quantified in 7
# tumors, ENSP00000357643 in 122 -- so the app read "protein detected in 7/122 tumors (6%)"
# for a gene detected in all of them. Silent and plausible, the shape this project keeps
# meeting.
#
# THE RULE IS PER-SAMPLE OR, NOT max(). The question the badge asks is "was this gene's
# product quantified in this tumor", and any isoform answering yes settles it. Taking the
# max over per-row counts is a different question (the best single isoform) and is wrong
# whenever two isoforms are quantified in non-identical sample sets -- measured on this
# file: OR exceeds max for 14 symbols, by up to 24 tumors.
#
# The returned vector has UNIQUE names by construction, and that is asserted here rather
# than at the call site: the defect was a lookup silently picking one of several rows, so
# the fix belongs where the vector is built, not everywhere it is read.

CPTAC_RAW_PATH <- function() file.path(ROOT, "data", "raw", "cptac", "cptac_proteomics.tsv.gz")

# The invariant that makes a by-name lookup safe, as a function rather than an inline
# stopifnot -- so a test can hand it a duplicated vector and watch it refuse. A guard
# buried inside a function whose correct path never trips it is a guard nothing checks:
# mutate it to if (FALSE) and every test still passes. This one is reachable from outside.
.cptac_unique <- function(cov) {
  if (anyDuplicated(names(cov)))
    stop(paste("cptac_coverage(): collapse left duplicate symbols --",
               sprintf("%d of %d names repeat.", sum(duplicated(names(cov))), length(cov)),
               "A by-name lookup would resolve to whichever row came first, which is the",
               "defect this function exists to remove"))
  cov
}

# list(coverage = named integer, tumors x symbol; n = tumors; rows/symbols/multi = shape).
# An absent file is a legitimate state (the raw tree is not required to run the app), and
# returns n = 0 so the caller can drop the badge -- every OTHER failure is loud.
cptac_coverage <- function(path = CPTAC_RAW_PATH()) {
  if (!file.exists(path))
    return(list(coverage = setNames(integer(0), character(0)),
                n = 0L, rows = 0L, symbols = 0L, multi = 0L))
  dt <- data.table::fread(path)
  gene <- as.character(dt[[1]])
  vals <- as.matrix(dt[, -1, with = FALSE])
  if (!nrow(vals) || !ncol(vals))
    stop(sprintf("cptac_coverage(): %s has no data (%d rows x %d sample columns)",
                 basename(path), nrow(vals), ncol(vals)))
  blank <- is.na(gene) | !nzchar(trimws(gene))
  if (any(blank))
    stop(sprintf(paste("cptac_coverage(): %d row(s) of %s carry no gene symbol -- a blank",
                       "key would collapse them all into one bogus feature"),
                 sum(blank), basename(path)))
  ok <- !is.na(vals)
  # rowsum() groups in C and returns one row per unique symbol, so the OR is
  # (count of isoform hits in this tumor) > 0. Names come back unique and sorted.
  per <- rowsum(ok * 1L, group = gene, reorder = TRUE)
  cov <- rowSums(per > 0L)
  storage.mode(cov) <- "integer"
  list(coverage = .cptac_unique(cov), n = ncol(vals), rows = nrow(vals),
       symbols = length(cov), multi = sum(table(gene) > 1L))
}

# The one-or-two-sentence hover explanation for the badge, declared beside the number it
# explains rather than in app.R. What a reader has to know before using this line is that
# it is a property of the ASSAY: a low count means the protein is hard to measure by mass
# spec, not that the patients lack it. The length cap is enforced by info_note_html().
CPTAC_BADGE_INFO <- paste(
  "Mass spec does not detect every protein in every sample, so this counts how measurable",
  "the gene's protein is here. It is not a patient result and says nothing about outcome.")
