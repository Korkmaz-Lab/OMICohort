# scan_lookup.R -- place a single Query-tab result within the genome-wide scan that is
# already on disk, WITHOUT recomputing anything and WITHOUT ever showing a q for a query
# the scan does not describe.
#
# The gap this closes: the Query tab exists to shop genes one at a time, and it reported
# a naked pooled p. That p is correct for the one hypothesis, but in a scan of ~1500 TFs
# a nominal p<0.05 is completely ordinary (ovarian OS: 733 of 1446 TFs clear it), so a
# reader over-reads it. The honest multiplicity anchor is the rank the TF ALREADY holds
# in the batch scan.
#
# Why a lookup and not a computed q: BH-FDR on a single hypothesis returns the input p
# unchanged, so a "corrected" q for one gene is the p wearing a different name -- a number
# that looks like a correction while doing nothing. So we do not compute; we look up the
# q the scan already assigned this TF.
#
# Why the match is STRICT: that stored q is only meaningful for a query run the SAME way
# the scan was. The scan is over VIPER activity, pools ALL the tissue's cohorts, at the
# tissue's declared horizon, and comes in two stratification flavours (marginal, and for
# breast a PAM50-adjusted one). Any difference -- a different score kind, a cohort subset,
# a different horizon, the other stratification -- means the stored q no longer describes
# this query, and we return "not applicable" rather than paste a stale number. This is the
# same disclose-don't-substitute rule as the small-k tau^2 flag and the PH "!" marker.
#
# Because the match is exact, the query's live pooled p equals the scan's stored p (verified
# to ~1e-15 by reproducing rows through get_survival), so the stored rank/q describes THIS
# query and no recomputation is needed.
#
# Pure functions, no I/O and no registry calls: the app resolves `scans`,
# `expected_cohorts` and `expected_horizon` and passes them in, so this file is testable
# against synthetic frames (tests/test_scan_lookup.R). app.R has no test harness, so any
# logic worth testing lives here and app.R only wires it -- the plots.R rule.

# NULL horizon (full follow-up) and a numeric horizon are different estimands, so they
# only match like-for-like. all.equal tolerates the numeric == numeric case.
.horizon_match <- function(a, b)
  (is.null(a) && is.null(b)) ||
  (!is.null(a) && !is.null(b) && isTRUE(all.equal(a, b)))

.ordinal <- function(n) {
  suf <- if (n %% 100 %in% 11:13) "th"
         else c("th","st","nd","rd","th","th","th","th","th","th")[(n %% 10) + 1]
  paste0(n, suf)
}

# Locate `feature` in the tissue's genome-wide scan, but only when the query's recipe
# matches the scan's exactly. Returns a list; $applicable is TRUE only for a clean hit.
# $status is one of: ok | kind | no_scan | cohorts | horizon | not_found.
#
#   scans             named list of the tissue's scan data.frames, keyed by basename
#                     ("mr_scan_os", "mr_scan_os_pam50adj", ...), as app.R loads them.
#   expected_cohorts  every cohort the scan pooled (cohorts_for(cancer_type, endpoint)).
#   expected_horizon  the horizon the scan was cut at, PER ENDPOINT
#                     (horizon_for(cancer_type, endpoint); NULL=full follow-up). The bare
#                     horizon_for(cancer_type) is wrong here: it collapses a tissue's
#                     per-endpoint taus with min(), which is a different estimand for every
#                     endpoint but one. See tests/test_scan_rank_wiring.R.
scan_rank_lookup <- function(feature, kind, cohorts, max_followup, adjust_strata,
                             endpoint, cancer_type, scans, expected_cohorts,
                             expected_horizon, fdr = 0.05, strat_var = NULL) {
  base <- list(applicable = FALSE, feature = feature,
               cancer_type = cancer_type, endpoint = endpoint)

  # The scan ranks VIPER activity; there is no expression/immune/cna scan to rank against.
  if (!identical(kind, "viper")) return(c(base, list(status = "kind")))

  # A stratified query belongs against the ADJUSTED scan, an unstratified one against the
  # marginal scan -- match each query to the scan run the SAME way, never across the two.
  # Tissues that never stratify (stratify_by=none -> adjust_strata FALSE) always resolve to
  # the marginal file, which is the only one they have.
  #
  # The suffix is DERIVED from the registry, not hardcoded. It was the literal "_pam50adj"
  # until 2026-08-27, which was silently correct only while breast was the sole stratified
  # tissue: the moment lgg declared stratify_by = idh (step 88), a stratified lgg query went
  # looking for mr_scan_os_pam50adj.csv -- a file that does not and will never exist -- and
  # reported "no scan" for a scan that was sitting right there as mr_scan_os_idhadj.csv.
  # This mirrors the suffix batch_regulators.R derives when WRITING the file, so the reader
  # and the writer cannot drift apart.
  strat <- isTRUE(adjust_strata)
  # strat_var is PASSED IN, not looked up. This file is a pure lookup/format layer with no
  # registry dependency -- deriving the name here by calling stratifier_for() would have
  # made it need R/registry.R sourced, which tests/test_scan_lookup.R (which loads this file
  # alone) caught immediately. Fail loud rather than defaulting: a caller that stratified but
  # did not say by what cannot have its scan file named, and silently answering "no scan"
  # would hide a wiring bug behind a normal-looking UI message.
  if (strat && (is.null(strat_var) || !nzchar(strat_var)))
    stop("scan_rank_lookup(): adjust_strata = TRUE but strat_var is missing, so the ",
         "adjusted scan's filename cannot be derived. Pass the stratifier the query was ",
         "fitted on (app.R: strata_vars_of(cancer_type)).")
  sv <- if (is.null(strat_var) || !nzchar(strat_var)) NA_character_ else strat_var
  base$strat_var <- sv
  # A MIXED pool has no single adjusted scan to rank against, and NA says so honestly.
  if (strat && is.na(sv))
    return(c(base, list(status = "no_scan", stratified = strat)))
  key <- paste0("mr_scan_", tolower(endpoint), if (strat) paste0("_", sv, "adj") else "")
  scan_df <- scans[[key]]
  if (is.null(scan_df) || !nrow(scan_df))
    return(c(base, list(status = "no_scan", stratified = strat)))

  # The scan pooled every endpoint-valid cohort of the tissue (expected_cohorts =
  # cohorts_for(cancer_type, endpoint)). Require the query to select all of them -- but
  # only that. An extra ticked cohort that is inert for this endpoint does not change the
  # pooled p, so it must not block the rank: SCANB is OS-only, so it drops out of every
  # DFS query AND out of the DFS scan alike, yet the UI leaves it ticked by default. A
  # genuine SUBSET -- deselecting a cohort the scan actually used -- does change the pool
  # and is still rejected here.
  if (!all(expected_cohorts %in% cohorts)) return(c(base, list(status = "cohorts")))
  # A horizon is part of the estimand, not a precision knob.
  if (!.horizon_match(max_followup, expected_horizon))
    return(c(base, list(status = "horizon")))

  i <- match(feature, scan_df$tf)
  if (is.na(i)) return(c(base, list(status = "not_found", stratified = strat)))

  rank <- match(feature, scan_df$tf[order(scan_df$p)])   # rank by p, don't trust file order
  row  <- scan_df[i, ]
  list(applicable = TRUE, status = "ok", feature = feature,
       cancer_type = cancer_type, endpoint = endpoint, stratified = strat,
       strat_var = sv,
       rank = rank, n = nrow(scan_df), q = row$q, p = row$p, HR = row$HR, k = row$k,
       clears_fdr = isTRUE(row$q < fdr),
       # Every breast row is k<=3 (tau^2 unidentified), so a breast rank rests on the
       # pooling rule; ovarian is k up to 13 and mostly identified. Surface it per row.
       pool_ci_identified = isTRUE(as.logical(row$pool_ci_identified)))
}

# The hover explanation for a rank that exists. The two unexplained terms in that sentence
# are `q` and which scan it is against -- both of which decide whether the rank means
# anything -- and neither fits on the line without doubling its length.
#
# It is attached to the `ok` state ONLY. The other five states are already complete
# sentences that say why there is no rank ("...so no rank is shown"), and one explanation
# stretched over six different messages would be wrong in five of them.
SCAN_RANK_INFO <- paste(
  "q is the p-value adjusted for testing every TF in the scan, so it says how surprising",
  "this rank is genome-wide. It is against the tissue's whole scan, not against your query.")

# Render a lookup result as one short line for the Query tab. Returns
# list(text, tone, info) with tone in {hit, miss, na} and info NULL unless there is a rank,
# or NULL when there is nothing worth saying.
format_scan_rank <- function(res) {
  # "PAM50-adjusted" was hardcoded here too -- see the note in scan_rank_lookup(). The
  # adjective now names whatever the registry actually declares, and falls back to a bare
  # "adjusted" if the variable is unknown, which is vague but never WRONG.
  .adj <- function(res) if (is.null(res$strat_var) || is.na(res$strat_var)) "adjusted"
                        else paste0(toupper(res$strat_var), "-adjusted")
  scan_lbl <- sprintf("%s %s %s scan", res$cancer_type, toupper(res$endpoint),
                      if (isTRUE(res$stratified)) .adj(res) else "marginal")
  if (identical(res$status, "ok")) {
    txt <- sprintf("%s ranks %s of %d TFs by p in the %s (VIPER activity); BH q = %s — %s",
                   res$feature, .ordinal(res$rank), res$n, scan_lbl, signif(res$q, 2),
                   if (isTRUE(res$clears_fdr)) "clears FDR<0.05." else "does not clear FDR<0.05.")
    if (!isTRUE(res$pool_ci_identified))
      txt <- paste0(txt, sprintf(
        " (τ² not identified at k=%s, so the ranking rests on the pooling rule.)", res$k))
    return(list(text = txt, tone = if (isTRUE(res$clears_fdr)) "hit" else "miss",
                info = SCAN_RANK_INFO))
  }
  msg <- switch(res$status,
    kind    = "Genome-wide rank applies to VIPER-activity queries only (the scan ranks VIPER activity), so it is not shown here.",
    no_scan = sprintf("No%s genome-wide scan exists for %s %s, so no rank is shown.",
                      if (isTRUE(res$stratified)) paste0(" ", .adj(res)) else "",
                      res$cancer_type, toupper(res$endpoint)),
    cohorts = "Genome-wide rank applies only to the full-cohort tissue scan; your cohort selection differs, so no rank is shown.",
    horizon = "Genome-wide rank applies only at the scan's follow-up horizon; your horizon differs, so no rank is shown.",
    not_found = sprintf("%s was not retained in the %s (too few events, or below VIPER's minimum regulon size in every cohort), so it has no genome-wide rank.",
                        res$feature, scan_lbl),
    NULL)
  if (is.null(msg)) return(NULL)
  list(text = msg, tone = "na", info = NULL)
}
