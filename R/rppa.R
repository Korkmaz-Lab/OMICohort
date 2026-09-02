# RPPA protein / phospho-protein survival -- a SEPARATE measurement layer.
#
# WHAT THIS IS NOT
# ----------------
# Not a validation cohort. RPPA exists for TCGA only, so a query here is always k = 1:
# no pooling, no I^2, no between-cohort heterogeneity. The meta-analysis machinery that
# the expression forest is built around contributes nothing, and the pooled diamond that
# is the forest's whole point must never be drawn here.
#
# The rows of an RPPA panel are ANTIBODIES measured on ONE set of patients, not cohorts
# measured on separate ones. Every row is the same people. That is the opposite of what
# a row means everywhere else in this tool, which is why the panel is drawn with a
# deliberately different grammar (see rppa_panel_plot() in R/plots.R) and why these
# results are never handed to forest_plot().
#
# WHY IT IS WORTH A LAYER
# -----------------------
# 57 of the 245 staged antibodies are phospho-specific (the union over the SIX staged
# tissues; BRCA alone carries 56), and 40 of the 205 mapped genes
# carry two or three antibodies -- Akt (Akt, Akt_pS473, Akt_pT308), Src (Src, Src_pY416,
# Src_pY527), EGFR, S6, GSK3, Rb, ERK, ER-alpha, HER2, mTOR. So this layer answers a
# question no other layer in the tool can: does TOTAL protein track survival, or only
# the activated form? Src is the clean case -- total, pY416 (activating) and pY527
# (inhibitory) are three different biological claims about one gene, on one cohort.
#
# WHY NOT IN config/cohorts.tsv
# -----------------------------
# cohorts_for(cancer_type, endpoint) is called in a dozen places with no role filter and
# every one of them means "the expression cohorts for this tissue". A registered RPPA
# cohort would appear as an extra forest row drawn from patients ALREADY on the plot --
# and for ovarian it would enter a k = 13 pooled diamond and double-count TCGA_OV. So
# RPPA carries its own clinical frame (data/processed/rppa_clinical.tsv, built by
# scripts/build_rppa.py from the curated Liu CDR) and hands it to get_survival() through
# its `clinical` argument. Nothing about the expression cohorts changes.
#
# The COHORT IDS are reused as-is (TCGA_BRCA, not TCGA_BRCA_RPPA) because the h5 naming
# convention is <kind>_<cohort>.h5 and get_feature(kind = "rppa") should keep working
# through validate_cohorts() unchanged. Reuse is safe precisely because RPPA never
# reaches cohorts_for().
#
# UNSTRATIFIED, DELIBERATELY
# --------------------------
# The expression engine stratifies TCGA_BRCA by PAM50. This layer does not: its clinical
# frame comes from the CDR and carries no PAM50, and 12 of its patients are not in
# clinical.db at all. get_survival() REFUSES adjust_strata = TRUE on an injected frame,
# so this is enforced rather than assumed. The consequence -- an RPPA HR is not on the
# same footing as an expression HR for the same gene -- is stated on the panel.
#
# AND THEY ARE NOT THE SAME PATIENTS EITHER (2026-08-26, step 82)
# ---------------------------------------------------------------
# The footing caveat above was the only one on the tab, and it is not the whole problem.
# The two layers are staged INDEPENDENTLY -- RPPA from the TCPA archives, expression from
# whichever array release each cohort was built on -- so neither contains the other. In
# four tissues that is a rounding difference and in two it is not. See
# rppa_expr_overlap(), which derives it from the files rather than stating it.

RPPA_MAP_PATH  <- function() file.path(PROC, "rppa_antibody_map.tsv")
RPPA_CLIN_PATH <- function() file.path(PROC, "rppa_clinical.tsv")
RPPA_H5 <- function(cohort) file.path(PROC, sprintf("rppa_%s.h5", tolower(cohort)))

RPPA_SOURCE <- "TCGA RPPA (MD Anderson), GDAC Firehose stddata__2016_01_28"

# A phospho-specific antibody: Akt_pS473, MAPK_pT202_Y204, Src_pY416.
RPPA_PHOSPHO_RE <- "_p[STY][0-9]"

rppa_is_phospho <- function(antibody) grepl(RPPA_PHOSPHO_RE, antibody)

# Which cohort carries RPPA for a tissue -- DERIVED from what is on disk, not listed.
# A hand-written tissue->cohort table would go stale the moment a tissue is staged or
# dropped, and would disagree with the h5 files silently.
rppa_cohort_for <- function(cancer_type, registry = COHORTS) {
  co <- registry$cohort[registry$cancer_type == cancer_type]
  co <- co[file.exists(vapply(co, RPPA_H5, character(1)))]
  if (length(co) > 1L)
    stop("rppa_cohort_for(): tissue '", cancer_type, "' has more than one RPPA cohort (",
         paste(co, collapse = ", "), ") -- a panel would silently show only one")
  if (!length(co)) NA_character_ else co
}

rppa_available <- function(cancer_type) !is.na(rppa_cohort_for(cancer_type))

# --- the antibody -> gene map -------------------------------------------------------

.RPPA_MAP_CACHE <- new.env(parent = emptyenv())

rppa_map <- function(path = RPPA_MAP_PATH()) {
  key <- path
  if (!is.null(.RPPA_MAP_CACHE[[key]])) return(.RPPA_MAP_CACHE[[key]])
  if (!file.exists(path)) return(NULL)
  m <- read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
  need <- c("antibody", "genes", "n_genes", "tissues")
  miss <- setdiff(need, names(m))
  if (length(miss))
    stop("rppa_map(): ", basename(path), " is missing column(s): ",
         paste(miss, collapse = ", "), " -- rebuild with scripts/build_rppa.py")
  if (!nrow(m)) stop("rppa_map(): ", basename(path), " has zero rows")
  if (anyDuplicated(m$antibody))
    stop("rppa_map(): duplicate antibody row(s): ",
         paste(unique(m$antibody[duplicated(m$antibody)]), collapse = ", "))
  m$gene_list   <- strsplit(m$genes, ",", fixed = TRUE)
  m$tissue_list <- strsplit(m$tissues, ",", fixed = TRUE)
  m$n_genes     <- as.integer(m$n_genes)
  # The count and the list are two encodings of one fact; if they disagree the file was
  # edited by hand rather than rebuilt, and everything downstream inherits the lie.
  if (!identical(m$n_genes, lengths(m$gene_list)))
    stop("rppa_map(): n_genes disagrees with the genes column -- rebuild, do not edit")
  .RPPA_MAP_CACHE[[key]] <- m
  m
}

# Every gene queryable in a cohort, sorted. This is the tab's vocabulary.
rppa_genes <- function(cohort, map = rppa_map()) {
  if (is.null(map) || is.na(cohort)) return(character(0))
  keep <- vapply(map$tissue_list, function(t) cohort %in% t, logical(1))
  sort(unique(unlist(map$gene_list[keep], use.names = FALSE)))
}

# The antibodies for one gene in one cohort, TOTAL PROTEIN FIRST then phospho forms.
# The ordering is the panel's argument: the unmodified protein is the reference the
# phospho rows are read against, so it belongs on top rather than wherever the
# alphabet puts it.
rppa_antibodies <- function(gene, cohort, map = rppa_map()) {
  if (is.null(map) || is.na(cohort) || !length(gene) || is.na(gene)) return(character(0))
  hit <- vapply(seq_len(nrow(map)), function(i)
    gene %in% map$gene_list[[i]] && cohort %in% map$tissue_list[[i]], logical(1))
  ab <- map$antibody[hit]
  ab[order(rppa_is_phospho(ab), ab)]
}

# Every gene an antibody measures, for the "does not resolve the isoform" caveat.
rppa_genes_of <- function(antibody, map = rppa_map()) {
  i <- match(antibody, map$antibody)
  if (is.na(i)) character(0) else map$gene_list[[i]]
}

# --- clinical -----------------------------------------------------------------------

.RPPA_CLIN_CACHE <- new.env(parent = emptyenv())

# Keyed by PATH, like rppa_map(). The first version cached under a single name and so
# returned the default frame for ANY path once it had been read -- which would have made
# every fail-loud guard below untestable, and silently ignored a caller pointing at a
# different file.
.rppa_clinical_raw <- function(path = RPPA_CLIN_PATH()) {
  if (!is.null(.RPPA_CLIN_CACHE[[path]])) return(.RPPA_CLIN_CACHE[[path]])
  if (!file.exists(path))
    stop("rppa clinical frame is missing: ", path,
         " -- run scripts/fetch_rppa.py then scripts/build_rppa.py")
  d <- read.delim(path, stringsAsFactors = FALSE, na.strings = "NA")
  need <- c("sample_id", "dataset_id", "os_time", "os_event",
            "dfs_time", "dfs_event", "dss_time", "dss_event")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop("rppa clinical frame is missing column(s): ", paste(miss, collapse = ", "))
  if (anyDuplicated(d$sample_id))
    stop("rppa clinical frame has duplicate sample_id(s) -- one row per patient is ",
         "assumed by every count the panel prints")
  .RPPA_CLIN_CACHE[[path]] <- d
  d
}

RPPA_ENDPOINTS <- c("OS", "DFS", "DSS")

# The frame get_survival(clinical = ) wants: sample_id, dataset_id, time, event.
rppa_clinical <- function(cohort, endpoint = "OS", path = RPPA_CLIN_PATH()) {
  if (!isTRUE(endpoint %in% RPPA_ENDPOINTS))
    stop("rppa_clinical(): unknown endpoint '", endpoint, "'. Known: ",
         paste(RPPA_ENDPOINTS, collapse = ", "))
  all <- .rppa_clinical_raw(path)
  d <- all[all$dataset_id == cohort, , drop = FALSE]
  if (!nrow(d))
    stop("rppa_clinical(): no rows for cohort '", cohort, "' -- staged cohorts are ",
         paste(sort(unique(all$dataset_id)), collapse = ", "))
  pre <- switch(endpoint, OS = "os", DFS = "dfs", DSS = "dss")
  d$time  <- d[[paste0(pre, "_time")]]
  d$event <- d[[paste0(pre, "_event")]]
  d[, c("sample_id", "dataset_id", "time", "event")]
}

# --- the panel ----------------------------------------------------------------------

RPPA_MIN_EVENTS <- 10

# One gene's antibodies in one tissue, each fitted independently on the SAME patients.
#
# Returns a list with $rows in display order. A row is either estimable (HR/ci/p) or
# skipped with a reason -- the same two shapes get_survival() itself produces, kept so
# the panel can show a named absence rather than dropping the antibody.
rppa_panel <- function(gene, cancer_type, endpoint = "OS") {
  cohort <- rppa_cohort_for(cancer_type)
  if (is.na(cohort))
    stop("rppa_panel(): no RPPA data staged for tissue '", cancer_type, "'")
  ab <- rppa_antibodies(gene, cohort)
  if (!length(ab))
    return(list(gene = gene, cohort = cohort, cancer_type = cancer_type,
                endpoint = endpoint, rows = list(), n = 0L, events = 0L,
                status = "absent",
                msg = sprintf("%s is not measured by any antibody on the RPPA panel in %s.",
                              gene, cohort)))
  clin <- rppa_clinical(cohort, endpoint)
  rows <- list()
  for (a in ab) {
    sc <- get_feature(a, cohort, kind = "rppa")
    r <- tryCatch(
      get_survival(sc, endpoint = endpoint, cohorts = cohort, clinical = clin,
                   adjust_strata = FALSE, min_events = RPPA_MIN_EVENTS),
      error = function(e) e)
    # The antibody has to ride on the result, not just in the row: km_plot() titles from
    # attr(res, "feature") and without it every RPPA curve is headed "score", which is the
    # one word that cannot distinguish Src_pY416 from Src_pY527.
    if (!inherits(r, "error")) r <- with_feature(r, a)
    if (inherits(r, "error")) {
      rows[[a]] <- list(antibody = a, genes = rppa_genes_of(a),
                        phospho = rppa_is_phospho(a),
                        skipped = TRUE, reason = conditionMessage(r), res = NULL)
      next
    }
    p1 <- r$per_cohort[[cohort]]
    if (is.null(p1) || isTRUE(p1$skipped)) {
      rows[[a]] <- list(antibody = a, genes = rppa_genes_of(a),
                        phospho = rppa_is_phospho(a), skipped = TRUE,
                        reason = if (is.null(p1)) "no estimable rows" else p1$reason,
                        n = if (is.null(p1)) 0L else p1$n,
                        events = if (is.null(p1)) 0L else p1$events, res = r)
      next
    }
    rows[[a]] <- list(
      antibody = a, genes = rppa_genes_of(a), phospho = rppa_is_phospho(a),
      skipped = FALSE, n = p1$n, events = p1$events,
      HR = p1$HR, lo = exp(p1$logHR - 1.96 * p1$se), hi = exp(p1$logHR + 1.96 * p1$se),
      p = p1$p, ph_violated = p1$ph_violated, res = r)
  }
  est <- Filter(function(x) isFALSE(x$skipped), rows)
  # Every row is the same patient set, so ONE n and ONE event count describe the whole
  # panel -- printing them per row would invite reading the rows as independent samples.
  # If they ever disagree the "same patients" claim is false and the panel must not
  # quietly average over it.
  #
  # It is not ALWAYS one set. 13 of the 208 ovarian antibodies carry missing values -- one
  # was not run on 415 of 425 patients -- so where two antibodies of a gene come from
  # different batches the usable n differs. Measured over ovarian's 37 multi-antibody
  # genes: 34 share n exactly, 3 differ by 10-20 patients. That is ordinary upstream
  # missingness, not a defect, so it is REPORTED (a range, and per-row n on the panel)
  # rather than warned about or averaged away.
  ns <- unique(vapply(est, `[[`, numeric(1), "n"))
  ev <- unique(vapply(est, `[[`, numeric(1), "events"))
  list(gene = gene, cohort = cohort, cancer_type = cancer_type, endpoint = endpoint,
       rows = unname(rows), n = if (length(ns) == 1L) as.integer(ns) else NA_integer_,
       events = if (length(ev) == 1L) as.integer(ev) else NA_integer_,
       n_range = if (length(ns)) as.integer(range(ns)) else integer(0),
       events_range = if (length(ev)) as.integer(range(ev)) else integer(0),
       n_varies = length(ns) > 1L, events_varies = length(ev) > 1L,
       n_estimable = length(est), status = if (length(est)) "ok" else "empty",
       msg = if (length(est)) NULL else
         sprintf("%s has %d antibody(ies) on the panel, none estimable at %s.",
                 gene, length(rows), endpoint))
}

# Phrased to stay true on the panels where n varies: "one cohort" holds even when an
# antibody was not run on part of it, whereas "the same patients" would not. The claim
# that has to survive intact is the second sentence.
# The third sentence is a KEY, not decoration. The marker fill has encoded phospho vs total
# protein since the panel was written, and until step 51 nothing on the page said so -- the
# first reader to look at a rendered panel read it as an inconsistency in the drawing, which
# is precisely what an unexplained encoding is. It is redundant with the "(P)" on the row
# label, deliberately: the fill is what makes the total/phospho contrast readable at a
# glance, and the label is what makes it readable at all.
#
# RPPA_CAPTION_LINES_MAX in R/plots.R budgets the height for this text. The two move
# together; tests/test_rppa_panel.R measures the wrap on a real device rather than trusting
# either.
RPPA_PANEL_CAPTION <- paste(
  "Rows are antibodies on one cohort, not independent cohorts.",
  "They are never combined, and there is no pooled estimate.",
  "Filled marker = phospho-specific antibody; open = total protein.")

# --- DepMap common-essentiality, on the panel's GENE (2026-08-20, step 51) -------------
# The flag transfers to this tab, and the chain is SHORTER here than on the expression tab,
# not longer. DepMap calls a gene common essential from a CRISPR knockout screen; a
# knockout removes the PROTEIN, and the protein is what this tab measures. The expression
# tab asks a reader to go mRNA -> assumed protein -> phenotype; here the measured molecule
# and the perturbed molecule are the same one.
#
# Two things it still is not:
#
#   * a patient result. Cell lines in a dish. Unchanged from the expression tab, and
#     DEPMAP_BADGE_INFO already says so.
#   * a statement about the PHOSPHO rows. "You need this protein to grow" is not "you need
#     it phosphorylated" -- a knockout is not a kinase-dead knock-in. This is why the badge
#     is drawn once for the panel, keyed on the gene, and never per row: a line sitting
#     beside mTOR_pS2448 would be making a claim DepMap did not test.
#
# Fire rate on this tab's own vocabulary, measured 2026-08-20 over all 205 mapped genes:
# 26 common essential (12.7%), 173 screened-not-essential, 6 never screened. Higher than
# the 3.0% seen across the TF vocabulary because an RPPA panel is deliberately built from
# core signalling proteins -- still rare enough that the flag carries information.
RPPA_DEPMAP_INFO <- paste(
  "Knockout removes the protein this tab measures, so the flag is about the gene, not one",
  "antibody. It says nothing about whether the phosphorylated form matters.")

# Every gene this tab can offer, across all five cohorts. The guard below needs a
# vocabulary and rppa_genes() is per-cohort; deriving the union here keeps the one
# definition of "a gene this tab knows" in the file that owns the map.
rppa_all_genes <- function(map = rppa_map()) sort(unique(unlist(map$gene_list)))

# list(text, tone, info) for the panel-level badge, or NULL. Delegates the STATE to
# depmap_note() -- the wording of the three states is defined once, in R/depmap.R, and this
# only replaces the hover explanation with the one that is true on this tab.
#
# The vocabulary check is a stop(), not a NULL, and it is the reason this wrapper exists at
# all. depmap_resolves() accepts any symbol-shaped string, so handing it an ANTIBODY --
# "mTOR_pS2448" -- returns the "not screened in this release" state rather than nothing:
# a confident sentence about a gene that does not exist, printed under a gene-first query.
# Silence would hide it; the app only ever passes a gene from rppa_genes(), so reaching
# this is a wiring defect and should say so.
rppa_depmap_note <- function(gene, lists, vocab = rppa_all_genes()) {
  if (length(gene) != 1L || is.na(gene) || !(gene %in% vocab))
    stop("rppa_depmap_note(): '", paste(gene, collapse = ", "), "' is not a gene in the ",
         "RPPA map. The badge is a statement about a GENE -- an antibody name reaching ",
         "DepMap would report 'not screened' for a gene that does not exist.")
  n <- depmap_note(gene, "rppa", lists)
  if (is.null(n)) return(NULL)
  n$info <- RPPA_DEPMAP_INFO
  n
}

# --- RPPA vs expression: is this the same patient set? (2026-08-26, step 82) -----------
#
# A reader who runs EGFR on the expression tab and then EGFR here is comparing two HRs and
# will assume they describe one cohort of patients. In four of the six staged tissues that
# assumption is close enough to true. In two it is badly false, and nothing on the tab said
# so until this step.
#
# Derived from the files, never transcribed -- measured 2026-08-26 at endpoint OS, over the
# patients the ENGINE would actually fit (finite time, event known, time > 0):
#
#   breast  872/875  99.7%      luad  350/354  98.9%      coad  345/346  99.7%
#   lusc    319/323  98.8%      ovarian 306/420 72.9%     gbm   127/231  55.0%
#
# GBM is the outlier the flag exists for: barely half of this panel's patients are in the
# GBM expression cohort at all, so "the RPPA result disagrees with the expression result"
# is, in GBM, only weakly a statement about the same disease in the same people.
#
# THE CUTOFF, AND WHY ITS VALUE DOES NOT MATTER. 0.90 is declared, not derived -- but the
# observed coverages are strongly bimodal (98.8-99.7% against 72.9% and 55.0%), so EVERY
# cutoff in (0.730, 0.988) classifies all six tissues identically. The constant is a label
# on a gap in the data, not a knife-edge, and tests/test_rppa_panel.R asserts that the gap
# is still there rather than trusting that it is -- if a future tissue lands mid-range the
# assertion fails and the cutoff has to be argued for on its merits instead.
RPPA_EXPR_OVERLAP_MIN <- 0.90

# How many of THIS PANEL's patients are also in the tissue's expression cohort of the same
# name. The denominator is the RPPA side on purpose: the question is whether an RPPA result
# can be read next to an expression one, and that is a property of the patients being fitted
# here. The endpoint matters -- rppa_clinical() drops rows with no usable time/event -- so
# it is a parameter rather than an assumption.
# `expr_h5` is a TEST SEAM, and it exists because of a guard that failed. The first version
# of the missing-matrix check below was tested by calling this with "lgg" -- which has no
# RPPA at all, so it tripped the stop() one line up and the assertion passed while the
# branch it named went untested. Mutation caught it (step 82). The two refusals are
# different claims and have to be reachable separately.
rppa_expr_overlap <- function(cancer_type, endpoint = "OS", expr_h5 = NULL) {
  cohort <- rppa_cohort_for(cancer_type)
  if (is.na(cohort))
    stop("rppa_expr_overlap(): no RPPA data staged for tissue '", cancer_type, "'")
  # The SAME eligibility rule get_survival() applies (R/survival_engine.R): finite time,
  # non-missing event, time > 0. Not a detail. rppa_clinical() keeps GBM's one patient at
  # time 0 and the engine drops them, so counting the raw frame printed "this panel's 232
  # patients" directly above a panel header reading "n = 231 patients" -- the note
  # contradicting the figure it annotates, which is the same defect as rounding 883/886 to
  # 100%. What this cannot promise is agreement with EVERY panel: an antibody with its own
  # missingness fits fewer, which is why the sentence names the COHORT and not the panel.
  cl <- rppa_clinical(cohort, endpoint)
  rs <- cl$sample_id[is.finite(cl$time) & !is.na(cl$event) & cl$time > 0]
  f  <- if (is.null(expr_h5)) file.path(PROC, sprintf("expr_%s.h5", tolower(cohort)))
        else expr_h5
  # Fail loud. A missing expression matrix here means the RPPA cohort has no same-named
  # expression cohort, and returning 0% would report a coverage crisis where the real
  # answer is that the comparison this note is about cannot be drawn at all.
  if (!file.exists(f))
    stop("rppa_expr_overlap(): ", basename(f), " does not exist, so the overlap between ",
         cohort, "'s RPPA and expression patients is not defined. This is a staging ",
         "question, not a coverage number -- do not report it as 0%.")
  es <- as.character(rhdf5::h5read(f, "/samples"))
  sh <- length(intersect(rs, es))
  list(cancer_type = cancer_type, cohort = cohort, endpoint = endpoint,
       n_rppa = length(rs), n_expr = length(es), n_shared = sh,
       frac = if (length(rs)) sh / length(rs) else NA_real_,
       low = if (length(rs)) sh / length(rs) < RPPA_EXPR_OVERLAP_MIN else NA)
}

# list(text, info, tone) in the shape rppa_depmap_note() uses. Returns a note in BOTH
# states rather than only the bad one: a panel that says nothing about coverage is
# indistinguishable from one that was never checked, and this tab's whole job is to be
# read alongside another tab's number.
rppa_overlap_note <- function(cancer_type, endpoint = "OS",
                              ov = rppa_expr_overlap(cancer_type, endpoint)) {
  # One decimal, ALWAYS. "%.0f" rendered breast's 883-of-886 as "100%", which is a number
  # the two counts beside it disprove -- a reader who checks the arithmetic finds the note
  # lying about the easiest thing on the panel to verify.
  pct <- sprintf("%.1f%%", 100 * ov$frac)
  if (isTRUE(ov$low))
    return(list(tone = "warn", info = RPPA_OVERLAP_INFO, text = sprintf(
      paste("Only %d of the %d %s patients measured here (%s) are also in that cohort's",
            "expression matrix, so this result and an expression result for the same gene",
            "are largely about DIFFERENT patients."),
      ov$n_shared, ov$n_rppa, ov$cohort, pct)))
  list(tone = "muted", info = RPPA_OVERLAP_INFO, text = sprintf(
    "%d of the %d %s patients measured here (%s) are also in that cohort's expression matrix.",
    ov$n_shared, ov$n_rppa, ov$cohort, pct))
}

RPPA_OVERLAP_INFO <- paste(
  "RPPA and expression are staged independently, so neither layer contains the other.",
  "Overlap is derived from the sample ids in the two matrices, not assumed.")

RPPA_PANEL_INFO <- paste(
  "RPPA measures protein with one antibody per row; phospho rows report a modified form.",
  "All rows share one patient set, so the rows cannot be pooled.")
