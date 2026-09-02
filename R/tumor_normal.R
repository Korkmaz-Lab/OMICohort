# Tumor vs matched normal-adjacent tissue -- the tool's one NON-SURVIVAL contrast.
#
# What it answers: is this feature higher in the tumor than in the same patient's own
# normal tissue? That is context for a gene, in the same family as the CPTAC
# detectability badge -- NOT a survival result. Nothing here produces an HR, and nothing
# here may reach a Cox fit (see the boundary note below).
#
# It costs no new data. TCGA's solid cohorts already ship their normal-adjacent samples
# inside the same expr_/viper_/immune_ h5s the survival path reads: 113 in TCGA_BRCA,
# 59 TCGA_LUAD, 51 TCGA_LUSC, 41 TCGA_COAD, 0 TCGA_OV. get_feature() returns them
# already -- it is .dedup_tcga(), at the survival boundary, that drops them.
#
# THE BOUNDARY, stated because the two paths read the same file with OPPOSITE rules.
# get_clinical() keeps ONE sample per patient (primary tumor preferred) because counting a
# patient twice is pseudo-replication in a Cox model. This file deliberately keeps two --
# a tumor AND a normal for the same patient -- which is the whole point of a paired
# contrast and would be a defect anywhere near survival. So: this file never calls
# get_survival(), never returns anything a survival function accepts, and its output
# carries a class of its own. tests/test_tumor_normal.R asserts that.
#
# WHY PAIRED, and why the paired difference is the number to read. The expr matrices are
# z-scored WITHIN cohort over every sample the cohort has -- normals included. So a
# normal's z encodes how many normals happen to be in the matrix: TCGA_BRCA's normal
# median sits near -1.8 partly because 113 of its 1226 samples are normal. Change that
# count and the number moves. The per-patient difference does not: each patient is their
# own control, so it is invariant to cohort composition. The two boxes are for looking;
# the paired difference is the quantity.
#
# WHY IT RE-STANDARDIZES, on the TUMOR samples. Only expr arrives z-scored; VIPER is raw
# NES (TCGA_BRCA/FOXM1: mean -0.51, sd 7.05), which is why get_survival() standardizes at
# fit time (survival_engine.R) so every HR in this tool is a per-SD HR. Left raw, this
# panel would report +1.87 for expression and +12.30 for VIPER activity -- two numbers in
# different units, side by side, looking comparable. So the same convention is applied
# here, and the reference population is the cohort's TUMOR samples: that makes one unit on
# this panel the SAME unit as the "per SD" in the forest above it. Standardizing over
# tumors+normals instead would inflate the SD by the very contrast being measured.
#
# Requires data_access.R (get_feature, SURV_COHORTS, .dedup rules) sourced first.

# TCGA barcode: TCGA-XX-XXXX(1:12)-<type 14:15><vial 16>. Sample type 01-09 is tumor,
# 10-19 normal, 20-29 control. Only "11" is SOLID TISSUE normal; "10" is blood-derived
# normal, which is a different tissue and must never be pooled with it -- comparing tumor
# expression against blood would be a tissue contrast wearing a tumor/normal label.
TN_TUMOR_TYPES  <- sprintf("%02d", 1:9)
TN_NORMAL_TYPE  <- "11"

# --- the fold change ------------------------------------------------------------------
#
# The panel's AXIS is per-SD (see the header), and stays that way. This is the companion
# NUMBER, because an SD is not interpretable on its own: it is scaled by how heterogeneous
# that cohort's tumors happen to be, so +2.2 SD is a different amount of biology for a
# tight gene than for a variable one. A ratio reads the same everywhere.
#
# It CANNOT be computed from the plotted values. Every TCGA expr matrix here is a z-score of
# DESeq2 VST, and a VST difference is not a log2 fold change -- measured against TPM on these
# exact pairs it runs 9-54% high and can flip sign near zero (scripts/stage_tn_foldchange.py
# records the numbers). So it is staged separately from GDC TPM and READ here, which is why
# .tn_fc() checks that the file it reads describes the same pairs this function just built.
TN_FC_KIND <- "expr"
.tn_fc_path <- function(cohort) file.path(PROC, sprintf("tnfc_%s.h5", tolower(cohort)))

# One gene's staged fold change, with its REASON when there is none -- four states, kept
# apart for the same reason $status is (a limitation must not read as "no change"):
#   ok         -- median paired log2(TPM+1) difference
#   ambiguous  -- the symbol has no single Ensembl gene, so no fold change was staged for it
#                 (45 of 20560 in breast; staging drops them rather than picking a paralog)
#   no_sidecar -- this cohort has no tnfc_ file
#   not_expr   -- VIPER NES and immune scores have no abundance underneath them to ratio
.tn_fc <- function(feature, cohort, kind, n_pair) {
  # Always BOTH fields, so a caller can never index a missing one into a silent NULL.
  no <- function(s) list(fc_status = s, log2fc = NA_real_)
  if (!identical(kind, TN_FC_KIND)) return(no("not_expr"))
  f <- .tn_fc_path(cohort)
  if (!file.exists(f)) return(no("no_sidecar"))
  idx <- match(feature, as.character(h5read(f, "/features")))
  if (is.na(idx)) return(no("ambiguous"))
  # The sidecar is a PRECOMPUTED statistic over a pairing this function rebuilds from the
  # expr h5 every call. Those two can only agree by construction, never by luck, so the one
  # thing that would make the printed number describe different patients than the drawn boxes
  # is asserted rather than assumed. Staging is cheap (15s); a silently stale sidecar is not.
  np <- h5readAttributes(f, "/")$n_pair
  if (!length(np) || as.integer(np) != as.integer(n_pair))
    stop(sprintf(paste("tumor_normal(): %s was staged over %s pairs but %s has %d here.",
                       "Re-run scripts/stage_tn_foldchange.py -- the fold change and the",
                       "plotted difference would describe different patients."),
                 basename(f), if (length(np)) as.character(np) else "an unrecorded number of",
                 cohort, n_pair))
  list(fc_status = "ok", log2fc = as.numeric(h5read(f, "/log2fc", index = list(idx))))
}

# Kinds this contrast is meaningful for. "cna" is EXCLUDED on purpose and not for lack of
# data: copy number is called against a matched normal reference, so a tumor-vs-normal CNA
# comparison is circular by construction. cna_*.h5 accordingly carries no normal samples,
# which means a silent path would report "no normal-adjacent samples" -- a true sentence
# that answers the wrong question, and exactly the quiet default this project bans.
TN_KINDS <- c("expr", "viper", "immune")

# Whether a set of sample ids is TCGA-barcoded, i.e. whether sample type is KNOWABLE.
# Derived from the ids, never from a list of cohort names: a TCGA cohort added to the
# registry tomorrow works with no edit here, and a non-TCGA one cannot be mislabelled.
.tn_is_tcga <- function(ids) length(ids) > 0L && all(grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z]", ids))

# The one sample per patient that the SURVIVAL path would have used, by the same rule as
# .dedup_tcga(): primary tumor first, then lowest vial letter. Reused rather than
# reinvented so the tumor box shows the value the Cox model actually saw for that patient.
.tn_pick <- function(ids) {
  if (!length(ids)) return(NA_character_)
  prio <- ifelse(substr(ids, 14, 15) == "01", 0L, 1L)
  ids[order(prio, substr(ids, 16, 16))][1]
}

# Tumor vs matched normal-adjacent, per cohort.
#
# Returns an object of class "tumor_normal": $feature, $kind, and $per_cohort, a named
# list with one entry per requested cohort carrying $status and, where estimable, the
# tumor/normal value vectors and the paired statistics.
#
# $status is FOUR distinct states, kept apart because collapsing them is how a limitation
# turns into a finding:
#   ok         -- has solid normals and at least one tumor/normal pair
#   no_normals -- TCGA cohort, zero solid-normal samples (TCGA_OV)
#   unlabelled -- sample ids carry no sample type, so tumor vs normal CANNOT be determined
#                 (METABRIC, SCAN-B, every GEO series). This is NOT "no normals": it is
#                 "this cohort cannot answer the question". clinical.db has no sample-type
#                 column, so the barcode is the only label that exists.
#   absent     -- the feature is not in this cohort's vocabulary, or the h5 is missing
tumor_normal <- function(feature, cohorts = SURV_COHORTS, kind = "expr", standardize = TRUE) {
  # Validated by hand rather than by match.arg() so the REASON survives. match.arg would
  # reject "cna" with "'arg' should be one of ..." -- true, and it teaches the caller
  # nothing about why copy number is excluded. (Written as match.arg first; the explicit
  # check below was dead code behind it until this was noticed.)
  if (length(kind) != 1L || !(kind %in% TN_KINDS))
    stop("tumor_normal(): no tumor/normal contrast for kind '", paste(kind, collapse = ", "),
         "'. Supported: ", paste(TN_KINDS, collapse = ", "), ". ",
         if (identical(kind, "cna"))
           paste("Copy number is called against a matched normal reference, so a",
                 "tumor-vs-normal CNA comparison is circular by construction.")
         else "")
  if (length(feature) != 1L || is.na(feature) || !nzchar(feature))
    stop("tumor_normal(): need exactly one feature name")
  cohorts <- validate_cohorts(cohorts)

  per <- lapply(cohorts, function(coh) {
    v <- suppressWarnings(get_feature(feature, coh, kind = kind))
    if (!length(v)) return(list(cohort = coh, status = "absent"))
    ids <- names(v)
    if (!.tn_is_tcga(ids)) return(list(cohort = coh, status = "unlabelled"))

    ty <- substr(ids, 14, 15)
    tum_ids <- ids[ty %in% TN_TUMOR_TYPES]
    nrm_ids <- ids[ty == TN_NORMAL_TYPE]
    # Per-SD, on the tumor reference (see header). Guarded on a positive finite SD for the
    # same reason get_survival() guards its own: a constant feature would otherwise divide
    # by zero and hand every downstream number a NaN that still plots.
    if (standardize && length(tum_ids)) {
      sdev <- stats::sd(v[tum_ids])
      if (is.finite(sdev) && sdev > 0) v <- (v - mean(v[tum_ids])) / sdev
    }
    # Anything that is neither is counted rather than dropped in silence -- a blood normal
    # or a control sample appearing here is a staging change this panel must not absorb.
    other_n <- length(ids) - length(tum_ids) - length(nrm_ids)

    base <- list(cohort = coh, n_tumor = length(tum_ids), n_normal = length(nrm_ids),
                 n_other = other_n, tumor = unname(v[tum_ids]), normal = unname(v[nrm_ids]))
    if (!length(nrm_ids)) return(c(base, list(status = "no_normals")))

    pt_t <- substr(tum_ids, 1, 12); pt_n <- substr(nrm_ids, 1, 12)
    both <- sort(intersect(pt_t, pt_n))
    if (!length(both)) return(c(base, list(status = "no_pairs", n_pair = 0L)))
    t_pick <- vapply(both, function(p) .tn_pick(tum_ids[pt_t == p]), character(1))
    n_pick <- vapply(both, function(p) .tn_pick(nrm_ids[pt_n == p]), character(1))
    d <- unname(v[t_pick] - v[n_pick])
    c(base, list(status = "ok", n_pair = length(both), patients = both,
                 diff = d, median_diff = stats::median(d),
                 p = stats::wilcox.test(d)$p.value),
      .tn_fc(feature, coh, kind, length(both)))
  })
  names(per) <- cohorts
  structure(list(feature = feature, kind = kind, per_cohort = per),
            class = "tumor_normal")
}

# The cohorts this contrast is estimable in, in the order given. The UI asks this before
# it draws anything: an empty answer means the panel is not shown at all, rather than a
# figure with no boxes in it.
tn_estimable <- function(tn) {
  stopifnot(inherits(tn, "tumor_normal"))
  names(Filter(function(x) identical(x$status, "ok"), tn$per_cohort))
}

# The sentence that keeps this panel from being read as an outcome result. Stated as a
# constant, not typed into app.R, so a test can assert the app still shows it.
TN_CAPTION <- paste("Tumor vs the same patient's own normal tissue. Context for the gene,",
                    "not a survival result. One unit is one SD of this cohort's tumors, the",
                    "same unit as the HR on the left.")

# Why each cohort that CANNOT be drawn is missing, ONE sentence per reason rather than per
# cohort: on breast, three separate lines saying the same thing about SCANB and METABRIC is
# noise that gets skipped, and the fact that must land -- that the two largest breast
# cohorts are silent here, rather than agreeing -- is the one most likely to be lost in it.
# Cohorts that can be drawn contribute nothing.
# Naming the cohorts is the informative form, but ovarian selects thirteen and a
# twelve-name list is a wall the reader skips. Past TN_NAME_MAX the count says it better.
TN_NAME_MAX <- 4L
.tn_join <- function(x) {
  if (length(x) == 1L) return(x)
  if (length(x) > TN_NAME_MAX) return(sprintf("%d of the selected cohorts", length(x)))
  paste(paste(utils::head(x, -1L), collapse = ", "), "and", utils::tail(x, 1L))
}
tn_notes <- function(tn) {
  stopifnot(inherits(tn, "tumor_normal"))
  st <- vapply(tn$per_cohort, `[[`, character(1), "status")
  co <- vapply(tn$per_cohort, `[[`, character(1), "cohort")
  # Singular and plural written out in full. The first version conjugated by suffix
  # ("carrie" + "s"/""), which happens to work for carries/carry and produces "doe not"
  # and "ha normal-adjacent samples" for everything else -- correct on the one cohort it
  # was eyeballed against, wrong on every list of two or more.
  say <- function(status, one, many) {
    who <- co[st == status]
    if (!length(who)) return(NULL)
    sprintf(if (length(who) == 1L) one else many, .tn_join(who))
  }
  out <- c(
    say("no_normals",
        "%s carries no normal-adjacent samples.",
        "%s carry no normal-adjacent samples."),
    say("no_pairs",
        "%s has normal-adjacent samples but none from a patient who also has a tumor, so no paired comparison is possible.",
        "%s have normal-adjacent samples but none from a patient who also has a tumor, so no paired comparison is possible."),
    say("unlabelled",
        "%s does not label tumor vs normal (sample ids are not TCGA barcodes), so it cannot be compared here.",
        "%s do not label tumor vs normal (sample ids are not TCGA barcodes), so they cannot be compared here."),
    say("absent",
        "%s does not carry this feature.",
        "%s do not carry this feature."))
  unknown <- setdiff(st, c("ok", "no_normals", "no_pairs", "unlabelled", "absent"))
  if (length(unknown))
    stop("tn_notes(): unknown status '", paste(unknown, collapse = "', '"), "'")
  out
}

# Why the fold change is missing when the panel is drawn but the number is not there.
#
# Separate from tn_notes(), which is about COHORTS: this is about the gene. Only "ambiguous"
# gets a sentence. "not_expr" needs none -- a VIPER activity panel has no abundance under it
# to take a ratio of, and nobody reads a NES box wondering where the fold change went --
# and "no_sidecar" cannot happen for a drawn cohort without the staging step having been
# skipped, which the n_pair assertion in .tn_fc() already turns into a hard error.
tn_fc_note <- function(tn) {
  stopifnot(inherits(tn, "tumor_normal"))
  est <- tn_estimable(tn)
  if (!length(est)) return(NULL)
  amb <- vapply(tn$per_cohort[est],
                function(x) identical(x$fc_status, "ambiguous"), logical(1))
  if (!any(amb)) return(NULL)
  sprintf(paste("No fold change for %s: more than one Ensembl gene carries this symbol, and",
                "staging drops those rather than attributing one gene's ratio to another's",
                "name. The per-SD difference above is unaffected."), tn$feature)
}
