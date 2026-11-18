# Cohort registry: single source of truth for which cohorts exist, their tissue,
# platform, role, per-cohort stratifier, and available endpoints. Read once from
# config/cohorts.tsv. Replaces the old hardcoded SURV_COHORTS + strata(pam50).
# Requires ROOT from data_access.R (source it first).

.registry_path <- function() file.path(ROOT, "config", "cohorts.tsv")

# The grammar of the `stratify_by` field (per-endpoint since step 102).
#
#   "none"                     no stratifier, whatever endpoint is asked for
#   "pam50"                    that stratifier, on EVERY endpoint the cohort declares
#   "OS:idh;DFS:none;DSS:idh"  one answer per endpoint, and it must cover every declared
#                              endpoint -- a partial map is a registry error, not a default
#
# WHY THE THIRD FORM EXISTS. The field was per-COHORT and a cohort can sit in more than one
# pool. Setting `TCGA_LGG = idh` for the OS pool therefore reached the DFS pool too, where its
# only partner GSE107850_LGG is `none` -- and that made lgg/DFS the project's only MIXED pool,
# a pooled HR averaging one stratified and one unstratified per-cohort estimate. Everything
# downstream had learned to DISCLOSE that mixture (engine print(), strata_k / strata_mixed in
# the CSV, "1 of 2" on screen, the forest's axis clause) and `strata_suffix()` refused to NAME
# it. Disclosure is the right answer to a mixture that has to exist. This one did not: it was
# an artefact of the field's granularity, and per-endpoint scoping removes it at the source
# rather than annotating it at every surface. See HANDOFF item (b).
#
# STEP 103 SUPERSEDED THE PARAGRAPH THAT STOOD HERE, and it is worth knowing why. It read:
# "GSE107850_LGG stays `none` on the evidence, not for convenience: it carries IDH calls, but
# 166 mut / 14 wt / 15 uncalled is not a stratum worth having, and stratifying on it would also
# drop the 15 uncalled patients." That judged the stratum by its MARGINAL DISTRIBUTION rather
# than by its effect on the estimand, which is measurable -- and measured, it is wrong three
# ways: (1) GSE43107_LGG is stratified on 74.4% IDH coverage while GSE107850_LGG, at 92.3%, was
# the only lgg cohort of six left unstratified, so the "drops the uncalled" objection applied
# harder to a cohort stratified without comment; (2) the 14 wildtype carry 10 events (71% vs 48%
# in mutants), and the strata() term moves that arm's HR beyond what dropping the uncalled does;
# (3) with BOTH arms stratified the two cohorts agree -- I2 65.6% -> 0.0% on SALL2, 76.1% -> 0.0%
# on RUNX1 -- and TCGA_LGG's PH violation disappears (p 0.00128 -> 0.229). The I2 was differing
# IDH composition, not biology: a marginal HR in a mixed population is a weighted average over
# two hazards and the weights differed by cohort, so the marginal pool averaged two quantities
# that were not the same parameter.
#
# So ALL SIX lgg cohorts are `idh` now and TCGA_LGG is a bare token again. The per-endpoint
# grammar below is deliberately KEPT and still tested: it is the right shape for a genuine
# per-endpoint difference, and having no user in the registry is the correct outcome when the
# inconsistency it was expressing gets resolved instead of encoded. The COST is disclosed rather
# than netted out: stratifying drops 18 IDH-uncalled patients (706 -> 688, 286 -> 273 events) and
# that missingness is NOT random -- 12 of the 15 uncalled GEO patients have events.
.parse_stratify_by <- function(s) {
  if (length(s) != 1 || is.na(s) || !nzchar(trimws(s))) return(NULL)
  parts <- trimws(strsplit(s, ";")[[1]])
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(NULL)
  if (!any(grepl(":", parts, fixed = TRUE))) return(list(scoped = FALSE, all = trimws(s)))
  if (!all(grepl(":", parts, fixed = TRUE)))
    stop("stratify_by mixes scoped and bare terms: '", s,
         "'. Scope every endpoint or none of them.")
  kv <- do.call(rbind, strsplit(parts, ":", fixed = TRUE))
  m  <- stats::setNames(trimws(kv[, 2]), trimws(kv[, 1]))
  if (anyDuplicated(names(m)))
    stop("stratify_by names an endpoint twice: '", s, "'")
  list(scoped = TRUE, map = m)
}

load_registry <- function(path = .registry_path()) {
  df <- read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
  need <- c("cohort","cancer_type","platform","role","stratify_by","endpoints","network_dir")
  stopifnot(all(need %in% names(df)))
  # read.delim leaves a trailing empty field as NA; normalize to "" so the
  # "contributed no network" case is one value, not two.
  df$network_dir[is.na(df$network_dir)] <- ""
  # A scoped `stratify_by` must cover EXACTLY the endpoints the same row declares. Anything
  # less is a partial map, and a partial map is where a silent default would live: ask for the
  # missing endpoint and something has to be assumed. Checked at LOAD so a bad edit fails on
  # the next command rather than inside a fit, and checked in both directions so a stratifier
  # scoped to an endpoint the cohort does not carry is caught too.
  for (i in seq_len(nrow(df))) {
    p <- .parse_stratify_by(df$stratify_by[i])
    if (is.null(p) || !p$scoped) next
    want <- trimws(strsplit(df$endpoints[i], ",")[[1]])
    if (!setequal(names(p$map), want))
      stop(sprintf(paste0("config/cohorts.tsv: %s scopes stratify_by to %s but declares ",
                          "endpoints %s. A scoped field must name every declared endpoint ",
                          "and no others."),
                   df$cohort[i], paste(sort(names(p$map)), collapse = "/"),
                   paste(sort(want), collapse = "/")))
  }
  df
}
COHORTS <- load_registry()

.split_field <- function(x) trimws(strsplit(x, ",")[[1]])

# Cohort ids for a tissue, optionally requiring a supported endpoint and/or a role.
cohorts_for <- function(cancer_type, endpoint = NULL, roles = NULL, registry = COHORTS) {
  d <- registry[registry$cancer_type == cancer_type, , drop = FALSE]
  if (!is.null(endpoint))
    d <- d[vapply(d$endpoints, function(e) endpoint %in% .split_field(e), logical(1)), , drop = FALSE]
  if (!is.null(roles))
    d <- d[vapply(d$role, function(r) any(.split_field(r) %in% roles), logical(1)), , drop = FALSE]
  d$cohort
}

# Fail loud on a cohort id that is not in the registry. Callers used to hand cohort
# vectors straight to a SQL "IN (...)" or an h5 path, where an unknown id contributed
# nothing and disappeared -- so a stale id (e.g. "TCGA" after it became "TCGA_BRCA")
# silently downgraded a two-cohort request to one, and a typo produced an empty
# result that read as "this cohort has no data". Returns the input so it can guard
# inline: cohorts <- validate_cohorts(cohorts).
validate_cohorts <- function(cohorts, registry = COHORTS) {
  if (!length(cohorts)) stop("no cohorts requested (empty cohort vector)")
  unknown <- setdiff(cohorts, registry$cohort)
  if (length(unknown))
    stop(sprintf("unknown cohort id(s): %s. Known cohorts: %s",
                 paste(unknown, collapse = ", "),
                 paste(registry$cohort, collapse = ", ")))
  cohorts
}


# The stratifier THIS cohort is fitted under FOR THIS ENDPOINT, or NA_character_ for none.
#
# `endpoint` is optional only because most cohorts answer the same way for all of them. A
# cohort whose field is scoped and a caller who did not say which endpoint is a question with
# no single answer, so it STOPS rather than picking one -- the alternative is a fit adjusted
# by a variable the caller never chose, which is the whole defect this change removes. Ask
# stratifiers_of() instead when the question really is "does this cohort ever stratify".
stratifier_for <- function(cohort, endpoint = NULL, registry = COHORTS) {
  i <- match(cohort, registry$cohort)
  p <- .parse_stratify_by(registry$stratify_by[i])
  if (is.null(p)) return(NA_character_)
  v <- if (!p$scoped) p$all else {
    if (is.null(endpoint) || !length(endpoint) || is.na(endpoint[1]) || !nzchar(endpoint[1]))
      stop(sprintf(paste0("stratifier_for('%s'): this cohort declares a PER-ENDPOINT ",
                          "stratifier (%s) and no endpoint was given. Pass one -- there is ",
                          "no single answer, and guessing one would fit a model the caller ",
                          "did not choose."),
                   cohort, registry$stratify_by[i]))
    if (!endpoint[1] %in% names(p$map))
      stop(sprintf(paste0("stratifier_for('%s', '%s'): the registry maps stratifiers for %s ",
                          "and says nothing about %s. A partial map is a registry error; fix ",
                          "config/cohorts.tsv."),
                   cohort, endpoint[1], paste(names(p$map), collapse = "/"), endpoint[1]))
    unname(p$map[[endpoint[1]]])
  }
  if (is.na(v) || v == "none") return(NA_character_)
  v
}

# Every stratifier a cohort uses on ANY endpoint. A DIFFERENT question from stratifier_for(),
# and separated because conflating the two is what a silent default would do: this one answers
# "is there a stratification control to offer at all", which is a UI question, while
# stratifier_for() answers "what is this fit conditioned on", which is an estimand question.
stratifiers_of <- function(cohort, registry = COHORTS) {
  p <- .parse_stratify_by(registry$stratify_by[match(cohort, registry$cohort)])
  if (is.null(p)) return(character(0))
  v <- if (p$scoped) unname(p$map) else p$all
  sort(unique(v[!is.na(v) & v != "none"]))
}

cancer_type_of <- function(cohort, registry = COHORTS)
  registry$cancer_type[match(cohort, registry$cohort)]

# The endpoints ONE cohort carries, split and cleaned. The inverse of cohorts_for()'s
# endpoint filter: that answers "which cohorts carry OS", this answers "what does SCANB
# carry", which is what the app has to say when a cohort is offered struck through.
#
# It exists because the same expression was written out three times in app.R --
# .split_field(COHORTS$endpoints[match(co, COHORTS$cohort)]) -- once per note, each with
# its own idea of how to treat a blank field. A cohort declaring no endpoint at all is a
# legitimate registry state (a network-only contributor), so it returns character(0)
# rather than NA or "", and every caller can then just ask for its length.
#
# Scalar only, and loudly so: registry$endpoints[match(v, ...)] on a vector yields one
# string per cohort, and .split_field takes [[1]] -- so a vectorised call would silently
# report the FIRST cohort's endpoints for all of them.
cohort_endpoints <- function(cohort, registry = COHORTS) {
  if (length(cohort) != 1)
    stop(sprintf("cohort_endpoints() is scalar; got %d cohort ids", length(cohort)))
  e <- .split_field(registry$endpoints[match(cohort, registry$cohort)])
  e[!is.na(e) & nzchar(e)]
}

# Declared follow-up horizon (months) for a cancer type's panel, or NULL for "observe
# to the end of each cohort's own follow-up".
#
# A tissue whose cohorts observe very unequal windows cannot be pooled at full
# follow-up without averaging different estimands: the meta-HR then mixes a 4-year
# association from one cohort with a 20-year one from another. Administrative censoring
# at a common tau makes every cohort answer the same question. It is NOT a subset --
# patients still at risk at tau are censored there, never dropped.
#
# The horizon is a property of a (TISSUE, ENDPOINT) pair, not of any one cohort -- and
# NOT of a tissue alone once a tissue reports two endpoints on differently-observed pools.
# So PANEL_HORIZONS is keyed cancer_type -> named vector of per-endpoint taus. The RULE is
# identical everywhere: tau = floor of the SHORTEST cohort window in that endpoint's pool
# (min over cohorts of the cohort's max observed time), fixed from the windows BEFORE
# looking at any HR, deliberately not the tau with the best p-value. A horizon has to be
# declarable, not shoppable (same reasoning as the app's fixed presets).
#
# ovarian OS = 48: GSE49997 stops observing anyone at 48.99 months, so 48 is the longest
#   horizon at which all 13 OS cohorts are fully observed. (tau = 36 gives smaller p for
#   both anchors and was rejected for exactly that reason.)
# luad OS = 73 / DFS = 99; lusc OS = 87 / DFS = 138: same rule, applied PER ENDPOINT --
#   because a subtype's OS and DFS pools observe genuinely different windows and are
#   different estimands, so one tau per subtype would force the DFS analysis onto the
#   shorter OS window for no reason (the OS-limiting cohort, GSE3141, is OS-only and is
#   not even in the DFS pool). Shortest windows, per pool: LUAD OS 73.3 (GSE3141_LUAD),
#   LUAD DFS 99.2 (GSE50081_LUAD), LUSC OS 87.3 (GSE3141_LUSC), LUSC DFS 138.5
#   (GSE8894_LUSC). Each tau is floored to its own pool's shortest, so every pooled HR --
#   whichever endpoint -- sits on a fully-observed estimand with no cohort truncated below
#   tau. The panels are frozen, so a short-window cohort cannot be dropped to lengthen tau.
# coad OS = 142 / DFS = 103: same per-endpoint rule. Shortest windows, per pool:
#   OS 142.55 (GSE17538), DFS 103.06 (GSE37892). The two pools are DISJOINT in a way the
#   other tissues' are not -- GSE14333 and GSE17538 share 131 patients and are deliberately
#   assigned to different endpoints so neither pool contains a patient twice (see
#   scripts/colon_panel_lib.R) -- which is precisely why colon's taus must be derived per
#   pool: the OS-limiting cohort is not in the DFS pool at all, and vice versa.
# gbm OS = 40: same rule, applied to the gbm/OS pool, which is k=11 since CGGA landed
#   (TCGA_GBM + 7 GEO arms + 3 CGGA arms; this line read "k=8 ... + 7 GEO arms" until
#   2026-08-31, describing the pool as it stood before step 71).
#   Shortest window is GSE4412_GPL96_GBM at 40.97 months, so 40 is the longest horizon at
#   which every gbm cohort is fully observed. Derived from the STAGED cohorts (2026-08-24,
#   step 71), not from the audit table, and fixed before any gbm HR was computed. Short in
#   absolute terms but not for the disease: median OS in TCGA_GBM is 12.5 months, so 40
#   months is past the third quartile of survival, unlike ovarian's tau = 48 which sits mid-
#   curve. gbm has ONE endpoint here (OS): its DFS pool is k=0 -- no GEO brain series
#   anywhere deposits recurrence time -- so there is no second tau to derive.
# lgg OS = 94 / DFS = 80: same per-endpoint rule. The two pools still share only TCGA_LGG
#   -- GSE43107 is OS-only, GSE107850 is DFS-only -- but OS is now k=5 (TCGA_LGG, GSE43107
#   and 3 CGGA arms) against DFS at k=2, since CGGA is OS-only. Shortest windows, per pool:
#   OS 94.59 (GSE43107_LGG), DFS 80.49 (GSE107850_LGG). Derived 2026-08-24 (step 73) from
#   the STAGED cohorts after the engine's time > 0 filter, and fixed before any lgg HR.
#
#   THE CGGA CAVEAT IS DISCHARGED, and this is what it said: both taus were provisional
#   because CGGA was admitted to the plan but not staged, and when its arms landed the pools
#   would grow and each tau would have to re-derive from the new pool -- BEFORE any published
#   HR, never after. CGGA has landed. All ten declared taus were re-derived from the live
#   clinical tables on 2026-08-31 (step 119) and every one matches, gbm's 40 and lgg's 94/80
#   included: CGGA's arms observe LONGER than the limiting cohorts (lgg 155-169 mo against
#   GSE43107's 94.59; gbm 118-146 against GSE4412_GPL96's 40.97), so they never touched the
#   floor. That is the rule working, and it is now ASSERTED rather than promised --
#   tests/test_panel_horizon.R derives every pair from PANEL_HORIZONS and requires equality
#   in both directions. Until step 119 it hand-listed four tissues and omitted these two.
# breast OS = 81 / DFS = 281, since release 0.1.8 (review item G19): the same per-endpoint
#   rule as every other tissue. Shortest windows, per pool: OS 81.28 (SCANB, against
#   TCGA_BRCA 282.71 and METABRIC 355.20), DFS 281.10 (TCGA_BRCA, against METABRIC 389.33).
#   Until 0.1.7 breast was the one tissue outside the rule, exempted because its anchor HRs
#   were frozen at full follow-up, and the reason stated for the exemption ("comparable
#   windows across its cohorts") was false at OS: a 4.37x spread, wider than gbm's 3.56x,
#   which has a declared tau. Step 116 measured what the rule does to the best-known breast
#   result, ESR1/OS: pooled HR 1.107 (p 0.015) at full follow-up against 0.999 at 81 months,
#   METABRIC moving from 1.080 to 0.716 and I2 from 0% to 89.5%. That is the estimand the
#   rule asks for, a common window, and the paper reports it; the frozen anchors are retired.
#
# horizon_for(cancer_type, endpoint) returns that endpoint's tau. Called without an
# endpoint it stays backward-compatible: a single-endpoint tissue (ovarian) returns its
# one horizon; a multi-endpoint tissue returns the most CONSERVATIVE (shortest) tau, so a
# bare call is never wrong, only cautious. Every scan/figure driver passes the endpoint
# and gets the exact one; the app now passes the SELECTED endpoint too (via tau_choices),
# so the bare form survives only as a defensive fallback, no longer a live app path.
PANEL_HORIZONS <- list(
  breast  = c(OS = 81, DFS = 281),
  ovarian = c(OS = 48),
  luad    = c(OS = 73, DFS = 99),
  lusc    = c(OS = 87, DFS = 138),
  coad    = c(OS = 142, DFS = 103),
  gbm     = c(OS = 40),
  lgg     = c(OS = 94, DFS = 80)
)

horizon_for <- function(cancer_type, endpoint = NULL) {
  h <- PANEL_HORIZONS[[cancer_type]]
  if (is.null(h) || !length(h)) return(NULL)          # undeclared tissue -> full follow-up
  if (is.null(endpoint))
    return(unname(if (length(h) == 1) h[[1]] else min(h)))
  v <- h[endpoint]
  if (length(v) != 1 || is.na(v)) return(NULL)         # endpoint not declared for this tissue
  unname(v)
}

# Endpoints any cohort of this tissue actually carries, in canonical order. Derived from
# the registry, never hand-listed.
endpoints_of <- function(cancer_type, registry = COHORTS)
  Filter(function(e) length(cohorts_for(cancer_type, e, registry = registry)) > 0,
         c("OS", "DFS", "DSS"))

# Which endpoint a tissue's analysis treats as PRIMARY -- the one its stored artifacts
# (mr_scan_<ep>.csv, figures, survtables) are cut at, and the one the app opens on and
# defaults its tau to. Since release 0.1.8 (review item G18) it is stated as a RULE, and
# tests/test_registry_endpoints.R derives all seven from it and fails if this constant
# disagrees:
#   1. If exactly one endpoint reaches POOL_MIN_K (5) cohorts, the conventional minimum for
#      random-effects meta-analysis, it is primary: OS in ovarian, gbm, lusc and lgg.
#   2. Otherwise an endpoint the TCGA Pan-Cancer Clinical Data Resource advises caution with
#      is set aside (Liu et al. 2018, Table 3: OS in TCGA-BRCA needs longer follow-up), and
#      the endpoint with more events in the tissue's TCGA cohort is chosen: DFS in breast,
#      luad and coad.
# HONEST HISTORY, so the rule is not mistaken for a pre-registration: these values were chosen
# one tissue at a time as each landed (breast DFS at step F after its first scans; coad on
# 2026-07-29 before any colon HR; lgg OS at step 73 on the CGGA ceiling) and the rule was
# written down afterwards, in 0.1.8, as the one statement that reproduces all seven. The
# paper states the rule and does not claim it was fixed in advance. Choosing a primary
# endpoint deletes nothing: every endpoint stays registered, scanned and queryable.
PRIMARY_ENDPOINT <- c(breast = "DFS", ovarian = "OS", luad = "DFS", lusc = "OS",
                      coad = "DFS", gbm = "OS", lgg = "OS")

# Human display label for a cancer type. The registry key is a lowercase slug (breast,
# ovarian, luad, lusc); the label is an EDITORIAL choice, not derivable from it -- full
# disease name, capitalised, with the TCGA-style abbreviation in parentheses where the
# slug itself IS the abbreviation (luad/lusc). Lives here beside PRIMARY_ENDPOINT so the
# app carries no cancer-type knowledge of its own. cancer_label_of() falls back to the
# raw slug for an unlisted tissue, so a newly added one still shows (labelled by slug),
# the same graceful degradation as the app's scan-label map.
CANCER_LABEL <- c(breast = "Breast", ovarian = "Ovarian",
                  luad = "Lung adenocarcinoma (LUAD)",
                  lusc = "Lung squamous cell carcinoma (LUSC)",
                  coad = "Colon adenocarcinoma (COAD)",
                  gbm  = "Glioblastoma (GBM)",
                  lgg  = "Lower-grade glioma (LGG)")
cancer_label_of <- function(cancer_type) {
  lab <- unname(CANCER_LABEL[cancer_type])
  ifelse(is.na(lab), cancer_type, lab)
}

# The primary endpoint for a tissue when it is declared AND its cohorts carry it, else the
# tissue's first available endpoint (so an undeclared tissue still resolves to something the
# UI can offer). NA only when the tissue has no endpoint at all.
primary_endpoint_of <- function(cancer_type) {
  eps <- endpoints_of(cancer_type)
  p   <- unname(PRIMARY_ENDPOINT[cancer_type])
  if (!is.na(p) && p %in% eps) p else if (length(eps)) eps[1] else NA_character_
}

# Resolve a REQUESTED endpoint against a tissue: keep it if the tissue carries it, else
# fall back to the tissue's declared primary.
#
# This is the rule the app's cohort selector runs on, and both of its inputs are states
# the Shiny UI actually passes. `endpoint_ui` is itself a renderUI, so on the first pass
# input$endpoint is NULL; and for one reactive beat after a tissue switch it still holds
# the PREVIOUS tissue's value -- DFS while the tissue is now ovarian, which carries none.
# Left unresolved, that one beat decides which cohort chips are live and pre-ticked, so
# the answer has to be a declared endpoint rather than whatever the widget happens to say.
#
# NULL, NA, character(0) and a foreign endpoint all resolve the same way, deliberately:
# every one of them means "the widget has no opinion this tissue can honour".
resolve_endpoint <- function(endpoint, cancer_type) {
  ok <- length(endpoint) == 1 && !is.na(endpoint) &&
        endpoint %in% endpoints_of(cancer_type)
  if (ok) endpoint else primary_endpoint_of(cancer_type)
}

# Regulon .rds paths for a tissue's own networks. Scoping this by cancer_type is
# what keeps a breast scan's n_targets annotation off the ovarian regulon — before
# the registry carried network_dir, callers globbed networks/*/ and silently picked
# up every tissue's network once ovarian gained one.
regulon_paths_for <- function(cancer_type, registry = COHORTS) {
  dirs <- registry$network_dir[registry$cancer_type == cancer_type &
                               nzchar(registry$network_dir)]
  p <- file.path(ROOT, "networks", dirs, sprintf("regulon_%s.rds", dirs))
  p[file.exists(p)]
}

# ---- the regulon drawer (release 0.1.11, item 4) -------------------------------------------
#
# What the drawer under a queried TF shows, as data; app.R only lays it out. The networks are
# exactly the ones behind the tissue's VIPER scores: breast's three cohorts are scored by
# metaVIPER with BOTH breast regulons (R/run_metaviper.R), every other tissue with its one TCGA
# network, and gse9891 (the ovarian replication network) scores nothing the app reads.
# tests/test_regulon_drawer.R proves that from the VIPER matrices themselves, cohort by cohort.
#
# Wording is part of the data: "inferred targets", "inferred activated / repressed". ARACNe
# infers a target from co-expression, and the sign is the direction of that co-expression.

# A reader's name for each network directory. Every network_dir in the registry must have one.
NETWORK_LABEL <- c(tcga_brca = "TCGA-BRCA", scanb = "SCAN-B", tcga_ov = "TCGA-OV",
                   tcga_luad = "TCGA-LUAD", tcga_lusc = "TCGA-LUSC", tcga_coad = "TCGA-COAD",
                   tcga_gbm = "TCGA-GBM", tcga_lgg = "TCGA-LGG")

REGULON_NOTE <- paste(
  "ARACNe infers these targets from co-expression across tumours; activated or repressed is",
  "the direction of that co-expression, not evidence of binding or of an experiment.")

# The network directories behind a tissue's VIPER scores (the ones on disk), in registry order.
regulon_networks_for <- function(cancer_type, registry = COHORTS)
  basename(dirname(regulon_paths_for(cancer_type, registry)))

# One TF's targets in one regulon, as a data frame ranked by |tfmode| x likelihood.
.regulon_targets <- function(x) {
  tm <- unname(x$tfmode); lk <- unname(x$likelihood)
  d <- data.frame(target = names(x$tfmode), tfmode = tm, likelihood = lk,
                  score = abs(tm) * lk, stringsAsFactors = FALSE)
  # A mode of exactly 0 has no direction; none of the eight regulons holds one (step 322), and
  # one would be a claim this drawer cannot word, so it stops.
  if (any(d$tfmode == 0)) stop(".regulon_targets(): a target with tfmode 0 has no direction")
  d$direction <- ifelse(d$tfmode > 0, "inferred activated", "inferred repressed")
  d[order(-d$score), , drop = FALSE]
}

# The drawer's data for `tf` in `cancer_type`, or NULL when `tf` is a regulator in none of the
# tissue's networks (a non-TF, a probe set, an empty selector). `absent_from` names a network of
# the tissue in which `tf` is not a regulator (breast can have one regulon and not the other).
regulon_of <- function(tf, cancer_type, registry = COHORTS, top = 10L) {
  if (length(tf) != 1L || is.na(tf) || !nzchar(tf)) return(NULL)
  dirs <- regulon_networks_for(cancer_type, registry)
  nets <- list(); absent <- character(0)
  for (d in dirs) {
    R <- readRDS(regulon_paths_for(cancer_type, registry)[dirs == d])
    if (is.null(R[[tf]])) { absent <- c(absent, d); next }
    tg <- .regulon_targets(R[[tf]])
    up <- tg[tg$tfmode > 0, , drop = FALSE]; dn <- tg[tg$tfmode < 0, , drop = FALSE]
    # The TFs that have `tf` among their own inferred targets, strongest first.
    rg <- do.call(rbind, lapply(setdiff(names(R), tf), function(t) {
      m <- R[[t]]$tfmode
      if (!(tf %in% names(m))) return(NULL)
      i <- match(tf, names(m)); v <- unname(m[i]); l <- unname(R[[t]]$likelihood[i])
      data.frame(tf = t, tfmode = v, likelihood = l, score = abs(v) * l,
                 direction = if (v > 0) "inferred activated" else if (v < 0) "inferred repressed" else stop("tfmode 0"),
                 stringsAsFactors = FALSE)
    }))
    if (is.null(rg)) rg <- data.frame(tf = character(0), tfmode = numeric(0), likelihood = numeric(0),
                                      score = numeric(0), direction = character(0))
    rg <- rg[order(-rg$score), , drop = FALSE]
    nets[[d]] <- list(dir = d, label = NETWORK_LABEL[[d]], n = nrow(tg), n_act = nrow(up),
                      n_rep = nrow(dn), top_act = utils::head(up, top), top_rep = utils::head(dn, top),
                      all = tg, regulators = utils::head(rg, top), n_regulators = nrow(rg))
  }
  if (!length(nets)) return(NULL)
  list(tf = tf, cancer_type = cancer_type, networks = nets, absent_from = absent)
}

regulon_summary_line <- function(net, tf)
  sprintf("Regulon of %s in the %s network: %d inferred targets, %d activated, %d repressed.",
          tf, net$label, net$n, net$n_act, net$n_rep)

# Every target of every network shown, one row each; `network` names the network on each row.
regulon_csv <- function(rg) {
  out <- do.call(rbind, lapply(rg$networks, function(net)
    data.frame(network = net$label, tf = rg$tf, target = net$all$target,
               direction = net$all$direction, tfmode = net$all$tfmode,
               likelihood = net$all$likelihood, score = net$all$score,
               stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}
regulon_csv_name <- function(rg)
  sprintf("regulon_%s_%s_%s.csv", gsub("[^A-Za-z0-9._-]+", "_", rg$tf), rg$cancer_type,
          paste(names(rg$networks), collapse = "_"))

# Every sentence the drawer prints, worded here so tests/test_regulon_drawer.R can hold them
# (app.R only lays them out). `integrated` is NULL for a tissue scored with one network.
regulon_sentences <- function(rg, top = 10L) {
  tf <- rg$tf
  dirs <- c(names(rg$networks), rg$absent_from)
  labs <- unname(NETWORK_LABEL[dirs])
  nets <- lapply(rg$networks, function(net) {
    n <- net$n_regulators
    # "other" (step 327, the author's wording): `tf` itself is never counted.
    reg <- if (n == 0) sprintf("No other TF in the %s network has %s among its inferred targets.", net$label, tf)
           else if (n == 1) sprintf("1 other TF in the %s network has %s among its inferred targets.", net$label, tf)
           else sprintf("%d other TFs in the %s network have %s among their inferred targets%s.", n, net$label, tf,
                        if (n > top) sprintf("; the %d strongest are shown", top) else "")
    list(summary = regulon_summary_line(net, tf), regulators = reg)
  })
  list(title = sprintf("Regulon of %s", tf),
       integrated = if (length(dirs) > 1)
         sprintf("%s VIPER scores combine the %s networks (metaVIPER), so each regulon is shown.",
                 cancer_label_of(rg$cancer_type),
                 paste(c(paste(labs[-length(labs)], collapse = ", "), labs[length(labs)]), collapse = " and ")),
       networks = nets,
       absent = vapply(rg$absent_from, function(d)
         sprintf("%s is not a regulator in the %s network.", tf, NETWORK_LABEL[[d]]), character(1),
         USE.NAMES = FALSE),
       note = REGULON_NOTE)
}
