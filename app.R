# Shiny front-end for the multi-cancer survival tool. Thin wrapper: no new science
# here, every number comes straight from get_feature()/get_survival()/the existing
# batch scan CSVs. Tabs:
#   - Single query:   live, on-demand single-feature lookup in ONE tissue (any gene/TF,
#                     any endpoint, any horizon).
#   - Multiple query: the same feature run in EVERY tissue at once, one panel each.
#                     Five INDEPENDENT get_survival() calls, never a merge -- see the
#                     block comment above the Multiple query tab.
#   - Browse: the batch scan / discovery-validation tables for the selected tissue,
#             with Step A/B annotations joined in where they exist (breast).
#
# Everything tissue-specific is resolved from the cohort registry (config/cohorts.tsv)
# — cohort list, endpoints, KM tabs, stratifier, feature vocabulary, scan tables. The
# app has NO hardcoded cohort or cancer-type knowledge: even the primary-endpoint judgement
# now lives in the registry (PRIMARY_ENDPOINT / primary_endpoint_of in R/registry.R), the
# single source of truth, so it cannot drift from the scan/figure drivers as it once did.
#
# Run:  shiny::runApp(".")   (from the project root, so the R/ and results/ relative
# paths below resolve — Shiny sets the working directory to the app's own directory).

options(shiny.autoload.r = FALSE)  # R/ is this project's plain function library, not
                                    # shiny-module support code - don't auto-source it all
                                    # (several files there are CLI scripts, not safe to source blindly)
suppressPackageStartupMessages({ library(shiny); library(DT) })
# Bootswatch "simplex" via shinythemes, if installed. Guarded so the app still runs
# (in stock Bootstrap) where the package is absent -- it is a look, not a dependency.
HAVE_SHINYTHEMES <- requireNamespace("shinythemes", quietly = TRUE)

source("R/data_access.R")
source("R/survival_engine.R")
source("R/plots.R")
source("R/scan_lookup.R")   # Query-tab genome-wide rank lookup (multiplicity context)
# BEFORE guide.R: GUIDE_CSS is built at source time from BRAND_* rather than substituted
# afterwards, so the constants have to exist first. It fails loudly here if they do not,
# which is the intent -- a stylesheet assembled from a missing colour is not a thing to
# discover in a browser.
source("R/brand.R")         # the two brand colours, measured from logo.png (step 112/113)
source("R/guide.R")         # Guide tab -- every figure on it DERIVED, never typed (step 109)
source("R/about.R")         # About tab -- provenance, licence, people (step 110)
source("R/tumor_normal.R")  # Single-query tumor vs normal-adjacent panel (NOT a survival result)
source("R/cptac.R")         # CPTAC-BRCA protein detectability (isoform-collapsed; NOT survival)
source("R/depmap.R")        # DepMap common-essentiality flag (cell lines; NOT survival)
source("R/rppa.R")          # RPPA protein/phospho survival -- SEPARATE layer, k=1, never pooled

# ---- tissue vocabulary, all resolved from the registry ---------------------------

CANCER_TYPES <- unique(COHORTS$cancer_type)

# The tissue's primary endpoint (an editorial judgement) and the endpoint list both come
# from the registry now — PRIMARY_ENDPOINT / primary_endpoint_of / endpoints_of live in
# R/registry.R so the scan drivers, the figure drivers and this app all read ONE copy.
# `default_endpoint_of` is kept as a local alias so the call sites below read unchanged.
cohorts_of <- function(ct) cohorts_for(ct)
default_endpoint_of <- primary_endpoint_of

# Does this tissue stratify at all, and BY WHAT? Read from config/cohorts.tsv, never
# assumed: breast is pam50, lgg is idh (since 2026-08-27), everything else is none.
# Drives whether the stratify checkbox is shown — offering it where every cohort is
# stratify_by=none would imply a control that silently does nothing.
# stratifiers_of(), not stratifier_for(): this answers "is there a control to OFFER",
# which spans endpoints, while stratifier_for() answers "what is THIS fit conditioned
# on" and since step 102 stops if asked without one. The engine does the per-endpoint
# resolution, so a tissue that stratifies on some endpoint keeps its checkbox and a
# pool where no cohort stratifies simply comes back with strata_k = 0.
# `ep` is not optional in spirit, only in signature. WITH an endpoint this resolves the pool
# that will actually be fitted -- cohorts_for(ct, ep), asked at ep -- which is the only
# answer that can honestly decide whether to OFFER a stratification control. WITHOUT one it
# unions across endpoints, which is right for the landing page's "which tissues ever
# stratify" summary and wrong everywhere a query exists.
#
# Step 102 shipped the union everywhere for one commit and the app immediately showed the
# defect: lgg/DFS offered "Stratify by IDH mutation status" on a pool where no cohort is
# stratified on DFS -- a control that silently does nothing, which is the exact thing the
# comment above output$strata_ui says not to build.
strata_vars_of <- function(ct, ep = NULL) {
  if (is.null(ep) || !length(ep) || !nzchar(ep[1]))
    return(unique(unlist(lapply(cohorts_of(ct), stratifiers_of))))
  co <- tryCatch(cohorts_for(ct, ep[1]), error = function(e) character(0))
  if (!length(co)) return(character(0))
  unique(Filter(Negate(is.na),
                vapply(co, stratifier_for, character(1), endpoint = ep[1])))
}

# Human names for the registry's stratifier VALUES.
#
# This map lives here and not in R/registry.R because it is a UI string and nothing outside
# the app wants it: the engine prints the raw column name (`strata(idh)`), which is exactly
# what a methods reader needs, and R/plots.R prints yes/no. CANCER_LABEL is in the registry
# because the engine and the plots both use it; this is not that.
#
# An unlisted stratifier falls back to its raw name rather than raising. That is a
# deliberate exception to fail-loud, and the reason is where the code runs: these are
# called inside renderUI, so an error blanks the panel. "Stratify by grade" is plain but
# never WRONG, whereas a blank control over a missing adjective is a real regression. The
# thing that MUST NOT happen — a label naming the wrong variable — is impossible here,
# because the name is derived from the same registry the engine reads.
STRATA_LABEL <- c(pam50 = "PAM50 subtype", idh = "IDH mutation status", stage = "stage")

# The ONE stratifier to name a tissue's adjusted scan file by, or NULL when there is not
# exactly one. Passed to scan_rank_lookup(), which refuses to guess a filename without it.
strata_var_for_scan <- function(ct, ep = NULL) {
  sv <- strata_vars_of(ct, ep)
  if (length(sv) == 1) sv else NULL
}

strata_label_of <- function(ct, ep = NULL) {
  sv <- strata_vars_of(ct, ep)
  if (!length(sv)) return(NULL)
  nm <- ifelse(sv %in% names(STRATA_LABEL), unname(STRATA_LABEL[sv]), sv)
  sprintf("Stratify by %s", paste(nm, collapse = " / "))
}

# Does the checkbox reach every cohort that will ACTUALLY BE POOLED? A control that applies
# to some cohorts and not others should say so on the control, rather than letting the pooled
# result carry the surprise — the engine discloses a mixed pool (step 85), but only after the
# user has already run it.
#
# Measured over the SELECTED ENDPOINT's pool, not the whole tissue, and that distinction is
# the difference between a useful note and a misleading one. lgg is 5 of 6 tissue-wide, so a
# tissue-wide note fires on lgg/OS and says a sixth cohort is unstratified — but that cohort
# is GSE107850_LGG, which is DFS-only and is not in the OS pool at all. The panel already
# shows it struck through as "not offered for OS". The honest statement is 5 of 5 for lgg/OS
# (silent) and 1 of 2 for lgg/DFS (worth saying). Falls back to the tissue when no endpoint
# is settled yet, which only happens on the first render.
strata_coverage_note <- function(ct, ep = NULL) {
  co <- if (!is.null(ep) && nzchar(ep)) tryCatch(cohorts_for(ct, ep),
                                                 error = function(e) cohorts_of(ct))
        else cohorts_of(ct)
  # Per ENDPOINT (step 102): the coverage note describes the pool actually fitted.
  v <- vapply(co, stratifier_for, character(1),
              endpoint = if (!is.null(ep) && nzchar(ep)) ep else NULL)
  if (!length(v) || all(is.na(v)) || !anyNA(v)) return(NULL)
  sprintf(paste0("Only %d of the %d pooled %s cohorts carry it; the rest are fitted ",
                 "unstratified, and the pooled result discloses the mix."),
          sum(!is.na(v)), length(v), cancer_label_of(ct))
}

# The Multiple query tab spans every tissue at once, so its note has to name them all.
# DERIVED, never a hand-written list: this sentence read "Breast (PAM50) only" until
# 2026-08-27, when lgg gained a stratifier and it silently became false. A hardcoded tissue
# list goes stale at exactly the moment the note matters.
mq_strata_note <- function() {
  hits <- Filter(function(ct) length(strata_vars_of(ct)) > 0, CANCER_TYPES)
  if (!length(hits)) return("No tissue declares a stratifier; this control does nothing.")
  # "<VAR> in <Tissue>", not "<Tissue> (<VAR>)": several cancer_label_of() values already
  # END in a parenthetical (e.g. "Lower-grade glioma (LGG)"), so the obvious phrasing renders
  # as "Lower-grade glioma (LGG) (IDH)".
  sprintf("%s; ignored by tissues with no stratifier.",
          paste(vapply(hits, function(ct)
                  sprintf("%s in %s", paste(toupper(strata_vars_of(ct)), collapse = "/"),
                          cancer_label_of(ct)),
                character(1)), collapse = ", "))
}

# ---- scan tables, discovered per tissue from results/<cancer_type>/ ---------------

# Human labels for the scan artifacts we know about. A CSV present in the tissue's
# results dir but absent here still shows up, labelled by filename — so a newly run
# scan appears in the app without a code change.
SCAN_LABELS <- c(
  mr_discval_os           = "OS - discovery->validation (full)",
  mr_discval_dfs          = "DFS - discovery->validation (full)",
  mr_discval_dfs_inverted = "DFS - discovery->validation, inverted split (sensitivity)",
  mr_discval_os_inverted  = "OS - discovery->validation, inverted split (sensitivity)",
  mr_scan_os              = "OS - marginal scan (all TFs)",
  mr_scan_dfs             = "DFS - marginal scan (all TFs)",
  mr_discval_os_idhadj    = "OS - discovery->validation, IDH-adjusted (lgg)",
  mr_scan_os_idhadj       = "OS - IDH-adjusted scan (all TFs)",
  mr_scan_os_pam50adj     = "OS - PAM50-adjusted scan (all TFs)",
  mr_scan_dfs_pam50adj    = "DFS - PAM50-adjusted scan (all TFs)",
  immune_mediation_summary = "Immune mediation test (Step D)",
  cptac_protein_validation = "CPTAC protein validation (Step E)",
  cna_mediation_summary   = "CNA mediation test (Step F)"
)
# Per-cohort detail tables and staging artifacts: present on disk, not useful as a
# standalone browser view (os_hits_annotated is joined into os_annotated instead).
SCAN_HIDE <- c("immune_mediation", "cna_mediation", "cohort_power_diagnostic",
               "phase2b_audit", "phase2b_duplicates", "phase2b_id_coverage",
               "os_hits_annotated")

.read_scan <- function(ct, name)
  read.csv(file.path(ROOT, "results", ct, paste0(name, ".csv")), stringsAsFactors = FALSE)

# Load every scan CSV for a tissue once at startup, keyed by basename.
.load_scans <- function(ct) {
  dir <- file.path(ROOT, "results", ct)
  if (!dir.exists(dir)) return(list())
  files <- setdiff(tools::file_path_sans_ext(basename(Sys.glob(file.path(dir, "*.csv")))),
                   SCAN_HIDE)
  files <- files[!grepl("^survtable_", files)]
  out <- setNames(lapply(files, function(f) tryCatch(.read_scan(ct, f), error = function(e) NULL)),
                  files)
  out <- out[!vapply(out, is.null, logical(1))]

  # Breast only: join the Step A/B annotation columns onto the OS discovery->validation
  # table, which is the curated headline view. Guarded on the annotation file existing,
  # so a tissue without Step A/B (ovarian) simply doesn't get this entry.
  ann_f <- file.path(dir, "os_hits_annotated.csv")
  if (file.exists(ann_f) && !is.null(out$mr_discval_os)) {
    ann <- read.csv(ann_f, stringsAsFactors = FALSE)
    keep <- intersect(c("tf", "chembl_target_id", "known_drugs", "druggable_bucket",
                        "top_pmids", "one_line_biology", "pam50_adj_HR", "pam50_adj_p",
                        "pam50_adj_q_genomewide", "pam50_adj_direction",
                        "pam50_subtype_independent"), names(ann))
    out$os_annotated <- merge(out$mr_discval_os, ann[keep], by = "tf", all.x = TRUE)
  }
  out
}
SCANS_BY_CT <- setNames(lapply(CANCER_TYPES, .load_scans), CANCER_TYPES)

scan_choices_of <- function(ct) {
  nm <- names(SCANS_BY_CT[[ct]])
  if (!length(nm)) return(character(0))
  lab <- ifelse(nm == "os_annotated",
                "OS - discovery->validation, Step A/B annotated",
                ifelse(nm %in% names(SCAN_LABELS), SCAN_LABELS[nm], nm))
  ord <- order(match(nm, c("os_annotated", names(SCAN_LABELS))), nm)
  setNames(nm[ord], lab[ord])
}

# ---- feature vocabularies, per tissue --------------------------------------------

# Union across the tissue's cohorts, per kind. File prefix == kind name exactly (see
# data_access.R's get_feature()); cohorts missing a given kind's h5 (e.g. the ovarian
# panel has no "cna") are silently skipped.
.feature_choices <- function(kind, cohorts) {
  feats <- character(0)
  for (coh in tolower(cohorts)) {
    f <- file.path(PROC, sprintf("%s_%s.h5", kind, coh))
    if (file.exists(f)) feats <- union(feats, as.character(rhdf5::h5read(f, "/features")))
  }
  sort(feats)
}
# Show the "Browse hits" (discovery) tab. Off while the focus is the plot engine.
SHOW_BROWSE <- FALSE

KINDS <- c("viper", "expr", "immune", "cna")
FEATURES_BY_CT <- setNames(lapply(CANCER_TYPES, function(ct)
  setNames(lapply(KINDS, .feature_choices, cohorts = cohorts_for(ct)), KINDS)), CANCER_TYPES)
KIND_LABEL <- c(viper = "VIPER", expr = "expr", immune = "immune", cna = "CNA")

# Kinds offered as a top-level query. Two of the four KINDS are deliberately NOT offered,
# for the same reason in two shapes: a kind belongs here only if the tool can say something
# about the answer, and neither of these has a scan, an anchor or a validation split behind
# it.
#
# "cna" was built as the COVARIATE for the Step F mediation test, not as a question in its
# own right. Standing alone it covers only TCGA_BRCA + METABRIC (no SCAN-B, no ovarian) and
# largely restates the cna_only HR the Step F badge already reports. The h5s stay --
# R/step_f_cna.R reads them, and get_feature(kind="cna") still works from the console.
#
# "immune" was withdrawn from the selectors on 2026-08-30 (step 108). It is 77 xCell/MCP/
# CIBERSORT cell-type scores on breast and 0 features on the other six tissues, and nothing
# ranks them: every scan in results/ ranks VIPER activity, so a queried immune score reaches
# the forest with no rank, no q, no discovery/validation split and no anchor. It was the one
# query in the app whose result could not be placed against anything the tool had computed.
# HIDDEN, not deleted, and deliberately so -- the immune_*.h5 files stay, get_feature(kind =
# "immune") still works from the console, tumour/normal still accepts it (TN_KINDS below),
# and the Step D immune-mediation badge is untouched and still renders. That badge is not
# this: it reports a pre-registered lineage score as a COVARIATE on a VIPER hit's HR, which
# is a statement about a result, where the selector offered a browse of scores with no
# result attached. Restoring the browse is this one vector; what it would need first is a
# reason to read a score, which the badge has and the selector did not.
QUERY_KINDS <- c("viper", "expr")

# ...and of those, only the ones that actually have data for this tissue: offering a score
# type for a tissue whose h5 for it was never staged would produce an empty selector and a
# confusing empty result. Both surviving kinds are present in all seven tissues, so this
# gate returns the full set today -- it stays because that is a fact about which files are
# staged, not a design guarantee, and the next tissue staged without one is what it is for.
kinds_of <- function(ct) {
  f <- FEATURES_BY_CT[[ct]]
  QUERY_KINDS[vapply(QUERY_KINDS, function(k) length(f[[k]]) > 0, logical(1))]
}
# All three Score-type radios take their labels through kind_choice_html() (2026-08-18,
# step 48), which attaches the VIPER hover explanation. All three deliberately, not just
# Single query: it is the same control with the same misleading label, and explaining it on
# one tab only is worse than not explaining it -- a reader who saw it once will read its
# absence on the next tab as "this one is different".
KIND_CHOICE_LABEL <- c(viper = "TF activity inferred from VIPER", expr = "mRNA expression",
                       immune = "Immune cell score", cna = "Copy number (CNA)")

# ---- Multiple query: the cross-tissue vocabulary ----------------------------------

# The tissues do NOT share a feature vocabulary, and they share less of one every time a
# tissue is added. Re-derived over all seven, 2026-08-26 (step 81): VIPER is 1598 symbols in
# the union but only 1189 present in EVERY tissue; expr is 58366 vs 20934. Both intersections
# fell as brain came in (1253 -> 1189, 21281 -> 20934) while both unions rose -- which is the
# arithmetic of the decision below, not a regression. The selector offers the UNION and a
# tissue that lacks the feature renders a panel saying so, rather than being dropped: a
# missing panel is indistinguishable from a null result, and silently showing six panels
# when seven were promised is exactly the quiet-default failure this project bans. Offering
# the INTERSECTION instead would have cost the selector 409 VIPER symbols outright.
FEATURES_ALL <- setNames(lapply(KINDS, function(k)
  sort(Reduce(union, lapply(CANCER_TYPES, function(ct) FEATURES_BY_CT[[ct]][[k]])))), KINDS)

# Score types offered on the Multiple query tab: any QUERY_KIND with data in at least one
# tissue. The threshold is "at least one" and not "all" for the reason in the paragraph
# above -- a tissue lacking the feature renders a panel saying so rather than being dropped.
# Until step 108 "immune" was the case that exercised this: it qualified on breast alone
# (77 features, 0 elsewhere), so choosing it yielded one estimable panel and six stated
# "not measured here" ones. With immune withdrawn, both survivors are present in all seven
# tissues, so this filter is currently a no-op over kinds and the per-tissue absence it
# handles is now reached per FEATURE rather than per kind.
MQ_KINDS <- QUERY_KINDS[vapply(QUERY_KINDS, function(k) length(FEATURES_ALL[[k]]) > 0, logical(1))]

# ---- CPTAC-BRCA proteomics (breast-only detectability badge) -----------------------

# CPTAC-BRCA proteomics (see scripts/stage_cptac.py): 122 tumors, NOT usable for
# survival (2 OS events, no DFS/DSS). Only used here as a protein-detectability
# check - was this gene's product actually quantified by mass spec, and in what
# fraction of tumors.
#
# The read, the isoform collapse and the reasons for both live in R/cptac.R, which
# app.R sources above. It used to be built inline here as
# setNames(rowSums(!is.na(vals)), genes), which indexed by a name that is NOT unique
# (707 symbols own several Ensembl-protein rows) and so reported the FIRST isoform's
# coverage as the gene's -- MKI67 read 7/122 when it is detected in all 122. Moved out
# because app.R has no test harness and the invariant that makes the lookup safe (unique
# names) has to be asserted somewhere a test can reach.
.cptac <- cptac_coverage()
CPTAC_COVERAGE <- .cptac$coverage
CPTAC_N <- .cptac$n

# ---- DepMap common essentiality (all five tissues) ---------------------------------

# Unlike the CPTAC badge above, this is not tissue-specific: it is a property of the GENE
# measured in cancer cell lines, so it applies wherever the gene can be queried. Read once
# at startup -- two small text files, 280 KB on disk, nothing per-request. See R/depmap.R
# for the three states, the kind gate and why the screened set is staged alongside the
# essentials list.
DEPMAP_LISTS <- depmap_lists()

DEFAULT_CT <- if ("breast" %in% CANCER_TYPES) "breast" else CANCER_TYPES[1]

# Everything the Guide tab states, derived ONCE at startup from the registry and the scans
# already in memory. No extra I/O, and no figure on that page is typed -- see R/guide.R's
# header for why a hand-written landing page would be the one document here with no
# staleness guard on it.
GUIDE_FACTS <- guide_facts(SCANS_BY_CT, FEATURES_BY_CT)
ABOUT_FACTS <- about_facts()

# Content-addressed, so the browser re-fetches the icon exactly when the icon changes. tools
# ships with R; no dependency added.
FAVICON_HREF <- paste0("favicon.png?v=",
                       substr(unname(tools::md5sum("www/favicon.png")), 1, 8))

# Tab titles are constants, not literals: navbarPage identifies a tab BY ITS TITLE, so the
# "Query selected TF ->" handler's updateNavbarPage() target has to be the same string as
# the tabPanel's. Two literals is how that silently stops navigating when one is renamed.
# NOT a nav item since step 111 (2026-08-30). The Guide is still navbarPage's FIRST
# tabPanel -- which is the only thing that makes it the landing page -- but its <li> is
# hidden in PAGE_CSS and the brand lockup carries the trip back. The string is still needed
# here and needed as a constant: the CSS selector that hides the item is keyed to it, and so
# is the brand handler's updateNavbarPage() target.
TAB_GUIDE  <- "Guide"
TAB_SINGLE <- "Single query"
TAB_MULTI  <- "Multiple query"
TAB_MGENE  <- "Multiple genes"
# Its own tab, deliberately. RPPA rows are antibodies on ONE cohort, not cohorts, so they
# must never reach the forest -- for ovarian an RPPA row would enter a k=13 pooled diamond
# built from TCGA_OV patients it shares. Separation is the guarantee; see R/rppa.R.
TAB_RPPA   <- "Protein (RPPA)"
# LAST, deliberately: provenance and people are what a reader comes back for, not what
# they open the tool to do. It is the Guide's counterpart on a different clock -- the
# Guide changes when the analysis changes, this changes when the project's circumstances do.
TAB_ABOUT  <- "About"
# Tissues with an rppa_<cohort>.h5 on disk. Resolved ONCE at startup from the files
# themselves, never listed: a hand-kept list is what goes stale when a tissue is staged.
RPPA_TISSUES <- Filter(rppa_available, CANCER_TYPES)
# Endpoints come from the curated CDR, which carries all three for every TCGA cohort.
RPPA_MAX_ANTIBODIES <- 8L   # dynamic KM outputs are pre-registered up to this many

# BRAND_PURPLE / BRAND_ORANGE / brand_tint() live in R/brand.R, sourced above. See that file
# for where the two hexes come from and for the one part of this palette that is deliberately
# NOT branded (the three .note-* severity colours below).
#
# WHICH COLOUR CARRIES THE BAR, and why it is not the louder half: the navbar keeps a light
# ground with a purple rule and purple type, rather than a saturated fill. The mark's white
# was made transparent (step 96) so the bar shows through it, and the mark is ~45% purple and
# ~24% orange -- a purple bar would swallow the helix, an orange bar would swallow the bars
# and the KM curve. A tinted ground is what lets the whole mark stay legible on it. Orange is
# spent where it is scarce and therefore useful: the ACTIVE tab.

# The left/right page inset + centred max-width, restored 2026-08-06. It was added
# 2026-07-28 as a global `.container-fluid` rule and was lost when app.R was later rolled
# back wholesale; the two-tab work did not remove it. Below 1400px the cap is inert and the
# rule is just the 5% inset; above it the content stops widening and centres, so the
# interpretation paragraphs do not sprawl to an unreadable line length on a wide monitor.
#
# Scoped to a CLASS rather than to `.container-fluid` globally. It was applied to the
# Single query tab alone while the Multiple query tab laid its panels two-across and wanted
# every pixel of width. Both halves of that have since gone: step 53 replaced the 2-up grid
# with one full-width card per tissue and caps the FIGURE at MQ_FOREST_W_PX instead of
# stretching it, so width is no longer what constrains a forest; and step 80 put all four
# tabs inside .narrow-page, because control panels that differ in width by 16px across tabs
# read as a rendering bug. The class stays scoped rather than global so the inset is a
# decision each tab makes, not something a Bootstrap container imposes.
# The one number PAGE_CSS cannot hold as a literal: .mq-forest's cap is MQ_FOREST_W_PX
# (R/plots.R), the same constant mq_forest_h_px() derives the height from. A second literal
# here is a width that can drift from the height computed against it, which is precisely
# how the RPPA panel came to reserve margin for a string it did not print.
PAGE_CSS <- "
/* Stacked RPPA KM curves. The heading carries the antibody because, unlike the Single
   query tab's KM tabs, there is no tab label to carry it -- and the curve's own title
   is inside the image, so it does not survive a screenshot crop. */
.rppa-km { margin-top: 18px; }
.rppa-km-head { margin: 0 0 2px 0; font-size: 15px; font-weight: 600; color: #333; }

.narrow-page { max-width: 1400px; margin: 0 auto; padding-left: 5%; padding-right: 5%; }

/* The four tabs' control panels must be the SAME WIDTH. They were not: measured 2026-08-25,
   Single query's .well came out 1111.5px and the other three 1095px. Cause is structural,
   not stylistic -- Single query is div(narrow-page, wellPanel(...)) while the other three
   are div(narrow-page, fluidPage(wellPanel(...))), and fluidPage() emits a .container-fluid
   whose 15px side padding narrows its contents by 30px. .narrow-page already owns the
   centring and the 5% inset, so that padding is pure duplication. Neutralised here rather
   than by removing the fluidPage() wrappers, which would be a layout change to three tabs
   to fix a margin on one. Asserted in tests/test_mq_layout.R. */
.narrow-page > .container-fluid { padding-left: 0; padding-right: 0; }

/* ...and once that was fixed a SECOND, smaller mismatch showed up underneath it: 1111.5px
   on Single query against 1125px on the other three, in the opposite direction. That one is
   not layout at all, it is the VERTICAL SCROLLBAR. Single query's page is taller than the
   viewport at rest (the badge lines), so it scrolls and loses ~15px of width; the other
   three at load are short enough not to. Reserving the gutter permanently makes the width
   independent of how much content a tab happens to have, which is the only way the four can
   be equal in every state rather than equal only after a query has been run on each. */
html { overflow-y: scroll; }

/* Cancer-type labels are long -- 'Lung squamous cell carcinoma (LUSC)' is 36 characters --
   and the Single query control is only ~154px wide (column(2) of a 12-column row), so the
   OPEN LIST was clipping every label except Breast and Ovarian. The Multiple genes tab did
   not show this only because its control sits in a column(3). Widening the column would
   have to take a unit from a neighbour; instead the dropdown is allowed to size to its
   content while the closed control keeps its column width. Scoped to the two cancer-type
   selectors, so no other dropdown's width behaviour changes. */
#cancer_type + .selectize-control .selectize-dropdown,
#mg_ct + .selectize-control .selectize-dropdown,
#rp_ct + .selectize-control .selectize-dropdown {
  width: auto !important;
  min-width: 100%;
  white-space: nowrap;
}
#cancer_type + .selectize-control .selectize-dropdown .option,
#mg_ct + .selectize-control .selectize-dropdown .option,
#rp_ct + .selectize-control .selectize-dropdown .option { white-space: nowrap; }

/* Brand lockup: the mark sits left of the wordmark, both vertically centred in the bar.
   navbar-brand ships a fixed line-height/padding for a text-only brand, so the image needs
   its own height rather than inheriting one. */
.navbar-brand { display: flex; align-items: center; }
.brand-wrap { display: inline-flex; align-items: center; gap: 9px; cursor: pointer; }
.brand-logo { height: 36px; width: auto; display: block; }
.brand-name { font-weight: 700; letter-spacing: 0.2px; font-size: 21px;
              color: __BRAND_PURPLE__; }
/* The brand IS the way back to the Guide (step 111): the Guide has no nav item, so the bar
   shows the four query tabs and About, and the lockup carries the return trip.
   Keyed to the tab's data-value, NOT to a position -- :nth-child would silently re-point at
   Single query the first time a tab is added -- and the value is substituted from TAB_GUIDE
   below, so the tab title is still stated exactly once in this file.
   Hiding the <a> is what collapses the item: in this navbar the padding lives on the anchor,
   not on the <li>. The :has() rule removes the empty <li> too where the browser supports it;
   either rule alone is sufficient, which is why both are here.
   SINGLE quotes around the value: PAGE_CSS is a double-quoted R string, and a double quote
   here does not escape -- it ends the string, and app.R stops parsing mid-selector. */
.navbar-nav > li:has(> a[data-value='__TAB_GUIDE__']) { display: none; }
.navbar-nav > li > a[data-value='__TAB_GUIDE__'] { display: none; }
/* actionLink brings its own link colour and underline; the brand has to keep the navbar's,
   in all three states, or the lockup changes colour on hover for no reason a reader can use. */
.navbar-brand a.brand-link, .navbar-brand a.brand-link:hover,
.navbar-brand a.brand-link:focus { color: __BRAND_PURPLE__; text-decoration: none; outline: none; }

/* The bar itself (step 112). Taller and in the brand's colours, which it was not: simplex
   ships a near-white bar with a red accent that belongs to the theme rather than to this
   tool, so the navbar was the one part of the page that looked like a default.
   Selectors are .navbar-default-prefixed to match simplex's OWN specificity -- overriding a
   theme with bare .navbar rules loses on specificity and needs !important to win, which then
   cannot itself be overridden. Matching the theme's shape is the way to stay overridable. */
.navbar-default { background-color: __BRAND_TINT__; border: 0;
                  border-bottom: 3px solid __BRAND_PURPLE__; min-height: 62px; }
.navbar-default .navbar-nav > li > a {
  color: __BRAND_PURPLE__; font-size: 15.5px; font-weight: 600; padding: 21px 17px;
  border-bottom: 3px solid transparent; margin-bottom: -3px; }
.navbar-default .navbar-nav > li > a:hover,
.navbar-default .navbar-nav > li > a:focus {
  color: __BRAND_PURPLE__; background-color: __BRAND_HOVER__;
  border-bottom-color: __BRAND_ORANGE__; }
/* The ACTIVE tab is the one place the orange is spent, and it is spent on the one thing a
   reader needs to read off the bar without thinking: where they are. */
.navbar-default .navbar-nav > .active > a,
.navbar-default .navbar-nav > .active > a:hover,
.navbar-default .navbar-nav > .active > a:focus {
  color: #fff; background-color: __BRAND_PURPLE__; border-bottom-color: __BRAND_ORANGE__; }
.navbar-brand { height: auto; padding-top: 10px; padding-bottom: 10px; }

/* Everything else the theme painted red (step 113). simplex ships #d9230f buttons and
   bootstrap ships #c7254e code chips; against a derived orange those read as a THIRD red
   that nearly matches and does not -- which looks like a mistake rather than a choice.
   Buttons take the purple, not the orange, and that is a contrast decision rather than a
   taste one: white on __BRAND_ORANGE__ is about 2.9:1, which fails at any text size, while
   white on __BRAND_PURPLE__ is about 12:1. Dark text on the orange would pass at ~7.1:1, but
   dark-on-orange is what a CAUTION control looks like, and .note-warn is amber a few inches
   below it.
   The three severity colours further down are NOT included here -- green, amber and red mean
   something to a reader before any word is read; see R/brand.R. (Written without the class
   prefix on purpose: tests/test_note_classes.R reads the DEFINED class set straight out of
   this stylesheet, and a class named in a comment would enter it as a phantom rule.)
   COMPOUND selectors (.btn.btn-primary, 0-2-0) rather than bare ones (.btn-primary, 0-1-0),
   and that is the fix for a real bug, not tidiness. simplex states these at 0-1-0 too, so a
   bare rule here wins only by coming LATER in the document -- and the order in which Shiny
   injects `theme` versus `header` is a property of the Shiny/htmltools build, not something
   this file controls. It happened to come out our way in one environment and red in another;
   at 0-2-0 the question stops being asked. Same reason `body code` carries an element more
   than bootstrap's bare `code`.
   AND `background-image: none` ON EVERY ONE OF THEM, which is the whole reason these buttons
   stayed red through two rounds of being called fixed. simplex does not paint its buttons
   with a colour;
   it paints them with `background-image: linear-gradient(#e72510, #d9230f 6%, #cb210e)`. A
   `background-color` override underneath a gradient changes a value nobody can see, so
   getComputedStyle reported our purple while every browser rendered the theme's red -- which
   is exactly how this was verified as working twice and was not: reading a computed property
   is not the same as looking at the pixels. A colour override on a
   themed control has to clear the image or it overrides nothing.
   Flat rather than a purple gradient of our own: the navbar and the active tab are flat, and
   a button with depth beside a bar without it is the odd one out. */
.btn.btn-default, .btn.btn-primary, .btn.btn-default:focus, .btn.btn-primary:focus {
  background-color: __BRAND_PURPLE__; background-image: none;
  border-color: __BRAND_PURPLE__; color: #fff; }
.btn.btn-default:hover, .btn.btn-primary:hover, .btn.btn-default:active,
.btn.btn-primary:active, .btn.btn-default.active,
.open > .dropdown-toggle.btn.btn-default {
  background-color: __BRAND_DARK__; background-image: none;
  border-color: __BRAND_DARK__; color: #fff; }
body code { color: __BRAND_PURPLE__; background-color: __BRAND_PALE__; }
a, a:focus { color: __BRAND_PURPLE__; }
a:hover { color: __BRAND_DARK__; }
/* Multiple query: one full-width card per tissue (step 53). The forest is capped inside
   it -- the card is what makes the column symmetric, the cap is what keeps the figure in
   the proportions it was designed at.
   The auto margins are not decoration (step 54). A capped block with no auto margins is
   flush LEFT, so the 900px forest sat against the card's left edge with 148px of dead
   space on the right while the KM curves and the table -- uncapped, so full width --
   filled it. Same card, two left edges, and the figure was the one that looked wrong.
   CENTRED rather than stretched to the card, which was the other way to square the two
   edges: the width is what mq_forest_h_px() derives the height FROM, so a fluid width
   against a height computed in R makes the aspect a function of the browser window. That
   is exactly the defect HANDOFF records against Single query's forest_slot, and it is not
   worth reintroducing here to save a margin. Cap the width, centre what is capped. */
.mq-card { border: 1px solid #ddd; border-radius: 4px; padding: 14px 16px; margin-bottom: 22px; }
.mq-card-head { margin-top: 0; margin-bottom: 2px; }
.mq-forest { max-width: __MQ_FOREST_W__px; margin-left: auto; margin-right: auto; }
.mq-more { display: inline-block; margin-top: 8px; font-size: 0.92em; }
.mq-detail { margin-top: 10px; border-top: 1px solid #eee; padding-top: 12px; }
.badge-block { max-width: 760px; margin-bottom: 6px; }
/* Export buttons are DEMOTED to an outline, and the selector has to be strong enough to
   stay that way. This is a 0-2-1 compound (a + .btn + .btn-export) and that number is
   load-bearing: it was `a.btn-export` (0-1-1) from 2026-08-06, which beat simplex's own
   `.btn-default` (0-1-0) with room to spare -- and then step 115 raised the brand button
   rule to `.btn.btn-default` (0-2-0) to fix the red-button bug and silently outranked it.
   Every download button in the app came out solid brand purple, the same weight as Run,
   which is the inverted hierarchy this class exists to prevent. Nothing failed; the page
   just quietly stopped saying which control was the primary one. Raising this to 0-2-1
   restores the gap and keeps it ahead of the brand rule wherever that rule goes next.
   `background` is the SHORTHAND on purpose -- it resets background-image, so the theme's
   gradient cannot show through underneath, which is the trap step 115 was about.
   The margins are what separate the buttons: the bar is right-aligned and wraps to two
   rows, so both a left and a bottom gap are needed or the buttons touch on both axes.
   text-decoration is cleared on every state because downloadButton renders an <a>, and
   Bootstrap underlines links on hover and focus -- these are buttons, not links. */
/* The one-line affordance that carries a figure's hover explanation. Muted and small on
   purpose: it is a way in for a reader who needs it, not a caption competing with the
   figure. `note-tuck` pulls a note up under the thing above it; this does the same, and is
   its own class rather than a reuse so that restyling notes cannot silently restyle help. */
.fig-help { margin: 2px 0 12px 0; font-size: 0.88em; color: #666; }

a.btn.btn-export, a.btn.btn-export:visited {
  background: #fff; color: #333; border: 1px solid #bbb; font-weight: normal;
  text-decoration: none; margin: 0 0 6px 6px; }
a.btn.btn-export:hover, a.btn.btn-export:focus, a.btn.btn-export:active {
  background: #f2f2f2; color: #111; border-color: #999; text-decoration: none; }

/* Cohort chips (2026-08-07). A checkboxGroupInput with the native box hidden and its label
   drawn as a pill: the INPUT BINDING stays Shiny's own, so input$cohorts is still a plain
   character vector and no call site changed -- only the paint is different.
   Selected/deselected are encoded by fill and border weight, not by hue: the palette here is
   red for the primary action and amber for warnings, and a fourth colour on 13 pills would
   compete with both. The id selectors (1,2,0) outrank Bootstrap's .checkbox-inline rules
   (0,2,0), so none of this needs !important. */
#cohorts .shiny-options-group { display: flex; flex-wrap: wrap; gap: 6px; }
#cohorts .checkbox-inline { margin: 0; padding: 0; }
#cohorts .checkbox-inline + .checkbox-inline { margin-left: 0; }
#cohorts .checkbox-inline input { position: absolute; opacity: 0; width: 0; height: 0; }
#cohorts .checkbox-inline span, .chip-dead {
  display: inline-block; padding: 3px 10px; border-radius: 11px; font-size: 12px;
  line-height: 1.5; border: 1px solid #ccc; background: #fff; color: #888; cursor: pointer; }
#cohorts .checkbox-inline span:hover { border-color: #999; color: #555; }
#cohorts .checkbox-inline input:checked + span {
  background: #ececec; border-color: #666; color: #222; }
#cohorts .checkbox-inline input:checked + span::before { content: '\\2713\\00a0'; }
#cohorts .checkbox-inline input:focus + span { box-shadow: 0 0 0 2px rgba(217,35,15,.35); }
.chip-dead { border-style: dashed; color: #aaa; text-decoration: line-through;
  cursor: not-allowed; }
.chip-row { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 6px; }

/* Secondary notes (2026-08-08). Every caveat, badge and disclosure line under a control or
   a plot. These were 18 inline style= strings that had drifted to FOUR sizes (80/85/88/90%)
   and three negative margin-tops (-8/-6/-4px) for one visual job, so two notes side by side
   on the same page did not match. The size is now declared ONCE here.

   85% is the plurality of what was there (8 of 18), so this is the value that moves fewest
   lines rather than a new preference. It is deliberately a single number: a note is either
   secondary text or it is not, and the app has no third rank of emphasis to encode.

   TONE is semantic, not decorative, and the four colours are the ones already in use across
   both tabs -- muted / good / warn / err. They must keep matching between tabs, which is why
   they are named here instead of retyped: the scan-rank badge is painted from the SAME three
   (good/warn/muted) on Single query and on every Multiple query panel, and a rank that
   clears FDR must not read as green on one tab and orange on the other.

   The two spacing modifiers exist because a note's relationship to what sits above it is
   genuinely different in two places: `tuck` pulls it under the widget it annotates (Shiny
   leaves a form-group gap that would otherwise read as a separate block), `gap` pushes it
   off a chip row that has no bottom margin. `stack` is for consecutive notes in one panel.
   Anything needing a fifth spacing is a layout problem, not a note problem.

   tests/test_note_classes.R pins this set, and asserts no inline font-size survives in the
   rest of app.R -- which is the only thing that stops the drift from starting again. */
.note { font-size: 85%; }
.note-muted { color: #888; }
.note-body  { color: #333; }
.note-good  { color: #2a6b2a; }
.note-warn  { color: #a05a00; }
.note-err   { color: #b00; }
.note-tuck  { margin-top: -6px; }
.note-gap   { margin-top: 6px; }
.note-stack { margin-bottom: 2px; }

/* The stratify control, tucked under the Endpoint select (2026-08-15). Bootstrap gives the
   select's .form-group a 15px bottom margin and .checkbox a 10px top margin, so stacking the
   two leaves a 25px gap that reads as a break between unrelated controls rather than as one
   column. The tuck closes it; the checkbox still clears the select's own baseline. */
.strata-tuck { margin-top: -12px; }

/* Hover explanations on a note (2026-08-18). The marker sits inline after the note text;
   the panel is absolutely positioned so opening it moves nothing on the page.
   THE TRIGGER IS THE ONLY PC-ONLY PART. `:hover` never fires on a touch device and there
   is no fallback -- that is a deliberate, recorded choice (this tool is PC-only), and
   swapping it for click means replacing the two selectors below with a toggled class.
   `:focus-within` is not decoration: it is what lets the marker be reached by keyboard,
   which is why the wrapper carries tabindex. */
.infowrap { position: relative; display: inline-block; margin-left: 5px; vertical-align: baseline; }
/* The marker's geometry is MEASURED, not eyeballed (2026-08-18). The first version was a
   13px box with line-height 12px and text-align centre, holding an italic sans `i`: the
   glyph landed left of centre (italic shifts it) and low (the line box did not match the
   content box), and the whole mark sat 1.55px BELOW the text's optical centre -- measured
   against a Range over one character of the note. inline-flex centres the glyph with no
   magic numbers, and the -0.14em is that 1.55px expressed against the note's 11.05px font
   so it survives a change of note size; it brings the offset to 0.15px. Width and height
   must stay equal or the circle becomes an ellipse. */
.infomark { display: inline-flex; align-items: center; justify-content: center;
            box-sizing: border-box; width: 14px; height: 14px; line-height: 1;
            font-size: 10px; font-weight: 700; font-style: normal;
            font-family: Georgia, serif;   /* falls back to any serif; the classic info `i` */
            border: 1px solid #aaa; border-radius: 50%; color: #777; cursor: help;
            vertical-align: middle; position: relative; top: -0.14em;
            -webkit-user-select: none; user-select: none; }
.infowrap:hover .infomark, .infowrap:focus-within .infomark { border-color: #333; color: #222; }
.infopop { display: none; position: absolute; left: -6px; top: 20px; z-index: 30;
           width: 300px; padding: 8px 10px; text-align: left; white-space: normal;
           color: #333; background: #fff; border: 1px solid #ccc; border-radius: 3px;
           box-shadow: 0 2px 8px rgba(0,0,0,0.15); }
.infowrap:hover .infopop, .infowrap:focus-within .infopop { display: block; }
.strata-tuck .checkbox { margin-top: 0; margin-bottom: 0; }

/* The Multiple-genes drawer (2026-08-16, step 45). The row shows what is compared ACROSS
   genes; clicking the arrow opens everything that describes ONE gene beneath it. Markup is
   built by mgene_drawer_html() in R/plots.R -- only the styling is here.

   The whole CELL is the click target, not the glyph inside it: a 10px arrow is a hard thing
   to hit, and DataTables gives the cell the .mg-control class through columnDefs, which is
   also what the delegated handler binds to. cursor + user-select say so on hover, since
   nothing else on this table is clickable. */
td.mg-control { cursor: pointer; text-align: center; color: #666;
                -webkit-user-select: none; user-select: none; }
td.mg-control:hover { color: #111; }
.mg-chev { display: inline-block; transition: transform 120ms ease; }
/* Rotated rather than swapped for a second character: one glyph that turns reads as the SAME
   control in two states, where a change of symbol reads as two controls. tr.mg-open is set
   by the click handler in mgene_dt_callback(). */
tr.mg-open .mg-chev { transform: rotate(90deg); }

.mg-drawer { padding: 10px 14px 12px 34px; background: #fafafa;
             border-left: 3px solid #ccc; }
.mg-drawer-head { font-weight: bold; margin-bottom: 4px; }
.mg-kvs:empty { display: none; }
.mg-kv { display: inline-block; margin-right: 16px; font-size: 85%; color: #333; }
.mg-kv b { font-weight: 600; margin-right: 4px; }
.mg-drawer .note { margin: 3px 0; }

/* The per-cohort table is survtable_display()'s, the same one the Single query tab draws --
   but it is plain HTML here rather than a DataTable, so it gets its borders explicitly
   instead of inheriting DT's. width:auto keeps it as wide as its numbers rather than
   stretching one cohort's row across the page. */
.mg-cohorts { margin-top: 8px; width: auto; border-collapse: collapse; font-size: 85%; }
.mg-cohorts th, .mg-cohorts td { padding: 3px 12px 3px 0; text-align: left;
                                 border-bottom: 1px solid #e5e5e5; white-space: nowrap; }
.mg-cohorts th { font-weight: 600; color: #555; }
/* The pooled row is a summary of the rows above it, not another cohort. On the Single query
   tab orderFixed pins it to the bottom; this table has no ordering, so the weight is what
   carries the distinction. */
.mg-cohorts tr.mg-pooled td { font-weight: 600; border-top: 1px solid #999; }
.mg-cap { margin-top: 6px; font-size: 80%; color: #888; }

/* Forest + tumor/normal panel, ONE row above the KM curves (2026-08-13).
   Flex rather than Bootstrap columns because the panel is a FIXED pixel size -- 422px from
   tn_size_in(), the same rule the PDF export uses -- and a 12-col grid would give it a
   percentage. At this page width col-sm-4 is 419px, three pixels under, which is not a
   rounding annoyance but a horizontal scrollbar on the whole page. The right column is
   pinned to the panel's own width and the forest takes what is left, so the row's outer
   edges land exactly on the KM's and the table's below it.
   min-width:0 on the forest because a flex item's default min-width is auto, which refuses
   to shrink below its content and pushes the panel off the right edge instead.
   wrap so a narrow window stacks them rather than squeezing the forest to nothing. */
.fp-row    { display: flex; gap: 24px; align-items: flex-start; flex-wrap: wrap; }
.fp-forest { flex: 1 1 520px; min-width: 0; }

/* The forest's box on the Single query tab (step 56). Its width is the BROWSER's -- the
   full page width with no tumour/normal panel, whatever is left beside one -- and it used
   to be given forest_height_in(k) * 96, the PDF's pixel HEIGHT, inside a 100%-width
   plotOutput. That fixes the wrong dimension: the height came from the design and the
   width from the window, so the figure came out at the WINDOW's aspect. Measured on
   ADNP/OS/ovarian (k=13, designed 7.5 x 8.9in, i.e. 720x854 px): 1237x854 with no panel
   and 665x854 with one -- the same query, 1.9x apart, neither of them the design.
   What that looks like is not bigger type. metafor's horizontal layout is scale-invariant
   here -- the interval region sat at 37-84% of the width at EVERY width tried -- while the
   text is drawn at a fixed size, so stretching only grows the gutters between the cohort
   column, the intervals and the HR column. Hence the cap: past ~900px the figure stops
   reading as a figure and starts reading as three columns adrift on a page.
   So: cap the width (the SAME constant .mq-forest uses -- one on-screen forest width for
   the tool, not one per tab), centre what is capped, and let CSS derive the height from
   whatever width the box actually gets. The aspect is inline because it is a function of k.
   contain:size is load-bearing, not a hint. Without it an auto height with an aspect-ratio
   resolves to the taller of the aspect height and the CONTENT -- and the content is the
   PREVIOUS query's image, which Shiny leaves in place while the new one renders. Measured:
   k=13 followed by k=3 kept 1068px instead of 461px, so the aspect silently reverted for
   exactly the queries that changed k. With it the height is a function of the width alone,
   which is the invariant this rule exists to state. */
.sq-forest { max-width: __MQ_FOREST_W__px; margin-left: auto; margin-right: auto;
             contain: size; }

/* Long feature symbols (2026-08-10). A feature name here is DATA, not a label: the offered
   vocabulary is 58366 distinct strings read out of the h5 files, 5837 of them '///'-joined probe
   sets, the longest 1466 characters, and the long ones contain no whitespace at all -- so
   nothing in them is a line-break opportunity and every default that relies on one fails.

   FOUR surfaces overflow, and they are the only four: measured in the running app at a
   1230px viewport with a 253-character symbol selected, by walking every element in the DOM
   and keeping the ones whose ink left their box -- not by reading the source and guessing
   which those would be. The guess would have found the first three and missed the fourth,
   which only exists after a query returns.

     .selectize-input .item        1835px laid out inside a 236px control, which scrolled
                                   the whole PAGE sideways (document scrollWidth 2300)
     .selectize-dropdown .option   already clipped at the control width, but with no marker
     prose (a validate() message, a .note)   1841px inside a 1080px column
     the survtable DT's <caption>  1841px, and it dragged the TABLE out with it

   The two selectize surfaces are FORM CONTROLS and must stay one line high, or choosing a
   long symbol re-flows the control row and moves the Run button: they truncate, with an
   ellipsis, and the whole symbol stays readable through the title= tooltip that
   FEATURE_ITEM_RENDER attaches. Prose WRAPS instead -- a message whose job is to say which
   feature was not found is useless with the feature cut off, and a paragraph has no fixed
   height to protect.

   The two wrap rules are NOT the same declaration, and the difference is the whole reason
   there are two. `break-word` breaks an over-long token but does not change an element's
   min-content width, so it is safe to inherit page-wide -- DT decides its column widths from
   min-content, and `anywhere` there would let every column collapse to one character. A
   <caption> is the case where min-content is exactly what has to give: the caption's
   min-content width sets the TABLE's width, so `break-word` leaves it 1841px and drags the
   table off the page. Hence `anywhere`, scoped to the one element where that is the point.
   tests/test_note_classes.R pins all four rules. */
.container-fluid { overflow-wrap: break-word; }
table caption { overflow-wrap: anywhere; }
.selectize-input .item, .selectize-dropdown .option {
  max-width: 100%; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
"

# gsub, not sub: the cap is now stated TWICE -- .mq-forest and .sq-forest -- and sub()
# replaces the first match only, so the second rule would ship the literal placeholder
# and cap nothing. A silent one, too: an unparseable declaration is dropped, the figure
# just goes back to full width. The test asserts no "__MQ_FOREST_W__" survives.
PAGE_CSS <- gsub("__MQ_FOREST_W__", MQ_FOREST_W_PX, PAGE_CSS, fixed = TRUE)
# Same reason, same failure mode: an unmatched placeholder leaves an unparseable selector,
# the browser drops the declaration, and the Guide's nav item quietly comes BACK -- next to
# a brand that also opens it. gsub because the value is used in two rules. Asserted in
# tests/test_guide_derived.R, which also refuses to let the title be typed here as a literal.
PAGE_CSS <- gsub("__TAB_GUIDE__", TAB_GUIDE, PAGE_CSS, fixed = TRUE)
# Same pattern again, four more placeholders. The two hexes come from BRAND_*; the two tints
# are MIXED from the purple rather than typed, so there is exactly one purple in this file and
# everything paler is a function of it.
for (.k in list(c("__BRAND_PURPLE__", BRAND_PURPLE), c("__BRAND_ORANGE__", BRAND_ORANGE),
                c("__BRAND_TINT__",  brand_tint(BRAND_PURPLE, 0.955)),
                c("__BRAND_HOVER__", brand_tint(BRAND_PURPLE, 0.88)),
                c("__BRAND_PALE__",  brand_tint(BRAND_PURPLE, 0.93)),
                c("__BRAND_DARK__",  brand_tint(BRAND_PURPLE, -0.28))))
  PAGE_CSS <- gsub(.k[[1]], .k[[2]], PAGE_CSS, fixed = TRUE)
rm(.k)
if (grepl("__BRAND_", PAGE_CSS, fixed = TRUE))
  stop("PAGE_CSS still carries a __BRAND_*__ placeholder after substitution -- the navbar ",
       "would ship an unparseable colour and fall back to the theme's default.")

# The chip above is ellipsised, so the whole symbol has to stay recoverable from the control
# itself -- the same rule the figure titles follow, where a cut is always marked and never
# silent. escape() is selectize's own HTML escaper: these strings come from the h5
# vocabularies rather than from a literal list, and "no symbol contains '<' today" is a fact
# about today's data, not a licence to interpolate raw. Defined ONCE and handed to both
# selectizeInputs, so the two chips cannot drift apart.
FEATURE_ITEM_RENDER <- I(paste0(
  "{item: function(item, escape) {",
  "return '<div class=\"item\" title=\"' + escape(item.label) + '\">'",
  " + escape(item.label) + '</div>';}}"))

# ---- UI ---------------------------------------------------------------------------

ui <- navbarPage(
  # Brand = mark + wordmark. The image is the GRAPHIC HALF of logo.png only (bars + KM
  # curve + helix); the source file also contains the words "OMICohort / KORKMAZ LAB"
  # underneath, and at navbar height that text renders as an unreadable smudge sitting
  # next to the real, legible wordmark beside it. The split point was derived from the
  # image rather than eyeballed -- the largest vertical run of blank rows (657 -> 724) is
  # the gap between the two halves. White is made transparent so the navbar colour shows
  # through. See scripts/make_logo.py; www/ is served by Shiny automatically.
  #
  # CLICKABLE since step 111, and load-bearing: the Guide has no nav item, so this is the
  # only route to it after the first navigation away. An actionLink rather than an onclick
  # that fakes a click on the hidden anchor -- the tab is switched by the same
  # updateNavbarPage() the "Query selected TF ->" button already uses, so there is one
  # navigation mechanism in this file instead of two, and it is one the tests can see.
  title = actionLink(
    "brand_home", class = "brand-link", title = "OMICohort \u2014 open the guide",
    label = tags$span(
      class = "brand-wrap",
      tags$img(src = "logo.png", class = "brand-logo", alt = "OMICohort"),
      tags$span(class = "brand-name", "OMICohort"))),
  id = "nav", windowTitle = "OMICohort",
  theme = if (HAVE_SHINYTHEMES) shinythemes::shinytheme("simplex") else NULL,
  # NOTE the ?v= on FAVICON_HREF. Browsers cache a favicon far more aggressively than any
  # other asset and will keep showing a stale one across restarts and hard reloads -- step 112
  # replaced the image and the old one went on rendering. The tag is the file's own md5, so it
  # changes exactly when the icon changes and never otherwise.
  #
  # www/favicon.png has existed since step 96 and was never linked, so nothing ever asked
  # for it: a browser requests /favicon.ico by default and this is a .png, which is exactly
  # the kind of asset that ships unused because no test can miss it. scripts/make_logo.py
  # writes it (the full square, wordmark included, at 64px); tests/test_about_derived.R now
  # holds the three together -- the file, the generator, and this link.
  header = tags$head(
    tags$link(rel = "icon", type = "image/png", href = FAVICON_HREF),
    tags$style(HTML(PAGE_CSS)), tags$style(HTML(GUIDE_CSS))),

  # Controls across the TOP, results beneath at full page width — the same shape as the
  # Multiple query tab (2026-08-06). The old sidebarLayout ran ~800px of controls beside a
  # ~570px forest, so it left a dead column, and it pinned the forest to the width-8 main
  # panel while the KM curves below it already had the whole page. Both plots now get the
  # full width; a forest's useful axis is the CI, a KM's is time, and both are horizontal.
  #
  # NOT a literal copy of the multi tab's bar: this tab has three controls that one does not
  # (cancer type, per-cohort selection, follow-up horizon) plus up to four Step D/E/F badges.
  # Those two are what made the sidebar tall, so they are handled rather than transplanted --
  # see the cohort selector's note in cohorts_ui, and the badge block below.
  # FIRST, so the app opens on it: navbarPage selects its first tabPanel by default, and a
  # reader who has never seen a forest plot should meet the explanation before the controls.
  # It is a static page -- no inputs, no reactives, nothing to recompute per session.
  # Its nav item is hidden (see PAGE_CSS): the panel has to EXIST and has to be first, but a
  # sixth item in the bar competed with the four tabs that actually do something. Do not
  # move it -- position is the landing rule, and the brand is the only other way in.
  tabPanel(TAB_GUIDE, div(class = "narrow-page", HTML(guide_html(GUIDE_FACTS)))),

  tabPanel(TAB_SINGLE,
   div(class = "narrow-page",
    wellPanel(
      fluidRow(
        # Display labels come from the registry (cancer_label_of); the VALUES stay the
        # lowercase slugs, so every downstream lookup is unchanged.
        column(2, selectInput("cancer_type", "Cancer type",
                              choices = setNames(CANCER_TYPES, cancer_label_of(CANCER_TYPES)),
                              selected = DEFAULT_CT)),
        column(2, uiOutput("kind_ui")),
        column(3, selectizeInput("feature", "Gene / TF symbol", choices = NULL, width = "100%",
                                 options = list(placeholder = "start typing a symbol...",
                                                render = FEATURE_ITEM_RENDER))),
        # The stratify control sits UNDER the endpoint (moved 2026-08-15 from its own
        # column in the row below). It is a model option -- it changes the Cox model, the
        # same way the endpoint and the horizon do -- so this is where its neighbours are.
        # Where it used to live, columns 7-8 of a 12-column row, it had a control to its
        # left and the export bar to its right and touched neither: a checkbox alone in a
        # gap reads as belonging to nothing, which is exactly how it read.
        column(2, uiOutput("endpoint_ui"),
                  div(class = "strata-tuck", uiOutput("strata_ui"))),
        # Fixed presets, NOT a free slider. The HR is genuinely sensitive to this
        # (ESR1/OS METABRIC: 0.95 at full follow-up, 0.59 at 60 months), so a
        # continuous control would let a user slide until the p-value looked good and
        # report that. A short list of pre-specified horizons keeps the choice
        # declarable rather than shoppable.
        # Choices come from tau_choices() in plots.R, not from a literal here: a tissue
        # with a declared panel horizon must DEFAULT to it FOR THE ENDPOINT being queried,
        # or the app contradicts the artifacts it serves (see tests/test_tau_choices.R).
        # The initial widget is built for DEFAULT_CT at its primary endpoint; the observer
        # below re-resolves it whenever the tissue OR endpoint changes.
        column(2, selectInput("max_fu", "Follow-up horizon",
                    choices  = tau_choices(DEFAULT_CT, default_endpoint_of(DEFAULT_CT))$choices,
                    selected = tau_choices(DEFAULT_CT, default_endpoint_of(DEFAULT_CT))$selected)),
        column(1, actionButton("run", "Run", class = "btn-primary",
                               style = "margin-top: 25px; width: 100%;"))
      ),
      fluidRow(
        column(8, uiOutput("cohorts_ui")),
        # class = "btn-export" demotes these to an outline (restored 2026-08-06 from the lost
        # 07-28 pass). Without it a downloadButton renders in the same fill as the actual
        # primary action -- the "Run" button -- inverting the hierarchy. downloadButton
        # renders as an <a>, so a compound selector beats the button fill on specificity
        # alone; no !important needed. THE NUMBER MOVED ONCE AND MUST NOT DRIFT AGAIN: this
        # said (0,1,1) vs the theme's (0,1,0) until step 115 raised the brand button rule to
        # (0,2,0) and took the demotion out without anything failing. See the rule itself in
        # PAGE_CSS, and tests/test_export_bar.R, which now asserts the outline rather than
        # trusting a comment about specificity to stay true.
        # ORDER MATCHES THE PAGE: forest, tumor/normal, KM, table -- the order the four
        # things appear going down the tab. A bar whose order differs from the figures'
        # makes the reader map button to figure by label, which is exactly the work the
        # ordering can do for free.
        #
        # ALL FOUR ARE uiOutputs (2026-08-15). Until then two were bare downloadButtons,
        # present from page load with no query behind them: fetching either before the first
        # Run returned HTTP 500, 140 bytes of HTML, and no file -- measured, not inferred.
        # That is the failed-download-with-no-message shape the KM button was made a
        # uiOutput to avoid in step 35; the same defect was simply sitting two buttons to its
        # left the whole time. One rule now covers the bar: a button exists only when the
        # thing it exports does. See output$dl_forest_ui and the three beside it.
        column(4, div(style = "margin-top: 25px; text-align: right;",
                      uiOutput("dl_forest_ui", inline = TRUE),
                      uiOutput("dl_tn_ui", inline = TRUE),
                      uiOutput("dl_km_ui", inline = TRUE),
                      uiOutput("dl_table_ui", inline = TRUE)))
      ),
      # In its own row, NOT a bare sibling of the one above. Bootstrap 3 columns are FLOATS, so
      # a bare block here flows AROUND them: with the taller chip selector the hint's first line
      # started beside the cohorts and its second line jumped back to the left margin, reading
      # as two unrelated fragments. A fluidRow clears the floats.
      fluidRow(column(12, uiOutput("ph_hint")))
    ),
    # The badge lines are multi-sentence PROSE, so they keep a capped measure instead of
    # running the full 1400px: a 200-character line is not readable, and these are the notes
    # most likely to be skipped when they are hard to read.
    #
    # Order is gene-context first (DepMap, CPTAC detectability -- properties of the gene,
    # true whatever this query returned), then the Step D/E/F badges, which are results
    # ABOUT this query. DepMap leads because it is the only one shown on all five tissues.
    div(class = "badge-block",
        uiOutput("depmap_badge"),
        uiOutput("cptac_badge"),
        uiOutput("cptac_validation_badge"),
        uiOutput("immune_mediation_badge"),
        uiOutput("cna_mediation_badge")),
    # Shape follows k via the same formula as the PDF export (plots.R), so the on-screen
    # forest and the downloaded one really are the same figure. The literal this replaced
    # (380px) still drew ovarian's 13 rows legibly -- the point is that two copies of the
    # sizing rule drifted apart, not that either one crashed. It took step 56 to make the
    # claim true rather than half true: a HEIGHT from forest_height_in(k) with the width
    # left to the browser gives the same rows at the window's aspect, not the design's.
    # The forest AND the tumor/normal panel, side by side (2026-08-13; it sat below the
    # table until then). The row is one output because the forest's width depends on
    # whether the panel is there -- two independent outputs cannot agree on that without
    # one of them reading the other's state anyway. The PLOTS stay separate outputs, so a
    # failure in either draws in its own box.
    uiOutput("forest_slot"),
    # A figure's help exists only when the figure does -- the SAME rule the export bar
    # follows, and for the same reason: an explanation of a plot that was never drawn is
    # an instruction for reading nothing. All of these are gated on the run counter below.
    uiOutput("sq_help_forest"),
    uiOutput("tau_note"),
    uiOutput("pooled_note"),
    uiOutput("scan_rank"),
    uiOutput("excluded_note"),
    # KM curves and the per-cohort table, full page width. A KM curve is read left-to-right
    # along time, so the width is the useful axis; the plot height still comes from
    # km_size_in() (synced to the PDF export), only width grows.
    fluidRow(column(12, uiOutput("sq_help_km"), uiOutput("km_tabs"),
                    uiOutput("sq_help_hr"), DTOutput("surv_table")))
   )
  ),

  # ---- Multiple query -------------------------------------------------------------
  #
  # ONE feature, every tissue, side by side -- and NOTHING is pooled across them. That is
  # a structural property, not a discipline: the server below makes one get_survival()
  # call per cancer_type with that tissue's own cohort list, and the engine's meta-pool
  # runs strictly within a cancer_type (the same firewall that keeps LUAD and LUSC apart).
  # Five independent survresult objects share no state, so there is no code path that
  # could combine them. Five panels, not four: LUAD and LUSC are separate cancer types by
  # design, and this view shows that split rather than hiding it under "lung".
  #
  # There is deliberately NO horizon control here. tau is declared PER (tissue, endpoint)
  # -- ovarian OS 48, luad OS 73/DFS 99, lusc OS 87/DFS 138, coad OS 142/DFS 103, breast
  # none -- so one global selector would run four of the five panels at another panel's
  # horizon and silently change the estimand. Each panel takes horizon_for(ct, endpoint);
  # the value travels on that panel's own x-axis (.forest_xlab) and note.
  tabPanel(TAB_MULTI,
    # narrow-page like the other three tabs. This was the one tab running edge-to-edge,
    # which is half of why it read as inconsistent beside them.
    div(class = "narrow-page",
     fluidPage(
      wellPanel(
        fluidRow(
          column(3, radioButtons("mq_kind", "Score type",
                                 choiceNames  = lapply(kind_choice_html(MQ_KINDS, KIND_CHOICE_LABEL), HTML),
                                 choiceValues = MQ_KINDS,
                                 selected = if ("viper" %in% MQ_KINDS) "viper" else MQ_KINDS[1])),
          column(3, selectizeInput("mq_feature", "Gene / TF symbol", choices = NULL,
                                   options = list(placeholder = "start typing a symbol...",
                                                  render = FEATURE_ITEM_RENDER))),
          # Endpoint is GLOBAL here -- the panels are only comparable if they answer the
          # same question. It opens on OS because that is the one endpoint every tissue
          # has a real pool for; the tissues' PRIMARY endpoints differ (breast DFS,
          # ovarian OS, luad DFS, lusc OS, coad DFS), which the note below states so the
          # tab cannot be mistaken for five primary analyses.
          column(2, selectInput("mq_endpoint", "Endpoint", c("OS", "DFS", "DSS"),
                                selected = "OS")),
          column(2, checkboxInput("mq_adjust_strata", "Stratify where available", value = TRUE),
                 tags$p(class = "note note-muted note-tuck", mq_strata_note())),
          column(2, actionButton("mq_run", "Run all panels", class = "btn-primary",
                                 style = "margin-top: 25px;"))
        ),
        # Prose only, NO numbers: the per-tissue horizons are stated on each panel, from
        # that panel's own result. A second copy of the tau values here would be a copy
        # that can drift from the one the query actually ran at.
        fluidRow(
          column(9,
            tags$p(class = "note note-body note-stack",
                   "Every tissue in the registry runs as its own panel, at its own declared ",
                   "follow-up horizon for the endpoint above. Nothing is pooled across panels. ",
                   "Each tissue's primary endpoint differs, so a single endpoint is not every ",
                   "panel's primary analysis — see the note on each panel.")),
          # The export appears only AFTER a run, for the same reason the panels do: a
          # download button with nothing behind it is a promise the tab cannot keep.
          column(3, style = "text-align:right;", uiOutput("mq_dl_ui")))
      ),
      # ONE forest help for the whole stack, not one per card: seven identical sentences
      # down a column is noise, and the panels are the same figure seven times. The KM and
      # HR help sit inside the per-tissue detail drawer instead, where they appear only for
      # a tissue a reader has actually opened.
      uiOutput("mq_help_forest"),
      uiOutput("mq_panels")
     ))
  ),

  # ---- Multiple genes ---------------------------------------------------------------
  #
  # The TRANSPOSE of the tab above: up to MGENE_MAX genes, ONE tissue. That swap is what
  # makes it the simpler of the two -- one tissue means one endpoint, one cohort list and
  # one declared horizon, so there is no five-horizons bookkeeping and no cross-panel
  # pooling question to answer. A row here is a loop over the Single query tab's own path:
  # its statistics are byte-identical to querying that gene singly.
  #
  # THE MULTIPLICITY ANSWER IS THE SAME ONE, AND IT IS NOT A BH OVER THIS LIST. Ten
  # hand-picked genes are not a family -- correcting within them would give a "corrected" q
  # that depends on which nine genes happened to be typed alongside. The anchor is the rank
  # each gene ALREADY holds in the genome-wide scan (R/scan_lookup.R), the same anchor the
  # Single query tab reports, and it is stated per gene in the note column.
  #
  # NO horizon control and NO cohort selector, unlike the Single query tab. Both are
  # deliberate: this tab runs cohorts_for(ct, ep) at horizon_for(ct, ep), which is exactly
  # the recipe the genome-wide scan was computed at, so the rank anchor above applies to
  # every row. A selector that let either drift would suppress the anchor on the tab whose
  # whole multiplicity story rests on it. The values are disclosed above the table, read
  # off the fitted results rather than off the controls.
  #
  # Tables only, no plots (2026-08-09). Ten forests is not a comparison a reader can hold;
  # the per-gene forest and KM curves are one click away on the Single query tab.
  tabPanel(TAB_MGENE,
    div(class = "narrow-page",
     fluidPage(
      wellPanel(
        fluidRow(
          column(3, selectInput("mg_ct", "Cancer type",
                                choices = setNames(CANCER_TYPES, cancer_label_of(CANCER_TYPES)),
                                selected = DEFAULT_CT)),
          column(3, uiOutput("mg_kind_ui")),
          column(2, uiOutput("mg_endpoint_ui")),
          column(2, div(style = "margin-top: 25px;", uiOutput("mg_strata_ui"))),
          column(2, actionButton("mg_run", "Run genes", class = "btn-primary",
                                 style = "margin-top: 25px; width: 100%;"))
        ),
        # One free-text box, not a multi-select: the list arrives PASTED -- a spreadsheet
        # column, a figure legend, a comma list out of a paper -- and mgene_resolve_symbols()
        # takes all of those forms. A ten-slot selectize would make the commonest input
        # (paste) the one gesture the control cannot accept.
        fluidRow(
          column(12, textAreaInput("mg_genes",
                     sprintf("Gene / TF symbols (up to %d)", MGENE_MAX), rows = 2,
                     width = "100%",
                     placeholder = "ESR1, MKI67, TP53 - commas, spaces, or one per line"))
        ),
        # The export sits WITH the controls, as on the Single query tab, not under the
        # table: it belongs to the query, and it has to be findable before a reader has
        # scrolled past ten rows.
        fluidRow(
          column(9, tags$p(class = "note note-body note-stack",
            "Each gene is analysed on its own over this tissue's cohorts for the endpoint ",
            "above, at that tissue and endpoint's declared horizon. Nothing is pooled across ",
            "genes, and no correction is applied across the list — a hand-picked set of ",
            "ten is not a multiple-testing family. The anchor is each gene's rank in the ",
            "genome-wide scan, reported per row.")),
          column(3, div(style = "text-align: right;", uiOutput("mg_dl_ui"))))
      ),
      uiOutput("mg_notes"),
      uiOutput("mg_help_hr"),
      uiOutput("mg_table_slot")
     ))
  ),

  tabPanel(TAB_RPPA,
    div(class = "narrow-page",
     fluidPage(
      wellPanel(
        fluidRow(
          column(3, selectInput("rp_ct", "Cancer type",
                                choices = setNames(RPPA_TISSUES, cancer_label_of(RPPA_TISSUES)),
                                selected = if (DEFAULT_CT %in% RPPA_TISSUES) DEFAULT_CT
                                           else RPPA_TISSUES[1])),
          column(3, uiOutput("rp_gene_ui")),
          column(2, selectInput("rp_endpoint", "Endpoint", choices = RPPA_ENDPOINTS,
                                selected = "OS")),
          # A Run button, like every other tab. Not decoration: the export bar invariant
          # (tests/test_export_bar.R) requires every downloadButton to sit behind a run
          # counter, because a live-reactive export can hand over a PDF of a query the
          # reader has already changed.
          column(2, actionButton("rp_run", "Run panel", class = "btn-primary",
                                 style = "margin-top: 25px; width: 100%;")),
          column(2, div(style = "text-align: right; margin-top: 25px;",
                        uiOutput("rp_dl_ui")))
        ),
        fluidRow(
          column(12, tags$p(class = "note note-body note-stack",
            "Protein and phospho-protein measured by reverse-phase protein array on TCGA ",
            "tumours. This is a SEPARATE layer, not another cohort: RPPA exists for TCGA ",
            "only, so every query here is a single cohort and nothing is pooled. Rows are ",
            "antibodies on the same patients \u2014 a phospho row reports the modified form of ",
            "the protein above it. HRs are unstratified and are not on the same footing as ",
            "this tissue's expression result for the same gene.")))
      ),
      uiOutput("rppa_overlap_note"),
      uiOutput("rppa_depmap_badge"),
      uiOutput("rppa_notes"),
      # No FOREST help here. The RPPA panel is one row per ANTIBODY with no pooled diamond
      # -- a deliberately different grammar (see rppa_panel_plot()) -- so the forest
      # explanation would be wrong on it rather than merely redundant.
      uiOutput("rppa_panel_slot"),
      uiOutput("rp_help_hr"),
      uiOutput("rppa_km_stack"),
      uiOutput("rp_help_km")
     ))
  ),

  # Browse is the DISCOVERY surface (batch scan / discovery-validation tables). Hidden
  # while the focus is the plot engine, but kept intact rather than deleted: the server
  # handlers below are untouched and flipping this constant restores the tab. The badges
  # do NOT depend on it -- they read SCANS_BY_CT directly.
  if (SHOW_BROWSE) tabPanel("Browse hits",
    fluidPage(
      fluidRow(
        column(8, uiOutput("browse_view_ui")),
        column(4, actionButton("query_selected", "Query selected TF ->", style = "margin-top: 25px;"))
      ),
      uiOutput("browse_empty"),
      DTOutput("browse_table")
    )
  ),

  # LAST. Placed after the Browse block rather than before it so About stays last whether or
  # not SHOW_BROWSE is on -- a conditional tabPanel evaluates to NULL, which navbarPage drops,
  # but the ORDER of the ones that survive is still the order they are written in.
  tabPanel(TAB_ABOUT, div(class = "narrow-page", HTML(about_html(ABOUT_FACTS))))
)

# .cohort_note() and .excluded_reason_line() -- the two sentences that stand in for a
# missing curve -- moved to R/plots.R on 2026-08-07 and are covered by
# tests/test_app_note_helpers.R, along with the endpoint-resolution rule that used to be
# inline in .endpoint_now() below (now registry.R's resolve_endpoint()). app.R is Shiny
# with no test harness, so pure logic left inline HERE is unverified; that rule was set in
# Phase 7c for .tau_arg/.tau_tag and applies to these three for the same reason.

# ---- server -------------------------------------------------------------------

server <- function(input, output, session) {

  # The brand lockup is the Guide's only nav item, in effect (step 111). If this observer
  # goes, the Guide is still the landing page and still unreachable afterwards -- a dead end
  # rather than a visible break, which is why tests/test_guide_derived.R asserts the pair
  # (hidden item + this handler) rather than either half.
  observeEvent(input$brand_home, updateNavbarPage(session, "nav", selected = TAB_GUIDE))

  # ---- hover help ON the figures (2026-08-30, step 117) -------------------------------
  # The explanations themselves are FIGURE_INFO in R/plots.R, declared beside the code that
  # draws each figure and length-capped by info_note_html(). Only the wiring is here.
  #
  # ONE rule, applied six times rather than six gates written six ways: a figure's help
  # exists exactly when the figure does. That is the same rule the export bar follows, and
  # it matters for the same reason -- an explanation of how to read a plot, sitting above a
  # tab that has never been run, is an instruction for reading nothing. The gate is the
  # tab's RUN COUNTER, not its result: a query that refused still drew no figure, and the
  # refusal message is the thing that belongs on screen at that moment.
  #
  # The MULTIPLE QUERY drawer and its table do NOT come through here. They are built inside
  # output$mq_panels, which already exists only for a tissue that produced a panel, so a
  # second gate there would be a gate on a gate.
  .fig_help <- function(which, counter) renderUI({
    n <- counter()
    if (is.null(n) || n == 0) return(NULL)
    HTML(figure_help_html(which))
  })
  output$sq_help_forest <- .fig_help("forest", reactive(input$run))
  output$sq_help_km     <- .fig_help("km",     reactive(input$run))
  output$sq_help_hr     <- .fig_help("hr",     reactive(input$run))
  output$mq_help_forest <- .fig_help("forest", reactive(input$mq_run))
  output$mg_help_hr     <- .fig_help("hr",     reactive(input$mg_run))
  output$rp_help_hr     <- .fig_help("hr",     reactive(input$rp_run))
  output$rp_help_km     <- .fig_help("km",     reactive(input$rp_run))

  ct <- reactive({ req(input$cancer_type); input$cancer_type })

  # The endpoint to resolve cohorts against. The RULE is registry.R's resolve_endpoint()
  # (tested); what stays here is only the reactive read of input$endpoint, which is the
  # part that cannot leave a Shiny session.
  .endpoint_now <- function(ct_) resolve_endpoint(input$endpoint, ct_)

  # --- tissue-driven controls ---
  output$kind_ui <- renderUI({
    k <- kinds_of(ct())
    radioButtons("kind", "Score type",
                 choiceNames  = lapply(kind_choice_html(k, KIND_CHOICE_LABEL), HTML),
                 choiceValues = k,
                 selected = if ("viper" %in% k) "viper" else k[1])
  })
  output$endpoint_ui <- renderUI({
    selectInput("endpoint", "Endpoint", endpoints_of(ct()),
                selected = default_endpoint_of(ct()))
  })
  # CHIP TOGGLES (2026-08-07). Every cohort is on screen at all times and a chip is a toggle:
  # selecting and dropping are the SAME gesture. The tag list this replaced was asymmetric --
  # dropping was an "x" on the tag, re-adding meant opening the dropdown and finding the name
  # again -- and that asymmetry is what made an accidental drop hard to undo.
  #
  # It is a checkboxGroupInput underneath, with the native box hidden by CSS. That is
  # deliberate: the input binding stays Shiny's own, so `input$cohorts` is still a plain
  # character vector and NO call site changed, exactly as when the checkbox column became a
  # tag list. Only the paint is different.
  #
  # ENDPOINT-AWARE: the choices are cohorts_for(ct, ENDPOINT), not cohorts_of(ct). A cohort
  # that does not carry the selected endpoint contributes nothing, and offering it PRE-TICKED
  # meant the user learned that only after pressing Run, from the excluded note. It is a
  # registry fact, known before the query -- 13 of the 15 (tissue, endpoint) pairs have one.
  #
  # Those cohorts are drawn STRUCK THROUGH IN PLACE rather than dropped from the list: a
  # cohort that simply vanishes is indistinguishable from one that never existed, which is the
  # quiet default this project keeps banning. They sit outside the checkbox group, so they
  # cannot be selected and cannot reach input$cohorts -- the chip is a statement, not a
  # control. The line beneath gives the reason and the count only; the names are already on
  # screen, and printing twelve of them in prose was the wall this presentation removes.
  output$cohorts_ui <- renderUI({
    ct_ <- ct()
    ep  <- .endpoint_now(ct_)
    all_co <- cohorts_of(ct_)
    co     <- cohorts_for(ct_, ep)
    dead   <- setdiff(all_co, co)
    tagList(
      checkboxGroupInput("cohorts",
                         sprintf("Cohorts (%d of %d carry %s)", length(co), length(all_co), ep),
                         choices = co, selected = co, inline = TRUE, width = "100%"),
      if (length(dead))
        div(class = "chip-row",
            lapply(dead, function(x)
              tags$span(class = "chip-dead",
                        title = sprintf("%s carries %s", x,
                                        paste(cohort_endpoints(x), collapse = ", ")), x))),
      if (length(dead))
        div(class = "note note-warn note-gap", .excluded_reason_line(dead, ep))
    )
  })
  # Only offer stratification where the tissue actually has a stratifier.
  output$strata_ui <- renderUI({
    # Per ENDPOINT (step 102). lgg is idh on OS and none on DFS, so the control and its
    # "no stratifier" message both have to follow the endpoint, not the tissue.
    sv <- strata_vars_of(ct(), input$endpoint)
    if (!length(sv)) return(tags$p(class = "note note-muted",
      sprintf("No stratifier for %s on %s - models are unstratified.",
              cancer_label_of(ct()), input$endpoint %||% "this endpoint")))
    # The label is DERIVED from the registry, not written here. It said "Stratify by PAM50
    # subtyping" for every tissue until 2026-08-27; the moment lgg gained stratify_by = idh
    # (step 88) that checkbox applied strata(idh) to glioma while naming PAM50 -- a control
    # that misreports the model it fits. The earlier fix in this spot removed a stray
    # sprintf() argument left over from a label that USED to be built from `sv`; deriving it
    # again is that idea restored, with the arity right.
    tagList(
      checkboxInput("adjust_strata", strata_label_of(ct(), input$endpoint), value = TRUE),
      { n <- strata_coverage_note(ct(), input$endpoint)
        if (is.null(n)) NULL else tags$p(class = "note note-muted note-tuck", n) },
      # What the strata CONTAIN, not just how far they reach. Derived in R/plots.R, where
      # it can be tested; the hover carries the WHO 2021 fact that no survival table can.
      { n <- idh_stratum_note(ct(), input$endpoint)
        if (is.null(n)) NULL else tags$p(class = "note note-warn note-tuck",
                                         HTML(info_note_html(n, IDH_STRATUM_INFO))) })
  })

  observeEvent(list(input$kind, input$cancer_type), {
    req(input$kind)
    updateSelectizeInput(session, "feature",
                         choices = FEATURES_BY_CT[[ct()]][[input$kind]], server = TRUE)
  }, ignoreNULL = FALSE)

  # Switching tissue OR endpoint re-resolves the horizon: the declared panel default for
  # the (tissue, endpoint) you moved TO, not whatever was selected before. Leaving the old
  # value in place would silently run e.g. a LUAD DFS query at the LUAD OS horizon (73 vs
  # 99), the same app-vs-artifact mismatch tau_choices exists to prevent. input$endpoint is
  # rendered dynamically (endpoint_ui), so req() guards the moment before it exists; the
  # endpoint also re-defaults when the tissue changes, which re-fires this observer.
  observeEvent(list(input$cancer_type, input$endpoint), {
    req(input$endpoint)
    tc <- tau_choices(ct(), input$endpoint)
    updateSelectInput(session, "max_fu", choices = tc$choices, selected = tc$selected)
  }, ignoreInit = TRUE)

  # --- badges (breast Step D/E/F results; absent for other tissues) ---
  .scans <- reactive(SCANS_BY_CT[[ct()]])
  .row_for <- function(df, tf) {
    if (is.null(df) || !("tf" %in% names(df))) return(NULL)
    i <- match(tf, df$tf)
    if (is.na(i)) return(NULL)
    df[i, ]
  }

  # Gene context, not a result: shown for every tissue, and it does NOT wait for a Run.
  # It answers a question about the gene the selector is showing, which is answerable
  # before any model is fit -- gating it on a result would imply it describes one.
  output$depmap_badge <- renderUI({
    # req() IS right here, unlike in res() below, and the difference is the output type.
    # A cancelled renderUI sends null and EMPTIES its div; it is renderPlot and renderDT
    # that keep the last value on screen. So emptying the selector clears these five
    # badges on its own -- checked on the page against the pre-step-57 source, after
    # waiting out the recalculating state, which is the step this measurement needs.
    req(input$feature, input$kind)
    n <- depmap_note(input$feature, input$kind, DEPMAP_LISTS)
    if (is.null(n)) return(NULL)
    tags$p(class = paste("note note-tuck", paste0("note-", n$tone)),
           HTML(info_note_html(n$text, n$info)))
  })

  output$cptac_badge <- renderUI({
    req(input$feature)
    if (CPTAC_N == 0 || ct() != "breast") return(NULL)
    n <- unname(CPTAC_COVERAGE[input$feature])
    # Both branches carry the same hover explanation: "not quantified" is exactly the case
    # where a reader is most likely to read absence of protein into absence of detection.
    if (is.na(n) || n == 0) {
      tags$p(class = "note note-muted note-tuck",
            HTML(info_note_html(
              sprintf("CPTAC-BRCA proteomics (independent MS cohort, n=%d): not quantified for this gene.",
                      CPTAC_N), CPTAC_BADGE_INFO)))
    } else {
      tags$p(class = "note note-good note-tuck",
            HTML(info_note_html(
              sprintf("CPTAC-BRCA proteomics: protein detected in %d/%d tumors (%.0f%%).",
                      n, CPTAC_N, 100 * n / CPTAC_N), CPTAC_BADGE_INFO)))
    }
  })

  output$cptac_validation_badge <- renderUI({
    req(input$feature)
    r <- .row_for(.scans()$cptac_protein_validation, input$feature)
    if (is.null(r) || is.na(r$p_viper)) return(NULL)
    sig <- r$p_viper < 0.05
    tags$p(class = paste("note note-tuck", if (sig) "note-good" else "note-muted"),
          sprintf("CPTAC protein validation (Step E): VIPER activity vs. measured protein r=%.2f, p=%.2g (n=%d).",
                  r$cor_viper_vs_protein, r$p_viper, r$n))
  })

  output$immune_mediation_badge <- renderUI({
    req(input$feature)
    r <- .row_for(.scans()$immune_mediation_summary, input$feature)
    if (is.null(r)) return(NULL)
    # Same rule as the CNA badge below, and now literally the same code: it used to be
    # written out here a second time (2026-08-18, step 48).
    mediated <- mediation_is_explained(r)
    tags$p(class = paste("note note-tuck", if (mediated) "note-warn" else "note-good"),
          sprintf("Immune mediation test (Step D): signal %s explained by %s infiltration (HR %.2f -> %.2f after adjusting, p=%.2g).",
                  if (mediated) "IS" else "NOT", r$immune_feature, r$HR_unadj, r$HR_adj, r$p_adj))
  })

  output$cna_mediation_badge <- renderUI({
    req(input$feature)
    r <- .row_for(.scans()$cna_mediation_summary, input$feature)
    if (is.null(r)) return(NULL)
    # The threshold rule and the wording both live in cna_mediation_note() (R/plots.R) now,
    # where a test can reach them: this is a VERDICT produced by three thresholds, and it
    # was previously written out here, untested, next to the sprintf that renders it.
    n <- cna_mediation_note(r)
    tags$p(class = paste("note note-tuck", if (n$mediated) "note-warn" else "note-good"),
           HTML(info_note_html(n$text, n$info)))
  })

  # --- the query ---
  res <- eventReactive(input$run, {
    # THE EMPTY SELECTOR IS NOT A req() CASE (step 57). A selectize can be emptied by
    # clicking into it and pressing backspace, so input$feature arriving here as "" is a
    # state the user can reach -- and req("") CANCELS, which on a Run button means the
    # click does NOTHING AT ALL: no result, no message, no console error, nothing that
    # distinguishes a refused query from a slow one. Measured on all four Run buttons
    # before this step; three were silent this way and only Multiple genes said why.
    # validate() carries the reason to every reader of res() instead, and those readers
    # already exist: the forest slot prints it, output$surv_table catches
    # `shiny.silent.error` and captions an empty table with it, and each KM box shows it.
    # input$endpoint and input$kind stay on req() -- they are selectInput/radioButtons
    # over fixed choices, so they have no empty state to reach, and req() is still the
    # right guard for the moment before a dynamically rendered control exists.
    req(input$endpoint)
    validate(need(nzchar(input$feature), "No gene or TF selected - pick a symbol above."))
    # NOT req(): the empty selection is reachable by removing the tags one at a time, and
    # req() CANCELS the render -- a cancelled forest KEEPS THE PREVIOUS QUERY'S PLOT, so the
    # last result would sit on screen with nothing selected under it. validate() replaces it
    # with the reason instead. (This outlived the All/None links that prompted it; the stale
    # render is the bug, the links were only the shortest route to it.)
    validate(need(length(input$cohorts) >= 1, "No cohorts selected - pick at least one."))
    score <- get_feature(input$feature, cohorts = input$cohorts, kind = input$kind)
    validate(need(length(score) > 0,
                  sprintf("'%s' not found in any selected cohort's %s data.", input$feature, input$kind)))
    # No stratify control for an unstratified tissue -> input$adjust_strata is NULL;
    # isTRUE() makes that read as FALSE rather than erroring.
    r <- get_survival(score, endpoint = input$endpoint, cohorts = input$cohorts,
                      adjust_strata = isTRUE(input$adjust_strata),
                      max_followup = .tau_arg(input$max_fu))
    validate(need(any(vapply(r$per_cohort, function(x) isFALSE(x$skipped), logical(1))),
                  "No selected cohort had enough events for this query."))
    with_feature(r, sprintf("%s (%s)", input$feature, KIND_LABEL[[input$kind]]))
  })

  # --- tumor vs normal-adjacent (context, not a survival result) ---
  #
  # Captured on the SAME run as res(), off input$run rather than live inputs, so the panel
  # always describes the query drawn above it. It follows the user's SELECTED cohorts: if
  # TCGA_BRCA is unticked the comparison genuinely is not available for that query, and
  # showing it anyway would report a cohort the result above does not include.
  #
  # The kind guard is the app side of tumor_normal()'s refusal to accept "cna". That
  # function stops rather than returning an empty result, because copy number called
  # against a matched normal has no tumor-vs-normal contrast to return -- so the caller
  # must not ask.
  #
  # NOTE, corrected 2026-08-14: this guard is currently UNREACHABLE from the UI, and the
  # earlier claim here that dropping it "turns a CNA query into an app-level error" was
  # wrong. QUERY_KINDS is a strict SUBSET of TN_KINDS -- viper and expr against TN_KINDS's
  # viper, expr and immune. It was set EQUALITY until step 108 withdrew immune from the
  # selectors; shrinking QUERY_KINDS moves in the safe direction, and the selector has
  # offered no copy-number query since CNA was removed as a query kind (see QUERY_KINDS) --
  # so kind is always in TN_KINDS by the time this runs. What the line actually guards is the day
  # QUERY_KINDS gets "cna" back: without it that restoration turns the panel into an app
  # error, at a call site nobody would think to look at. It stays for that reason, and the
  # step-42 note recording a live "CNA query" is corrected in BUILD_LOG -- that state was
  # reached by setting the input out of band, not by anything a user can click.
  tn_res <- eventReactive(input$run, {
    req(input$feature, input$kind)
    if (!(input$kind %in% TN_KINDS)) return(NULL)
    if (!length(input$cohorts)) return(NULL)
    # NOT wrapped in tryCatch(). The two returns above are the DELIBERATE absences, and
    # they are the only ones: a tryCatch here would additionally swallow every unforeseen
    # error into the same silent empty panel, which is indistinguishable from "this query
    # has no comparison to show" -- the quiet default this project bans. An error here
    # belongs on screen, in this slot, where it can be seen and fixed.
    tumor_normal(input$feature, cohorts = input$cohorts, kind = input$kind)
  })

  # What the row should LOOK like -- and the ONLY place a tumor/normal failure is caught.
  #
  # The forest now shares a container with the panel, so anything this reactive throws would
  # take the survival result down with it. That is the boundary R/tumor_normal.R exists to
  # keep: a context panel must never be able to remove an outcome figure from the screen. So
  # layout is decided defensively here, while output$tn_panel below still calls tn_res()
  # WITHOUT a catch -- on failure this reserves the column and the real error draws inside
  # the panel's own box, which is loud, in the right place, and costs the forest nothing.
  # `slot` and `dl` are SEPARATE fields, and the error branch is the only place they differ:
  # it reserves the column (slot) while offering no button (dl). "Give this failure somewhere
  # visible to draw" and "there is a figure here worth exporting" are different questions, and
  # collapsing them into one flag would put a download button in front of a panel that raised
  # -- the failed-download-with-no-message shape output$dl_km_ui exists to avoid.
  tn_layout <- reactive({
    tryCatch({
      tn <- tn_res()
      if (is.null(tn)) return(list(slot = FALSE, dl = FALSE, k = 0L, notes = NULL))
      est <- tn_estimable(tn)
      # No estimable cohort -> no panel AT ALL, not an empty one. Every cohort in this query
      # being unable to answer is not a finding about the gene, and a figure with no boxes in
      # it would read as one. The reasons still print, so the absence is never silent.
      list(slot = length(est) > 0L, dl = length(est) > 0L, k = max(1L, length(est)),
           notes = c(tn_notes(tn), tn_fc_note(tn)))
    }, error = function(e) list(slot = TRUE, dl = FALSE, k = 1L, notes = NULL))
  })

  output$forest_slot <- renderUI({
    r  <- tryCatch(res(), error = function(e) NULL)
    k  <- if (is.null(r)) 3 else max(1, .n_estimable(r))
    lay <- tn_layout()

    # ONE box, built once and placed in whichever branch runs (step 56). Not two calls: the
    # branches exist because the WIDTH differs -- full page, or what is left beside the
    # tumour/normal panel -- and the whole point of the fix is that the figure's proportions
    # must not. Two plotOutput() calls would be two places for that to stop being true, and
    # the difference would show up as a squashed forest in one branch only, which reads as a
    # data problem rather than a layout one. app.R now names "forest_plot" exactly once.
    #
    # The aspect goes inline because it is a function of k; everything about it that is NOT
    # a function of k -- the cap, the centring, contain:size -- is in .sq-forest, where the
    # comment explains why the last of those is load-bearing.
    box <- div(class = "sq-forest",
               style = sprintf("aspect-ratio: %g / %g;", FOREST_WIDTH_IN, forest_height_in(k)),
               plotOutput("forest_plot", height = "100%"))

    # No panel: the forest goes back to the FULL width rather than sitting in a two-thirds
    # column beside an empty one. Ovarian and every CNA query take this branch, so it is not
    # an edge case, and a permanently reserved blank column would also give the "why it is
    # missing" note top billing beside the forest title that an absence has not earned.
    if (!lay$slot)
      return(tagList(
        box,
        if (length(lay$notes))
          tags$p(class = "note note-muted note-tuck", paste(lay$notes, collapse = " "))))

    s <- tn_size_in(lay$k)
    w <- round(s[["width"]] * 96)
    div(class = "fp-row",
        div(class = "fp-forest", box),
        div(style = sprintf("flex: 0 0 %dpx;", w),
            plotOutput("tn_panel", width = sprintf("%dpx", w),
                       height = sprintf("%dpx", round(s[["height"]] * 96))),
            tags$p(class = "note note-muted note-tuck",
                   paste(c(TN_CAPTION, lay$notes), collapse = " "))))
  })

  # res = 96, NOT Shiny's default 72 (2026-08-15).
  #
  # The box above is sized `tn_size_in() * 96` -- inches converted to pixels at 96 per inch.
  # Shiny's default renderPlot then reads those pixels back at 72 per inch, so the device it
  # opens is 5.86in wide where the export is 4.4in. Text is sized in POINTS, so the same
  # figure gets 4.4/5.86 = 75% of its export text size on screen: the two stat lines above the
  # boxes are drawn at cex 0.62, which lands at 7px, and at 7px they are ink rather than
  # numbers. That is the reported symptom -- "the statistics are only in the download" -- and
  # they were never missing, only unreadable. Matching the two numbers makes the on-screen
  # panel the export, which is what tn_size_in()'s header already claims it is.
  #
  # ONLY THIS OUTPUT, and the boundary is not arbitrary: I checked all five renderPlots in
  # this file. The forest and the KMs take their HEIGHT from the same inch rule but their
  # WIDTH from the layout (100% of a flex or grid column), so no device size of theirs is the
  # export's and there is no promise here to keep. tn_panel is the only plot whose width AND
  # height both come from a sizing function, so it is the only one that can be the same
  # figure in both places.
  output$tn_panel <- renderPlot({
    tn <- req(tn_res())
    validate(need(length(tn_estimable(tn)) >= 1, "Nothing to plot."))
    tn_plot(tn)
  }, res = 96)

  # OFFERED ONLY WHEN THE PANEL IS ACTUALLY ON SCREEN, and off the SAME reactive that decides
  # whether to draw it -- not a second reading of tn_res(), so button and figure cannot drift
  # into disagreeing about whether there is anything to export. Same two-gate form as
  # output$dl_km_ui, and gated on the RUN COUNTER first for the reason recorded there: before
  # the first click the eventReactive raises a silent condition, which tn_layout()'s tryCatch
  # turns into the error branch. That branch already answers dl = FALSE, so the counter is
  # belt and braces -- but it means the button's absence before a run does not depend on how
  # req()'s condition happens to be classified.
  output$dl_tn_ui <- renderUI({
    if (is.null(input$run) || input$run == 0) return(NULL)
    if (!isTRUE(tn_layout()$dl)) return(NULL)
    downloadButton("dl_tn", "Tumor vs normal (PDF)", class = "btn-export")
  })

  # NAMED ENTIRELY FROM THE RESULT -- tn_export_name() reads feature, kind and cohorts off
  # tn_res() and nothing off a live input. That is the rule mq_res() states and the reason it
  # carries feature_id: inputs can be changed WITHOUT pressing Run, so a name built from them
  # can label a file with a query the file does not contain, and a filename outlives the
  # session that made it. The other three Single query buttons still read input$feature live;
  # they predate the rule and are not touched here.
  #
  # The name is tn_export_name()'s and not a sprintf here, for the same reason mg_dl_table
  # calls mgene_export_name(): app.R has no test harness, so a rule written inline in this
  # file can only ever be checked by grepping source.
  output$dl_tn <- downloadHandler(
    filename = function() tn_export_name(tn_res()),
    content  = function(file) tn_plot(tn_res(), file = file)
  )

  # (.n_estimable() -- the per-query count of estimable cohorts, not the tissue's count --
  # moved to R/plots.R on 2026-08-07 alongside mq_panel_notes(), which needs it.)

  output$forest_plot <- renderPlot({
    req(res())
    validate(need(.n_estimable(res()) >= 1, "Nothing to plot."))
    forest_plot(res())
  })

  # k is a per-query quantity, not a tissue constant: a regulator below VIPER's
  # minimum regulon size drops out of individual cohorts, so the ovarian panel can
  # pool 12 cohorts for one TF and fewer for the next. State it rather than let the
  # user assume the cohort count they ticked is the k they got.
  output$pooled_note <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL)
    req(r, r$pooled)
    tags$p(class = "note note-body",
      sprintf("Pooled over k=%d of %d selected cohorts (%s, %s).",
              r$pooled$k, length(input$cohorts), ct(), r$endpoint))
  })

  # Where this query's TF sits in the genome-wide scan already on disk -- the multiplicity
  # anchor the Query tab used to lack. Captured on the SAME run as res() so the rank
  # describes exactly the result shown above it, and computed by lookup, never recompute:
  # the strict match guarantees the query's live p equals the scan's stored p, so the
  # stored rank/q is this query's. On any recipe mismatch it says so rather than paste a
  # stale q (scan_rank_lookup / format_scan_rank live in R/scan_lookup.R -- app.R untested).
  scan_rank <- eventReactive(input$run, {
    # req(), NOT the validate() that res() and mq_res() carry (step 57). This is a NOTE,
    # not a result: with no feature selected, scan_rank_lookup("") answers honestly that
    # nothing by that name was retained in the scan, and format_scan_rank() renders the
    # answer as "<feature> was not retained in the breast DFS PAM50-adjusted scan" -- a
    # sentence with NO SUBJECT, which reads as a finding about the previous gene. Measured
    # on the page, and it survived the first cut of this step. A third copy of "no gene
    # selected" is not the fix either: the forest slot and the table caption already print
    # it, so the right answer here is silence, and req() is how a renderUI says nothing.
    req(nzchar(input$feature))
    ep <- input$endpoint; ct_ <- ct()
    # EFFECTIVE stratification, mirroring get_survival's own gate: stratify only if the
    # tissue actually has a stratifier. Not raw input$adjust_strata -- Shiny keeps a
    # removed input's last value, so switching from breast (checkbox on) to ovarian (no
    # checkbox) leaves input$adjust_strata==TRUE stale. The engine ignores it for ovarian
    # (stratify_by=none), so the lookup must too, or it hunts for an adjusted scan the
    # unstratified query never used.
    strat_eff <- isTRUE(input$adjust_strata) && length(strata_vars_of(ct_, ep)) > 0
    scan_rank_lookup(
      feature = input$feature, kind = input$kind, cohorts = input$cohorts,
      max_followup = .tau_arg(input$max_fu), adjust_strata = strat_eff,
      strat_var = strata_var_for_scan(ct_, ep),
      endpoint = ep, cancer_type = ct_, scans = SCANS_BY_CT[[ct_]],
      # horizon_for(ct_, ep), NOT the bare horizon_for(ct_): since horizons became
      # per-endpoint (2026-07-27) the bare call returns min() over the tissue's taus, which
      # is the OTHER endpoint's window for luad/DFS (73 vs the scan's 99) and coad/OS (103
      # vs 142). That silently suppressed the rank on two of the seven scans -- including
      # LUAD's primary endpoint -- because the lookup correctly refused a horizon it could
      # not match. It failed SAFE, which is why nothing caught it. tests/test_scan_rank_wiring.R
      # pins this against each scan's own recorded max_fu.
      expected_cohorts = cohorts_for(ct_, ep), expected_horizon = horizon_for(ct_, ep))
  })
  output$scan_rank <- renderUI({
    req(input$run)
    r <- tryCatch(scan_rank(), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    f <- format_scan_rank(r)
    if (is.null(f)) return(NULL)
    # format_scan_rank()'s tone, mapped to the shared note tones rather than to hex here.
    # The Multiple query panels paint the same three from the same source, which is the
    # point of naming them in PAGE_CSS: a rank that clears FDR cannot end up green on one
    # tab and orange on the other.
    tone <- switch(f$tone, hit = "note-good", miss = "note-warn", na = "note-muted")
    # The marker appears only when format_scan_rank() supplied an explanation, which is the
    # `ok` state alone -- the five "no rank, because..." messages already say why, and
    # deciding that HERE would mean app.R re-deriving which state it is looking at.
    tags$p(class = paste("note note-tuck", tone),
           if (is.null(f$info)) f$text else HTML(info_note_html(f$text, f$info)))
  })

  # Progressive disclosure: the horizon control is only worth reaching for when an
  # effect actually moved during follow-up, so say so exactly then and stay silent
  # otherwise. Suppressed once a horizon is set — at that point the user has acted, and
  # a hint that kept nagging would be an invitation to keep sliding.
  output$ph_hint <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL)
    req(r, is.null(r$max_followup))
    bad <- Filter(function(x) isTRUE(x$ph_violated), r$per_cohort)
    if (!length(bad)) return(NULL)
    tags$p(class = "note note-warn note-tuck",
      sprintf("Proportional hazards rejected in %s. The HR there is a time-average of an effect that changed during follow-up; a shorter horizon may describe it better.",
              paste(names(bad), collapse = ", ")))
  })

  # A truncated result must never be readable as a full-follow-up one, in the app as
  # much as in the exported figures.
  output$tau_note <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL)
    req(r, r$max_followup)
    tags$p(class = "note note-body",
      sprintf("Follow-up truncated at %g months: patients still at risk are censored there, not excluded. HRs below describe that window only.",
              r$max_followup))
  })

  output$excluded_note <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL)
    req(r)
    usable <- names(Filter(function(x) isFALSE(x$skipped), r$per_cohort))
    missing <- setdiff(input$cohorts, usable)
    if (!length(missing)) return(NULL)
    tags$p(class = "note note-warn",
          paste(vapply(missing, .cohort_note, character(1), r = r), collapse = " "))
  })

  # KM tabs are rebuilt per tissue. The output slots themselves are registered ONCE
  # below for every cohort in the registry — registering them inside a renderUI would
  # re-create observers on every tissue switch and leak them.
  output$km_tabs <- renderUI({
    # Endpoint-aware for the same reason the selector is: with the selector no longer
    # OFFERING a cohort that cannot carry the endpoint, a KM tab for that same cohort sitting
    # below it -- reading "SCANB carries no DFS data" -- would be the two controls
    # contradicting each other about what is in the query.
    co <- cohorts_for(ct(), .endpoint_now(ct()))
    do.call(tabsetPanel, c(id = "km_tabset",
      # Height from km_size_in(), never a literal: the exported PDF uses the same function,
      # so the vertical proportions -- crucially the curve panel vs the number-at-risk
      # panel -- cannot drift between screen and file. 96 px/in is the CSS reference pixel
      # a plotOutput height is measured in.
      #
      # WIDTH is deliberately NOT synced: on screen it is fluid (100% of the container,
      # ~813px at a desktop width) because this is a browser-first tool, while the export
      # is a fixed 7in. So the two are the same height and the same layout, not the same
      # aspect ratio. Base graphics re-solve the layout per device, and the risk-table
      # alignment was checked at both aspects.
      lapply(co, function(x) tabPanel(x, plotOutput(paste0("km_", x),
        height = sprintf("%dpx", round(km_size_in()[["height"]] * 96)))))))
  })

  # A container that is hidden or not yet laid out reports its size as 0, and shiny:::startPNG
  # opens the device and calls plot.new() BEFORE the expression below is evaluated. A
  # zero-width (or zero-height) device fails there with
  #     Warning: Error in graphics::plot.new: figure margins too large
  # naming this output. No guard inside the render expression can prevent it -- res() has not
  # been touched yet when it fires. So floor the size Shiny asks the device for.
  #
  # Measured, not assumed: startPNG throws that error if and only if a dimension is exactly 0.
  # One pixel is fine. The fallback is therefore not a "big enough" guess -- it is the design
  # size from km_size_in(), the honest stand-in for "the browser has not told us yet", and it
  # is replaced by the real width the moment the client reports one.
  #
  # WIDTH stays fluid: the function returns the client's own value whenever it has one, so the
  # deliberate screen-vs-export width difference documented on output$km_tabs is preserved.
  .km_dim <- function(co, which) function() {
    v <- session$clientData[[sprintf("output_km_%s_%s", co, which)]]
    if (is.null(v) || !is.numeric(v) || !is.finite(v) || v < 1)
      round(km_size_in()[[which]] * 96)
    else v
  }

  for (coh in COHORTS$cohort) {
    local({
      co <- coh
      output[[paste0("km_", co)]] <- renderPlot(
        {
          r <- res()
          validate(need(!is.null(r$per_cohort[[co]]) && isFALSE(r$per_cohort[[co]]$skipped),
                        .cohort_note(co, r)))
          km_plot(r, co)
        },
        width  = .km_dim(co, "width"),
        height = .km_dim(co, "height"))
    })
  }

  # Display view, not the export: feature/endpoint/max_fu are already in the sidebar, so
  # on screen they are constant columns of noise. Values arrive pre-formatted from
  # survtable_display() (HR, CI and the PH "!" are strings), so no formatRound here.
  output$surv_table <- renderDT({
    # An explicit EMPTY TABLE for every refused query, not req() and not validate().
    # req() cancels the render, and a DT that is not re-rendered KEEPS THE PREVIOUS QUERY'S
    # ROWS -- deselecting every chip and pressing Run left the last query's three rows sitting
    # under a forest slot that correctly read "No cohorts selected". validate() is no better
    # here: renderDT does not surface a validation message, it blanks the cells and leaves the
    # empty row skeletons behind, which reads as a result with missing numbers. A real empty
    # table carrying the reason is the only form that cannot be misread. (Same failure as the
    # stale DT on the Multiple query tab, step 16, fixed the same way: leave nothing on screen
    # that can go stale.)
    #
    # Catching `shiny.silent.error` covers EVERY refusal res() can raise -- the empty cohort
    # selection, "'X' not found in any selected cohort's expr data", "no selected cohort had
    # enough events" -- rather than only the one that exposed the bug. A REAL error is not
    # caught and still shows in Shiny's red block: this must not quieten a genuine failure.
    d <- tryCatch(survtable_display(res()),
                  shiny.silent.error = function(e) conditionMessage(e))
    if (is.character(d)) {
      # A bare req() carries no message. That is the pre-run state, where there is nothing to
      # say and nothing stale to clear, so the slot stays empty exactly as it did before.
      req(nzchar(d))
      return(datatable(data.frame(Cohort = character(0)), rownames = FALSE, caption = d,
                       options = list(dom = "t")))
    }
    datatable(d, rownames = FALSE, caption = attr(d, "caption"),
              options = c(list(pageLength = 10, dom = "tip"),
                          survtable_dt_options(d)))
  })

  # .feature_tag() on the FEATURE only, and nothing else in these names. The other three
  # fragments are registry values, and that was checked rather than assumed: all 41 cohort +
  # cancer-type ids and all 3 endpoints are already `[A-Za-z0-9._-]`, longest 14 characters
  # (derived from config/cohorts.tsv, 2026-08-10). Sanitising them would be a no-op dressed up
  # as a guard. Same division mq_export_dir() makes.
  #
  # These three are the LAST of the five handlers step 35 measured, and the mildest -- which is
  # why they were left when mq_export_dir() and mq_dl_all were fixed in step 36. What they
  # actually buy, split the way step 35 split it:
  #   * the "/" in 5837 of 58366 offered features: TIDINESS. Measured end to end -- Chrome
  #     already saves "C2///CFB" as "C2___CFB", the vocabulary has ZERO collisions under that
  #     mapping, and the true symbol is inside the file. Nothing was broken; it is now
  #     consistent with the other three tabs instead of relying on browser behaviour that
  #     nothing promises.
  #   * the 353 names over NAME_MAX after that substitution: CORRECTNESS. (353 is what step
  #     35 OBSERVED failing, not a vocabulary count -- see .feature_tag()'s note, which records
  #     363 for the vocabulary and why the two are different quantities.) Those did not land
  #     at all -- served 200, no file in ~/Downloads, reproduced away from the app against a
  #     10-line static server. FEATURE_NAME_BUDGET is what fixes that, and it is the half of
  #     .feature_tag() that is not cosmetic.
  # THE BAR'S ONE RULE, applied to the forest and the table too (2026-08-15).
  #
  # Both were bare downloadButtons, on screen from page load. Before the first Run, res() is
  # an eventReactive that has not fired, so the handler raised and Shiny answered the fetch
  # with HTTP 500 and 140 bytes of HTML -- no file, and nothing on the page to say why.
  # Measured on a fresh session, not reasoned about. It is the same defect step 35 found on
  # the KM button, and the fix is the same two gates:
  #
  #   1. the RUN COUNTER, first, for the reason mq_dl_ui records -- before the first click an
  #      eventReactive raises a SILENT condition, which cancels the whole renderUI rather than
  #      returning from it, so a tryCatch alone cannot answer the pre-run question.
  #   2. the same condition the OUTPUT BESIDE IT uses, so button and figure cannot disagree
  #      about whether there is anything to export.
  #
  # tryCatch AFTER the counter, and it is not the swallowed-error kind: a validate() inside
  # res() is a condition, so without this the message prints INSIDE THE EXPORT BAR, once per
  # button -- which is what dl_km_ui did (measured: unticking every cohort put "No cohorts
  # selected - pick at least one." in the bar and wrapped Table onto its own line). The reason
  # belongs in the forest slot, where it is already printed; a bar has no room for prose and
  # saying it three more times is the same sentence four times over.
  output$dl_forest_ui <- renderUI({
    if (is.null(input$run) || input$run == 0) return(NULL)
    r <- tryCatch(res(), error = function(e) NULL)
    if (is.null(r) || .n_estimable(r) < 1) return(NULL)
    downloadButton("dl_forest", "Forest (PDF)", class = "btn-export")
  })

  # The CSV, not the figure, so the gate is weaker on purpose: survtable() has a row for every
  # selected cohort INCLUDING the ones that did not estimate, and that row -- the reason a
  # cohort produced nothing -- is worth exporting even when no forest can be drawn. So this
  # asks only that the query returned, not that anything was estimable.
  output$dl_table_ui <- renderUI({
    if (is.null(input$run) || input$run == 0) return(NULL)
    r <- tryCatch(res(), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    downloadButton("dl_table", "Table (CSV)", class = "btn-export")
  })

  output$dl_forest <- downloadHandler(
    filename = function() sprintf("forest_%s_%s_%s%s.pdf", ct(), .feature_tag(input$feature),
                                  input$endpoint, .tau_tag(res())),
    content  = function(file) forest_plot(res(), file = file)
  )
  output$dl_table <- downloadHandler(
    filename = function() sprintf("survtable_%s_%s_%s%s.csv", ct(), .feature_tag(input$feature),
                                  input$endpoint, .tau_tag(res())),
    content  = function(file) write.csv(survtable(res()), file, row.names = FALSE)
  )
  # OFFERED ONLY WHEN THE SELECTED TAB HAS A CURVE TO EXPORT.
  #
  # output$km_tabs builds a tab for every cohort in cohorts_for(ct, endpoint) -- every cohort
  # CARRYING the endpoint, not every cohort that estimated -- because a tab reading
  # "TCGA_OV: no usable data for this query." is the disclosure, and dropping those tabs would
  # hide why a cohort is missing. The plot output guards that case with validate(need(...)).
  # This download did not, so km_plot() ran on a cohort with no fit, stopped, and Shiny
  # returned HTTP 500 with nothing anywhere on the page (measured 2026-08-10, BUILD_LOG step
  # 35: reachable by pure clicking, and the DEFAULT-selected tab is the tissue's first cohort,
  # which for ovarian/expr is one of the four with no fit).
  #
  # The same shape mg_dl_ui was written to avoid on the Multiple genes tab: a button offered
  # after a run that resolved nothing hands the browser a failed download with no message.
  # Same two-gate form, and the same condition the plot beside it uses -- not a second reading
  # of it, so the two cannot drift.
  #
  # NULL rather than a disabled anchor, for two reasons. A Bootstrap `.disabled` <a> only
  # stops navigating because of a `pointer-events` rule, which is styling standing in for a
  # guard. And the tab immediately below already states why there is no curve, so a note here
  # would be the same sentence twice (the step-31 rule).
  # Gated on the RUN COUNTER, not on tryCatch(res()), for the reason mq_dl_ui records: before
  # the first click an eventReactive raises a silent condition that cancels the whole renderUI.
  # tryCatch added 2026-08-15, for the reason recorded above output$dl_forest_ui: the bare
  # res() here printed a validate() message inside the export bar.
  output$dl_km_ui <- renderUI({
    if (is.null(input$run) || input$run == 0) return(NULL)
    r <- tryCatch(res(), error = function(e) NULL); co <- input$km_tabset
    if (is.null(r) || is.null(co) || is.null(r$per_cohort[[co]])) return(NULL)
    if (!isFALSE(r$per_cohort[[co]]$skipped)) return(NULL)
    downloadButton("dl_km", "KM, selected cohort (PDF)", class = "btn-export")
  })

  # The KM the user is actually looking at -- the SELECTED tab, not a fixed cohort. The
  # cohort belongs in the filename for the same reason tau does (.tau_tag): two exports of
  # one feature/endpoint differ only by cohort, and a name that omits it silently
  # overwrites. input$km_tabset is the tabsetPanel id set in the km_tabs renderUI.
  output$dl_km <- downloadHandler(
    filename = function() sprintf("km_%s_%s_%s_%s%s.pdf", ct(), .feature_tag(input$feature),
                                  input$km_tabset, input$endpoint, .tau_tag(res())),
    content  = function(file) km_plot(res(), input$km_tabset, file = file)
  )

  # --- Multiple query -------------------------------------------------------------

  observeEvent(input$mq_kind, {
    req(input$mq_kind)
    updateSelectizeInput(session, "mq_feature", choices = FEATURES_ALL[[input$mq_kind]],
                         server = TRUE)
  }, ignoreNULL = FALSE)

  # One entry per cancer type, and every entry is a COMPLETE account of what happened to
  # that panel -- including the failures. A panel that quietly vanished would read as "no
  # effect" when it actually means "not measured" or "the call raised". Four statuses:
  #   ok     - a survresult with at least one estimable cohort
  #   absent - structural: the tissue has no cohort for this endpoint, or does not measure
  #            this feature at all. Stated, never silently dropped.
  #   empty  - measured, but no cohort cleared min_events
  #   error  - the call raised. The condition message is shown verbatim on the panel
  #            rather than swallowed into a generic "no data", which is the exact
  #            silent-acceptance pattern every shipped bug in this project came from.
  # The tryCatch is per-panel ON PURPOSE: one tissue raising must not blank the other four,
  # but it must be visible on its own panel.
  mq_res <- eventReactive(input$mq_run, {
    # Same rule as res() above, and this tab is where it was WORST (step 57). The panels
    # area is a renderUI, so the cancelled render did not merely leave the pre-run line
    # standing -- it REMOVED it, and the tab measured 0px of content under the controls.
    # The ZIP button appeared anyway, because it was gated on the run counter alone, and
    # its href served HTTP 500 with an empty <p> for a message. See output$mq_dl_ui for
    # the other half of that fix.
    req(input$mq_endpoint, input$mq_kind)
    validate(need(nzchar(input$mq_feature),
                  "No gene or TF selected - pick a symbol above, then press Run all panels."))
    feat <- input$mq_feature; kind <- input$mq_kind; ep <- input$mq_endpoint
    strat <- isTRUE(input$mq_adjust_strata)
    withProgress(message = sprintf("%s / %s over %d panels", feat, ep, length(CANCER_TYPES)),
                 value = 0, {
      setNames(lapply(CANCER_TYPES, function(ct) {
        incProgress(1 / length(CANCER_TYPES), detail = cancer_label_of(ct))
        co  <- cohorts_for(ct, ep)
        tau <- horizon_for(ct, ep)   # NULL => full follow-up, disclosed in the panel note
        # feature/kind ride along so a panel entry is a COMPLETE account on its own: the
        # export builds filenames and the CSV's feature column from the list, not from the
        # inputs, so a download cannot be labelled with a feature the run did not use.
        # Both forms of the feature ride along so a panel entry is a COMPLETE account on
        # its own: `feature` is the labelled form the result carries into the CSV,
        # `feature_id` the bare symbol that can go in a filename. The export reads them
        # from here rather than from the live inputs, which can be changed without
        # pressing Run -- a name built from input$mq_feature would label the file with a
        # query the zip does not contain.
        e <- list(ct = ct, feature = sprintf("%s (%s)", feat, KIND_LABEL[[kind]]),
                  feature_id = feat, kind = kind, endpoint = ep, tau = tau,
                  n_sel = length(co), status = "ok", msg = NULL, res = NULL)
        if (!length(co)) {
          e$status <- "absent"
          e$msg <- sprintf("No cohort in this panel carries %s.", ep)
          return(e)
        }
        if (!(feat %in% FEATURES_BY_CT[[ct]][[kind]])) {
          e$status <- "absent"
          e$msg <- sprintf("'%s' is not measured in this panel's %s data (%d %s cohort%s).",
                           feat, KIND_LABEL[[kind]], length(co), ep,
                           if (length(co) == 1) "" else "s")
          return(e)
        }
        r <- tryCatch({
          score <- get_feature(feat, cohorts = co, kind = kind)
          # feat IS in this tissue's vocabulary (checked immediately above), so an empty
          # score here is an inconsistency between the vocabulary and the matrices, not a
          # legitimate "no data" outcome. Raise it onto the panel.
          if (!length(score))
            stop(sprintf("'%s' is in the %s vocabulary for %s but returned no values",
                         feat, kind, ct))
          with_feature(get_survival(score, endpoint = ep, cohorts = co,
                                    adjust_strata = strat, max_followup = tau),
                       sprintf("%s (%s)", feat, KIND_LABEL[[kind]]))
        }, error = function(err) err)
        if (inherits(r, "error")) {
          e$status <- "error"; e$msg <- conditionMessage(r); return(e)
        }
        e$res <- r
        if (.n_estimable(r) < 1) {
          e$status <- "empty"
          e$msg <- sprintf("No %s cohort here had enough events for this query.", ep)
        } else {
          # Where this panel's TF sits in THIS TISSUE's genome-wide scan. Resolved here and
          # carried on the entry rather than looked up at render time, so the sentence that
          # reaches the screen and the one the ZIP writes into panel_note are the same
          # object -- the step-29 rule. mq_panel_notes() only formats it.
          #
          # The recipe matches the scan's by construction on this tab: co IS
          # cohorts_for(ct, ep) and tau IS horizon_for(ct, ep), the two things the Single
          # query tab has to check because the user can change them. What can still differ
          # is stratification, so use the EFFECTIVE value the engine used, not the raw
          # checkbox: mq_adjust_strata defaults TRUE but only applies where the registry
          # declares a stratifier -- breast (pam50) and, since 2026-08-27, lgg (idh) -- so
          # the tissues with none must resolve to their marginal scan. The condition below
          # is what enforces that; this comment names no fixed tissue count on purpose.
          e$rank <- scan_rank_lookup(
            feature = feat, kind = kind, cohorts = co, max_followup = tau,
            adjust_strata = strat && length(strata_vars_of(ct, ep)) > 0,
            strat_var = strata_var_for_scan(ct, ep),
            endpoint = ep, cancer_type = ct, scans = SCANS_BY_CT[[ct]],
            expected_cohorts = co, expected_horizon = tau)
        }
        e
      }), CANCER_TYPES)
    })
  })

  # ONE ruler for the whole run (step 55). A reactive rather than a call inside each
  # panel's renderPlot: five panels asking the same pure function the same question five
  # times is not just waste, it is five chances for a future edit to make one of them ask
  # a different question -- and a panel on its own axis is invisible next to four that
  # share one. Shiny caches this per run; every panel reads the one value.
  #
  # NULL below two estimable panels, which is exactly what forest_plot() takes to mean
  # "your own axis", so the one-panel run keeps the figure it has always had.
  mq_xrange <- reactive(mq_forest_xrange(mq_res()))

  # --- cross-panel export ---
  # TWO gates, the same pair output$dl_forest_ui carries on the Single query tab, and for
  # the reasons recorded there:
  #   1. the RUN COUNTER, because before the first click an eventReactive raises a silent
  #      condition that CANCELS the whole renderUI rather than returning from it, so a
  #      tryCatch alone cannot answer the pre-run question (the same trap that once made
  #      mq_panels come up blank instead of showing its pre-run line);
  #   2. tryCatch(mq_res()), because a run that REFUSED must not leave a button behind.
  # This was the last export bar in the app still on the counter alone (step 57). Pressing
  # Run with an emptied feature selector produced a button whose href served HTTP 500 and
  # an error page with an EMPTY message -- measured, not inferred -- because the filename
  # function calls mq_res() and mq_res() had refused. The reason belongs in the panels
  # area below, where it is now printed; a bar has no room for prose.
  output$mq_dl_ui <- renderUI({
    if (is.null(input$mq_run) || input$mq_run == 0) return(NULL)
    if (is.null(tryCatch(mq_res(), error = function(e) NULL))) return(NULL)
    downloadButton("mq_dl_all", "All panels (ZIP)", class = "btn-export")
  })

  # One zip: the long-format table of EVERY panel (including the ones that produced
  # nothing) plus a forest per estimated panel and a KM per estimated cohort.
  #
  # The filename and the contents come from mq_res(), never from the live inputs. The
  # selectors can be changed without pressing Run, so a name built from input$mq_feature
  # would label the file with a query the zip does not contain.
  output$mq_dl_all <- downloadHandler(
    filename = function() {
      R <- mq_res()
      # .feature_tag(), not the raw symbol: the ZIP'S OWN name has the same exposure its
      # entries do. Sanitising only inside mq_export_dir() would build a good archive that
      # then does not land -- measured 2026-08-10, a name over NAME_MAX is refused by the
      # browser with no file and no message, one step further out than the failure inside.
      sprintf("panels_%s_%s.zip", .feature_tag(R[[1]]$feature_id), R[[1]]$endpoint)
    },
    content = function(file) mq_export_zip(mq_res(), file)
  )

  output$mq_panels <- renderUI({
    # Gate on the button's own counter, NOT on tryCatch(mq_res()): before the first click
    # an eventReactive raises a silent condition that cancels the whole renderUI, so the
    # tryCatch branch never ran and the tab came up blank instead of showing this line.
    if (is.null(input$mq_run) || input$mq_run == 0)
      return(tags$p(style = "color:#888;",
        sprintf("Pick a feature and press 'Run all panels' - one panel per cancer type (%s).",
                paste(cancer_label_of(CANCER_TYPES), collapse = ", "))))
    R <- mq_res()
    cards <- lapply(CANCER_TYPES, function(ct) {
      e <- R[[ct]]
      # Same height rule as the single-query tab and the PDF export: from the number of
      # cohorts actually estimated in THIS panel, so a k=13 ovarian forest and a k=3
      # breast one are both drawn at their own proportions instead of a shared literal.
      ok <- identical(e$status, "ok")
      k  <- if (ok) max(1, .n_estimable(e$res)) else 1
      div(class = "mq-card",
        div(NULL,
            h4(cancer_label_of(ct), class = "mq-card-head"),
            # A panel with no forest reserves no forest-sized space: "auto" lets the one-line
            # reason sit under the heading instead of floating at the top of an empty box
            # the size of a plot that was never drawn.
            # Capped, not stretched to the card -- see MQ_FOREST_W_PX in R/plots.R.
            div(class = "mq-forest",
                plotOutput(paste0("mq_forest_", ct),
                           height = if (ok) sprintf("%dpx", mq_forest_h_px(k)) else "auto")),
            uiOutput(paste0("mq_note_", ct)),
            # The KM curves and the per-cohort table now live behind a per-tissue
            # toggle, CLOSED by default. Measured at the five panels of the day, forest +
            # KM tabset + DataTable came to 41,100px of page; the forests are a small
            # fraction of that, so collapsing the detail is what actually removes the
            # clutter, and it leaves the page reading as one forest per tissue down the
            # column -- which is the comparison the tab exists to make. At the seven panels
            # of 2026-08-26 the collapsed page measures 6,416px (BUILD_LOG step 81), so the
            # argument got stronger with the tissue count, not weaker.
            #
            # It is a server-side toggle, NOT a hidden div, and that is load-bearing.
            # A plotOutput created inside a display:none container does not render: the
            # RPPA KM curves sat in `recalculating` for 28s with R at 0.0% CPU when
            # suspended, and drew into a 0-pixel device when unsuspended (step 50). The
            # closed state here removes the outputs FROM THE DOM instead, which is the
            # same reason the table is omitted rather than req()-guarded below.
            if (ok) actionLink(paste0("mq_more_", ct), MQ_DETAIL_SHOW, class = "mq-more"),
            uiOutput(paste0("mq_detail_", ct))))
    })
    # One card per row, full width. No grid: see the note on MQ_FOREST_W_PX in R/plots.R
    # for why the tissues cannot be laid out two-up. The count is read from the registry
    # and has gone 5 -> 7 without touching this line; the argument never depended on it.
    do.call(tagList, cards)
  })
  # This tab is hidden at page load, and Shiny SUSPENDS outputs that are hidden when they
  # are created -- so mq_panels never computed its pre-run line and the tab opened as a
  # blank page under the controls, which reads as broken rather than as "not run yet".
  outputOptions(output, "mq_panels", suspendWhenHidden = FALSE)

  # Registered ONCE per cancer type, at server start -- not inside the renderUI above.
  # Building output slots inside a renderUI re-creates their observers on every re-render
  # and leaks them (the same reason the KM slots are registered in a loop up there).
  for (mct in CANCER_TYPES) {
    local({
      ct_ <- mct

      output[[paste0("mq_forest_", ct_)]] <- renderPlot({
        e <- mq_res()[[ct_]]
        validate(need(identical(e$status, "ok"),
                      if (identical(e$status, "error"))
                        "Query failed in this panel - see the message below."
                      else e$msg))
        forest_plot(e$res, xrange = mq_xrange())
      })

      # The SENTENCES come from mq_panel_notes() in R/plots.R; what is left here is styling.
      # That split is deliberate: the ZIP export pastes the same sentences into the CSV's
      # panel_note column, and a panel note written twice would be free to drift from the
      # one the reader saw on screen. Severity travels as the vector's names.
      output[[paste0("mq_note_", ct_)]] <- renderUI({
        e <- mq_res()[[ct_]]
        # mq_panel_notes() names each sentence by severity; this maps that name to the
        # shared note classes in PAGE_CSS. The SAME three the Single query tab's scan_rank
        # badge resolves to (note-good / note-warn / neutral), which is the point of naming
        # them once: a rank that clears FDR cannot read green on one tab and orange on the
        # other. An unrecognised severity falls through to the neutral tone rather than
        # erroring -- a sentence with no colour is still readable, a blank panel is not.
        cls <- function(sev) paste("note note-stack",
                                   switch(sev, warn = "note-warn", hit = "note-good",
                                          "note-body"))

        if (identical(e$status, "error"))
          return(tags$p(class = "note note-stack note-err",
                        tags$b("Query failed in this panel: "), e$msg))
        # absent/empty: the reason is already printed in the plot slot, exactly where the
        # forest would have been, which is where a reader looking for the missing panel
        # looks. Repeating it here was two copies of one sentence.
        if (!identical(e$status, "ok")) return(NULL)

        n <- mq_panel_notes(e)
        do.call(tagList, lapply(seq_along(n), function(i)
          tags$p(class = cls(names(n)[i]), n[[i]])))
      })

      # The detail slot. Returns NULL when closed, so the KM outputs and the DataTable are
      # ABSENT FROM THE DOM rather than hidden -- a plotOutput inside a display:none
      # container does not render (step 50's RPPA KM bug), and a DT that is merely
      # req()-cancelled keeps the PREVIOUS query's rows on screen.
      #
      # Both this and the label below read mq_detail_open() on the same counter, so they
      # cannot disagree about whether the panel is open.
      output[[paste0("mq_detail_", ct_)]] <- renderUI({
        if (!mq_detail_open(input[[paste0("mq_more_", ct_)]])) return(NULL)
        e <- mq_res()[[ct_]]
        if (!identical(e$status, "ok")) return(NULL)
        div(class = "mq-detail",
            # One tab per cohort THIS RUN actually estimated. The tab list is built here
            # rather than in its own renderUI because it depends on the result the layout
            # already holds; the plot slots themselves are registered once at server start
            # (below), never inside a renderUI, or their observers would be re-created on
            # every run and leak.
            #
            # Nothing is refit -- km_plot() draws from r$km, the survfit get_survival()
            # already computed. With the detail closed by default, the page draws ZERO
            # curves until a reader asks for one, where it used to draw one set per panel.
            HTML(figure_help_html("km")),
            do.call(tabsetPanel, c(
              id = paste0("mq_km_tabset_", ct_),
              lapply(names(Filter(function(x) isFALSE(x$skipped), e$res$per_cohort)),
                     function(co) tabPanel(co, plotOutput(paste0("mq_km_", co),
                       height = sprintf("%dpx", round(km_size_in()[["height"]] * 96))))))),
            HTML(figure_help_html("hr")),
            DTOutput(paste0("mq_table_", ct_)))
      })

      # Label follows the same counter. A re-rendered actionLink reports 0, which
      # mq_detail_open() reads as closed -- the same answer the fresh element already shows,
      # so a new run cannot leave "Hide ..." sitting above nothing.
      observeEvent(input[[paste0("mq_more_", ct_)]], {
        updateActionButton(session, paste0("mq_more_", ct_),
          label = if (mq_detail_open(input[[paste0("mq_more_", ct_)]])) MQ_DETAIL_HIDE
                  else MQ_DETAIL_SHOW)
      }, ignoreInit = TRUE)

      output[[paste0("mq_table_", ct_)]] <- renderDT({
        e <- mq_res()[[ct_]]
        req(identical(e$status, "ok"))
        d <- survtable_display(e$res)
        datatable(d, rownames = FALSE, caption = attr(d, "caption"),
                  options = c(list(pageLength = 5, dom = "tp", scrollX = TRUE),
                              survtable_dt_options(d)))
      })

      # One KM slot per cohort in this panel. The ids are mq_km_<cohort>: a cohort belongs
      # to exactly one cancer type so the name is unique, and the mq_ prefix keeps these
      # off the single-query tab's km_<cohort> slots -- an output id binds to exactly ONE
      # slot, so sharing the ids would have made the two tabs fight over one plot.
      for (mco in cohorts_of(ct_)) {
        local({
          co_ <- mco
          output[[paste0("mq_km_", co_)]] <- renderPlot({
            e <- mq_res()[[ct_]]
            # validate(), never req(): a cancelled render KEEPS THE PREVIOUS QUERY'S PLOT
            # on screen, which is exactly how a stale table once sat under a "not measured
            # here" message. validate() replaces the plot with the reason instead.
            validate(need(identical(e$status, "ok"), "Not estimated in this panel."))
            r <- e$res
            validate(need(!is.null(r$per_cohort[[co_]]) && isFALSE(r$per_cohort[[co_]]$skipped),
                          .cohort_note(co_, r)))
            km_plot(r, co_)
          })
        })
      }
    })
  }

  # --- Multiple genes ---------------------------------------------------------------
  #
  # WIRING ONLY. Every judgement this tab makes lives in R/plots.R and R/scan_lookup.R and
  # is tested there (test_multigene_input / _table / _display): which symbols resolve, what
  # is said once above the table versus per gene, which columns the table keeps. What is
  # left here is the reactive plumbing, which is the part that cannot leave a Shiny session
  # -- the same split set in Phase 7c for .tau_arg/.tau_tag.

  mg_ct <- reactive({ req(input$mg_ct); input$mg_ct })

  # Tissue-driven controls, the same three the Single query tab renders and for the same
  # reasons: only the score kinds this tissue has data for, only the endpoints its cohorts
  # carry, and a stratify checkbox only where a stratifier exists.
  output$mg_kind_ui <- renderUI({
    k <- kinds_of(mg_ct())
    radioButtons("mg_kind", "Score type",
                 choiceNames  = lapply(kind_choice_html(k, KIND_CHOICE_LABEL), HTML),
                 choiceValues = k,
                 selected = if ("viper" %in% k) "viper" else k[1])
  })
  output$mg_endpoint_ui <- renderUI({
    selectInput("mg_endpoint", "Endpoint", endpoints_of(mg_ct()),
                selected = default_endpoint_of(mg_ct()))
  })
  output$mg_strata_ui <- renderUI({
    sv <- strata_vars_of(mg_ct(), input$mg_endpoint)
    if (!length(sv)) return(tags$p(class = "note note-muted",
      sprintf("No stratifier for %s on %s - models are unstratified.",
              cancer_label_of(mg_ct()), input$mg_endpoint %||% "this endpoint")))
    # Derived from the registry -- see the note above output$strata_ui's copy of this line.
    tagList(
      checkboxInput("mg_adjust_strata", strata_label_of(mg_ct(), input$mg_endpoint), value = TRUE),
      { n <- strata_coverage_note(mg_ct(), input$mg_endpoint)
        if (is.null(n)) NULL else tags$p(class = "note note-muted note-tuck", n) },
      { n <- idh_stratum_note(mg_ct(), input$mg_endpoint)
        if (is.null(n)) NULL else tags$p(class = "note note-warn note-tuck",
                                         HTML(info_note_html(n, IDH_STRATUM_INFO))) })
  })

  # One entry per RESOLVED gene, every entry a complete account of what happened to it --
  # the same four statuses and the same per-gene tryCatch as mq_res(), for the same reason:
  # one gene raising must not blank the other nine, but it must be visible on its own row.
  # The resolver's output rides along beside the entries because the symbols that produced
  # NO entry (unrecognised, over the cap) are exactly the ones with no row to appear on.
  mg_res <- eventReactive(input$mg_run, {
    req(input$mg_ct, input$mg_kind, input$mg_endpoint)
    ct_ <- input$mg_ct; kind <- input$mg_kind; ep <- input$mg_endpoint
    # EFFECTIVE stratification, not the raw checkbox: Shiny keeps a removed input's last
    # value, so switching from breast (checkbox on) to ovarian (no checkbox) leaves
    # input$mg_adjust_strata TRUE and stale. The engine ignores it for a tissue with no
    # stratifier, so the scan lookup must too -- see scan_rank's note on the Single tab.
    strat <- isTRUE(input$mg_adjust_strata) && length(strata_vars_of(ct_, input$mg_endpoint)) > 0
    rs  <- mgene_resolve_symbols(input$mg_genes, FEATURES_BY_CT[[ct_]][[kind]])
    co  <- cohorts_for(ct_, ep)
    tau <- horizon_for(ct_, ep)          # NULL => full follow-up, disclosed in the header
    entries <- withProgress(
      message = sprintf("%d gene%s, %s %s", length(rs$ok),
                        if (length(rs$ok) == 1) "" else "s", cancer_label_of(ct_), ep),
      value = 0, {
        lapply(rs$ok, function(g) {
          incProgress(1 / max(1L, length(rs$ok)), detail = g)
          # depmap is attached HERE, at construction, not beside e$rank below: it is a
          # property of the symbol, so it holds for the "absent", "empty" and "error"
          # entries too -- each of which returns early. e$rank is genuinely conditional
          # (a rank describes a fit), this is not.
          e <- list(ct = ct_, feature = sprintf("%s (%s)", g, KIND_LABEL[[kind]]),
                    feature_id = g, kind = kind, endpoint = ep, tau = tau,
                    n_sel = length(co), status = "ok", msg = NULL, rank = NULL, res = NULL,
                    depmap = depmap_note(g, kind, DEPMAP_LISTS))
          if (!length(co)) {
            e$status <- "absent"
            e$msg <- sprintf("No cohort of this tissue carries %s.", ep)
            return(e)
          }
          r <- tryCatch({
            score <- get_feature(g, cohorts = co, kind = kind)
            if (!length(score)) NULL else
              with_feature(get_survival(score, endpoint = ep, cohorts = co,
                                        adjust_strata = strat, max_followup = tau),
                           e$feature)
          }, error = function(err) err)
          if (inherits(r, "error")) {
            e$status <- "error"; e$msg <- conditionMessage(r); return(e)
          }
          # An empty score is a LEGITIMATE outcome here, not an inconsistency, and the
          # message says which fact it is. FEATURES_BY_CT is built over cohorts_for(ct) --
          # every cohort of the tissue -- while the query runs over cohorts_for(ct, ep).
          # So a gene carried only by an endpoint-inert cohort (SCANB is OS-only) is
          # genuinely in this tissue's vocabulary and genuinely absent from this endpoint's
          # matrices. The Multiple query tab raises "in the vocabulary but returned no
          # values" for the same case, which is true but names the wrong cause.
          if (is.null(r)) {
            e$status <- "absent"
            e$msg <- sprintf(
              "'%s' is in this tissue's %s feature list but is not measured in any of the %d cohort%s that carry %s.",
              g, KIND_LABEL[[kind]], length(co), if (length(co) == 1) "" else "s", ep)
            return(e)
          }
          e$res <- r
          if (.n_estimable(r) < 1) {
            e$status <- "empty"
            e$msg <- sprintf("No %s cohort here had enough events for '%s'.", ep, g)
          } else {
            # The recipe matches the scan's BY CONSTRUCTION on this tab: co IS
            # cohorts_for(ct, ep) and tau IS horizon_for(ct, ep), because this tab offers no
            # control over either. That is the whole reason it offers none.
            e$rank <- scan_rank_lookup(
              feature = g, kind = kind, cohorts = co, max_followup = tau,
              adjust_strata = strat, strat_var = strata_var_for_scan(ct_, ep),
              endpoint = ep, cancer_type = ct_,
              scans = SCANS_BY_CT[[ct_]], expected_cohorts = co, expected_horizon = tau)
          }
          e
        })
      })
    list(resolved = rs, entries = entries)
  })

  # The header: what is true of the whole query, said once. Gated on the button's own
  # counter, NOT on tryCatch(mg_res()) -- before the first click an eventReactive raises a
  # silent condition that cancels the whole renderUI, which is how the mq tab once came up
  # blank instead of showing its pre-run line.
  output$mg_notes <- renderUI({
    if (is.null(input$mg_run) || input$mg_run == 0)
      return(tags$p(class = "note note-muted",
        sprintf("Type up to %d gene or TF symbols and press 'Run genes'.", MGENE_MAX)))
    R <- mg_res()
    # Zero entries is not an error state -- a typo'd list, or mRNA symbols pasted against a
    # VIPER vocabulary, resolves to nothing. It is the case where silence is worst, so the
    # resolver's own sentences are printed on their own. mgene_header_notes() cannot run
    # here: with no entry there is no tissue, endpoint or horizon to describe.
    n <- if (length(R$entries)) mgene_header_notes(R$entries, R$resolved)
         else c(warn = "No gene was analysed.", mgene_resolved_notes(R$resolved))
    cls <- function(sev) paste("note note-stack",
                               switch(sev, warn = "note-warn", hit = "note-good",
                                      "note-body"))
    do.call(tagList, lapply(seq_along(n), function(i)
      tags$p(class = cls(names(n)[i]), n[[i]])))
  })

  # The table EXISTS ONLY WHEN THERE IS ONE. Guarding renderDT with req() instead is not
  # enough: req() cancels the render, and a DT that is not re-rendered KEEPS THE PREVIOUS
  # QUERY'S TABLE on screen -- which is how a "not measured here" panel once sat directly
  # above a live table of another feature's numbers. Omitting the output from the layout
  # removes it from the DOM instead, and a DOM node that is gone cannot go stale.
  output$mg_table_slot <- renderUI({
    if (is.null(input$mg_run) || input$mg_run == 0) return(NULL)
    if (!length(mg_res()$entries)) return(NULL)
    DTOutput("mg_table")
  })

  output$mg_table <- renderDT({
    R <- mg_res()
    d <- mgene_display(R$entries)
    # dom = "t": no pager and no search box. The page length is the cap, so every row the
    # user asked for is on screen at once -- a paged table would hide genes behind a
    # control, which is the same silent shortfall the resolver's total accounting prevents
    # one layer up.
    #
    # ordering = FALSE, and it is NOT cosmetic. HR / p / I2 are FORMATTED STRINGS (they carry
    # "" for a gene that produced nothing, which no numeric column can), so DataTables sorts
    # them LEXICALLY. Observed on this table: sorting ascending by p put FOXM1 at p=4.48e-07
    # -- the strongest gene in the query, 27th of 1501 in the genome-wide scan -- BELOW JUN at
    # p=0.844, because the string "4.48e-07" sorts after "0.8". A reader clicking `p` to find
    # the top gene gets the reverse of the answer, with nothing on screen to contradict it.
    #
    # Disabling it rather than sorting numerically is also the honest design here: the row
    # order is the comparison the user typed, and ranking ten hand-picked genes by nominal p
    # is exactly the reading this tab's multiplicity note exists to discourage. The ordering
    # that IS meaningful -- each gene's rank among ~1500 in the genome-wide scan -- is in the
    # note column already.
    #
    # NOTE: survtable_display() formats HR / p / log-rank p / PH p the same way, and the
    # Single query and Multiple query tabs leave their tables SORTABLE -- so they got the
    # other answer, `orderData` against hidden rank columns (survtable_dt_options(), step
    # 40, 2026-08-11). The two tabs differ in what a ROW is, which is why they differ here:
    # there a row is a cohort and the ordering among cohorts is a real question, so the fix
    # is to make sorting correct; here a row is a gene the user typed and the order is the
    # comparison itself, so the fix is to not offer it.
    #
    # selection = "none": DT selects rows on click by default, and nothing here consumes a
    # selection. A row that highlights and then does nothing implies an affordance the tab
    # does not have -- the same objection the stratify checkbox's note makes to offering a
    # control for a tissue that ignores it. The DRAWER is not that affordance: it is opened
    # from one cell (td.mg-control), not from the row, so the rest of the row stays inert.
    #
    # EVERY ARGUMENT BELOW IS DERIVED FROM `d`, none of them written out here (2026-08-16,
    # step 45). options / escape / callback are all sets of COLUMN INDICES, and app.R has no
    # test harness -- an index written as a literal is wrong silently, the table still draws,
    # and the arrow expands whichever cell the off-by-one landed on. The three helpers live
    # in R/plots.R next to the function that builds the columns, where the tests can reach
    # them, which is the rule survtable_dt_options() was extracted under in Phase 7c.
    #
    # escape is a whitelist of two columns and MUST NOT become FALSE: the `gene` column holds
    # h5 feature strings and this tab's input is symbols the user pasted. See mgene_dt_escape().
    datatable(d, rownames = FALSE, selection = "none",
              escape   = mgene_dt_escape(d),
              callback = JS(mgene_dt_callback(d)),
              options  = mgene_dt_options(d))
  })

  # --- export ---
  #
  # The button EXISTS ONLY WHEN THERE IS SOMETHING TO EXPORT, and the second gate is not
  # cosmetic: mgene_survtable() stops on an empty entry list (by design -- a CSV of no genes
  # is not a thing), so a button shown after a run that resolved nothing would hand the
  # browser a failed download with no message anywhere. The zero-resolved case is the common
  # one it guards -- a typo'd list, or mRNA symbols pasted against a VIPER vocabulary -- and
  # it is exactly the case where mg_notes is already saying "No gene was analysed." Same two
  # gates as mg_table_slot, for the same reason.
  output$mg_dl_ui <- renderUI({
    if (is.null(input$mg_run) || input$mg_run == 0) return(NULL)
    if (!length(mg_res()$entries)) return(NULL)
    downloadButton("mg_dl_table", "Table (CSV)", class = "btn-export")
  })

  # The long table: one row per gene PER COHORT plus a pooled row, which is the per-cohort
  # detail the on-screen table deliberately does not show. Both the name and the contents
  # come from mg_res(), never from the live inputs.
  output$mg_dl_table <- downloadHandler(
    filename = function() mgene_export_name(mg_res()$entries),
    content  = function(file)
      write.csv(mgene_survtable(mg_res()$entries, mg_res()$resolved), file, row.names = FALSE)
  )

  # EVERY renderUI on this tab, in one place, AFTER they are all defined.
  #
  # Shiny suspends an output that is HIDDEN WHEN IT IS CREATED, and every tab but the first is
  # hidden at page load. This tab is the first in the app to be both a non-default tab and
  # built from renderUI -- the Single query tab renders the same three controls this way but
  # is the one visible at load, and the Multiple query tab is hidden but uses static inputs --
  # so the trap is new here, and it cost the tab three of its five controls on the first
  # render: no score type, no endpoint, no stratify checkbox, just a tissue selector and a
  # Run button that could not run. Without input$mg_endpoint, req() in mg_res cancels and the
  # button does nothing AT ALL, with no error anywhere.
  #
  # Found by LOOKING at the rendered page, with the whole suite green -- nothing in R/ can see
  # a suspended output and app.R has no harness. tests/test_multigene_wiring.R derives this
  # list from the source and fails if any renderUI here is left off it.
  for (o in c("mg_kind_ui", "mg_endpoint_ui", "mg_strata_ui", "mg_notes", "mg_table_slot",
              "mg_dl_ui"))
    outputOptions(output, o, suspendWhenHidden = FALSE)

  # --- browse ---
  output$browse_view_ui <- renderUI({
    ch <- scan_choices_of(ct())
    if (!length(ch)) return(NULL)
    selectInput("browse_view", "Table", choices = ch, width = "100%")
  })
  output$browse_empty <- renderUI({
    if (length(scan_choices_of(ct()))) return(NULL)
    tags$p(style = "color:#a05a00;",
      sprintf("No scan results yet for %s. Run: Rscript R/discovery_validation.R %s",
              ct(), ct()))
  })
  output$browse_table <- renderDT({
    req(input$browse_view)
    df <- .scans()[[input$browse_view]]
    req(!is.null(df))
    ord <- if ("q" %in% names(df)) list(list(match("q", names(df)) - 1, "asc")) else list()
    datatable(df, filter = "top", rownames = FALSE, selection = "single",
              options = list(pageLength = 25, scrollX = TRUE, order = ord))
  })

  # ---- Protein (RPPA) tab -----------------------------------------------------------
  #
  # Self-contained on purpose: nothing here calls cohorts_for(), forest_plot(), or
  # survtable(), and nothing on the other tabs reads these reactives. That isolation IS
  # the guarantee that an antibody row can never turn up as a forest row.

  rp_ct <- reactive({ req(input$rp_ct); input$rp_ct })

  output$rp_gene_ui <- renderUI({
    co <- rppa_cohort_for(rp_ct())
    g  <- rppa_genes(co)
    # Empty is a real state (a tissue staged with no map yet), not an impossible one, so
    # it gets a disabled control naming the cause rather than an empty dropdown.
    if (!length(g))
      return(selectInput("rp_gene", "Gene", choices = character(0)))
    selectizeInput("rp_gene", sprintf("Gene (%d on the panel)", length(g)),
                   choices = g, selected = if ("AKT1" %in% g) "AKT1" else g[1])
  })

  rppa_res <- eventReactive(input$rp_run, {
    # Same rule as res() and mq_res() (step 57). TWO ways to reach the empty gene here:
    # clearing the selectize, and a tissue whose RPPA panel carries no genes at all, which
    # output$rp_gene_ui renders as choices = character(0). Both made Run panel a button
    # that did nothing, on a tab whose pre-run state is also empty -- so the click was
    # indistinguishable from not having clicked.
    req(input$rp_endpoint)
    validate(need(nzchar(input$rp_gene), "No gene selected - pick one from the panel above."))
    g <- input$rp_gene
    # A stale input$rp_gene survives a tissue switch for one reactive flush, so a gene
    # valid in breast can arrive here labelled ovarian. Checking membership rather than
    # letting rppa_panel() return "absent" keeps a transient from rendering as a finding.
    if (!(g %in% rppa_genes(rppa_cohort_for(rp_ct())))) return(NULL)
    rppa_panel(g, rp_ct(), input$rp_endpoint)
  })

  # Gene context, not a result -- the same contract as output$depmap_badge on the Single
  # query tab, and deliberately the same behaviour: keyed on the SELECTED gene, not on
  # rppa_res(), so it does not wait for a Run. Gating it on the panel would imply it
  # describes the panel, and it does not: it describes the gene, and it is equally true
  # before any model is fit. R/rppa.R owns the phospho caveat that rides with it.
  output$rppa_overlap_note <- renderUI({
    req(input$rp_ct, input$rp_endpoint)
    n <- rppa_overlap_note(input$rp_ct, input$rp_endpoint)
    tags$p(class = paste("note note-stack", paste0("note-", n$tone)),
           HTML(info_note_html(n$text, n$info)))
  })

  output$rppa_depmap_badge <- renderUI({
    req(input$rp_gene)
    n <- rppa_depmap_note(input$rp_gene, DEPMAP_LISTS)
    if (is.null(n)) return(NULL)
    tags$p(class = paste("note note-tuck", paste0("note-", n$tone)),
           HTML(info_note_html(n$text, n$info)))
  })

  output$rppa_notes <- renderUI({
    p <- rppa_res()
    if (is.null(p)) return(NULL)
    if (identical(p$status, "ok") && !isTRUE(p$n_varies)) return(NULL)
    tags$p(class = "note note-warn note-stack",
           p$msg %||% "This panel's rows do not share one patient set.")
  })

  output$rppa_panel_slot <- renderUI({
    # tryCatch, not a bare rppa_res() (step 57): a refusal message belongs to
    # output$rppa_notes, and printing it here as well would be the same sentence twice --
    # the second copy inside a box the size of a panel that was never drawn.
    p <- tryCatch(rppa_res(), error = function(e) NULL)
    if (is.null(p) || !identical(p$status, "ok")) return(NULL)
    plotOutput("rppa_panel", height = sprintf("%dpx",
      round(rppa_panel_size_in(length(p$rows))[["height"]] * 96)))
  })

  output$rppa_panel <- renderPlot({
    p <- rppa_res(); req(!is.null(p), identical(p$status, "ok"))
    rppa_panel_plot(p)
  }, res = 96)

  output$rp_dl_ui <- renderUI({
    if (is.null(input$rp_run) || input$rp_run == 0) return(NULL)
    # tryCatch for the reason recorded above output$dl_forest_ui: a bare rppa_res() prints
    # a validate() message INSIDE THE EXPORT BAR, once per button (step 57).
    p <- tryCatch(rppa_res(), error = function(e) NULL)
    if (is.null(p) || !identical(p$status, "ok")) return(NULL)
    # btn-export like every other download in the app (step 114). It was the one that was
    # not, so with the buttons now on the brand purple it would have been the only download
    # rendering as a primary action -- the exact hierarchy inversion the class exists to fix.
    downloadButton("dl_rppa", "Download panel (PDF)", class = "btn-export")
  })
  output$dl_rppa <- downloadHandler(
    filename = function() rppa_export_name(rppa_res()),
    content  = function(file) rppa_panel_plot(rppa_res(), file = file)
  )

  # One KM per ESTIMABLE antibody, STACKED rather than in a tabset.
  #
  # Two reasons, and the second is why it is not a workaround. First: a plotOutput that is
  # hidden at creation inside a renderUI-built tabsetPanel never renders. Suspended, it
  # stays in `recalculating` forever because switching the tab does not resume it (checked
  # on the page: 28s, R at 0.0% CPU); unsuspended, it renders eagerly into a 0-pixel-wide
  # container and km_plot()'s graphics::layout() dies with "invalid graphics state". Only
  # the first tab ever works. (The Single query tab's KM tabset has the same shape and the
  # same behaviour -- see BUILD_LOG; that one is NOT this step's to fix.)
  #
  # Second, and the reason this is the right shape anyway: the whole point of this tab is
  # COMPARING an antibody against the other forms of the same protein. Putting the curves
  # behind tabs hides exactly the comparison the panel above just set up. There are at
  # most four.
  .rppa_est <- reactive({
    # tryCatch, same reason as output$rppa_panel_slot (step 57): output$rppa_km_stack reads
    # this bare, so a refusal raised here printed the SAME SENTENCE a second time, once in
    # the notes and once where the curves would be. A refusal means no estimable rows.
    p <- tryCatch(rppa_res(), error = function(e) NULL)
    if (is.null(p)) return(list())
    head(Filter(function(r) isFALSE(r$skipped), p$rows), RPPA_MAX_ANTIBODIES)
  })

  output$rppa_km_stack <- renderUI({
    est <- .rppa_est()
    if (!length(est)) return(NULL)
    tagList(lapply(seq_along(est), function(i)
      div(class = "rppa-km",
          tags$h4(est[[i]]$antibody, class = "rppa-km-head"),
          plotOutput(paste0("rppa_km_", i),
                     height = sprintf("%dpx", round(km_size_in()[["height"]] * 96))))))
  })

  # Outputs are registered for every SLOT up front, not per render: a renderPlot created
  # inside another reactive is re-registered on each flush and Shiny warns, and the tab
  # count changes with the gene. RPPA_MAX_ANTIBODIES bounds it (the busiest gene on the
  # panel carries 4).
  for (.i in seq_len(RPPA_MAX_ANTIBODIES)) {
    local({
      i <- .i
      output[[paste0("rppa_km_", i)]] <- renderPlot({
        est <- .rppa_est(); req(length(est) >= i)
        km_plot(est[[i]]$res, est[[i]]$res$cohorts[1] %||% rppa_cohort_for(rp_ct()))
      }, res = 96)
    })
  }

  # UNSUSPEND. Same failure the mg_ loop above exists for, and the RPPA tab is its second
  # instance: an output on a tab that is HIDDEN AT PAGE LOAD is suspended by Shiny and
  # never computes until something forces a re-render, so the control or the plot simply
  # is not there. Verified by looking at the rendered page -- the second KM tab sat in
  # `recalculating` for 28s with the R process at 0.0% CPU, which is what a suspended
  # output looks like, not a slow one.
  #
  # renderUI OUTPUTS ONLY. Unsuspending a renderPlot here is actively harmful and was
  # tried first: an eager render happens while the container is still 0 pixels wide, and
  # km_plot()'s graphics::layout() dies with "invalid graphics state" -- the plot then
  # stays broken because nothing invalidates it when the tab is finally shown. A plot
  # must stay suspended so that its first render happens at its real size.
  for (o in c("rp_gene_ui", "rp_dl_ui", "rppa_overlap_note", "rppa_depmap_badge", "rppa_notes",
              "rppa_panel_slot", "rppa_km_stack"))
    outputOptions(output, o, suspendWhenHidden = FALSE)

  observeEvent(input$query_selected, {
    df <- .scans()[[input$browse_view]]
    sel <- input$browse_table_rows_selected
    req(length(sel) == 1)
    updateRadioButtons(session, "kind", selected = "viper")
    updateSelectizeInput(session, "feature", choices = FEATURES_BY_CT[[ct()]]$viper,
                         selected = df$tf[sel], server = TRUE)
    updateNavbarPage(session, "nav", selected = TAB_SINGLE)
  })
}

shinyApp(ui, server)
