# R/guide.R -- the Guide tab: what OMICohort is, what ARACNe -> VIPER buys, and how to read
# what comes back. Added 2026-08-30 (step 109).
#
# WHY THIS FILE EXISTS RATHER THAN A BLOCK OF HTML IN app.R. Two reasons, and the second is
# the one that will still matter in a year.
#
#   1. app.R has no test harness. The house rule is already written at R/plots.R's
#      KIND_INFO: "Each explanation is declared BESIDE the text it explains, never in app.R
#      -- an explanation that drifts from the sentence it annotates is exactly the failure a
#      reader cannot see." A whole page of explanation deserves the rule more, not less.
#
#   2. EVERY FIGURE ON THE PAGE IS DERIVED, NONE IS TYPED. The four markdown documents have
#      tests/test_doc_integrity.R holding them to the files; a Shiny tab is invisible to it.
#      A hand-written landing page would therefore be the one document in this project with
#      no staleness guard on it at all -- and it would be the first thing a new reader sees.
#      So guide_facts() reads the registry and the scans the app has already loaded, and
#      guide_html() may only interpolate what guide_facts() returned. Nothing here is
#      allowed to state a number the files do not currently say.
#
# The prose is the exception to that and is deliberately so: sentences about what a hazard
# ratio MEANS cannot be derived, they can only be written and read. What is guarded is that
# no sentence carries a stale count beside it.
#
# Returns character HTML, not shiny tags, for the reason R/plots.R gives at .html_escape():
# R/ is a plain function library and a web dependency has no business in it. app.R wraps the
# result in HTML() exactly once.


# ---- 1. the facts -------------------------------------------------------------------

# Pick the scan whose estimand the app actually reports for (tissue, endpoint). When a
# tissue declares a stratifier the ADJUSTED scan is that estimand -- breast reports the
# PAM50-adjusted scan, lgg the IDH-adjusted one -- and the marginal sibling is kept only for
# comparison. Derived from the names on disk rather than listed, so a tissue that gains or
# loses an adjusted arm moves this page with it.
.guide_scan_name <- function(ep, scans) {
  base <- paste0("mr_scan_", tolower(ep))
  cand <- grep(paste0("^", base), names(scans), value = TRUE)
  if (!length(cand)) return(NULL)
  adj <- setdiff(cand, base)
  if (length(adj)) adj[[1]] else base
}

# Same, for the discovery -> validation table, with one hard exclusion: `_inverted` is a
# SENSITIVITY arm (the split run the other way round) and is circular by construction, so it
# must never be read as a result. Naming it here rather than filtering it downstream keeps
# the exclusion where a reader of this file will see it.
.guide_discval_name <- function(ep, scans) {
  cand <- grep(paste0("^mr_discval_", tolower(ep)), names(scans), value = TRUE)
  cand <- cand[!grepl("_inverted$", cand)]
  if (length(cand)) cand[[1]] else NULL
}

# The FDR threshold, declared ONCE. It is both counted with (n_fdr) and printed in the
# caption under the table; two copies of it is how a page ends up counting at one threshold
# and claiming another.
GUIDE_FDR <- 0.05

.guide_col <- function(d, nm) if (!is.null(d) && nm %in% names(d)) d[[nm]] else NULL

# The endpoints a tissue was actually SCANNED at, read off the scan names -- not
# endpoints_of(), which returns every endpoint the registry declares. Those are not the same
# set and the difference is exactly the thing a reader would be misled by: breast and lgg
# both register DSS at k=1, below POOL_MIN_K, so it is declared and never scanned. Listing
# "OS, DFS, DSS" against every tissue would promise seven analyses that do not exist.
.guide_scanned_eps <- function(scans, primary) {
  ep <- sub("^mr_scan_([a-z]+).*$", "\\1", grep("^mr_scan_", names(scans), value = TRUE))
  ep <- toupper(unique(ep))
  c(intersect(primary, ep), sort(setdiff(ep, primary)))    # primary first, then the rest
}

# One row of the data table. Every field is NULL-safe: a tissue staged but not yet scanned
# must render as a tissue with no scan, never as an error and never as a silent absence.
.guide_tissue <- function(ct, scans, registry) {
  co   <- cohorts_for(ct, registry = registry)
  eps  <- endpoints_of(ct, registry = registry)
  prim <- primary_endpoint_of(ct)

  strat <- unique(unlist(lapply(co, function(x) stratifiers_of(x, registry = registry))))
  strat <- setdiff(strat, c("none", NA, ""))

  sn <- .guide_scan_name(prim, scans)
  s  <- if (!is.null(sn)) scans[[sn]] else NULL
  q  <- .guide_col(s, "q")
  k  <- .guide_col(s, "k"); n <- .guide_col(s, "n"); ev <- .guide_col(s, "events")

  dn <- .guide_discval_name(prim, scans)
  d  <- if (!is.null(dn)) scans[[dn]] else NULL

  list(ct = ct, label = cancer_label_of(ct), n_cohorts = length(co),
       endpoints = eps, scanned_eps = .guide_scanned_eps(scans, prim),
       primary = prim, strata = strat,
       # horizon_for() returns NULL for a tissue that declares no tau (breast runs full
       # follow-up). is.na(NULL) is logical(0) and `if (logical(0))` is an ERROR rather than
       # FALSE, so it is normalised HERE, once, instead of at each place that reads it.
       horizon = { h <- tryCatch(horizon_for(ct, prim), error = function(e) NULL)
                   if (!length(h)) NA_real_ else as.numeric(h[[1]]) },
       scan = sn, n_tf = if (is.null(s)) NA_integer_ else nrow(s),
       k = if (is.null(k)) NA_integer_ else stats::median(k, na.rm = TRUE),
       n = if (is.null(n)) NA_integer_ else max(n, na.rm = TRUE),
       events = if (is.null(ev)) NA_integer_ else max(ev, na.rm = TRUE),
       n_fdr = if (is.null(q)) NA_integer_ else sum(q < GUIDE_FDR, na.rm = TRUE),
       discovered = if (is.null(d)) NA_integer_ else sum(d$discovered, na.rm = TRUE),
       replicated = if (is.null(d)) NA_integer_ else sum(d$replicated, na.rm = TRUE))
}

# scans_by_ct  -- app.R's SCANS_BY_CT (already read at startup; this adds no I/O)
# features_by_ct -- app.R's FEATURES_BY_CT, for the TF vocabulary the selectors offer
guide_facts <- function(scans_by_ct, features_by_ct = NULL, registry = COHORTS) {
  cts <- sort(unique(registry$cancer_type))
  ti  <- lapply(cts, function(ct)
    .guide_tissue(ct, if (is.null(scans_by_ct[[ct]])) list() else scans_by_ct[[ct]], registry))
  names(ti) <- cts

  # Networks IN USE, resolved through the registry rather than by listing networks/: the
  # directory also holds gse9891, the second ovarian network that was evaluated and
  # deliberately NOT committed. Counting the directory would put it on the page as though
  # the tool queried it.
  # FAIL LOUD IF A DECLARED NETWORK IS NOT ON DISK. Two things used to conspire here, and
  # the result was a false sentence rendered with no error (found 2026-08-31, step 119, and
  # confirmed by moving one file aside: the page went from 8 to 7 in silence):
  #
  #   - regulon_paths_for() ends in `p[file.exists(p)]`, so a missing regulon is DROPPED, not
  #     reported. That is right for its other caller -- batch_regulators() needs the ones that
  #     are there -- and wrong here, where the length IS the claim.
  #   - this line then wrapped it in tryCatch(..., error = character(0)), which turns "I could
  #     not tell" into "zero" and hands it to .g_int() as a fact.
  #
  # A deployment that ships code and data but not networks/ (50 MB of regulons that the app
  # never reads -- VIPER scores are precomputed into viper_*.h5) would have published a landing
  # page reading "0 tissue-specific ARACNe networks covering 1,598 transcription factors".
  # The registry is passed through as well: guide_facts() takes a `registry` argument and this
  # call ignored it, so the count answered to the global COHORTS whatever it was given.
  want <- unique(registry$network_dir[nzchar(registry$network_dir)])
  nets <- unique(unlist(lapply(cts, regulon_paths_for, registry = registry)))
  if (length(nets) < length(want)) {
    have <- basename(dirname(nets))
    stop(sprintf(paste0(
      "the registry declares %d ARACNe network(s) but only %d regulon file(s) are on disk; ",
      "missing: %s.\nThe Guide states this count as fact, so it must not be quietly short. ",
      "Expected at networks/<dir>/regulon_<dir>.rds under %s.\nIf this is a deployment that ",
      "deliberately omits networks/, ship the regulons (about 50 MB) -- the count and the ",
      "candidate-regulator range on the landing page both come from that directory."),
      length(want), length(nets),
      paste(setdiff(want, have), collapse = ", "), ROOT))
  }

  tf <- if (is.null(features_by_ct)) NA_integer_ else
    length(unique(unlist(lapply(features_by_ct, function(f) f[["viper"]]))))

  # The candidate-regulator list each ARACNe run was seeded with, read off the network
  # directories the registry actually points at -- NOT by listing networks/, for the same
  # reason `nets` is not: gse9891 has one too and is deliberately not in use.
  # networks/tcga has no tf_list_*.txt of its own: that network was reused from the earlier
  # TCGA-BRCA work and its list lives outside this repository. So this range covers the
  # networks whose list is HERE, which is what a figure on a page should be able to prove.
  tfl <- lapply(unlist(lapply(unique(dirname(nets)), list.files,
                              pattern = "^tf_list_.*[.]txt$", full.names = TRUE)),
                function(f) readLines(f, warn = FALSE))
  reg <- if (length(tfl)) range(lengths(tfl)) else c(NA_integer_, NA_integer_)

  list(n_cohorts = nrow(registry), n_types = length(cts), tissues = ti,
       n_networks = length(nets), n_tf = tf,
       n_reg_min = reg[[1]], n_reg_max = reg[[2]],
       patients = sum(vapply(ti, function(t) if (is.na(t$n)) 0L else as.integer(t$n), integer(1))),
       events = sum(vapply(ti, function(t) if (is.na(t$events)) 0L else as.integer(t$events), integer(1))))
}


# ---- 2. the page --------------------------------------------------------------------

.g_int <- function(x) if (is.na(x)) "n/a" else formatC(as.integer(x), format = "d", big.mark = ",")
.g_esc <- function(x) .html_escape(x)
.g_join <- function(x) paste(.g_esc(x), collapse = ", ")

# THE ENDPOINT COLUMN CARRIES THE PRIMARY MARK, rather than a column of its own. Until step
# 142 this table had nine columns and two of them were about endpoints: one listing what was
# scanned and one naming which is primary. The primary endpoint is always one of the scanned
# ones, so the second column repeated a value already on the row and cost a ninth of the width
# to do it. Bolding it inside the list says the same thing in one column.
#
# tests/test_guide_derived.R section 3 reads this row for the presence of each SCANNED endpoint
# and the absence of each declared-but-unscanned one, with fixed = TRUE. Wrapping the primary
# in <b> does not disturb that: the substring is still there. Do not switch the join to
# something that puts markup BETWEEN the letters of an endpoint name.
.guide_eps <- function(t) {
  if (!length(t$scanned_eps)) return("none")
  paste(vapply(t$scanned_eps, function(e)
    if (identical(e, t$primary)) paste0("<b>", .g_esc(e), "</b>") else .g_esc(e),
    character(1)), collapse = ", ")
}

.guide_table <- function(f) {
  rows <- vapply(f$tissues, function(t) sprintf(
    paste0("<tr><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td>",
           "<td>%s</td><td>%s / %s</td><td>%s</td><td>%s</td></tr>"),
    .g_esc(t$label), .g_int(t$n_cohorts), .guide_eps(t),
    if (length(t$strata)) paste0("<code>", .g_join(t$strata), "</code>") else "none",
    if (is.na(t$horizon)) "full" else paste0(.g_int(t$horizon), " mo"),
    .g_int(t$k), .g_int(t$n),
    .g_int(t$n_fdr), if (is.na(t$replicated)) "n/a" else
      sprintf("%s of %s", .g_int(t$replicated), .g_int(t$discovered))), character(1))

  paste0(
    # Eight columns still do not fit a narrow window. Scrolling the TABLE rather than the page
    # is the rule the rest of the app follows for wide content.
    '<div class="guide-wrap"><table class="guide-tbl"><thead><tr>',
    '<th>Tissue</th><th>Cohorts</th><th>Endpoints</th><th>Stratified by</th>',
    '<th>Horizon &tau;</th><th>Cohorts pooled / patients</th>',
    # THE THRESHOLD IS INTERPOLATED, NEVER TYPED. Section 6 of tests/test_guide_derived.R
    # allows exactly one occurrence of the number in non-comment source: the GUIDE_FDR
    # declaration. A header that spelled it out would be a second copy, and a page that counts
    # at one threshold while its heading claims another is unreadable from the outside.
    '<th>Pass q &lt; ', format(GUIDE_FDR), '</th><th>Replicated</th>',
    '</tr></thead><tbody>', paste(rows, collapse = ""), '</tbody></table></div>',

    # WHAT THE COLUMNS MEAN, as a list rather than the single long sentence this used to be.
    # The old paragraph ran five clauses together separated by semicolons and defined only two
    # of the nine headings; "Clear FDR" in particular was a heading a reader could not decode
    # and a definition that assumed they already knew what a false discovery rate controls.
    # A reader meets this table before anything else on the page, so a column they cannot read
    # is a number they will either ignore or misread.
    '<p class="guide-fine">What the columns mean:</p>',
    '<ul class="guide-list guide-fine">',
    '<li><b>Cohorts</b>: how many separate patient series this tissue has in the tool. They ',
    'are always analysed one at a time and then combined, never merged into a single pool.</li>',
    '<li><b>Endpoints</b>: which survival endpoints were scanned across all regulators. The ',
    '<b>primary</b> one is in bold; it is the endpoint the figures in this row describe. An ',
    'endpoint appears here only if it was actually scanned, so a tissue whose third endpoint ',
    'never reached the minimum pool size does not list it.</li>',
    '<li><b>Stratified by</b>: a clinical variable each cohort&rsquo;s model was given its own ',
    'baseline risk for, so that patients are compared within a group rather than across ',
    'groups. &ldquo;none&rdquo; means no stratifier was declared for that tissue.</li>',
    '<li><b>Horizon &tau;</b>: follow-up beyond this many months is not used, which keeps a ',
    'few very long-followed patients from carrying a result. &ldquo;full&rdquo; means the ',
    'whole of the available follow-up is used.</li>',
    '<li><b>Cohorts pooled / patients</b>: how many cohorts and how many patients actually ',
    'entered the combined estimate for the primary endpoint. This can be fewer than ',
    '<b>Cohorts</b>, because a cohort joins only if it recorded that endpoint.</li>',
    '<li><b>Pass q &lt; ', format(GUIDE_FDR), '</b>: how many regulators cleared that ',
    'threshold in the scan for this row. The q value is a Benjamini-Hochberg false discovery ',
    'rate, which is a statement about the whole list rather than about any one regulator: of ',
    'the regulators called here, roughly ', format(100 * GUIDE_FDR), ' in every hundred are ',
    'expected to be there by chance. Where a stratifier is declared, the stratified scan is ',
    'the one counted, since that is the estimate the app reports.</li>',
    '<li><b>Replicated</b>: regulators found in one set of cohorts and then tested in a ',
    'separate set that took no part in finding them. It reads as <i>held up of found</i>. ',
    'The zeros are real results and not missing data: in those tissues nothing survived the ',
    'second look.</li>',
    '</ul>',
    '<p class="guide-fine">Every figure in this table is computed when the app starts, from ',
    'the same registry and result tables the analysis itself was run on. None of it is typed ',
    'into this page, so it cannot disagree with the analysis it describes. For what these ',
    'numbers do not license, see <i>What a hit does not mean</i>.</p>')
}

guide_html <- function(f) {
  paste0(
'<div class="guide">',

# The lockup, not the word. logo_full.png is mark + wordmark + "KORKMAZ LAB", trimmed to its
# own ink by scripts/make_logo.py, so this is the one place the whole thing is shown at a size
# it was drawn for. It stands IN the <h2> rather than above it: the image is the heading, so
# printing "OMICohort" as text beside it would say the name twice on the same line -- exactly
# the duplication the navbar mark was cropped to avoid (see scripts/make_logo.py).
'<div class="guide-hero">',
'<h2 class="guide-h1"><img src="logo_full.png" class="guide-logo" alt="OMICohort"></h2>',
'<p class="guide-lead">A survival-analysis tool for asking whether a molecular score predicts ',
'outcome across <b>', .g_int(f$n_cohorts), ' patient cohorts</b> in <b>', .g_int(f$n_types),
' cancer types</b>, covering <b>', .g_int(f$patients), ' patients</b> and <b>', .g_int(f$events),
' recorded events</b>. Each cohort is fitted on its own and the results combined by meta-analysis, never ',
'merged into one pile. Unlike the public KM plotters, the score you query does not have to be a ',
'gene&rsquo;s expression.</p>',
'</div>',

# ---------------------------------------------------------------- the distinguishing part
'<h3 class="guide-h2">What makes this different: regulon activity, not expression</h3>',

'<p>Most survival tools can answer one question: <i>does this gene&rsquo;s mRNA level predict ',
'outcome?</i> That question has a known weakness when the gene is a <b>transcription factor</b>. ',
'A TF does its work as a protein, and how much of that protein is active is set by things mRNA ',
'does not see, such as phosphorylation, cofactor availability, whether it is in the nucleus at all. ',
'A TF can be flat at the mRNA level and still be driving a tumour, and it can be highly ',
'expressed and doing nothing.</p>',

'<p>OMICohort asks a different question: <b>is this TF&rsquo;s regulatory program switched on?</b> ',
'That is measured from the behaviour of the genes the TF controls, in two steps.</p>',

'<div class="guide-steps">',
'<div class="guide-step"><span class="guide-num">1</span><div>',
'<b>ARACNe</b> builds the network. From tumour expression <i>within each tissue</i>, it infers ',
'which genes each TF regulates, using mutual information and removing indirect edges. The output ',
'is a <b>regulon</b> per TF: its inferred target set, each target carrying a confidence and a ',
'direction. This is run per tissue, not borrowed; a breast regulon and a glioma regulon for ',
'the same TF are different sets, because the TF does different work in each.',
' The candidate regulators are the <i>same list</i> for every tissue: the human transcription ',
'factors of the <a href="https://humantfs.ccbr.utoronto.ca/">Human Transcription Factors</a> ',
'database (v1.01; Lambert <i>et al.</i>, <i>Cell</i> 2018), intersected per network with the ',
'genes that survive expression filtering, which leaves <b>', .g_int(f$n_reg_min), ' to ',
.g_int(f$n_reg_max), '</b> regulators per network. Sharing the list is what makes a TF ',
'comparable across tissues: the regulons differ because the tissue does, not because the ',
'candidate set did.',
'</div></div>',
'<div class="guide-step"><span class="guide-num">2</span><div>',
'<b>VIPER</b> scores the activity. For one patient, it asks whether that TF&rsquo;s targets move ',
'together in the direction the regulon predicts; activated targets up, repressed targets ',
'down. Strong coordinated movement means the program is on. The result is a normalised enrichment ',
'score (NES) per TF <i>per patient</i>, which is what the survival model then uses.',
'</div></div>',
'</div>',

'<p class="guide-callout"><b>The consequence to keep in mind:</b> VIPER activity and the ',
'TF&rsquo;s own mRNA can point in opposite directions, and when they do, that is a finding rather ',
'than a contradiction. They are measuring different things. Both are queryable here (the ',
'<i>Score type</i> selector) so you can see when they disagree.</p>',

'<p>This is the part of the tool with no public equivalent: <b>', .g_int(f$n_networks),
' tissue-specific ARACNe networks</b> covering <b>', .g_int(f$n_tf), ' transcription factors</b>, ',
'each TF&rsquo;s activity tested against survival in every cohort that carries the endpoint, then ',
'pooled. <b>', .g_int(f$n_types), ' cancer types were scanned genome-wide across ',
'the full regulator list</b>, with a discovery set and a held-out validation set.</p>',

# ---------------------------------------------------------------- the data
'<h3 class="guide-h2">The data</h3>',
'<p>Public cohorts only, with clinical follow-up: TCGA, METABRIC, SCAN-B, CGGA and curated GEO ',
'series. A cohort is included for an endpoint only if it actually carries that endpoint, ',
'which is why the pool size <i>k</i> changes when you switch between OS and DFS.</p>',
.guide_table(f),

# ---------------------------------------------------------------- reading the output
'<h3 class="guide-h2">Reading what comes back</h3>',

'<h4 class="guide-h3">The hazard ratio</h4>',
'<p>Every result is a <b>hazard ratio per one standard deviation</b> of the score, within that ',
'cohort. HR&nbsp;=&nbsp;1.25 means: patients one SD higher on this score have 25% higher ',
'instantaneous risk of the event at any given moment, compared with otherwise-similar patients.</p>',
'<ul class="guide-list">',
'<li><b>HR &gt; 1</b>: higher score, worse outcome. <b>HR &lt; 1</b>: protective.</li>',
'<li>It is <b>not</b> a fold change, and not "25% more patients died". It is a rate ratio.</li>',
'<li><b>Read the confidence interval, not the point estimate.</b> HR 1.9 [0.8, 2.4] and ',
'HR 1.2 [1.15, 1.25] are not "1.9 is stronger" but the first is compatible with no effect ',
'and the second is not.</li>',
'<li>Per-SD units are what make cohorts comparable. A raw score unit means something different ',
'in each cohort; an SD does not.</li>',
'</ul>',

'<h4 class="guide-h3">The forest plot</h4>',
'<p>One row per cohort, so you can see the evidence before it is combined. The box is the ',
'cohort&rsquo;s HR and its area reflects how much that cohort contributes; the horizontal line is ',
'the CI; the <b>diamond at the bottom is the pooled estimate</b>.</p>',
'<ul class="guide-list">',
'<li><b>Pooling is random-effects Restricted Maximum Likelihood (REML)</b>, not fixed-effect. A fixed-effect pool assumes every ',
'cohort is estimating the identical number and treats disagreement between them as noise. Cohorts ',
'here differ in platform, treatment era and case mix, so they are not estimating the identical ',
'number, and a random-effects pool says so by widening the interval.</li>',
'<li><b>I&sup2; is the disagreement/heterogeneity</b> between cohorts, as a percentage of the total variation. ',
'Near 0% means they tell the same story. High I&sup2; with a significant diamond is the case to ',
'be careful with, the average may not describe any actual cohort.</li>',
'<li><b>A row that does not clear significance is not a failed replication.</b> Small cohorts have ',
'wide intervals. The question is whether the effects point the same way.</li>',
'</ul>',

'<h4 class="guide-h3">The Kaplan-Meier plot</h4>',
'<p>The KM curves split patients at the <b>median</b> of the score and show the observed survival ',
'of the two halves. It is the most intuitive panel and the easiest to over-read, so:</p>',
'<ul class="guide-list">',
'<li><b>The KM is not the model.</b> The HR is fitted on the <i>continuous</i> score, using every ',
'patient&rsquo;s actual value. The KM throws that away and dichotomises at the median, purely to ',
'draw a picture. They can disagree in emphasis; the HR is the result.</li>',
'<li>Curves are drawn per cohort, never on merged patients; pooling happens at the level of ',
'effect estimates, never by concatenating datasets.</li>',
'<li>Vertical ticks are censored patients (last known alive, or lost to follow-up).</li>',
'<li>Separation late in a curve where few patients remain is the least reliable part of it.</li>',
'</ul>',

'<h4 class="guide-h3">Three things the panel tells you that are easy to miss</h4>',
'<ul class="guide-list">',
'<li><b>The horizon &tau;.</b> Each tissue declares a follow-up horizon and patients are ',
'administratively censored there. Cohorts follow patients for very different lengths of time, and ',
'without a common horizon a pooled HR partly reflects who was followed longest rather than who ',
'did worse. The horizon in force is printed under the forest.</li>',
'<li><b>Stratification.</b> Where a tissue declares a stratifier (breast on PAM50 subtype, glioma ',
'on IDH status), the Cox model gives each stratum its <i>own baseline hazard</i> while sharing one ',
'coefficient. It removes the confound of subtype composition differing between cohorts. It ',
'assumes the effect is the same size in every stratum, which is an assumption, and it is stated ',
'rather than hidden: the forest axis names it, and the table reports how many cohorts were ',
'actually stratified.</li>',
'<li><b>The rank and q-value.</b> On a VIPER query the panel also reports where that TF sits in ',
'the genome-wide scan for the same tissue and endpoint, with a Benjamini-Hochberg q. This is ',
'the multiplicity context: a p of 0.01 means something different when it is one test than when it ',
'is one of ', .g_int(f$n_tf), '. The rank is looked up from the scan on disk, and only shown when ',
'the query&rsquo;s recipe matches the scan&rsquo;s exactly.</li>',
'</ul>',

# ---------------------------------------------------------------- caveats
'<h3 class="guide-h2">What a hit does <i>not</i> mean</h3>',
'<ul class="guide-list">',
'<li><b>Prognostic is not causal, and not predictive of treatment benefit.</b> Everything here is ',
'an association with outcome in patients who were treated in whatever way their cohort was ',
'treated.</li>',
'<li><b>Hit counts are inflated by correlation.</b> Regulon activities are strongly correlated ',
'with one another, so "N regulators clear FDR" is not N independent findings. Where the tool ',
'reports a replicated set, read its block structure: correlated programs travel together.</li>',
'<li><b>A zero in the replicated column is usually power, not absence.</b> The split fixes which ',
'cohorts discover and which validate, so a discovery arm with few events finds little regardless ',
'of what is there.</li>',
'<li><b>Proportional hazards is an assumption and it is tested.</b> Each fit is checked, and the ',
'panel warns when a cohort fails; the scans carry the count. A failed check means the single HR ',
'is a time-average of an effect that changed during follow-up, not that the effect is ',
'absent.</li>',
'<li><b>Regulons are inferred, not measured.</b> ARACNe edges are statistical, from one cohort&rsquo;s ',
'expression. A regulon is a good aggregate readout, not a validated target list.</li>',
'</ul>',

# ---------------------------------------------------------------- the tabs
'<h3 class="guide-h2">The tabs</h3>',
'<ul class="guide-list">',
'<li><b>Single query:</b> One score, one tissue, one endpoint: forest across cohorts, ',
'per-cohort KM, tumour-vs-normal context, and the genome-wide rank.</li>',
'<li><b>Multiple query:</b> One score across <i>all</i> ', .g_int(f$n_types),
' cancer types at once, as independent panels. Each panel keeps its own cohorts and its own ',
'horizon, and they are <b>never pooled with each other</b>: a pan-cancer diamond over different ',
'diseases would be an average of things that are not the same quantity.</li>',
'<li><b>Multiple genes:</b> Several genes at once in one tissue, for comparison at a glance.</li>',
'<li><b>Protein (RPPA) / Reverse-Phase Protein Array</b>: Antibody-based protein ',
'and phospho-protein levels, measured directly rather than inferred. Kept ',
'deliberately separate: RPPA exists for one cohort per tissue, so those results are single-cohort ',
'and must never enter a pooled diamond.</li>',
'</ul>',
# Two closing lines were dropped here on 2026-08-30 (step 114), both at the user's call, and
# both worth recording because each was carrying something:
#   - "This page has no tab of its own..." pointed at the brand as the way back. The
#     affordance now rests entirely on the lockup's title= tooltip, which the Guide's test
#     asserts in its place -- weaker, and deliberately so.
#   - "N of N tissues are at that level today" was vacuous the moment the ratio hit 7 of 7,
#     which it has. n_scanned came off guide_facts() with it rather than being left derived
#     for nobody.

# ---------------------------------------------------------------- how it is built
# MOVED HERE FROM R/about.R on 2026-08-30 (step 118), at the user's call, and it belongs here
# for a reason worth keeping: on About it ended by telling the reader that the Guide explains
# each step -- a section whose closing move is to send you elsewhere for the rest of itself is
# a section on the wrong page. That sentence is gone; this IS the elsewhere. What it gained
# instead is the reciprocal pointer at the end, so the two pages now refer to each other
# rather than each carrying half of both subjects.
'<h3 class="guide-h2">How it was built</h3>',
'<p>Tissue-specific ARACNe networks, VIPER regulon activity per patient, per-cohort Cox ',
'models and random-effects meta-analysis: the steps above, in the order the pipeline ',
'runs them. The engine is R, with <code>survival</code> for the models, <code>metafor</code> ',
'for the pooling, and Shiny for this interface.</p>',
# NAMED WITHOUT ITS PATH since step 142. This paragraph said "the cohort registry
# (config/cohorts.tsv)", which is true, and which means nothing to a reader who has a browser
# and not a checkout. The POINT of the sentence survives without the filename: what matters to
# a visitor is that one record governs the whole app and that every page derives from it, not
# where that record sits on somebody's disk.
'<p>One registry ships inside the tool and is the single record of which cohorts, endpoints, ',
'horizons and stratifiers exist. Every page of this app, including every figure on this one, ',
'reads what it states from that registry when the app starts rather than repeating it, which ',
'is why nothing here can quietly disagree with what was actually analysed. If a cohort were ',
'added or an endpoint dropped, this page would say so on the next start without anyone ',
'editing it.</p>',
'<p class="guide-fine">Where the data came from and on whose terms, the licence, how to ',
'cite the tool and who built it are on the <b>About</b> tab.</p>',

'</div>')
}


# ---- 3. the page's own styling ------------------------------------------------------
# Kept here rather than in app.R's PAGE_CSS so the guide is one file -- prose, figures and
# presentation together. Every class is prefixed `guide-`, which is what keeps it clear of
# the `.note*` set that tests/test_note_classes.R pins by name.
GUIDE_CSS <- local({
  # local(), not a top-level loop: tests/test_source_side_effects.R requires every top-level
  # expression in R/ to be an assignment or a load, so that source()ing a file cannot RUN
  # anything -- the rule that came out of step 18, where sourcing drivers overwrote 26
  # artifacts. An assignment whose right-hand side does the work satisfies it honestly rather
  # than by hiding the work somewhere the check does not look.
  css <- "
/* No max-width since step 112: the Guide and About pages now run to the same edges as the
   four query tabs, which is what .narrow-page (1400px, 5% inset) already sets for all six.
   The 62rem cap was there for MEASURE -- prose at this width is a long line to track back
   from -- and that cost is real and accepted: matching the other tabs was the ask, and a
   landing page that is visibly narrower than every other tab reads as a different site. */
.guide { max-width: none; margin: 0 auto; padding: 0.5rem 0 3rem 0; line-height: 1.62; }
.guide-logo { height: 165px; width: auto; display: block; margin-left: -6px; }
/* Hero: lockup on the left, the lead paragraph beside it rather than under it (step 113).
   Under it, the logo left a wide empty band to its right on a page that is now full width,
   and the first sentence started 190px down. Scoped to `.guide-hero` so the base .guide-h1
   and .guide-lead margins keep working for the About page, which has no hero. It wraps on a
   narrow window because the row is flex-wrap, not because a breakpoint was written for it. */
.guide-hero { display: flex; flex-wrap: wrap; align-items: center; gap: 0.4rem 2rem;
              margin: 0.4rem 0 1.8rem 0; }
.guide-hero .guide-h1 { margin: 0; flex: 0 0 auto; }
.guide-hero .guide-lead { margin: 0; flex: 1 1 26rem; }
.guide p, .guide li { font-size: 108%; }
.guide-h1 { margin: 0.6rem 0 0.2rem 0; font-weight: 700; }
.guide-lead { font-size: 123% !important; color: #333; margin-bottom: 1.6rem; }
.guide-h2 { margin: 2.1rem 0 0.6rem 0; padding-bottom: 0.3rem; font-weight: 700;
            border-bottom: 2px solid #e6e6e6; }
/* MEASURED IN THE BROWSER, not reasoned about. Every percentage in this block resolves
   against the same 13px inherited base, so `.guide p` at 108% is 14.04px -- and `.guide-h3`
   at its original 104% was 13.52px, a sub-heading rendering SMALLER than the paragraph
   beneath it. That was true on the Guide from step 110 and nothing could see it: these
   percentages LOOK as though they were relative to the body text, and are not. It surfaced
   only when step 118 added a fourth level, which landed at 13px -- underneath the third as
   well. The ladder now reads as one: h2 23 / h3 15.86 / h4 14.56 / body 14.04 / fine 12.61,
   read off the live page rather than derived from the percentages, since the base is the
   thing that was wrong. Anything changed here has to be looked at in a browser again; no
   test in this suite can compute a font-size. */
.guide-h3 { margin: 1.4rem 0 0.35rem 0; font-weight: 700; font-size: 122%; }
/* A fourth level, added in step 118 for the About page: its source families sit inside
   `Licence and data terms` under an h4 of their own, and at .guide-h3 a single source name
   would carry the same weight as the heading introducing the whole list of them. Defined
   here and not in an About stylesheet because R/about.R deliberately declares none -- two
   stylesheets for two pages of the same document is how they start to look like different
   sites. tests/test_guide_derived.R's both-directions class check therefore reads BOTH
   rendered pages; until step 118 it read only the Guide, which was safe exactly as long as
   About used no class the Guide did not. */
.guide-h4 { margin: 1.1rem 0 0.3rem 0; font-weight: 700; font-size: 112%; color: #444; }
/* The institutional marks on About. Their aspect ratios differ by roughly a factor of five,
   the lab mark is square, the university's is a wide lockup, so the HEIGHT is declared per
   image in ABOUT_LOGOS and only the row is styled here; one shared height would render one of
   them at several times the visual weight of another. align-items:center puts a square and a
   wide lockup on a shared optical centre rather than a shared top edge. */
.guide-marks { display: flex; flex-wrap: wrap; align-items: center; gap: 0.6rem 2.2rem;
               margin: 1.2rem 0 0.4rem 0; }
.guide-mark { width: auto; display: block; }
.guide-list { margin: 0.4rem 0 0.9rem 0; }
.guide-list li { margin-bottom: 0.42rem; }
.guide-fine { font-size: 97% !important; color: #666; margin-top: 0.7rem; }
.guide-callout { border-left: 3px solid __BRAND_ORANGE__; background: __BRAND_WARM__;
                 padding: 0.7rem 0.9rem; margin: 1rem 0; }
.guide-steps { margin: 0.9rem 0 1.2rem 0; }
.guide-step { display: flex; gap: 0.85rem; align-items: flex-start; margin-bottom: 0.8rem; }
.guide-num { flex: 0 0 1.7rem; height: 1.7rem; line-height: 1.7rem; text-align: center;
             border-radius: 50%; background: __BRAND_PURPLE__; color: #fff; font-weight: 700; }
.guide-tbl { width: 100%; border-collapse: collapse; margin: 0.8rem 0; font-size: 99%; }
.guide-tbl th { text-align: left; border-bottom: 2px solid #ccc; padding: 0.4rem 0.55rem;
                white-space: nowrap; }
.guide-tbl td { border-bottom: 1px solid #eee; padding: 0.36rem 0.55rem; }
.guide-tbl tr:hover td { background: #fafafa; }
.guide-wrap { overflow-x: auto; }
"
  # Same substitution pass app.R runs over PAGE_CSS, and for the same reason: the colours are
  # stated once, in R/brand.R, and every stylesheet takes them from there. R/brand.R is sourced
  # BEFORE this file (see app.R) so the constants exist; the check below turns a missed one
  # into an error here rather than a declaration the browser drops in silence.
  for (k in list(c("__BRAND_PURPLE__", BRAND_PURPLE), c("__BRAND_ORANGE__", BRAND_ORANGE),
                 c("__BRAND_WARM__",   brand_tint(BRAND_ORANGE, 0.94))))
    css <- gsub(k[[1]], k[[2]], css, fixed = TRUE)
  if (grepl("__BRAND_", css, fixed = TRUE))
    stop("GUIDE_CSS still carries a __BRAND_*__ placeholder after substitution.")
  css
})
