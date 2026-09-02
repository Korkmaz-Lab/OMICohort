# Plotting for the survival tool: KM curves (per cohort) and a meta-analysis forest
# plot, both driven straight off a `survresult` object from get_survival(). Base
# graphics + metafor::forest only (no survminer/ggplot dependency).
#
# A survresult is plot-shape-stable: metaVIPER changes the NUMBERS inside it, never
# its structure, so these plots work identically before/after the SCAN-B network.

# PDF devices in descending order of font handling. The base pdf() device references the
# base-14 fonts BY NAME instead of embedding the font programs: a viewer that ships
# Helvetica substitutes it and the figure looks perfect, so the defect is invisible on
# the machine that made it, while anywhere else every glyph vanishes and only the curves
# survive. Journals generally require embedded fonts, so pdf() is the last resort, not
# the default.
#
# Selected by RUNTIME PROBE, not by capabilities(): on this machine capabilities("cairo")
# reports TRUE while cairo_pdf() dies on a missing libXrender (no XQuartz). A flag that
# lies is worse than no flag, so the only trustworthy test is opening the device.
.PDF_DEVICES <- list(
  cairo = function(f, w, h) cairo_pdf(f, width = w, height = h),
  quartz = function(f, w, h) quartz(type = "pdf", file = f, width = w, height = h),
  base  = function(f, w, h) pdf(f, width = w, height = h)
)

.open_pdf <- function(file, width, height) {
  before <- grDevices::dev.cur()
  for (open in .PDF_DEVICES) {
    opened <- tryCatch({ suppressWarnings(open(file, width, height)); TRUE },
                       error = function(e) FALSE)
    if (opened && !identical(grDevices::dev.cur(), before)) return(invisible(TRUE))
  }
  stop("no usable PDF device")
}

.open_dev <- function(file, width = 7, height = 5) {
  if (is.null(file)) return(FALSE)
  ext <- tolower(tools::file_ext(file))
  if (ext == "pdf") .open_pdf(file, width, height)
  else png(file, width = width * 120, height = height * 120, res = 120)
  TRUE
}

# Full y-axis wording for an endpoint code (avoids the redundant "DFS survival ...").
.endpoint_label <- function(ep) {
  switch(ep,
         OS  = "Overall survival probability",
         DFS = "Disease-free survival probability",
         DSS = "Disease-specific survival probability",
         sprintf("%s survival probability", ep))
}

# Corner annotation for a KM panel. Split out from km_plot() so the wording is testable
# without rendering a device.
#
# The KM curves themselves are non-parametric and need no PH assumption -- they are
# right either way. The annotation is what can mislead: under a PH violation the Cox
# HR is an average of an effect that changed over follow-up, and the log-rank additionally
# LOSES POWER because early and late differences cancel in its sum. ESR1/METABRIC shows
# both failures at once (HR 0.95 p=0.30, log-rank p=0.056) on curves that visibly cross.
# So when PH is rejected, say so on the figure rather than letting the number stand alone.
#
# A truncated fit gets a second note, for the same reason: at tau=60 METABRIC/ESR1 is
# HR 0.59 where full follow-up gives 0.95. Those are different estimands, so a figure
# showing one without naming its horizon is wrong rather than merely incomplete.
.km_label <- function(r, max_followup = NULL) {
  lab <- sprintf("log-rank p = %.2g\nHR/SD = %.2f  (p = %.2g)",
                 r$logrank_p, r$HR, r$p)
  if (!is.null(max_followup))
    lab <- paste0(lab, sprintf("\nfollow-up truncated at %g mo", max_followup))
  if (r$strata_adjusted)
    lab <- paste0(lab, sprintf("\n%s-stratified", r$strata_var))
  if (isTRUE(r$ph_violated) && !is.na(r$ph_p))
    lab <- paste0(lab, sprintf("\n! PH rejected (p = %.2g)\n   HR above is a time-average", r$ph_p))
  lab
}

# Which corner the stats annotation goes in. A KM curve starts at 1.0 and decays, so it
# vacates the top-right only if it actually falls; when survival stays high -- a
# truncated horizon, or an indolent cohort -- the top-right is where the curves are and
# the text lands on them. Pick the half the curves have left empty.
.km_corner <- function(km) {
  lowest <- min(km$surv, na.rm = TRUE)
  if (is.finite(lowest) && lowest > 0.5) "bottomright" else "topright"
}

# Declared geometry for a KM figure, in inches, shared by the screen and the file export.
#
# Exported for the same reason as forest_height_in(): app.R hardcoded a 420px plotOutput
# while km_plot() took .open_dev()'s generic 7x5 defaults, so one figure had two shapes
# maintained in two places and nothing flagged it when they diverged. The browser view is
# the reference -- it is what this project is built around and what is confirmed to look
# right -- so 4.5in is the export conforming to the screen (420px / 96 = 4.4in), not the
# other way round.
#
# Landscape by construction: a KM curve is read left-to-right along time.
# Height allows for the number-at-risk panel beneath the curves (~22% of the figure).
KM_SIZE_IN <- c(width = 7, height = 5.5)
km_size_in <- function() KM_SIZE_IN

# The median split, computed ONCE. The curves, the legend's per-arm n/events, and the
# number-at-risk table must all describe the same two groups; three call sites each
# recomputing "who is high" is precisely how the KM and the Cox drifted into describing
# different patients (Phase 12). r$data is the engine's fitted sample (d2), so this split
# is the same one survfit() was given.
.km_arms <- function(r) {
  d <- r$data
  med <- stats::median(d$score)
  list(d = d, med = med, hi = d$score >= med)
}

# Number still under observation at each of `times`, per arm.
#
# "At risk at t" = follow-up has reached t without the event or censoring having happened
# first, i.e. time >= t. Counted straight off the fitted rows rather than parsed out of
# summary(survfit), which returns a strata-blocked frame that is easy to mis-slice --
# tests/test_km_risk_table.R cross-checks this against survfit's own n.risk.
.km_risk <- function(r, times) {
  a <- .km_arms(r)
  n_at <- function(tt) vapply(times, function(t) sum(tt >= t), integer(1))
  list(high = n_at(a$d$time[a$hi]), low = n_at(a$d$time[!a$hi]))
}

KM_COL <- c(low = "#2c7fb8", high = "#d95f0e")

# KM curves (median-split high vs low) for ONE cohort, with the log-rank p and a
# number-at-risk panel beneath.
#
# risk_table = FALSE restores the pre-Phase-16 single-panel figure. Kept because the risk
# panel costs vertical space, and a caller embedding the curve somewhere small may want the
# curves alone -- not because the counts are optional information.
km_plot <- function(res, cohort, file = NULL, risk_table = TRUE) {
  r <- res$per_cohort[[cohort]]
  if (is.null(r) || isTRUE(r$skipped) || is.null(r$km))
    stop(sprintf("no KM fit for cohort '%s'", cohort))
  s <- km_size_in()
  opened <- .open_dev(file, width = s[["width"]], height = s[["height"]])
  # Capture BEFORE layout()/par() are touched, and restore before closing -- par() with no
  # device open opens the default one (Phase 14). layout() must be undone too: it is device
  # state, so leaving it set would tile whatever the caller draws next. That matters for the
  # app, which renders into a device it owns and reuses.
  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    suppressWarnings(graphics::par(op))
    if (opened) grDevices::dev.off()
  })

  feat <- attr(res, "feature") %||% "score"
  # ONE median split, shared by the curves, the legend and the risk panel.
  a <- .km_arms(r); d <- a$d; med <- a$med; hi <- a$hi
  n_hi <- sum(hi); n_lo <- sum(!hi)
  ev_hi <- sum(d$event[hi] == 1); ev_lo <- sum(d$event[!hi] == 1)

  if (risk_table) {
    # Curves and risk panel stacked. Left/right margins are IDENTICAL between the two so
    # the panels share a horizontal coordinate system -- that is what makes a count sit
    # under its tick. The x axis stays on the curve panel; only "Months" moves down.
    graphics::layout(matrix(c(1L, 2L), nrow = 2L), heights = c(3.6, 1))
    graphics::par(mar = c(2.6, 4.1, 3.0, 2.1))
  }
  # Same treatment as the forest, and it has to be set BEFORE plot() for the same reason:
  # mar takes effect at the next plot.new(). The cohort and endpoint are the `tail` here --
  # a KM of the wrong cohort is not a smaller error than a KM of the wrong feature, and
  # clipping ate that end first. cex.main is 1.2 by default, so the fit is measured at 1.2.
  tl <- .title_lines(feat, sprintf("  |  %s  (%s)", cohort, res$endpoint),
                     w = .title_width(), cex = 1.2)
  plot(r$km, col = unname(KM_COL[c("low", "high")]), lwd = 2,
       xlab = if (risk_table) "" else "Months", ylab = .endpoint_label(res$endpoint),
       mark.time = TRUE,
       main = paste(tl, collapse = "\n"))
  # Split is at the cohort MEDIAN of the (within-cohort z-scored) score: top 50% = "high".
  legend("bottomleft", bty = "n", lwd = 2, col = c("#d95f0e", "#2c7fb8"),
         title = sprintf("split at median (z = %.2f)", med), title.adj = 0,
         legend = c(sprintf("high  n=%d, ev=%d", n_hi, ev_hi),
                    sprintf("low   n=%d, ev=%d", n_lo, ev_lo)))
  # One legend entry PER LINE. legend() sizes its box by entry count and ignores
  # embedded newlines, so a single \n-joined string is measured as one line high: from
  # "bottomright" the box is anchored at the plot floor and every line after the first
  # is drawn below it and clipped away. (Invisible while the annotation always sat in
  # "topright", where the overflow grew down into empty space.)
  legend(.km_corner(r$km), bty = "n", text.col = "grey20",
         legend = strsplit(.km_label(r, res$max_followup), "\n")[[1]])

  if (risk_table) {
    # Tick positions are taken from the curve panel that was just drawn, not recomputed, so
    # a count cannot land between ticks. Likewise the x range comes from par("usr") -- which
    # already includes R's 4% axis padding -- and is replayed with xaxs="i" so the second
    # panel maps time to x exactly as the first did.
    ticks <- graphics::axTicks(1)
    usr_x <- graphics::par("usr")[1:2]
    rk <- .km_risk(r, ticks)

    graphics::par(mar = c(2.7, 4.1, 0.4, 2.1))
    graphics::plot.new()
    graphics::plot.window(xlim = usr_x, ylim = c(0.3, 2.9), xaxs = "i")
    graphics::mtext("Number at risk", side = 3, adj = 0, line = -0.6,
                    cex = 0.8, col = "grey20", font = 2)
    for (i in seq_along(ticks)) {
      graphics::text(ticks[i], 2.0, rk$high[i], cex = 0.85, col = KM_COL[["high"]])
      graphics::text(ticks[i], 1.1, rk$low[i],  cex = 0.85, col = KM_COL[["low"]])
    }
    # Row labels in the left margin, aligned to the same rows. las = 1 keeps them
    # horizontal; the ylab column is already reserved by the shared left margin.
    graphics::mtext("high", side = 2, at = 2.0, las = 1, line = 0.4,
                    cex = 0.85, col = KM_COL[["high"]])
    graphics::mtext("low",  side = 2, at = 1.1, las = 1, line = 0.4,
                    cex = 0.85, col = KM_COL[["low"]])
    graphics::mtext("Months", side = 1, line = 1.3, cex = 1)
  }
  invisible(res)
}

`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a

# p-value formatting for every table and figure in the tool, defined ONCE.
#
# It was defined twice as a local, inside survtable_display() and inside mgene_display(),
# and the two copies had already diverged: one was vectorized via ifelse(), the other was
# scalar and NULL-safe. Neither was wrong for its own caller, which is why nothing caught
# it -- and tn_plot() would have been the third. This version does both jobs: NULL and
# length-0 return "" rather than the character(0) that silently drops a sprintf field.
fmt_p <- function(v) {
  if (is.null(v) || !length(v)) return("")
  ifelse(is.na(v), "", formatC(v, format = "g", digits = 3))
}

# --- figure titles that do not lie about being cut ----------------------------------------
#
# Longest prefixes of `s` that each fit within `w` inches, in order. Bisection rather than a
# character-at-a-time walk because strwidth() is a device call.
#
# strwrap() is NOT usable here and that is the whole reason this exists: it breaks on
# WHITESPACE, and the names that overflow are `///`-joined probe sets with none
# ("C2///CFB", "FAM90A9P///FAM90A10P///..."). Handed one, strwrap() returns it whole.
.chunk_to_width <- function(s, w, cex, font) {
  out <- character(0)
  while (nchar(s) > 0L) {
    n <- nchar(s)
    if (graphics::strwidth(s, units = "inches", cex = cex, font = font) <= w) {
      out <- c(out, s); break
    }
    lo <- 1L; hi <- n
    while (lo < hi) {
      mid <- (lo + hi + 1L) %/% 2L
      if (graphics::strwidth(substr(s, 1L, mid), units = "inches", cex = cex, font = font) <= w)
        lo <- mid else hi <- mid - 1L
    }
    out <- c(out, substr(s, 1L, lo)); s <- substr(s, lo + 1L, n)
  }
  out
}

# A plot title as the lines to draw, wrapped to the device and truncated ALOUD if it still
# does not fit. Returns a character vector, one element per line.
#
# What it is fixing, measured 2026-08-10 (BUILD_LOG step 38): the title was a single
# sprintf() handed to title()/plot(), and a name wider than the device was CLIPPED AT BOTH
# ENDS with no sign that anything had been cut. `pdftotext` on a delivered forest recovered a
# 74-character middle slice of a 266-character title -- it does not begin where the symbol
# begins and it has lost the trailing "(expr) vs OS", so what remains reads like a complete
# name for a different feature. A wrong label is worse than a missing one.
#
# `tail` is the part that MUST survive: the kind and endpoint (and, on a KM, the cohort).
# That is what says what the figure IS, and clipping took it first, being at the end. It is
# never cut -- the marker is inserted before it.
#
# max_lines = 2 is derived, not picked -- and the DEFINITION is written down here, because the
# 2026-08-10 version of this paragraph gave figures ("97.0%", "98.6%", "626") with no statement
# of which title, geometry or cex they were measured at, which is exactly why nobody could tell
# they had gone stale when brain doubled the vocabulary. Re-derived 2026-08-27 over the union of
# every feature the app offers (58366 DISTINCT strings across seven tissues x two query kinds;
# .feature_tag()'s note says 59964 because it counts a symbol once per kind), on the KM at its
# export geometry -- 7x5.5in, mar c(2.6, 4.1, 3.0, 2.1), cex.main 1.2, so w = 6.60in -- with the
# longest cohort id in the tail: the title fits on ONE line for 96.4% and within TWO for 98.7%.
# A third line would buy 0.5% and cost every long-titled figure a row of plot area. The
# remaining 743 get "[+<n> chars]", which is the same rule .feature_tag() applies to a filename:
# state the cut. The forest is roomier and reaches the same conclusion with margin (7.5in,
# cex 1.1, w = 7.10in: 98.1% / 99.1%, 522 marked, a third line +0.3%). tests/test_feature_
# vocabulary.R re-derives all of it and fails if any figure moves. The number is there so the reader knows the SCALE of what is
# missing -- "[+11 chars]" and "[+1352 chars]" are very different figures -- and the whole
# symbol is in the exported CSV's `feature` column, which is where it has always been.
.title_lines <- function(feat, tail, w, cex, font = 2L, max_lines = 2L) {
  full <- paste0(feat, tail)
  if (graphics::strwidth(full, units = "inches", cex = cex, font = font) <= w) return(full)
  ch <- .chunk_to_width(full, w, cex, font)
  if (length(ch) <= max_lines) return(ch)
  head_lines <- ch[seq_len(max_lines - 1L)]
  shown <- sum(nchar(head_lines))
  cut <- nchar(full) - shown
  # The marker states how many characters were dropped, and its own width depends on that
  # number, which changes how many characters fit beside it. Iterate to the fixed point
  # rather than guessing: at most a few passes, since the count only shrinks.
  repeat {
    mark <- sprintf("[+%d chars]%s", cut, tail)
    room <- w - graphics::strwidth(mark, units = "inches", cex = cex, font = font)
    if (room <= 0) return(c(head_lines, mark))
    keep <- .chunk_to_width(substr(full, shown + 1L, nchar(full)), room, cex, font)[1L]
    new_cut <- nchar(full) - shown - nchar(keep)
    if (new_cut == cut) return(c(head_lines, paste0(keep, mark)))
    cut <- new_cut
  }
}

# The widest a centred main title can be and still land on the page.
#
# NOT par("fin"). title(main=) centres on the PLOT region, not the figure region, so a title
# wider than the plot region overhangs equally into BOTH side margins -- and the narrower of
# the two is what stops it. Measured 2026-08-10 by drawing a string of known width with
# distinctive end characters and reading the PDF back: at the KM's geometry (7in, mar left
# 4.1 / right 2.1) a title sized to par("fin") comes back MISSING ITS LAST CHARACTER, and one
# sized to this comes back whole. That is how this function got written -- the first draft
# used par("fin") - 0.1 and produced a KM title reading "... (OS" with the ")" clipped off,
# which is the original defect in miniature and would have shipped looking like a fix.
.title_width <- function() {
  mai <- graphics::par("mai")
  graphics::par("pin")[1L] + 2 * min(mai[2L], mai[4L])
}

# NOTE there is deliberately NO top-margin adjustment for a wrapped title, and no growing of
# the device height either.
#
# The height cannot be grown: how many lines a title takes is a property of the device it is
# drawn on, and forest_height_in() is called by app.R to size a plotOutput BEFORE any device
# exists -- the screen device is far wider than the 7.5in export, so the same title wraps
# differently in the two places.
#
# The margin does not need to be. This step first added a .title_mar() that widened mar[3],
# then measured it: the default top margin is 0.60in on both figures, and two lines cost
# 0.44in on the forest (cex 1.1) and 0.48in on the KM (cex 1.2). The adjustment bought
# headroom nothing needed and took ~0.24in off the plot region to do it. Deleted. What is
# left in its place is an ASSERTION, in tests/test_forest_label_fit.R and test_km_size.R,
# that two lines still fit the untouched margin -- so raising max_lines or the title cex
# fails loudly instead of silently clipping the top line off the page.

# x-axis label for the forest. Split out like .km_label() so the truncation disclosure
# is testable without rendering a device. Truncation belongs on the AXIS, not the title:
# it qualifies what the HR scale means, and it applies to every row including the pooled
# diamond (a common tau is the only reason the rows are poolable in the first place).
# WHY THIS IS THE ONLY SURFACE, and why the pooled label below is NOT a second one.
#
# A "! MIXED: 1 of 2 fitted with strata(idh)" line under the diamond was written and then
# REMOVED, because it would be a FOURTH label line and forest_height_in() buys room for a
# third: its allowance is `if (k >= 2 && k < POOL_MIN_K) 1.3`, keyed on k alone. Making the
# height depend on the label's actual line count means changing a signature that app.R reads
# for a CSS aspect-ratio (app.R:1226) and that four test files pin as 2.4 + 0.5 * k -- a
# layout contract, rewritten for a case with zero artifacts on disk (no mixed pool is
# scannable, step 93, and lgg/DFS has no panel, step 100). The x-axis states the mixture
# exactly, print() warns in full sentences and survtable() carries pool_strata_k, so the
# fact is not lost; only the second rendering of it is.
#
# What the figure must say about stratification, or NULL when there is nothing to say.
#
# The forest was the one surface that did not disclose it. print() warns, survtable() carries
# `pool_strata_k`, survtable_display() renders "2 of 3" -- but a PDF travels ALONE, further from
# its context than any of those, and lgg's forest showed a pooled n of 1050 where the marginal
# pool has 1106 with nothing on the page to say why. Step 100 put the first IDH-stratified
# figures in results/; this is the annotation they needed.
#
# NULL at strata_k == 0, deliberately, and for the same reason .strat_col() renders "" there:
# an all-unstratified pool is one model applied throughout, and "strata: 0 of 5" would read as
# a defect rather than as the ordinary case. So this changes no unstratified figure -- today
# that is every tissue except breast (pam50) and lgg (idh).
#
# NULL also when `pooled` is absent or the field predates it: an old survresult cannot be
# asked, and inventing "unstratified" for it would be a claim rather than a reading.
.forest_strata_phrase <- function(res) {
  p <- res$pooled
  if (is.null(p) || is.null(p$strata_k) || is.na(p$strata_k) || p$strata_k < 1L) return(NULL)
  v <- paste(p$strata_vars, collapse = "/")
  # A mixed pool must never render as a plain "strata(v)": that would claim of the whole
  # diamond what is true of part of it, which is exactly the misnaming step 93 refused for
  # scan FILENAMES. Same rule here -- the count travels with the name or the name does not go.
  if (isTRUE(p$strata_mixed)) sprintf("strata(%s) in %d of %d", v, p$strata_k, p$k)
  else sprintf("strata(%s)", v)
}

.forest_xlab <- function(res) {
  x <- sprintf("HR per SD  (%s", res$endpoint)
  if (!is.null(res$max_followup))
    x <- paste0(x, sprintf(", truncated at %g mo", res$max_followup))
  # The stratification sits with the horizon on purpose. This parenthesis is the figure's
  # ESTIMAND line -- what was measured, over what follow-up, under what model -- and a
  # strata() term belongs to that list exactly as "truncated at 94 mo" does. It is also the
  # only annotation that reaches BOTH branches of forest_plot(): the k=1 branch draws no
  # pooled label at all, so anything written only into the summary line would silently skip
  # the single-cohort figures the app can produce.
  s <- .forest_strata_phrase(res)
  if (!is.null(s)) x <- paste0(x, ", ", s)
  paste0(x, ")")
}

# Forest plot of per-cohort per-SD HRs + the REML pooled diamond. Reconstructs the
# rma from the per-cohort logHR/se stored on the survresult (same inputs the engine
# pooled), so the diamond matches res$pooled exactly.
# Device height (inches) for a forest of k estimable cohorts: title + header + k study
# rows + gap + summary + axis.
#
# Exported because the SCREEN needs it too. app.R used to hardcode 380px (~3.9in) while
# the PDF export used this formula (8.9in at ovarian's k=13), so one figure had two
# shapes maintained in two places. A 13-row forest does still draw legibly at 380px --
# this is about the two devices not drifting, not about a crash.
#
# The k<POOL_MIN_K term buys vertical room for the summary label's THIRD line (the tau^2
# caveat -- see .pooled_mlab()). A text line is a fixed number of INCHES while the summary
# slot is a fixed number of ROWS, so the squeeze is worst where the plot is shortest: at
# k=2 the three-line block is drawn through the rule above the pooled row. Deliberately
# NOT monotone in k -- k=4 carries the caveat and k=5 does not, so the k=4 figure is the
# taller one.
forest_height_in <- function(k) 2.4 + 0.5 * k + if (k >= 2 && k < POOL_MIN_K) 1.3 else 0

# The width forest_height_in() is a height FOR. It was a literal 7.5 inside forest_plot()
# and nowhere else, which was fine while only the PDF cared. It stopped being fine when the
# Multiple query tab began sizing its on-screen forests (step 53): a height in inches is
# meaningless without the width it was designed against, and a second literal 7.5 in app.R
# is a number that can drift from this one. Same rule as RPPA_HR_HEADER -- one constant,
# measured and drawn.
FOREST_WIDTH_IN <- 7.5

# The forest's two column headers. Constants because step 55 added a THIRD place that has
# to say the same thing -- .forest_axis() hands metafor a throwaway panel to lay out, and a
# layout is only transferable if the panel it was computed for carried the same headers.
# Same rule as RPPA_HR_HEADER, which now reads from this one rather than repeating it: step
# 51 matched the RPPA panel's header to the forest's on purpose, and two literals holding
# equal strings record that as a coincidence instead of as the decision it was.
FOREST_SLAB_HEADER <- "Cohort"
FOREST_HR_HEADER   <- "HR [95% CI]"

# --- Multiple query: the stacked-panel layout (2026-08-21, step 53) ---------------------
# The panels were a 2-up grid and are now one full-width card per tissue. The grid could
# not be made to work here, for two reasons that are both structural rather than cosmetic:
#
#   * five tissues into rows of two is 2, 2, 1 -- the last card is half-width and alone.
#     Five will not divide by two, and the tissue count is read from the registry.
#   * cards in one row are never the same height. Forest height is derived per panel from
#     THAT panel's k, so breast (k=3) draws at 499px beside ovarian (k=13) at 854px. That
#     rule is correct and shared with the single-query tab and the PDF export; a grid is
#     simply the one layout in which being correct looks like a defect.
#
# The forest is capped rather than stretched to the card. A full-width card is what makes
# the column symmetric, but a forest drawn at 1240x499 is a 2.5:1 box holding a figure
# designed at 7.5:5.2 -- the CI region turns into empty space and the labels drift away
# from the plot. So the CARD is full width and the FIGURE keeps its designed proportions.
#
# Since step 56 this caps the SINGLE query forest too (.sq-forest in app.R's PAGE_CSS), so
# the MQ_ in the name is now under-specific: it is one on-screen forest width for the tool.
# Left as it is on purpose -- the name reaches app.R, two tests, NAMESPACE and a man page,
# and renaming an exported symbol to tidy a prefix is churn that this comment can prevent
# instead. What must NOT happen is a second 900 appearing next to it for the other tab.
MQ_FOREST_W_PX <- 900L

# On-screen height for a k-row forest at MQ_FOREST_W_PX, in the design's own proportions.
# Not forest_height_in(k) * 96: that is the PDF's pixel height at the PDF's width, and it
# only happened to look right while the card was ~562px wide.
mq_forest_h_px <- function(k)
  as.integer(round(MQ_FOREST_W_PX * forest_height_in(k) / FOREST_WIDTH_IN))

# Whether a panel's detail (KM curves + per-cohort table) is open, from the toggle's own
# click counter. PARITY, deliberately, and not a reactiveValues flag beside it: re-rendering
# the panels re-creates the actionLink and Shiny reports the fresh element's counter as 0,
# so a separate flag would survive a re-render that the counter does not and the label would
# end up disagreeing with the content. One counter, read the same way in both places, cannot
# disagree with itself. n = 0 (never clicked, or just re-created) is CLOSED, which is also
# the state a new run should land in.
#
# length(n) == 1L carries the NULL case on its own -- length(NULL) is 0 -- so there is no
# is.null() guard here. There was one; mutating it away changed nothing, which is how it
# was found. An unreachable guard reads as protection and is not.
mq_detail_open <- function(n)
  length(n) == 1L && !is.na(n) && as.integer(n) %% 2L == 1L

# The toggle's two labels. In R/ rather than app.R for the reason every other user-facing
# string in this project is: app.R has no test harness.
MQ_DETAIL_SHOW <- "Show cohort KM curves and per-cohort table"
MQ_DETAIL_HIDE <- "Hide cohort KM curves and per-cohort table"

# The superscripts are PLOTMATH, and that distinction is the whole safety argument.
#
# Never a literal "I<sup>2</sup>" or "tau<sup>2</sup>" character. R's pdf() device defaults to
# Helvetica in a single-byte encoding that has neither the Greek tau nor superscript-two, and
# it drops both SILENTLY -- no warning, no error, no condition of any kind. That is not
# hypothetical: the caveat was once written with those glyphs and every k<5 forest between
# that edit and 2026-07-29 rendered ".... not identified! (k=4 < 5)", found only by looking at
# the figures. Re-verified 2026-08-07 -- a literal superscript-two still comes out of base
# pdf() as "I..".
#
# Plotmath does not take that path at all. `I^2` is CODE, not a string, so the superscript is
# TYPESET from an ordinary "2" and the Greek comes from the base-14 Symbol font; both appear in
# the PDF's text layer where the literals produced nothing. It is also why this stays inside
# tests/test_pdf_fonts.R's rule that every string LITERAL in this file is ASCII -- the rule is
# untouched, because no glyph here is a literal.
#
# A diamond drawn from an unidentified tau^2 must not look like any other diamond. The
# figure is the artifact most likely to be read alone, so the caveat travels ON it --
# same reason the horizon is on the x-axis and the PH verdict is in the KM corner.
.pooled_mlab <- function(n, ev, i2, p, k) {
  lines <- list(sprintf("Pooled (REML):  n=%d, ev=%d", n, ev),
                bquote(I^2 * .(sprintf(" = %.0f%%,  p = %.2g", i2, p))))
  if (k < POOL_MIN_K)
    lines <- c(lines, list(bquote("! " * tau^2 *
                                  .(sprintf(" not identified (k=%d < %d)", k, POOL_MIN_K)))))
  lines
}

# Stack the label's lines into one plotmath expression, LEFT-ALIGNED.
#
# atop() is the only way to get several plotmath lines into a single label ("\n" is a string
# feature and does not survive an expression), but it CENTRES its lines on each other, which
# left the shorter line indented under the longer one and out of line with every cohort label
# in the column. So each line is padded on the right with an invisible phantom() sized to the
# difference; centring a padded line puts its visible part flush left.
#
# The pad is measured in INCHES, before metafor draws: strwidth() in user units needs a plot
# that does not exist yet, while the inches ruler works on a bare device -- and inches is the
# honest unit here anyway, since that is what the collision this label already survived was
# measured in. Widths are compared as a RATIO, so the device's cex cancels out.
#
# The phantom is built from a prefix of the widest line rather than from repeated spaces:
# real text gives real font metrics, so the match is as exact as the font allows.
#
# Every line is wrapped in textstyle(), which is NOT cosmetic. Three lines need a NESTED
# atop, and plotmath shrinks each nesting level -- so the k<POOL_MIN_K label came out with
# its first two lines visibly smaller than the caveat under them, as if the pooled n/ev were
# a footnote to the warning. textstyle() pins every line to the base size, which also keeps
# the widths measured above equal to the widths finally drawn.
.stack_left <- function(lines) {
  if (length(lines) == 1L) return(as.expression(lines[[1]]))
  w <- vapply(lines, function(e) graphics::strwidth(as.expression(e), units = "inches"),
              numeric(1))
  ref <- lines[[which.max(w)]]
  ref <- if (is.character(ref)) ref else paste(deparse(ref), collapse = "")
  padded <- lapply(seq_along(lines), function(i) {
    d <- max(w) - w[i]
    e <- if (d <= 0) lines[[i]] else {
      ks <- seq_len(nchar(ref))
      pw <- vapply(ks, function(k) graphics::strwidth(substr(ref, 1, k), units = "inches"),
                   numeric(1))
      bquote(.(lines[[i]]) * phantom(.(substr(ref, 1, ks[which.min(abs(pw - d))]))))
    }
    bquote(textstyle(.(e)))
  })
  as.expression(Reduce(function(a, b) bquote(atop(.(a), .(b))), padded))
}

# --- the shared x-axis across Multiple query panels (2026-08-21, step 55) ---------------
#
# Every panel used to get its own axis, because metafor::forest() computes one from the
# CIs it is handed and forest_plot() passed neither alim= nor xlim=. For ONE panel that is
# right and is still what happens (xrange = NULL below). For the Multiple query tab, whose
# whole purpose is one gene across five tissues, it means five panels stacked down a column
# -- which read as one figure -- are five different rulers. Measured for FOXM1/DFS: ovarian
# drew [0.80, 1.05] and luad [0.00, 5.00], so the same horizontal position was HR 0.9 in one
# panel and HR 2.5 in the next, with nothing on screen saying so.
#
# This is NOT the hazard the step-52 note recorded. That note said the ranges "happen to
# land close together, so nothing looks wrong today". They do not: for ESR1/OS the narrowest
# panel spans 20% of the widest one's axis. The claim was wrong when it was written, and the
# numbers above are from the live engine.

# The pooled fit, in ONE place. .forest_bounds() has to know where the diamond will land
# before anything is drawn, and forest_plot() then draws it -- two rma() calls that must
# agree or the axis is computed for a diamond the figure does not contain. They cannot
# disagree if there is only one call. (Deliberately NOT .pool_rma() from the engine: that
# one retries with a damped Fisher step, this one does not, and the range has to describe
# the fit that is actually DRAWN. The divergence is pre-existing and recorded in HANDOFF.)
.forest_rma <- function(yi, vi) rma(yi = yi, vi = vi, method = "REML")

# Everything one forest DRAWS, as a range on the HR scale: every estimable cohort's CI, and
# the pooled diamond's when there are two or more. The pooled row is included because it is
# the row the tab exists to compare -- an axis that clipped the diamond would be arrowing
# the one estimate a reader came for.
#
# qnorm(0.975), not a literal 1.96, because that is the multiplier metafor actually draws
# with -- 1.959964 -- and this range has to describe the picture exactly.
#
# Note which way round that cuts, because the first version of this comment had it backwards
# and the test caught it: 1.96 is the LARGER of the two, so a range built from it errs a hair
# WIDE and can never clip. The reason to match is therefore not a safety margin. It is that a
# range which is only approximately the figure's is a number that cannot be checked against
# the figure -- and checking it against the figure is exactly what the invariant below does.
.forest_bounds <- function(res) {
  est <- Filter(function(x) isFALSE(x$skipped) && !is.null(x$logHR), res$per_cohort)
  if (!length(est)) return(NULL)
  yi <- vapply(est, `[[`, numeric(1), "logHR")
  se <- vapply(est, `[[`, numeric(1), "se")
  z  <- stats::qnorm(0.975)
  lo <- yi - z * se
  hi <- yi + z * se
  if (length(est) >= 2) {
    # If this fails, forest_plot() fails on the same panel for the same reason and draws
    # nothing -- so the diamond that is not there cannot need room. The cohorts still do.
    m <- tryCatch(.forest_rma(yi, se^2), error = function(e) NULL)
    if (!is.null(m)) { lo <- c(lo, m$ci.lb); hi <- c(hi, m$ci.ub) }
  }
  exp(c(min(lo), max(hi)))
}

# The union of what every panel of one Multiple query run draws, or NULL.
#
# NULL below two estimable panels is the honest answer, not a shortcut: with one panel there
# is no cross-panel comparison to protect, and forcing an axis would differ from metafor's
# own choice for no reason a reader could see. NULL is also what a caller passes to get
# today's behaviour, so the two meet.
mq_forest_xrange <- function(panels) {
  b <- Filter(Negate(is.null), lapply(panels, function(e)
    if (is.null(e$res)) NULL else .forest_bounds(e$res)))
  if (length(b) < 2) return(NULL)
  c(min(vapply(b, `[`, numeric(1), 1L)), max(vapply(b, `[`, numeric(1), 2L)))
}

# The axis for a shared range: metafor's OWN layout for it.
#
# The alternative was to compute alim/xlim/at here -- and that would have invented a second
# axis convention for this tool, one that five panels use and every other forest does not.
# Instead the range is handed to metafor on a throwaway device and the (alim, xlim, at) it
# chooses is read back and reused. The rule applied is metafor's, exactly as on every other
# forest in the app; the only thing that changed is the DATA the rule is applied to.
#
# Three measured facts make this safe, and none of them is obvious:
#
#   * the triple is DEVICE-INDEPENDENT. alim/xlim/at are user coordinates computed from the
#     range alone -- identical at 7.5x5, 9.375x6 and 4x9 -- so the on-screen panel and the
#     exported PDF get the same ruler, which is the entire point.
#   * a nested pdf()/dev.off() RESTORES the outer device and leaves it drawable, so this can
#     be called from inside forest_plot() with the real device already open.
#   * xlim cannot be skipped. metafor honours alim= but computes xlim= from the DATA, so a
#     shared alim wider than one panel's own CIs draws an axis that runs off the plot region
#     -- observed as a forest whose axis carried a single tick label and no line.
#
# The containment stop() is unreachable with metafor 4.x: 400 random ranges over [0.05, 40]
# all rounded OUTWARD. It is here because the failure it catches is silent -- a clipped alim
# turns real intervals into arrows -- and because "metafor rounds outward" is an assumption
# about someone else's code, not a fact about ours.
.forest_axis <- function(xr) {
  # The caller's device is restored BY NAME, not by trusting dev.off() to fall back to it.
  # dev.off() makes "the next device in the list" current, which is the caller's only when
  # the caller's is the only one open -- and under Shiny it is not: renderPlot holds its own
  # device while this runs. Measured with two devices open, dev.off() on the nested pdf left
  # device 2 current when the caller was device 3, so metafor drew into the wrong device and
  # every panel came back as "figure margins too large" or "attempt to select less than one
  # element in get1index". The first version of the test for this opened ONE device and so
  # could not see it.
  cur <- grDevices::dev.cur()
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f, width = FOREST_WIDTH_IN, height = 5)
  on.exit({
    grDevices::dev.off()
    if (cur != 1L) grDevices::dev.set(cur)   # 1 is the null device: nothing to go back to
    unlink(f)
  })
  op <- graphics::par(mar = c(5, 4, 3, 2))
  z <- metafor::forest(x = xr, ci.lb = c(xr[1], xr[1]), ci.ub = c(xr[2], xr[2]),
                       refline = 1, slab = c("a", "b"),
                       header = c("Cohort", FOREST_HR_HEADER))
  graphics::par(op)
  if (!(z$alim[1] <= xr[1] && z$alim[2] >= xr[2]))
    stop(sprintf(paste("metafor chose an axis [%g, %g] that does not contain the shared",
                       "range [%g, %g]; intervals would be drawn as arrows on the panel",
                       "that defined the range"), z$alim[1], z$alim[2], xr[1], xr[2]))
  list(alim = z$alim, xlim = z$xlim, at = z$at)
}

# Fail loud on a range that cannot be an HR axis. Every one of these would otherwise reach
# metafor as something it half-accepts: a reversed pair silently flips the axis, a zero or
# negative bound is not on the HR scale at all, and NA propagates into xlim and takes the
# whole figure out with a message about graphics parameters.
.check_xrange <- function(xr) {
  if (length(xr) != 2L || !is.numeric(xr) || any(!is.finite(xr)))
    stop("forest_plot(xrange=): need two finite numbers, got '",
         paste(utils::capture.output(str(xr)), collapse = " "), "'")
  if (any(xr <= 0))
    stop(sprintf("forest_plot(xrange=): HR bounds must be positive, got [%g, %g]", xr[1], xr[2]))
  if (!(xr[1] < xr[2]))
    stop(sprintf("forest_plot(xrange=): need xrange[1] < xrange[2], got [%g, %g]", xr[1], xr[2]))
  invisible(xr)
}

# xrange = NULL is today's figure, unchanged: metafor is called with no alim/xlim/at at all,
# not with NULL ones. That distinction is load-bearing -- forest.rma rejects alim = NULL with
# "Argument 'alim' must be of length 2", so the arguments have to be ABSENT, which is why the
# call below goes through do.call() and a list that is empty in the single-panel case.
forest_plot <- function(res, file = NULL, xrange = NULL) {
  .require_survresult(res, "forest_plot")
  est <- Filter(function(x) isFALSE(x$skipped) && !is.null(x$logHR), res$per_cohort)
  if (length(est) < 1) stop("no estimable cohort to plot")
  yi <- vapply(est, `[[`, numeric(1), "logHR")
  vi <- vapply(est, `[[`, numeric(1), "se")^2
  n  <- vapply(est, `[[`, numeric(1), "n")
  ev <- vapply(est, `[[`, numeric(1), "events")
  slab <- sprintf("%-8s  n=%d, ev=%d", names(est), n, ev)

  h <- forest_height_in(length(est))
  opened <- .open_dev(file, width = FOREST_WIDTH_IN, height = h)
  op <- graphics::par(mar = c(5, 4, 3, 2))
  feat <- attr(res, "feature") %||% "score"
  # BEFORE metafor draws, not after. par(mar=) takes effect at the next plot.new(), so
  # growing the top margin once the rows are on the page moves nothing and the extra title
  # line runs off the device -- which is the bug this whole step is about, reintroduced one
  # layer up. Measuring here is possible because strwidth(units = "inches") needs a device
  # and font metrics, not a plot region, and the device is already open.
  tl <- .title_lines(feat, sprintf(" vs %s", res$endpoint), w = .title_width(), cex = 1.1)
  # ONE on.exit, and par() before dev.off() -- the order is load-bearing. R runs on.exit
  # expressions in registration order, so two separate calls closed the device first and
  # then restored par() with nothing open. par() is not a passive query: setting it OPENS
  # the default device, which under Rscript is Rplots.pdf in the working directory. That
  # leaked a stray file and left its device open on every forest_plot(file = ...) call --
  # which is why deleting Rplots.pdf never stuck, the next driver run remade it.
  #
  # Restoring par first is correct in both modes: when we opened the device the restore is
  # moot (par dies with it), and when the caller owns the device (file = NULL, the app's
  # renderPlot) the restore is required and the device must be left open.
  on.exit({ graphics::par(op); if (opened) grDevices::dev.off() })
  xlab <- .forest_xlab(res)

  # EMPTY when xrange is NULL, so metafor is called exactly as it was before step 55 and
  # the single-panel figure is unchanged -- see the note on forest_plot()'s signature.
  ax   <- if (is.null(xrange)) list() else .forest_axis(.check_xrange(xrange))
  head <- c(FOREST_SLAB_HEADER, FOREST_HR_HEADER)

  if (length(est) >= 2) {
    m <- .forest_rma(yi, vi)
    mlines <- .pooled_mlab(sum(n), sum(ev), m$I2, m$pval, length(est))
    mlab   <- .stack_left(mlines)
    z <- do.call(metafor::forest, c(list(m, transf = exp, refline = 1, slab = slab,
                                         xlab = xlab, header = head, mlab = mlab), ax))
    z$mlab <- mlab; z$mlab_lines <- mlines
  } else {
    # Single cohort: no pooling — plot the lone estimate with its normal CI.
    z <- do.call(metafor::forest, c(list(x = yi, vi = vi, transf = exp, refline = 1,
                                         slab = slab, xlab = xlab, header = head), ax))
    z$mlab <- NULL; z$mlab_lines <- NULL
  }
  graphics::title(main = paste(tl, collapse = "\n"), cex.main = 1.1)
  # metafor::forest returns its resolved layout (xlim, alim, cex, ...) invisibly, and it is
  # carried out on the result for the same reason metafor bothers to return it: the label
  # column's width is a FRACTION of the device and the labels' width is in INCHES, so
  # whether they collide is a property of the drawn figure and cannot be checked from the
  # strings alone. tests/test_forest_label_fit.R reads this to measure the clearance
  # between the pooled label and the pooled diamond on a real render. Nothing in the app
  # consumes it; the plot is still drawn for its side effect and `res` is unchanged.
  #
  # The label that was ACTUALLY drawn rides along, and that detail is load-bearing: the
  # first version of that test rebuilt the label itself and so measured a string the figure
  # need not contain -- it passed unchanged against a deliberately reverted one-line label.
  # A geometry check has to read the geometry's own inputs or it is not a check.
  invisible(structure(res, forest_layout = z))
}

# --- follow-up horizon plumbing for the app ---------------------------------
# A selectInput value is a string, so the UI choice has to be converted before it can
# reach get_survival(). The "full" sentinel / "" / NULL / NA mean "full follow-up" and
# must become NULL -- the engine's no-op -- while anything unparseable raises rather than
# quietly falling back to full follow-up, which would show untruncated numbers under a UI
# claiming a horizon. (Same fail-loud principle as validate_cohorts().)
# "" is still accepted: it was the original full-follow-up value, and a bookmarked URL or
# a saved session can still carry it.
.tau_arg <- function(x) {
  if (is.null(x) || !length(x) || (length(x) == 1 && is.na(x))) return(NULL)
  if (is.character(x) && (!nzchar(x) || identical(x, TAU_FULL))) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v) || v <= 0) stop(sprintf("invalid follow-up horizon: '%s'", x))
  v
}

# Filename fragment for a truncated result. Without it, exporting one feature/endpoint
# at full follow-up and again at tau=60 yields two files with the SAME name -- the
# second silently overwrites the first, and the survivor cannot be told apart.
.tau_tag <- function(res) {
  if (is.null(res$max_followup)) "" else sprintf("_tau%g", res$max_followup)
}

# Any string, made safe to sit in one path segment. The single definition the whole app
# shares: mgene_export_name() sanitises a gene LIST with it, .feature_tag() one symbol,
# and the two must not be able to disagree about what "safe" means (the reason
# .mgene_query_fields() was extracted rather than copied, step 33 stage 4).
.safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# How many characters of a filename ONE feature symbol may claim. Not a filesystem limit --
# the longest wrapper this app builds around a feature is 41 bytes (measured over the
# registry: `km_<ct>_<feature>_<cohort>_<endpoint>_tau<n>.pdf` at its longest ct, cohort and
# endpoint), so NAME_MAX 255 would allow 214. 80 is a readability budget, the same one
# MGENE_NAME_BUDGET applies to a joined gene list; a 214-character filename is legal and
# unusable. Kept as its own constant because the two budget DIFFERENT things -- a list that
# can shed members versus a single symbol that cannot.
FEATURE_NAME_BUDGET <- 80L

# One feature symbol as a filename fragment.
#
# This exists because of what the app's own vocabularies actually contain, measured against
# these matrices (2026-08-09, re-derived 2026-08-10, re-derived again over seven tissues
# 2026-08-26, figures re-checked 2026-08-27, re-derived again 2026-08-30):
# of the 59964 features the selectors offer across seven tissues x two query kinds -- 58366
# DISTINCT strings, since a symbol offered under both
# viper and expr is counted twice in 59964 and once in 58366, and every figure in this paragraph
# is over the distinct set -- **5837 contain "/"** -- probe sets mapping to several genes, e.g.
# "C2///CFB" -- and 363 are longer than 215 bytes, the longest 1466. The longest feature that
# does NOT contain a "/" is 40 characters, so a sanitiser that also caps length covers both
# problems and neither can occur without the other. Adding brain moved the offer from 45956
# entries to 60041 and changed exactly one of these figures (5836 -> 5837): the design constants
# this paragraph justifies are unaffected, which is the reason to re-derive rather than assume.
# Withdrawing immune from the selectors (step 108) moved the offer again, 60041 -> 59964 and the
# distinct set 58443 -> 58366, and changed NONE of the figures above: the 77 cell-type scores it
# took away are short and slash-free, so they were never in any of the sets being counted here.
#
# What went wrong without it, measured end to end on 2026-08-10 (BUILD_LOG step 35): the
# Multiple query ZIP put such a name through file.path(), so pdf() was handed a path whose
# dirname does not exist, forest_plot() died as "no usable PDF device", zip() never ran, and
# the user got NO FILE AND NO MESSAGE -- the button simply did nothing. Not a cosmetic defect.
# The Single query buttons failed one step further out and just as quietly (step 37 uses this
# too): the server served 200 and the BROWSER dropped the 353 over-NAME_MAX names on the floor,
# no file in ~/Downloads at all. (353 is what that EXPERIMENT observed, not a count of the
# vocabulary -- do not read it against the 363 above and conclude one of them drifted. The
# vocabulary figure is 363 and has not moved since 2026-08-10; 353 is how many of those names
# were exercised and failed to land. Checked 2026-08-27, when 353 turned out to be reproducible
# under no reading of the vocabulary at all, which is the point: it never was one.) Two different layers refusing, one cap fixing both.
#
# Over budget the name is truncated and SAYS SO: "__trunc<n>", n being the length of the
# SANITISED name (not the raw symbol -- "///" collapses to "_", so the two differ, and the
# sanitised length is the one the collision counts below were measured on). The marker is
# the point. A bare clipped prefix NAMES A FEATURE THAT WAS NOT RUN,
# which is what mgene_export_name() refuses when it emits "7genes" rather than a clipped list;
# "__trunc253" cannot be read as a symbol by anyone, and the full name is inside the file on
# every row anyway. The length is there because it is a free discriminator: across the whole
# vocabulary the prefix alone collides 19 times among the 708 names this truncates, the prefix
# plus length 5 times. Truncation cannot be made injective -- that is pigeonhole, not a bug.
#
# Do NOT read that residual as harmless on the grounds that a browser will disambiguate a
# repeated name. It may not: on 2026-08-10, in the course of this very work, Chrome silently
# OVERWROTE a same-named CSV in ~/Downloads (1902 bytes replaced by 6) rather than adding a
# suffix, while other files in that same folder carry `-1`/`-2` suffixes from other browsers.
# The behaviour is not something to lean on. What makes the residual acceptable is not the
# browser: it is that the authoritative statement of which feature a file came from is INSIDE
# the file, never in its name -- the `feature` column of survtable() / mq_survtable(), which
# carries the symbol whole and unsanitised.
#
# NOT the figure titles, which is what this comment claimed until 2026-08-10 and was wrong for
# exactly the names that need it. Once step 37 made these exports reachable at all, the PDFs
# could be read for the first time: for the 253-character feature the forest title recovered as
# a 74-character MIDDLE SLICE (`pdftotext` on the delivered file) and the KM's as 63 -- the
# device clips a centred title at the page box, so both ends went, including the `(expr) vs OS`
# suffix. Step 38 fixed that (.title_lines(): the title now wraps to the device and marks any
# residual cut), so a figure DOES name its feature again -- but only within two lines plus a
# `[+n chars]` marker, which for the 743 longest features of 58366 is still not the whole
# symbol. So the sentence above stands as written: the authoritative statement of which
# feature a file came from is the `feature` column INSIDE it, not its name and not its title.
.feature_tag <- function(x) {
  s <- .safe_name(as.character(x)[1])
  if (nchar(s) <= FEATURE_NAME_BUDGET) return(s)
  sprintf("%s__trunc%d", substr(s, 1, FEATURE_NAME_BUDGET), nchar(s))
}

# Follow-up-horizon choices for the app's selector, for one cancer type.
#
# Returns list(choices = <named character vector for selectInput>, selected = <default>).
# Values are strings because that is what Shiny inputs carry; "" means full follow-up.
# Every value round-trips through .tau_arg().
#
# A tissue with a DECLARED panel horizon (horizon_for()) defaults to it, rather than to
# full follow-up. Without this the app contradicts its own artifacts: since Phase 10 every
# stored ovarian number is fit at tau=48, so an ovarian query defaulting to full follow-up
# returns a different HR than the CSV the same app serves in its Browse tab.
#
# The stock presets stay reachable, including full follow-up -- the declared horizon is a
# default, not a cage. It is labelled "(panel default)" so the selection reads as a
# standard rather than as something arbitrary left in the box.
# Full follow-up carries an explicit sentinel value, NOT "". Shiny renders this as a
# selectize widget and selectize DROPS options whose value is the empty string -- it
# reserves "" for its own placeholder state. While "" was also the default that was
# harmless (the empty state meant the same thing), but the moment a tissue defaults to a
# real horizon the option vanishes from the list and full follow-up becomes unreachable
# for that tissue. Found by driving the running app; no R-level assertion could see it.
TAU_FULL <- "full"
.TAU_PRESETS <- c("60", "120", "180")

tau_choices <- function(cancer_type, endpoint = NULL) {
  # Endpoint-AWARE: when the caller knows the query endpoint (the app always does), the
  # default is THAT endpoint's declared horizon, so the selector matches the pooled artifact
  # the app serves for it -- LUAD DFS @ tau=99, LUSC OS @ tau=87, not the conservative min.
  # With endpoint = NULL it falls back to the bare (shortest per-endpoint) horizon, so any
  # legacy endpoint-agnostic caller is unchanged. horizon_for() returns NULL for an endpoint
  # a tissue does not declare, which correctly lands that combination on full follow-up.
  tau <- horizon_for(cancer_type, endpoint)
  vals <- c(TAU_FULL, .TAU_PRESETS)
  labs <- c("Full follow-up", sprintf("%s months", .TAU_PRESETS))
  if (!is.null(tau)) {
    tv <- format(tau, trim = TRUE)
    keep <- vals != tv                      # never offer the same horizon twice
    vals <- c(vals[keep], tv)
    labs <- c(labs[keep], sprintf("%s months (panel default)", tv))
    # numeric order, with the non-numeric sentinel pinned first
    ord <- order(suppressWarnings(as.numeric(vals)), na.last = FALSE)
    vals <- vals[ord]; labs <- labs[ord]
  }
  list(choices = setNames(vals, labs),
       selected = if (is.null(tau)) TAU_FULL else format(tau, trim = TRUE))
}

# Tag a survresult with the feature name so plot titles/axes read nicely:
#   res <- with_feature(get_survival(get_feature("ESR1", kind="viper")), "ESR1 (VIPER)")
with_feature <- function(res, feature) { attr(res, "feature") <- feature; res }

# Tidy per-cohort + pooled HR/CI table for a survresult - the exact numbers behind the
# forest (so the "TCGA crosses 1, pooled doesn't" story is documented, not eyeballed).
# Per-cohort CI is the Cox Wald CI (logHR +/- 1.96*se); pooled row is the REML diamond.
# The survival figures take a `survresult` and nothing else.
#
# Added 2026-08-11 with R/tumor_normal.R, whose header claims its output can never reach a
# survival function. It could not, but only by ACCIDENT: a tumor_normal object has a
# $per_cohort too, so survtable() walked it and died on "arguments imply differing number
# of rows" and forest_plot() on "no estimable cohort" -- incidental shape errors that name
# nothing, and that a small change to either object's shape could turn into a garbage
# table instead of a stop. A claimed boundary that rests on an accident is not a boundary.
.require_survresult <- function(res, fn) {
  if (!inherits(res, "survresult"))
    stop(fn, "(): expected a survresult from get_survival(), got ",
         paste(class(res), collapse = "/"),
         ". Survival tables and figures describe outcome; this object does not.")
  invisible(TRUE)
}

survtable <- function(res) {
  .require_survresult(res, "survtable")
  feat <- attr(res, "feature") %||% "score"
  rows <- list()
  for (coh in names(res$per_cohort)) {
    r <- res$per_cohort[[coh]]
    if (isTRUE(r$skipped) || is.null(r$logHR)) {
      rows[[coh]] <- data.frame(feature = feat, endpoint = res$endpoint, cohort = coh,
        n = r$n %||% NA, events = r$events %||% NA, HR = NA, ci_lb = NA, ci_ub = NA,
        p = NA, logrank_p = NA, ph_p = NA_real_, ph_violated = NA,
        strata_adj = NA,
        pool_ci_identified = NA, pool_strata_k = NA_integer_,
        HR_fe = NA_real_, p_fe = NA_real_,
        stringsAsFactors = FALSE)
      next
    }
    rows[[coh]] <- data.frame(feature = feat, endpoint = res$endpoint, cohort = coh,
      n = r$n, events = r$events, HR = r$HR,
      ci_lb = exp(r$logHR - 1.96 * r$se), ci_ub = exp(r$logHR + 1.96 * r$se),
      p = r$p, logrank_p = r$logrank_p,
      # The PH verdict travels WITH the HR. It was previously reported only by print()
      # and the KM annotation, so it disappeared the moment a result was tabulated or
      # exported -- leaving a time-averaged HR looking identical to a constant one in
      # the artifact most likely to be read later.
      ph_p = r$ph_p, ph_violated = r$ph_violated,
      strata_adj = r$strata_adjusted,
      # Pooling-rule disclosure is a property of the POOLED row only; per-cohort rows are
      # single Cox fits with no between-study variance in them. NA here, never a value
      # carried down from the pooled row -- the same reason ph_p is NA on the pooled row.
      pool_ci_identified = NA, pool_strata_k = NA_integer_,
      HR_fe = NA_real_, p_fe = NA_real_,
      stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  if (!is.null(res$pooled)) {
    p <- res$pooled
    tab <- rbind(tab, data.frame(feature = feat, endpoint = res$endpoint,
      cohort = sprintf("POOLED (k=%d, I2=%s%%)", p$k,
                       ifelse(is.na(p$I2), "NA", sprintf("%.0f", p$I2))),
      n = NA, events = NA, HR = p$HR, ci_lb = p$ci_lb, ci_ub = p$ci_ub,
      # No PH test exists for a pooled estimate -- it is a weighted average of per-cohort
      # coefficients, not a fit with residuals. NA, never a value carried up from a row.
      p = p$p, logrank_p = NA, ph_p = NA_real_, ph_violated = NA,
      strata_adj = NA,
      # FALSE = the CI on this row assumes a between-study variance that k cohorts cannot
      # identify (k < POOL_MIN_K). HR_fe/p_fe are the fixed-effect fit on the same inputs,
      # so a reader can see whether the finding survives the other pooling rule. They are
      # equal to HR/p exactly when tau^2 = 0, i.e. when the choice was moot.
      pool_ci_identified = isTRUE(p$tau2_identified),
      # How many of this pool's k cohorts carried a strata() term. Pooled-row-only, the same
      # shape as pool_ci_identified: a per-cohort row is one fit and already says so in
      # strata_adj. Read it against the k in the cohort label -- anything between 1 and k-1
      # means the pooled HR averages stratified and unstratified estimates. NA_integer_ when
      # the result predates the field (a hand-built pooled list in a test), never 0, because
      # "no cohort was stratified" and "unknown" are different statements.
      pool_strata_k = if (is.null(p$strata_k)) NA_integer_ else as.integer(p$strata_k),
      # Explicit NULL guards, not %||% -- that helper tests nzchar() and is string-shaped,
      # so it misreads a numeric (and would error on a length-0 numeric).
      HR_fe = if (is.null(p$HR_fe)) NA_real_ else p$HR_fe,
      p_fe  = if (is.null(p$p_fe))  NA_real_ else p$p_fe,
      stringsAsFactors = FALSE))
  }
  # Follow-up horizon, NA when untruncated. Constant down the column by construction --
  # it is a property of the whole result, not of a row -- and it earns the space anyway:
  # this table gets written to results/survtable_*.csv, where it outlives the session
  # that produced it. Without the column a tau=60 export is permanently
  # indistinguishable from a full-follow-up one, and those are different estimands.
  tab$max_fu <- if (is.null(res$max_followup)) NA_real_ else res$max_followup
  rownames(tab) <- NULL
  tab
}

# The `stratified` column of survtable_display(). Per-cohort rows are one fit each and read
# yes/no. The POOLED row is not a fit -- it has no strata term of its own, which is why
# survtable() leaves its strata_adj at NA -- but it does have a COMPOSITION: how many of the
# k cohorts behind the diamond carried one. Rendering that as "2 of 3" puts the mixture in
# the on-screen table, which is the one place a reader meets the pooled HR without also
# reading print() or the CSV. Blank when no cohort in the pool was stratified, because then
# there is no mixture to warn about and "0 of 3" would read as a defect rather than a design.
.strat_col <- function(tab, res) {
  out <- ifelse(is.na(tab$strata_adj), "", ifelse(tab$strata_adj, "yes", "no"))
  pl  <- grepl("^POOLED", tab$cohort)
  p   <- res$pooled
  if (any(pl) && !is.null(p) && !is.null(p$strata_k) && p$strata_k > 0)
    out[pl] <- sprintf("%d of %d", p$strata_k, p$k)
  out
}

# On-screen view of survtable(). Same rows, fewer columns, readable headers.
#
# The export and the display want opposite things. survtable() repeats feature /
# endpoint / max_fu on every row because the CSV lands in results/ and has to be
# self-describing once nothing else is attached to it. The app already shows all three
# in the sidebar, so on screen they are three constant columns of noise -- and
# ci_lb/ci_ub are two columns saying one thing. Dropping them here is what makes the
# PH column fit without the table becoming unreadable.
#
# tau is the one dropped constant that still has to reach the reader, so it moves to a
# caption attribute rather than disappearing: a truncated table of HRs that does not
# say so is the same failure the plot layer already guards against.
survtable_display <- function(res) {
  tab <- survtable(res)
  d <- data.frame(
    cohort   = tab$cohort,
    n        = tab$n,
    events   = tab$events,
    HR       = ifelse(is.na(tab$HR), "", sprintf("%.2f", tab$HR)),
    `95% CI` = ifelse(is.na(tab$ci_lb), "",
                      sprintf("%.2f - %.2f", tab$ci_lb, tab$ci_ub)),
    p        = fmt_p(tab$p),
    `log-rank p` = fmt_p(tab$logrank_p),
    # "!" carries the warning from print() and the KM annotation into the table. Without
    # it a time-averaged HR reads exactly like a constant one.
    `PH p`   = ifelse(is.na(tab$ph_p), "",
                      paste0(fmt_p(tab$ph_p),
                             ifelse(isTRUE_vec(tab$ph_violated), "  !", ""))),
    stratified = .strat_col(tab, res),
    check.names = FALSE, stringsAsFactors = FALSE)
  # HIDDEN SORT KEYS. Four of the columns above are FORMATTED STRINGS -- they have to be,
  # because "" is what a skipped cohort's HR is and no numeric column can hold that, and
  # because the PH column carries a "!" -- and DataTables therefore sorts them LEXICALLY.
  # Observed on the Multiple-genes table in step 33 and identical here: ascending by `p`
  # puts 4.48e-07 BELOW 0.844, because the string "4.48e-07" sorts after "0.8". A reader
  # clicking `p` to find the strongest cohort gets the reverse of the answer with nothing
  # on screen to contradict it.
  #
  # The key is a RANK, not the value, and not a sentinel. A sentinel has to be "above
  # everything in this column", which p-values make easy (2) and HRs do not -- HR is
  # positive and unbounded, so any constant is a guess about the data. Rank sidesteps
  # that: `na.last = TRUE` gives the blanks the largest ranks, so they leave the numbers
  # in one block at one end instead of interleaving. Inf is not usable here for the
  # separate reason that it does not survive JSON -- jsonlite writes it as the STRING
  # "Inf", which would turn the column back into a lexical sort.
  keyed <- list(`HR` = tab$HR, `p` = tab$p,
                `log-rank p` = tab$logrank_p, `PH p` = tab$ph_p)
  knames <- paste0(SURVTABLE_SORT_PREFIX, seq_along(keyed))
  for (i in seq_along(keyed))
    d[[knames[i]]] <- rank(keyed[[i]], na.last = TRUE, ties.method = "first")
  # The POOLED row is a SUMMARY of the rows above it, not a fourteenth cohort, and its
  # meaning is carried by its position. Sorting by any column drops it among the cohorts,
  # where it reads as one of them. This column pins it back to the bottom via orderFixed;
  # it is derived from the label survtable() writes, never from row position.
  d[[SURVTABLE_POOLED_KEY]] <- as.integer(grepl("^POOLED", tab$cohort))
  attr(d, "sort_keys") <- setNames(knames, names(keyed))
  if (!is.null(res$max_followup))
    attr(d, "caption") <- sprintf(
      "Follow-up truncated at %g months; patients still at risk are censored there. '!' = proportional hazards rejected (p < %g).",
      res$max_followup, PH_ALPHA)
  d
}

# Column names the app hides. Prefixed rather than listed so a fifth keyed column cannot
# be added without survtable_dt_options() seeing it.
SURVTABLE_SORT_PREFIX <- ".sortkey"
SURVTABLE_POOLED_KEY  <- ".sortkey_pooled"

# The one display column that gets NO ordering at all. "0.82 - 1.05" is a RANGE: there is
# no single quantity to sort it by, and three readers would expect three different ones
# (lower bound, upper bound, width). Sorting it lexically is wrong; picking one of the
# three silently and calling it "sorted by CI" is also wrong. HR is the column that
# orders this information, and it is already next to it.
SURVTABLE_NO_ORDER <- "95% CI"

# The DataTables `options` fragment that makes survtable_display()'s table sortable
# CORRECTLY. Lives here rather than inline in app.R for the reason set in Phase 7c: app.R
# is Shiny with no test harness, and the whole content of this function is a set of COLUMN
# INDICES that are wrong silently -- an off-by-one points `p` at the log-rank key and the
# table still sorts, just by the wrong number.
#
# Indices are 0-based over the data frame's columns, which is what `targets` means when
# the table is drawn with rownames = FALSE. Both call sites pass rownames = FALSE; the
# assertion for that is in tests/test_survtable_display.R, because it cannot be seen here.
survtable_dt_options <- function(d) {
  keys <- attr(d, "sort_keys")
  if (is.null(keys) || !length(keys))
    stop("survtable_dt_options(): this table carries no sort keys. It did not come from ",
         "survtable_display(), and sorting it would sort formatted strings lexically.")
  idx <- function(nm) {
    j <- match(nm, names(d))
    if (is.na(j)) stop("survtable_dt_options(): no column '", nm, "' in the display table")
    j - 1L
  }
  hidden <- c(unname(keys), SURVTABLE_POOLED_KEY)
  list(
    columnDefs = c(
      lapply(hidden, function(k)
        list(targets = idx(k), visible = FALSE, searchable = FALSE)),
      lapply(seq_along(keys), function(i)
        list(targets = idx(names(keys)[i]), orderData = idx(keys[[i]]))),
      list(list(targets = idx(SURVTABLE_NO_ORDER), orderable = FALSE))),
    # `pre`, not `post`: pre-ordering is applied BEFORE the user's, so the cohort rows
    # (0) always precede the pooled row (1) whichever column is clicked. `post` only
    # breaks ties and would leave the pooled row floating.
    orderFixed = list(pre = list(idx(SURVTABLE_POOLED_KEY), "asc")))
}

# isTRUE() is scalar-only; this is the vectorised form for a logical column that may
# hold NA (skipped cohorts and the pooled row both do).
isTRUE_vec <- function(x) !is.na(x) & x

# --- why a cohort is not in this query: the app's two absence notes ----------
# Both are pure text over the registry and one result, and both live here rather than in
# app.R for the reason set in Phase 7c for .tau_arg/.tau_tag: app.R is Shiny with no test
# harness, so anything inline there is unverified. These two are the sentences that stand
# in for a missing curve, and a wrong one is READ AS A FACT ABOUT THE DATA -- there is
# nothing else on screen to contradict it. Unexported: the app sources R/ directly.

# Why a requested cohort has no curve for this query. An endpoint the cohort simply does
# not carry is a structural fact worth naming explicitly (SCAN-B is OS-only) rather than
# a generic "no data" message.
#
# The three branches are ordered by how much they know, most-informed first: the registry
# fact outranks the engine's reason because "carries no DFS data" explains WHY there were
# too few events, while "0 events" alone invites the reader to conclude the biology is
# null. The last line is the honest fallback for a cohort the result never mentions.
.cohort_note <- function(co, r, registry = COHORTS) {
  eps <- cohort_endpoints(co, registry)
  if (length(eps) && !(r$endpoint %in% eps))
    return(sprintf("%s carries no %s data (endpoints: %s).", co, r$endpoint,
                   paste(eps, collapse = ", ")))
  ent <- r$per_cohort[[co]]
  if (!is.null(ent) && isTRUE(ent$skipped)) return(sprintf("%s: %s.", co, ent$reason))
  sprintf("%s: no usable data for this query.", co)
}

# The line under the cohort chips. It states the REASON and a COUNT only -- the cohort
# names are already on screen as struck-through chips, and repeating twelve of them in
# prose was the wall this presentation exists to remove. Grouped by what each one does
# carry, read from the registry, never hand-listed.
.excluded_reason_line <- function(dead, endpoint, registry = COHORTS) {
  carries <- vapply(dead, function(co) {
    e <- cohort_endpoints(co, registry)
    if (!length(e)) "none" else paste(e, collapse = "/")
  }, character(1))
  n <- table(carries)
  parts <- vapply(names(n), function(k) {
    v <- as.integer(n[[k]])
    if (identical(k, "none")) sprintf("%d declare no endpoint", v)
    else sprintf("%d carr%s %s only", v, if (v == 1) "ies" else "y", k)
  }, character(1))
  # ASCII dash, not an em dash: the project copy of this file is under a non-ASCII string
  # literal ban (tests/test_pdf_fonts.R) and the two copies must stay character-identical.
  sprintf("Struck through: not offered for %s - %s.", endpoint,
          paste(parts, collapse = ", "))
}

# %||% is defined above (used here too, before its file position is irrelevant at call time).

# --- the Multiple query tab's panels: disclosures, table, export -------------
# A `panels` list is what app.R's mq_res() returns: one entry per cancer type, each a
# COMPLETE account of what happened to that panel, including the failures. The fields used
# here are ct / feature / endpoint / tau / n_sel / status / msg / res, and `status` is one
# of "ok", "absent" (structural: no cohort carries the endpoint, or the feature is not
# measured), "empty" (measured, nothing cleared min_events) or "error".

# Estimable cohorts in a result -- a per-QUERY count, not a tissue constant, because a
# small regulon drops out of individual cohorts.
.n_estimable <- function(r)
  sum(vapply(r$per_cohort, function(x) isFALSE(x$skipped), logical(1)))

# ---- disclosures shared by the two multi-result tabs --------------------------------
#
# The Multiple query tab (one feature, five tissues) and the Multiple-gene tab (up to ten
# genes, one tissue) have to say three IDENTICAL things about any one result: that a lone
# cohort is not a meta-analysis, which of the tissue's cohorts are not in it and why, and
# where the feature sits in the genome-wide scan. Factored out when the second tab arrived
# (step 33) rather than copied, for the reason mq_panel_notes' own comment already gives
# about the tau^2 caveat -- "a second copy written for the export would be free to drift
# from the one the reader saw". Here the drift would be worse than cosmetic: the rank
# sentence carries an FDR verdict, and one tab calling a gene a hit while the other calls
# it a miss is a contradiction the reader cannot resolve.
#
# What is NOT shared is the pooling sentence, because its wording is genuinely different
# ("nothing is pooled with any other PANEL" vs "with any other GENE") and a single
# parameterised string would say neither well.

MSG_SINGLE_COHORT <-
  "Single cohort: this is one Cox fit, not a meta-analysis - no heterogeneity is estimable."

# The tissue's cohorts that this result did NOT estimate, each with its own reason.
# Returns a `warn`-named string, or character(0) when every cohort made it in -- so a
# caller can c() it unconditionally.
.missing_cohort_note <- function(e) {
  r <- e$res
  usable  <- names(Filter(function(x) isFALSE(x$skipped), r$per_cohort))
  missing <- setdiff(cohorts_for(e$ct, e$endpoint), usable)
  if (!length(missing)) return(character(0))
  c(warn = paste(vapply(missing, .cohort_note, character(1), r = r), collapse = " "))
}

# The multiplicity anchor for one result, as a severity-named sentence.
#
# `suppress` names the lookup statuses that must NOT be printed per result because on the
# calling tab they are a property of the QUERY and identical on every row -- printing one
# would repeat a single sentence N times on screen and in N CSV rows. The caller states
# those once instead.
#
#   Multiple query   suppresses "kind" only. It varies one TISSUE, so no_scan and
#                    not_found are genuinely per-panel facts (luad has no OS scan).
#   Multiple-gene    suppresses every status except ok/not_found: the tissue, endpoint,
#                    cohorts and horizon are all fixed across the query there, so each of
#                    those verdicts is the same on all ten rows.
#
# Severity names map onto the same three tones both tabs already use, so one badge does not
# change colour depending on which tab it is read from.
.rank_note <- function(e, suppress = "kind") {
  if (is.null(e$rank) || e$rank$status %in% suppress) return(character(0))
  f <- format_scan_rank(e$rank)
  if (is.null(f)) return(character(0))
  setNames(f$text, switch(f$tone, hit = "hit", miss = "warn", "info"))
}

# Everything one panel has to disclose, as sentences. Returned as a character vector NAMED
# by severity ("info" / "warn") so a caller can style it without knowing what the sentences
# say; the export pastes them into one cell and the app wraps each in its own <p>.
#
# ONE source, two consumers, and that is the point. These sentences are not decoration --
# the endpoint is global across panels so they stay comparable, which means for most
# tissues it is NOT that tissue's primary analysis, and a table of five HRs that does not
# say so invites exactly the reading the panel note exists to prevent. An export that
# dropped them would be a lossier artifact than the screen it came from; a second copy
# written for the export would be free to drift from the one the reader saw.
mq_panel_notes <- function(e) {
  if (!identical(e$status, "ok")) return(c(warn = e$msg))
  r <- e$res
  # NOT `%||%`: that operator tests nzchar(), a string test, and errors on a length-0
  # numeric. An explicit is.null test is the thing meant here.
  k <- if (!is.null(r$pooled) && !is.null(r$pooled$k)) r$pooled$k else .n_estimable(r)
  out <- c(info = sprintf(
    "Pooled over k=%d of %d %s cohorts in this panel. Nothing is pooled with any other panel.",
    k, e$n_sel, e$endpoint))
  # The horizon is read off the RESULT, never off the selector: what is stated is what was
  # fit. Both branches are explicit -- an undeclared horizon has to SAY "full follow-up",
  # because silence would read as a truncated fit.
  out <- c(out, info = if (is.null(r$max_followup))
    sprintf("Full follow-up: no panel horizon is declared for %s %s.",
            cancer_label_of(e$ct), e$endpoint)
  else
    sprintf("Horizon %g months (this panel's declared default for %s). Patients still at risk are censored there, not excluded.",
            r$max_followup, e$endpoint))
  # k=1 draws a lone estimate with no diamond, so the forest carries no pooling caveat of
  # its own. (k=2..4 does: forest_plot writes the tau^2 line into the summary label, so
  # repeating it here would be a second copy that can drift.)
  if (k < 2)
    out <- c(out, warn = MSG_SINGLE_COHORT)
  pep <- primary_endpoint_of(e$ct)
  if (!identical(pep, e$endpoint))
    out <- c(out, warn = sprintf("%s's primary endpoint is %s; %s is a secondary analysis here.",
                                 cancer_label_of(e$ct), pep, e$endpoint))
  out <- c(out, .missing_cohort_note(e))
  # The multiplicity anchor, when the caller resolved one (e$rank, a scan_rank_lookup
  # result). Severity is the lookup's own verdict, which is why this adds a THIRD name: a
  # rank that clears FDR is not a caveat and must not be styled as one, and it is not
  # neutral either. hit / warn / info map onto the app's green / orange / grey.
  #
  # `kind` is dropped deliberately: it says the query is not VIPER activity, which is a
  # property of the QUERY and identical on every panel, so printing it per panel would
  # repeat one sentence across all of them. Every other status is panel-specific and kept.
  out <- c(out, .rank_note(e))
  out
}

# One row of survtable()'s shape for a panel that produced no result at all. Its column
# names are not asserted here because rbind() below enforces them the moment ANY panel did
# produce a table -- and tests/test_mq_export.R pins the two shapes against each other for
# the case where none did.
.mq_null_row <- function(e) data.frame(
  feature = e$feature, endpoint = e$endpoint, cohort = NA_character_,
  n = NA_integer_, events = NA_integer_, HR = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
  p = NA_real_, logrank_p = NA_real_, ph_p = NA_real_, ph_violated = NA, strata_adj = NA,
  pool_ci_identified = NA, pool_strata_k = NA_integer_, HR_fe = NA_real_, p_fe = NA_real_,
  max_fu = if (is.null(e$tau)) NA_real_ else as.numeric(e$tau),
  stringsAsFactors = FALSE)

# Every panel of a Multiple query as ONE long table -- including the panels that produced
# nothing.
#
# That last clause is the whole design. A panel that quietly vanished from the file would
# read as "no effect" when it actually means "not measured here" or "the call raised", and
# a reader with only the CSV has nothing on screen to correct them. So `absent`, `empty`
# and `error` panels each contribute a row carrying their status and their reason, and the
# row count is the reader's guarantee that all five panels are accounted for.
#
# `cancer_type` is what makes the rows long rather than merged. Nothing is pooled across
# panels here any more than on screen; the tissues sit in one file, not in one model.
mq_survtable <- function(panels) {
  stopifnot(length(panels) > 0)
  do.call(rbind, lapply(panels, function(e) {
    tab <- if (is.null(e$res)) .mq_null_row(e) else survtable(e$res)
    cbind(data.frame(cancer_type = e$ct, status = e$status,
                     # Declared, and machine-readable rather than only in the prose:
                     # comparing it to `endpoint` is how a reader sees which panels are
                     # reporting a secondary analysis.
                     primary_endpoint = primary_endpoint_of(e$ct),
                     panel_note = paste(mq_panel_notes(e), collapse = " "),
                     stringsAsFactors = FALSE),
          tab, stringsAsFactors = FALSE)
  }))
}

# Write the whole Multiple query to `dir`: the long table, plus one forest per estimated
# panel and one KM per estimated cohort. Returns the paths written, relative to `dir`.
#
# Filenames reuse the Single query tab's convention exactly -- including .tau_tag(), which
# exists because two exports of one feature/endpoint that differ only by horizon otherwise
# collide and the survivor cannot be told apart. Here the horizon differs BETWEEN PANELS by
# construction (each tissue has its own declared tau), so the tag is doing real work in
# every zip, not guarding a rare case.
#
# Panels with nothing estimable get no figure and are not silently absent from the export:
# they are in the CSV, with their status and their reason.
mq_export_dir <- function(panels, dir) {
  stopifnot(dir.exists(dir), length(panels) > 0)
  # Read the feature off the PANELS, not from a caller-supplied string: every panel of one
  # run answers the same query, so disagreement here means the list was assembled from two
  # runs and the zip would be labelled with one of them.
  feature <- unique(vapply(panels, function(e) e$feature_id, character(1)))
  if (length(feature) != 1)
    stop(sprintf("panels disagree about the feature (%s) -- this is not one query",
                 paste(feature, collapse = ", ")))
  # SANITISED before it reaches file.path(). Until 2026-08-10 the raw symbol went straight in,
  # and for the 5837 features containing "/" that made pdf() open a path under a directory
  # that does not exist: the whole export died as "no usable PDF device" with nothing on
  # screen. See .feature_tag(). The tag is used for FILENAMES ONLY: mq_survtable()'s `feature`
  # column still carries the true symbol whole, and that CSV is in this same zip, so within an
  # archive the sanitising costs the name's fidelity and nothing else. (An earlier version of
  # this comment also named the forest and KM titles as carrying it. They do not for the names
  # that need it -- the device clips a long centred title at the page box. Measured 2026-08-10,
  # step 37, once the export was reachable at all. .feature_tag() records the numbers.)
  ftag <- .feature_tag(feature)
  # The SAME shared axis the screen used, from the same function on the same panels (step
  # 55). Not because the number is expensive, but because a zip whose forests are on five
  # private rulers while the tab's are on one shared ruler is a file that disagrees with
  # the screen it came from -- the failure this whole export is guarded against.
  xr <- mq_forest_xrange(panels)
  out <- "survtable_all_panels.csv"
  utils::write.csv(mq_survtable(panels), file.path(dir, out), row.names = FALSE)
  for (e in panels) {
    if (is.null(e$res) || .n_estimable(e$res) < 1) next
    tag <- .tau_tag(e$res)
    f <- sprintf("forest_%s_%s_%s%s.pdf", e$ct, ftag, e$endpoint, tag)
    forest_plot(e$res, file = file.path(dir, f), xrange = xr)
    out <- c(out, f)
    for (co in names(Filter(function(x) isFALSE(x$skipped), e$res$per_cohort))) {
      f <- sprintf("km_%s_%s_%s_%s%s.pdf", e$ct, ftag, co, e$endpoint, tag)
      km_plot(e$res, co, file = file.path(dir, f))
      out <- c(out, f)
    }
  }
  out
}

# The whole Multiple query as one zip at `zipfile`.
#
# In R/ rather than in app.R's downloadHandler for the reason the rest of this file exists:
# app.R has no test harness. The two fiddly parts are both here where they can be checked --
# zip() stores the paths it is HANDED, so it has to run from inside the staging directory or
# every entry unpacks under a tmp/... tree; and a non-zero status leaves a truncated or
# absent archive, which a download would otherwise hand over looking exactly like a good one.
mq_export_zip <- function(panels, zipfile) {
  zipfile <- normalizePath(zipfile, mustWork = FALSE)
  d <- file.path(tempdir(), sprintf("mq_export_%s", as.integer(Sys.time())))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  files <- mq_export_dir(panels, d)
  owd <- setwd(d); on.exit(setwd(owd), add = TRUE)
  st <- utils::zip(zipfile, files, flags = "-q9X")
  if (!identical(st, 0L) || !file.exists(zipfile))
    stop(sprintf("zip failed (status %s) while writing %d files", st, length(files)))
  files
}

# %||% is defined above (used here too, before its file position is irrelevant at call time).

# ---- Multiple-gene query: input resolution -----------------------------------------
#
# The Multiple query tab runs ONE feature across five tissues. This is its TRANSPOSE: up
# to MGENE_MAX features in ONE tissue, one result each, never pooled with one another.
#
# This function is the input boundary, and it exists because the failure here is SILENT.
# The list is typed or pasted by hand, so some of it will be typos, retired aliases, or
# symbols this tissue's matrices simply do not carry. A resolver that quietly returns the
# 7 it recognised out of 10 produces a table that looks complete: the 3 that vanished
# appear nowhere, and a reader comparing genes has nothing on screen to correct them. That
# is the quiet-default failure this project bans, and it is worse here than in most places
# because the reader chose those genes deliberately and will assume all of them were run.
#
# So nothing is dropped quietly. Every distinct token the caller passes comes back in
# exactly one of `ok` / `unknown` / `over_cap`, and the caller is expected to state the
# non-empty ones. Returning a structure rather than stop()-ing is deliberate and matches
# get_feature(): a mistyped symbol is a legitimate user-input outcome, not a caller bug.
# The hard fail stays at get_survival(), the analytical boundary.
#
# Pure, and takes the vocabulary as an ARGUMENT rather than reaching for FEATURES_BY_CT --
# the scan_lookup.R rule. app.R owns the vocabulary; the logic that can be tested against
# it lives here (tests/test_multigene_input.R).
#
#   raw     what the user typed: one string or several. Split on commas, semicolons and
#           any whitespace, so a pasted spreadsheet column, "A, B, C" and "A B C" all work.
#   vocab   the symbols this tissue actually carries for this score type.
#   max_n   the cap. Set to 10 because that is the size the tab was scoped to; a larger
#           list is a genome-wide scan and belongs in batch_regulators(), which is built
#           for it and reads each cohort's matrix once instead of once per gene.
#
# Returns list(ok, unknown, over_cap, n_input):
#   ok        canonical symbols to run -- input order, deduplicated, at most max_n
#   unknown   tokens matching nothing in vocab (see the ambiguity note below)
#   over_cap  recognised symbols BEYOND max_n, named rather than silently truncated
#   n_input   distinct symbols the caller asked for, so a caller can say "10 of 14"
#             without re-deriving it and getting a different number
MGENE_MAX <- 10L

mgene_resolve_symbols <- function(raw, vocab, max_n = MGENE_MAX) {
  if (is.null(raw)) raw <- character(0)
  # Fail loud on the shapes that are caller bugs rather than user typos. An unnamed factor
  # or a numeric column read from a file would otherwise tokenise into its integer codes
  # and come back as a list of plausible-looking "unknown symbols".
  if (!is.character(raw))
    stop("`raw` must be character (what the user typed); got ",
         paste(class(raw), collapse = "/"))
  if (!is.character(vocab))
    stop("`vocab` must be a character vector of the symbols this tissue carries; got ",
         paste(class(vocab), collapse = "/"))
  if (!is.numeric(max_n) || length(max_n) != 1 || !is.finite(max_n) || max_n < 1)
    stop("max_n must be a single positive number of symbols")

  tok <- unlist(strsplit(raw, "[,;[:space:]]+"))
  tok <- tok[!is.na(tok) & nzchar(tok)]
  if (!length(tok))
    return(list(ok = character(0), unknown = character(0),
                over_cap = character(0), n_input = 0L))

  # Exact match wins. Only if there is none do we fold case, and only when the fold is
  # UNAMBIGUOUS: if two vocabulary entries differ from each other by case alone, picking
  # either one is a guess, so the token goes to `unknown` and the caller says so. Guessing
  # here would silently analyse a different gene than the one asked for -- the one outcome
  # a survival table must never produce.
  .canon <- function(t) {
    if (t %in% vocab) return(t)
    h <- which(toupper(vocab) == toupper(t))
    if (length(h) == 1L) return(vocab[h])
    NA_character_
  }
  canon <- vapply(tok, .canon, character(1), USE.NAMES = FALSE)

  # Deduplicate on the CANONICAL name, not the typed token, so "ESR1" and "esr1" in one
  # list are one gene and not two identical rows.
  ok_all  <- unique(canon[!is.na(canon)])
  unknown <- unique(tok[is.na(canon)])
  list(ok       = utils::head(ok_all, max_n),
       unknown  = unknown,
       over_cap = if (length(ok_all) > max_n) ok_all[-seq_len(max_n)] else character(0),
       n_input  = length(ok_all) + length(unknown))
}

# ---- Multiple-gene query: disclosures and the result table --------------------------

# Rank-lookup statuses that are a property of the QUERY on this tab rather than of the
# gene. The tissue, endpoint, cohort set, horizon and score kind are all FIXED across a
# Multiple-gene run, so each of these verdicts is the identical sentence on all ten rows.
# mgene_header_notes() states them once instead.
#
# Only two statuses genuinely vary gene to gene and stay on the row: `ok` (the rank itself)
# and `not_found` (this gene was dropped from the scan, which is a fact about the gene).
#
# This is the step-31 reasoning transposed. That step suppressed `kind` across five panels
# because it described the query, not the panel; here the query fixes four more of the six
# statuses, because it is the tissue that is held constant instead of the feature.
MGENE_RANK_HEADER_STATUSES <- c("kind", "no_scan", "cohorts", "horizon")

# What ONE gene's result has to disclose, as severity-named sentences -- the per-row half
# of the split described above.
#
# `count_line` drops the leading pooled-count sentence. It is FALSE for the on-screen table
# and TRUE for the CSV, and that is not a style preference: on screen the `cohorts` column
# states k of n as a number and the header states the independence across genes, so the
# sentence is the same two facts a third time. In the CSV neither of those is present -- a
# row there is read with nothing else attached to it -- so it stays. Same reasoning as
# survtable_display() dropping the three constant columns survtable() keeps for the export.
mgene_row_notes <- function(e, count_line = TRUE) {
  if (!identical(e$status, "ok")) return(c(warn = e$msg))
  r <- e$res
  k <- if (!is.null(r$pooled) && !is.null(r$pooled$k)) r$pooled$k else .n_estimable(r)
  # Wording deliberately differs from the Multiple query tab's: "any other GENE" is the
  # misreading available here, and it is the one that matters, because ten genes in one
  # tissue look far more poolable than one gene in five tissues does.
  out <- if (!count_line) character(0) else c(info = sprintf(
    "Pooled over k=%d of %d %s cohorts. Nothing is pooled with any other gene.",
    k, e$n_sel, e$endpoint))
  if (k < 2) out <- c(out, warn = MSG_SINGLE_COHORT)
  out <- c(out, .missing_cohort_note(e))
  out <- c(out, .rank_note(e, suppress = MGENE_RANK_HEADER_STATUSES))
  out
}

# The three fields that describe the QUERY rather than a gene, with the one-run check that
# is what makes them meaningful.
#
# One run means one tissue, one endpoint and one score kind. Disagreement means the list was
# assembled from more than one run, and anything said about "the query" would then describe
# only part of it. Same check mq_export_dir() makes on the feature, for the same reason.
#
# Extracted from mgene_header_notes() (step 33 stage 4, 2026-08-09) rather than copied,
# because the second caller is the export FILENAME: a name is a claim about the whole file,
# so the two must raise on exactly the same input. A private copy of the loop would let the
# header refuse a two-run list while the download quietly labelled it with the first entry's
# tissue -- which is the failure a filename can least afford, since it outlives the session.
.mgene_query_fields <- function(entries) {
  if (!length(entries)) stop(".mgene_query_fields(): no entries")
  f <- c("ct", "endpoint", "kind")
  setNames(lapply(f, function(k) {
    v <- unique(vapply(entries, function(e) as.character(e[[k]])[1], character(1)))
    if (length(v) != 1)
      stop(sprintf("mgene entries disagree on %s (%s) -- assembled from more than one run",
                   k, paste(v, collapse = ", ")))
    v
  }), f)
}

# What the WHOLE query has to disclose, stated once above the table.
#
# `resolved` is an mgene_resolve_symbols() result, and passing it is what makes the
# disclosure complete: symbols that were unrecognised or over the cap HAVE no row, which is
# precisely the problem they represent, so the only place they can be reported is here.
mgene_header_notes <- function(entries, resolved = NULL) {
  if (!length(entries)) stop("mgene_header_notes(): no entries")
  .mgene_query_fields(entries)
  taus <- unique(vapply(entries, function(e)
    if (is.null(e$tau)) NA_real_ else as.numeric(e$tau), numeric(1)))
  if (length(taus) != 1)
    stop(sprintf("mgene entries disagree on the horizon (%s) -- more than one run",
                 paste(taus, collapse = ", ")))

  e1 <- entries[[1]]; ct <- e1$ct; ep <- e1$endpoint
  out <- c(info = sprintf(
    "%d gene%s, %s %s. Each gene is meta-analysed independently over this tissue's cohorts; nothing is pooled across genes.",
    length(entries), if (length(entries) == 1) "" else "s", cancer_label_of(ct), ep))

  # The horizon is read off a RESULT wherever there is one: what is stated is what was fit,
  # not what the selector held. e$tau is the fallback for a run where nothing estimated at
  # all, and it is the same number by construction.
  fit <- Filter(function(e) !is.null(e$res), entries)
  mfu <- if (length(fit)) fit[[1]]$res$max_followup else e1$tau
  out <- c(out, info = if (is.null(mfu))
    sprintf("Full follow-up: no horizon is declared for %s %s.", cancer_label_of(ct), ep)
  else
    sprintf("Horizon %g months (the declared default for %s %s). Patients still at risk are censored there, not excluded.",
            mfu, cancer_label_of(ct), ep))

  pep <- primary_endpoint_of(ct)
  if (!identical(pep, ep))
    out <- c(out, warn = sprintf("%s's primary endpoint is %s; %s is a secondary analysis here.",
                                 cancer_label_of(ct), pep, ep))

  # The query-level multiplicity verdict, once. Read off the first entry carrying one: by
  # MGENE_RANK_HEADER_STATUSES' definition it is the same verdict on every entry that has it.
  rk <- Filter(function(e) !is.null(e$rank) &&
                 e$rank$status %in% MGENE_RANK_HEADER_STATUSES, entries)
  if (length(rk)) {
    f <- format_scan_rank(rk[[1]]$rank)
    if (!is.null(f))
      out <- c(out, setNames(f$text, switch(f$tone, hit = "hit", miss = "warn", "info")))
  }

  out <- c(out, mgene_resolved_notes(resolved))
  out
}

# What the user asked for that produced NO ROW AT ALL, as severity-named sentences.
#
# Split out of mgene_header_notes() (step 33 stage 3, 2026-08-09) because it is the only thing that CAN
# be said when nothing resolved. That case is not exotic -- a typo'd list, or an mRNA symbol
# pasted against a VIPER vocabulary, resolves to zero genes -- and it is the case where
# silence is worst: an empty page after pressing Run reads either as a broken tool or, worse,
# as "none of your genes came out". There are no entries to describe, so the rest of the
# header cannot be built; these two sentences still can.
mgene_resolved_notes <- function(resolved) {
  if (is.null(resolved)) return(character(0))
  out <- character(0)
  if (length(resolved$unknown))
    out <- c(out, warn = sprintf(
      "Not recognised in this tissue's feature list, NOT analysed: %s.",
      paste(resolved$unknown, collapse = ", ")))
  # The cap is read back off the resolver's own output rather than from MGENE_MAX, so a
  # caller that passed a different max_n cannot be told the wrong limit.
  if (length(resolved$over_cap))
    out <- c(out, warn = sprintf(
      "Over the %d-gene limit, NOT analysed: %s.",
      length(resolved$ok), paste(resolved$over_cap, collapse = ", ")))
  out
}

# The whole Multiple-gene query as ONE long table -- including the genes that produced
# nothing, for the reason mq_survtable() gives: a gene that quietly vanished would read as
# "no effect" when it actually means "not measured here" or "the call raised".
#
# `query_note` repeats on every row and that is deliberate. The CSV lands in results/ with
# nothing attached to it, and the step-29 rule is that the file must not say less than the
# screen: the horizon, the independence-across-genes statement and the unrecognised symbols
# are all on screen above the table, so they have to be IN the file, and a constant column
# is the only place a long table can put them. It also means a reader who filters the file
# down to one gene still has the disclosures attached.
mgene_survtable <- function(entries, resolved = NULL) {
  stopifnot(length(entries) > 0)
  qn <- paste(mgene_header_notes(entries, resolved), collapse = " ")
  do.call(rbind, lapply(entries, function(e) {
    tab <- if (is.null(e$res)) .mq_null_row(e) else survtable(e$res)
    cbind(data.frame(cancer_type = e$ct, status = e$status,
                     primary_endpoint = primary_endpoint_of(e$ct),
                     query_note = qn,
                     gene_note = paste(mgene_row_notes(e), collapse = " "),
                     stringsAsFactors = FALSE),
          tab, stringsAsFactors = FALSE)
  }))
}

# How many characters of the filename the gene list may claim. Not a filesystem limit --
# NAME_MAX is 255 bytes -- but a budget that leaves the name readable in a Downloads folder
# and leaves headroom for the tissue/kind/endpoint prefix and the browser's " (1)" suffix.
MGENE_NAME_BUDGET <- 80L

# The download's FILENAME, read off the ENTRIES and never off the live selectors.
#
# Same rule mq_dl_all states: the controls can be changed without pressing Run, so a name
# built from input$mg_ct would label the file with a query it does not contain.
# .mgene_query_fields() also raises on a list assembled from two runs, which is the one case
# where a name taken from the first entry would be a confident lie -- and a filename outlives
# the session that produced it, so it is the worst place in the app to be confidently wrong.
#
# No .tau_tag() here, unlike every other export in the app, and no stratification flag. Both
# are deliberate. This tab has no horizon control -- tau is horizon_for(ct, ep), so it cannot
# vary once the tissue and endpoint are in the name -- and stratification, which CAN vary,
# is carried per row by survtable()'s own `strata_adj` column, as the horizon is by `max_fu`.
# The tag exists to stop two DIFFERENT estimands sharing one name; here the difference is
# inside the file on every row, where a reader will find it without having read the name.
#
# THE SYMBOLS ARE SANITISED, and that is not defensive boilerplate for a case that cannot
# happen. Measured 2026-08-09 against this repo's own matrices: ovarian mRNA offers 32756
# features of which 5828 CONTAIN "/" -- probe sets mapping to several genes, e.g.
# "ABCB6///ATG9A" -- and the longest, from GSE13876, is 1466 characters. Both are reachable
# here: mgene_resolve_symbols() splits on commas, semicolons and whitespace but NOT on "/",
# so such a name resolves as one symbol, gets analysed, and arrives intact at the export.
# Unsanitised, ten of them are a path separator inside a Content-Disposition filename and a
# name several times over NAME_MAX.
#
# Over budget the gene segment becomes "<n>genes" rather than a truncated list, for the same
# reason mgene_resolve_symbols() refuses an ambiguous case-fold: a clipped "ABCB6___ATG9A-TP"
# NAMES A GENE THAT WAS NOT RUN. "7genes" claims less and is true. Nothing is lost either
# way -- the real symbols are in the file's own `feature` column, one row group each.
#
# NOTE this was one of FIVE handlers building a filename by sprintf'ing a feature symbol into
# it, and the only one that sanitised. The other four were fixed in two steps, in the order
# step 35's measurement ranked them: mq_dl_all / mq_export_dir on 2026-08-10 (step 36, a silent
# no-op, the worst severity), then dl_forest / dl_km / dl_table on the Single query tab (step
# 37). Every filename in the app now goes through .safe_name(), here or via .feature_tag().
mgene_export_name <- function(entries) {
  q <- .mgene_query_fields(entries)
  ids <- .safe_name(vapply(entries, function(e) as.character(e$feature_id)[1], character(1)))
  genes <- paste(ids, collapse = "-")
  if (nchar(genes) > MGENE_NAME_BUDGET) genes <- sprintf("%dgenes", length(ids))
  sprintf("mgenes_%s_%s_%s_%s.csv", .safe_name(q$ct), .safe_name(q$kind),
          .safe_name(q$endpoint), genes)
}

# How many characters of the on-screen `gene` cell one symbol may claim.
#
# 40 is DERIVED, not chosen: the longest feature in the whole 58366-string vocabulary that
# does not contain a "/" is 40 characters (the same measurement FEATURE_NAME_BUDGET's note
# records). So this budget never touches a symbol that names a single gene -- it can only
# ever shorten a "///"-joined probe set, which is exactly the class that has no whitespace
# and therefore no line-break opportunity.
#
# Why the cell and not the CSS: DataTables measures its column widths from the rendered
# content and then writes them as inline styles, so no `overflow-wrap` rule can reach them
# (measured 2026-08-10 -- the mgene table came out 2405px wide inside a 1050px scroller, the
# gene column alone claiming 1855px, and applying `overflow-wrap: anywhere` to the cells,
# to the wrapper, or page-wide changed the layout by exactly zero pixels). Every OTHER
# column was then off-screen: you could not see the HR without scrolling 1.8k pixels
# sideways. The width has to stop being 1855px at the source, which is here.
#
# mgene_survtable() -- the CSV -- is untouched and carries every symbol whole. That is the
# same display/export split survtable_display() makes against survtable().
MGENE_CELL_BUDGET <- 40L

# One symbol as a table cell. Truncated names SAY SO, the same rule .feature_tag() follows
# for filenames and .title_lines() for figure titles: a bare clipped prefix NAMES A FEATURE
# THAT WAS NOT RUN, and on this tab that prefix would sit in a row of HRs belonging to
# something else. `[+n chars]` cannot be misread as a symbol.
#
# The marker is inside the budget, not added to it. Capping the symbol at 40 and then
# appending the marker gives a 52-character cell -- the exact mistake that shipped a clipped
# forest axis label on 2026-08-10 before it was measured (BUILD_LOG step 38). The loop is a
# fixed point because the marker's own width depends on how many characters it announces.
.gene_cell <- function(g, budget = MGENE_CELL_BUDGET) {
  g <- as.character(g)
  vapply(g, function(s) {
    if (is.na(s) || nchar(s) <= budget) return(s)
    keep <- budget
    repeat {
      mark <- sprintf("[+%d chars]", nchar(s) - keep)
      new_keep <- budget - nchar(mark)
      if (new_keep <= 0L) return(mark)
      if (new_keep == keep) break
      keep <- new_keep
    }
    paste0(substr(s, 1L, keep), sprintf("[+%d chars]", nchar(s) - keep))
  }, character(1), USE.NAMES = FALSE)
}

# Every SCALAR one gene's row can state, computed once and named once.
#
# Split out of mgene_display() on 2026-08-16 (step 45), when the table gained an expanding
# per-gene drawer. That is what makes the row/drawer boundary a DECISION rather than a
# layout accident: the values are built here, MGENE_ROW_FIELDS below names the ones the
# COLLAPSED row shows, and every field it does not name falls to the drawer -- by
# complement, not by a second list. So moving a field between the two is an edit to one
# vector, no value is recomputed, no format string is duplicated, and a field cannot be
# moved OUT of the row into nothing, because there is no third place for it to land.
# tests/test_mgene_drawer.R pins that partition.
#
# A LIST, not a character vector: `n` and `events` stay numeric so DataTables right-aligns
# them and so the display keeps saying they are counts. Coercing everything to string here
# would make the partition tidier and the table worse.
#
# `n` and `events` are SUMMED over the cohorts this gene actually estimated. Same quantity
# forest_plot() prints on the pooled diamond, derived the same way (sum over the unskipped
# per-cohort fits), and it is emphatically NOT constant down the column: FEATURES_BY_CT is
# the UNION over a tissue's cohorts, so a gene can be in the tissue's vocabulary and still be
# missing from one cohort's matrix -- that gene is then pooled over fewer patients than the
# row above it, and a reader comparing the two HRs has to be able to see that. Note this
# survives the k-count too: two genes can BOTH read "3 of 4" and be pooled over different
# patients, because a cohort can carry the gene for only part of its own series. That is why
# `n` stayed in the row when the drawer was added and `I2` did not.
#
# The cohort count comes from the per-cohort list rather than from pooled$k. Those are the
# same number for anything get_survival() produced (metafor is handed exactly the unskipped
# fits), but the per-cohort list is the primary record and `n`/`events` are summed over it,
# so taking all three from one place means the row cannot disagree with itself.
.mgene_fields <- function(e) {
  ok  <- identical(e$status, "ok") && !is.null(e$res)
  est <- if (ok) Filter(function(x) isFALSE(x$skipped) && !is.null(x$logHR),
                        e$res$per_cohort) else list()
  po  <- if (ok) e$res$pooled else NULL
  nv  <- sum(vapply(est, function(x) isTRUE(x$ph_violated), logical(1)))
  list(
    gene     = .gene_cell(e$feature_id),
    cohorts  = sprintf("%d of %d", length(est), e$n_sel),
    n        = if (length(est)) sum(vapply(est, `[[`, numeric(1), "n")) else NA_real_,
    events   = if (length(est)) sum(vapply(est, `[[`, numeric(1), "events")) else NA_real_,
    HR       = if (is.null(po)) "" else sprintf("%.2f", po$HR),
    `95% CI` = if (is.null(po)) "" else sprintf("%.2f - %.2f", po$ci_lb, po$ci_ub),
    p        = if (is.null(po)) "" else fmt_p(po$p),
    # Blank at k=1 rather than "0%": there is no between-study variance to report from one
    # cohort, and a printed 0 would read as "these cohorts agree perfectly".
    I2       = if (is.null(po) || is.na(po$I2)) "" else sprintf("%.0f%%", po$I2),
    PH       = if (nv == 0) "" else sprintf("%d of %d  !", nv, length(est)))
}

# WHICH of those the collapsed row shows, in the order it shows them. THE ONE LINE TO EDIT
# to move a field up or down; everything else about the split is derived from it.
#
# The default split follows one rule: the row keeps what is compared ACROSS genes -- which is
# the whole purpose of this tab -- and the drawer keeps what describes ONE gene's internals.
# `I2` is agreement among that gene's own cohorts, so it goes down. `n`/`events` differ gene
# to gene and are the weight behind the HR beside them, so they stay up. `PH` stays up for
# the reason survtable_display() carries the "!": a pooled estimate has NO
# proportional-hazards test -- it is a weighted average of coefficients, not a fit with
# residuals -- so on a table of pooled rows that verdict has nowhere else to go, and hiding a
# WARNING behind a click is not the same act as hiding a DETAIL behind one. What went down
# is the prose note, which was pinned at 40% of the table width and was the clutter the
# drawer was asked for.
#
# `gene` is required to be here and the test says so: a row whose identifying symbol is in
# the drawer cannot be read at all until it is opened.
MGENE_ROW_FIELDS <- c("gene", "cohorts", "n", "events", "HR", "95% CI", "p", "PH")

# The two structural columns, named rather than positional for the reason
# SURVTABLE_SORT_PREFIX is: both are addressed by INDEX in JavaScript, where being wrong is
# silent -- an off-by-one hands the drawer opener the `p` column and it renders a p-value as
# the panel body.
MGENE_CHEVRON_COL <- " "        # header stays blank; the cell carries the disclosure arrow
MGENE_DRAWER_COL  <- ".drawer"  # the pre-rendered panel, hidden, read by the click handler

# The ON-SCREEN table: one row per gene, in the order the genes were typed, with the drawer
# for each gene carried alongside it.
#
# The drawer is built HERE, into a hidden column, rather than fetched when the row is opened.
# Two reasons, and neither is performance. First, every value it shows is already in
# `entries` -- the fits ran on Run -- so a round trip would re-derive what is in hand. Second,
# a drawer fetched later reads the reactive LATER: the same class of defect as a download
# handler reading input$feature instead of the result it was built from, where the panel that
# opens describes a query the row above it is not from. Rendering both from one `entries`
# means the row and its drawer cannot disagree.
mgene_display <- function(entries) {
  stopifnot(length(entries) > 0)
  miss <- setdiff(MGENE_ROW_FIELDS, names(.mgene_fields(entries[[1]])))
  if (length(miss))
    stop("MGENE_ROW_FIELDS names fields .mgene_fields() does not build: ",
         paste(miss, collapse = ", "))
  rows <- lapply(entries, function(e) {
    f <- .mgene_fields(e)
    d <- data.frame(c(setNames(list(MGENE_CHEVRON_CELL), MGENE_CHEVRON_COL),
                      f[MGENE_ROW_FIELDS],
                      setNames(list(mgene_drawer_html(e)), MGENE_DRAWER_COL)),
                    check.names = FALSE, stringsAsFactors = FALSE)
    d
  })
  d <- do.call(rbind, rows)
  rownames(d) <- NULL
  d
}

# The disclosure arrow. A literal character rather than a CSS ::before, so that what the
# cell MEANS survives with styling off -- and so the click target is the cell (td.mg-control)
# rather than a pseudo-element no handler can bind to.
MGENE_CHEVRON_CELL <- '<span class="mg-chev">&#9656;</span>'

# HTML-escape, written here rather than taken from htmltools::htmlEscape.
#
# Not invented-here: R/plots.R is library code that ships inside pkg/OMICohort, whose
# Imports are the engine's (survival, metafor, rhdf5, RSQLite, DBI + base). htmltools is a
# web dependency, and adding one to the STATISTICS package so that one string can gain five
# entity substitutions is the wrong trade -- the package would then need a web stack to
# install for a reason that has nothing to do with what it computes.
#
# `&` FIRST and the order is load-bearing: escaping "<" to "&lt;" and then escaping "&" would
# double-encode it into "&amp;lt;" and print the entity instead of the character.
.html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;",  x, fixed = TRUE)
  x <- gsub("<", "&lt;",   x, fixed = TRUE)
  x <- gsub(">", "&gt;",   x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

# ---- hover explanations on a note (2026-08-18) ---------------------------------------
# A note on screen says WHAT was measured; some of them need a sentence saying what the
# number means before a reader can use it. That sentence goes behind a small marker and
# appears on hover, so the page stays as short as it is now.
#
# HOVER IS THE TRIGGER, AND IT IS A CHOICE WITH A COST. A hover panel does not exist on a
# touch device -- there is no fallback, it simply never opens -- and it cannot be selected
# or copied. That is acceptable here only because this tool is stated to be PC-only. If
# that ever changes, the swap is CSS: the `:hover` selector in PAGE_CSS becomes a class
# that a click toggles, and nothing in this file moves. Everything else -- the marker, the
# text, the markup, where it sits -- is already independent of how it is opened.
#
# THE LENGTH CAP IS ENFORCED, NOT REQUESTED. The brief was one or two sentences, and a
# tooltip that grows into a paragraph is both unreadable on hover (it vanishes when the
# pointer moves) and the beginning of documentation living in a popup. A convention cannot
# hold that -- the next edit that adds "just one more clause" looks exactly like the code
# around it -- so it is a stop().
INFO_MAX_SENTENCES <- 2L
INFO_MAX_CHARS     <- 220L

# Sentence count from terminal punctuation followed by a space or the end of the string.
# Deliberately crude, and it OVER-counts: "e.g." mid-sentence reads as a boundary. That is
# the acceptable direction for a guard -- it can only ask for shorter text, never let a
# longer one through -- and the fix when it fires on real prose is to rephrase, not to teach
# this function about abbreviations. tests/test_info_hover.R pins the over-count.
.info_sentences <- function(x)
  length(regmatches(x, gregexpr("[.!?](\\s|$)", x))[[1]])

info_note_html <- function(text, info) {
  if (length(text) != 1 || length(info) != 1)
    stop("info_note_html() is scalar; got ", length(text), " text and ", length(info), " info")
  if (!nzchar(trimws(info)))
    stop("info_note_html(): empty explanation -- drop the marker instead of showing a blank one")
  n <- nchar(info)
  if (n > INFO_MAX_CHARS)
    stop(sprintf(paste("info_note_html(): explanation is %d characters (limit %d). It is a",
                       "hover tooltip, not documentation -- shorten it or put it in README."),
                 n, INFO_MAX_CHARS))
  ns <- .info_sentences(info)
  if (ns > INFO_MAX_SENTENCES)
    stop(sprintf("info_note_html(): explanation is %d sentences (limit %d)", ns, INFO_MAX_SENTENCES))
  sprintf(paste0('%s<span class="infowrap" tabindex="0">',
                 '<span class="infomark">i</span>',
                 '<span class="infopop">%s</span></span>'),
          .html_escape(text), .html_escape(info))
}

# ---- the notes that carry a hover explanation (2026-08-18, step 48) -------------------
# Each explanation is declared BESIDE the text it explains, never in app.R: app.R has no
# test harness, and an explanation that drifts from the sentence it annotates is exactly the
# failure a reader cannot see. info_note_html() enforces the length of every one of them.

# The score-type radio. "TF activity" is the one label on the page that is actively
# misleading if the reader has not met VIPER -- it is read as "the TF's expression", and
# then the VIPER and mRNA tabs look like they disagree with each other.
KIND_INFO <- c(viper = paste(
  "VIPER scores a TF from the behaviour of its target genes, not from its own mRNA.",
  "So it can move opposite to the gene's own expression, and a TF with flat mRNA can still score high."))

# Radio labels as HTML, one per kind, with the marker attached to the kinds that have an
# explanation. Returns a character vector; app.R wraps each in HTML(). EVERY label is
# escaped, including the ones with no marker -- a mix of escaped and raw strings behind one
# HTML() call is how a page ends up with a stray &amp; in it.
kind_choice_html <- function(kinds, labels, infos = KIND_INFO) {
  miss <- setdiff(kinds, names(labels))
  if (length(miss))
    stop("kind_choice_html(): no label for kind(s): ", paste(miss, collapse = ", "))
  vapply(kinds, function(k)
    if (k %in% names(infos)) info_note_html(labels[[k]], infos[[k]]) else .html_escape(labels[[k]]),
    character(1), USE.NAMES = FALSE)
}
# ---- how to read the FIGURES themselves (2026-08-30, step 117) -----------------------
# The Guide tab explains all of this at length. This is the other half of that item: a
# reader with a forest in front of them should not have to leave the result to find out what
# the diamond is. The Guide carries the NARRATIVE; these three carry the DEFINITIONS, at the
# figure, where the question is actually asked.
#
# Three, because there are three things on the page a reader can misread independently, and
# each of the three has a specific misreading worth spending a tooltip on:
#
#   forest -- that a wide interval crossing 1 is a small effect rather than an unresolved one;
#   km     -- that the curves ARE the model. They are not: the split is at the median and the
#             HR printed beside them is the continuous Cox on every patient, so the log-rank
#             p and the HR are answering two different questions about the same cohort. This
#             is the one that has actually confused people, and it is why the KM note names
#             both numbers rather than explaining the curve;
#   hr     -- that HR is a probability, or a risk over the whole study. It is per SD of the
#             score, and it is an average over follow-up -- which is exactly the assumption
#             results/rmst_sensitivity.csv was written to test, and it is stated here rather
#             than only in a document a web reader never opens.
#
# Every one goes through info_note_html(), so the 2-sentence / 220-character cap is enforced
# on them the same as on every other explanation in this file. That cap is why these say one
# thing each instead of three: a tooltip is not the Guide.
FIGURE_HELP_LABEL <- "How to read this"

FIGURE_INFO <- c(
  forest = paste(
    "One row per cohort: the square is that cohort's hazard ratio and the line its 95% CI;",
    "the diamond is the pooled estimate. An interval crossing 1 is unresolved, not small."),
  km = paste(
    "The curves split each cohort at its median score, but the HR quoted with them is the",
    "continuous Cox over every patient. The log-rank p tests the split, not that HR."),
  hr = paste(
    "Per one SD of the score, so 1.20 is 20% more hazard per SD and above 1 is worse.",
    "It averages the whole follow-up, so where proportional hazards is rejected it is not",
    "one constant ratio.")
)

# The line that carries the marker. `which` is checked against the declared set rather than
# defaulted, so a typo is a stop() at render time and not a figure that silently loses its
# help -- the failure mode of the whole mechanism is a marker that is simply absent, which
# looks identical to a page that never had one.
figure_help_html <- function(which, infos = FIGURE_INFO, label = FIGURE_HELP_LABEL) {
  if (length(which) != 1L || !(which %in% names(infos)))
    stop("figure_help_html(): unknown figure '", paste(which, collapse = ", "),
         "'. Declared: ", paste(names(infos), collapse = ", "))
  sprintf('<p class="fig-help">%s</p>', info_note_html(label, infos[[which]]))
}

# ---- the IDH stratum: what strata(idh) does and does not fix (2026-08-28, step 98) ----
# Stratifying a glioma pool on IDH gives the wildtype patients their own baseline hazard.
# That is the right correction for the confound, and it is NOT a claim that the stratum is
# one disease. In TCGA-LGG the 94 IDH-wildtype patients split almost exactly along the
# published methylation classes: 68 Classic-/Mesenchymal-like, which run HR 15.3 against the
# IDH-mutant patients, and 26 PA-like, which run 1.86 and are not distinguishable from them.
# 65 of the 94 meet a WHO 2021 glioblastoma criterion on the two of three markers TCGA
# publishes (TERT promoter, +7/-10) -- a LOWER bound, EGFR amplification not being in that
# table. So what a reader needs before reading a glioma HR is not "some patients are
# wildtype", it is how much of the EVIDENCE they are: a quarter of the lgg/OS pool's
# patients and nearly half of its events.
#
# ON SCREEN goes the arithmetic, derived per pool. In the HOVER goes the WHO fact, which no
# survival table can supply. Full derivation: SUMMARY.md / README.md, both pinned by
# tests/test_lgg_idh_caveat.R.
IDH_STRATUM_INFO <- paste(
  "WHO 2021 calls an IDH-wildtype diffuse glioma with TERT mutation or +7/-10 a glioblastoma.",
  "strata(idh) gives that group its own baseline hazard; it does not make it one disease.")

# Pure: counts a frame, never reads one. `clin` is a get_clinical() result (its `time` and
# `event` already resolved to the endpoint); `tau` is the horizon, or NULL for full
# follow-up. The eligibility filter mirrors the engine's (survival_engine.R:183) EXCEPT for
# `is.finite(score)`, which is per-TF and unknowable before a query runs -- so this counts
# the endpoint's ELIGIBLE patients and a given TF's fit may drop a few more. Said out loud
# because step 96 quoted a cohort's count where the artifact carried the fitted one.
idh_pool_composition <- function(clin, tau = NULL) {
  miss <- setdiff(c("idh", "time", "event"), names(clin))
  if (length(miss))
    stop("idh_pool_composition(): clinical frame is missing ", paste(miss, collapse = ", "))
  d  <- clin[is.finite(clin$time) & !is.na(clin$event) & clin$time > 0, ]
  ev <- as.logical(d$event) & (if (is.null(tau)) TRUE else d$time <= tau)
  g  <- ifelse(is.na(d$idh), "uncalled", d$idh)
  bad <- setdiff(unique(g), c("mut", "wt", "uncalled"))
  if (length(bad))                       # the staged vocabulary is two tokens; anything
    stop("idh_pool_composition(): unexpected idh value(s): ",   # else is a staging bug
         paste(sort(bad), collapse = ", "))
  list(n = nrow(d), events = sum(ev),
       n_wt  = sum(g == "wt"),       ev_wt  = sum(ev[g == "wt"]),
       n_mut = sum(g == "mut"),      ev_mut = sum(ev[g == "mut"]),
       n_unc = sum(g == "uncalled"), ev_unc = sum(ev[g == "uncalled"]))
}

# The sentence, for one (tissue, endpoint). NULL unless the tissue's DECLARED stratifier is
# idh -- so it cannot appear on breast, and it appears on any tissue that later declares idh
# without an edit here. Resolves the endpoint the way strata_coverage_note() does, falling
# back to the whole tissue on the first render, when input$endpoint is still NULL.
#
# This is the first thing in the plotting layer to read clinical data. It has to: the note
# has to be true before Run, and the app otherwise learns patient counts only from a fitted
# survresult, which does not exist yet. The read is one clinical.db query (~10 ms).
idh_stratum_note <- function(ct, ep = NULL) {
  co <- if (!is.null(ep) && nzchar(ep))
          tryCatch(cohorts_for(ct, ep), error = function(e) cohorts_for(ct))
        else cohorts_for(ct)
  if (!length(co)) return(NULL)
  # The stratifier is read off the registry here rather than taken from app.R's
  # strata_vars_of(): app.R has no test harness, and this decision is the whole gate.
  epr <- if (!is.null(ep) && nzchar(ep)) ep else "OS"
  # Asked FOR THIS ENDPOINT (step 102). lgg is idh on OS/DSS and none on DFS, so a
  # tissue-level answer would put an IDH note under a DFS pool fitted without strata(idh).
  sv <- unique(Filter(Negate(is.na),
                      vapply(co, stratifier_for, character(1), endpoint = epr)))
  if (!identical(unname(sv), "idh")) return(NULL)
  k <- idh_pool_composition(get_clinical(cohorts = co, endpoint = epr),
                            horizon_for(ct, epr))
  if (!k$n_wt || !k$events) return(NULL)
  sprintf("IDH-wildtype is %s of the %s patients here (%.0f%%) but %s of the %s events (%.0f%%).",
          format(k$n_wt, big.mark = ","), format(k$n, big.mark = ","),
          100 * k$n_wt / k$n,
          format(k$ev_wt, big.mark = ","), format(k$events, big.mark = ","),
          100 * k$ev_wt / k$events)
}

# The breast-only Step F badge. What moves into the hover is the METHOD; what stays on
# screen is the VERDICT (IS / NOT) and both hazard ratios, because a verdict behind a hover
# is a verdict most readers never see.
CNA_MEDIATION_INFO <- paste(
  "Compares the hazard ratio before and after adjusting for the gene's own copy number,",
  "to test whether the signal is just amplification. Breast only (Step F).")

# The mediation RULE moved here from app.R unchanged, and it is the reason this is a
# function rather than a sprintf at the call site: it is a threshold rule producing a
# verdict, and it was previously untestable. Requiring a genuine unadjusted association
# (p_unadj < 0.05) before saying "explained away" is deliberate -- at genome-wide scale many
# near-null TFs show large attenuation_pct / p_adj swings that are noise on an already-null
# effect, not mediation (HANDOFF.md, Step F genome-wide note).
# THE RULE ITSELF, shared. It was written out TWICE -- once for the CNA badge and once,
# character for character, for the immune one (Step D) -- which is how two copies of a
# threshold drift apart without either looking wrong. Both badges call this now, so a change
# to what "explained away" means cannot apply to one mediator and not the other.
# NA in any input is "not established", never a silent TRUE.
MEDIATION_ATTENUATION_MIN <- 50
MEDIATION_ALPHA <- 0.05
mediation_is_explained <- function(r) {
  need <- c("attenuation_pct", "p_adj", "p_unadj")
  miss <- setdiff(need, names(r))
  if (length(miss))
    stop("mediation_is_explained(): row is missing ", paste(miss, collapse = ", "),
         " -- it did not come from a *_mediation_summary.csv")
  !is.na(r$attenuation_pct) && r$attenuation_pct >= MEDIATION_ATTENUATION_MIN &&
    !is.na(r$p_adj) && r$p_adj > MEDIATION_ALPHA &&
    !is.na(r$p_unadj) && r$p_unadj < MEDIATION_ALPHA
}

cna_mediation_note <- function(r) {
  need <- c("HR_unadj", "HR_adj", "p_cna_only", "HR_cna_only")
  miss <- setdiff(need, names(r))
  if (length(miss))
    stop("cna_mediation_note(): row is missing ", paste(miss, collapse = ", "),
         " -- it did not come from cna_mediation_summary.csv")
  mediated <- mediation_is_explained(r)
  txt <- sprintf("Copy number: signal %s explained by own-locus CNA (HR %.2f -> %.2f after adjusting).",
                 if (mediated) "IS" else "NOT", r$HR_unadj, r$HR_adj)
  if (!is.na(r$p_cna_only) && r$p_cna_only < 0.05)
    txt <- paste(txt, sprintf("CNA alone also predicts OS (HR %.2f, p %.2g).",
                              r$HR_cna_only, r$p_cna_only))
  list(text = txt, mediated = mediated, info = CNA_MEDIATION_INFO)
}

# THIS IS NOT OPTIONAL DEFENCE. The vocabulary these symbols come from is not a list of tidy
# gene names: ovarian mRNA offers 32756 features of which 5828 contain "/" and the longest is
# 1466 characters (the measurement mgene_export_name()'s note records). Feature strings reach
# this function from the h5 column names, they are interpolated into markup, and DataTables
# is told not to escape the column that carries it -- so anything angle-bracketed in a
# probe-set name would be parsed as markup rather than shown. Every interpolated value below
# goes through .html_escape() for that reason, including ones that "cannot" contain a
# bracket today.
#
# An EMPTY value is drawn as an em dash rather than left blank. Blank is what `I2` is at k=1
# and what every field of a gene that produced no fit is, and a label followed by nothing
# reads as a panel that failed to fill rather than as a quantity that does not exist. The
# field is never dropped for being empty: a missing row is the silent loss this whole file
# is arranged against.
.mgene_kv <- function(label, value) {
  v <- if (length(value) != 1 || is.na(value) || !nzchar(as.character(value)))
    "&mdash;" else .html_escape(value)
  sprintf('<span class="mg-kv"><b>%s</b> %s</span>', .html_escape(label), v)
}

# One gene's per-cohort table, as the drawer's body.
#
# survtable_display() rather than a private formatter, and that is the point of the drawer:
# what opens under a gene here is the SAME table the Single query tab shows for that gene,
# built by the same function with the same "!" convention and the same pooled row. A second
# formatter would drift from it silently -- two views of one fit, rounding or labelling
# differently, with nothing on screen to say which is authoritative.
#
# The hidden sort keys come off by PREFIX, the way survtable_dt_options() addresses them. A
# fifth keyed column added to survtable_display() would otherwise appear here as a bare
# integer column called `.sortkey5` in the middle of the statistics.
.mgene_cohort_table <- function(res) {
  d <- survtable_display(res)
  cap <- attr(d, "caption")
  d <- d[, !startsWith(names(d), SURVTABLE_SORT_PREFIX), drop = FALSE]
  head <- paste(sprintf("<th>%s</th>", .html_escape(names(d))), collapse = "")
  body <- paste(vapply(seq_len(nrow(d)), function(i) {
    cells <- paste(sprintf("<td>%s</td>",
                           .html_escape(vapply(d[i, ], function(v)
                             if (is.na(v)) "" else as.character(v), character(1)))),
                   collapse = "")
    # The pooled row is a SUMMARY of the rows above it, not another cohort. On the Single
    # query tab its position is pinned by orderFixed; this table has no ordering at all, so
    # the distinction has to be carried by the row's own class instead.
    sprintf('<tr class="%s">%s</tr>',
            if (grepl("^POOLED", d$cohort[i])) "mg-pooled" else "", cells)
  }, character(1)), collapse = "")
  sprintf('<table class="mg-cohorts"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>%s',
          head, body,
          if (is.null(cap)) "" else sprintf('<p class="mg-cap">%s</p>', .html_escape(cap)))
}

# ONE gene's drawer, as an HTML string.
#
# Returns a string rather than a Shiny tag list because it is carried inside a DataTables
# CELL: the click handler hands `row.data()[j]` to `row.child()`, which takes markup. That
# also keeps this function testable without Shiny loaded, the same reason
# survtable_dt_options() returns a plain list of indices.
#
# WHAT IS IN HERE IS EVERYTHING THE ROW DOES NOT SAY, and it is derived as the COMPLEMENT of
# MGENE_ROW_FIELDS rather than listed. A field moved out of the row therefore appears here
# without this function being touched, and -- the case that matters -- a field dropped from
# MGENE_ROW_FIELDS without being dropped from .mgene_fields() shows up in the drawer instead
# of vanishing from the app.
#
# The FULL symbol heads the panel, untruncated. .gene_cell() clips the row's copy to 40
# characters because DataTables measures its column widths from the rendered content (the
# 2405px table that measurement records), but the drawer is a block element with a
# line-break opportunity, so this is where a "///"-joined probe set can finally be read whole.
#
# THE COUNT LINE COMES BACK ON. mgene_row_notes(count_line = FALSE) is what the old note
# column used, because the `cohorts` column stated k of n two cells away and the header
# stated the independence across genes -- the sentence was the same two facts a third time.
# In an opened drawer it is not: the panel is read as a unit about one gene, the way a CSV
# row is, so it takes the same TRUE the export takes.
#
# A gene that produced no fit still gets a drawer, and it says so. The alternative -- an
# empty panel -- is the failure this project keeps naming: the arrow would open onto nothing
# and read as a broken tool rather than as "this gene was not measured here".
mgene_drawer_html <- function(e) {
  f <- .mgene_fields(e)
  extra <- setdiff(names(f), MGENE_ROW_FIELDS)
  kv <- if (length(extra))
    paste(vapply(extra, function(k) .mgene_kv(k, f[[k]]), character(1)), collapse = "")
  else ""
  # Severity names come from mgene_row_notes(); the classes are PAGE_CSS's, the same ones
  # output$mg_notes paints the header sentences with, so one severity looks the same
  # wherever it is said.
  nt <- mgene_row_notes(e, count_line = TRUE)
  notes <- paste(vapply(seq_along(nt), function(i) {
    sev <- names(nt)[i]
    cls <- switch(if (is.null(sev) || is.na(sev)) "" else sev,
                  warn = "note-warn", hit = "note-good", "note-body")
    sprintf('<p class="note %s">%s</p>', cls, .html_escape(nt[[i]]))
  }, character(1)), collapse = "")
  body <- if (is.null(e$res)) "" else .mgene_cohort_table(e$res)
  sprintf(paste0('<div class="mg-drawer"><div class="mg-drawer-head">%s</div>',
                 '<div class="mg-kvs">%s</div>%s%s%s</div>'),
          .html_escape(e$feature), kv, .mgene_depmap_html(e), notes, body)
}

# The DepMap common-essentiality line for one gene's drawer, or "" when there is none.
#
# It sits between the key/value block and the notes because of what it IS: gene context,
# true of the symbol whatever this query returned, where everything below it describes THIS
# run. That is the same order the Single query tab's badge stack uses.
#
# It is NOT routed through mgene_row_notes(), and that is deliberate rather than
# convenient. Those notes are also what the `gene_note` column of the CSV export carries
# and what output$mg_notes prints above the table, so adding a line there would change an
# export nobody asked to change, and would repeat the flag once per gene above a table
# that already shows it per gene. It also could not carry the hover: those strings are
# .html_escape()d by the loop above, which would print the marker's markup as text.
#
# e$depmap is attached by app.R (depmap_note()), exactly as e$rank is attached from
# scan_rank_lookup(). R/plots.R therefore keeps its pure-rendering contract and does not
# read data/raw or source R/depmap.R -- tests/test_mgene_drawer.R loads this file alone.
.mgene_depmap_html <- function(e) {
  d <- e$depmap
  if (is.null(d)) return("")
  need <- c("text", "tone", "info")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop(".mgene_depmap_html(): e$depmap is missing ", paste(miss, collapse = ", "),
         " -- it did not come from depmap_note()")
  sprintf('<p class="note note-%s">%s</p>', d$tone, info_note_html(d$text, d$info))
}

# Column position of one of this table's structural columns, 0-BASED for JavaScript.
#
# Fail-loud and shared by the two consumers below for the reason survtable_dt_options()'s
# header gives: the whole content of that function and this one is a set of COLUMN INDICES,
# and an index that is wrong is wrong SILENTLY -- the table still draws, the arrow still
# clicks, and what expands is whichever cell the off-by-one landed on.
.mgene_col_index <- function(d, nm) {
  j <- match(nm, names(d))
  if (is.na(j))
    stop("mgene table has no column '", nm, "' -- it did not come from mgene_display()")
  j - 1L
}

# The DataTables `options` fragment for the Multiple-genes table.
#
# `ordering = FALSE` is kept from step 33 and is NOT cosmetic: HR / p / I2 are formatted
# strings, DataTables would sort them lexically, and the row order here IS the comparison the
# user typed. The drawer does not change that argument -- it adds a second reason, since a
# child row is attached to its parent and re-sorting would carry open panels around the table.
#
# `dom = "t"` for the reason it always was: the page length is the cap, so every gene the user
# asked for is on screen at once rather than behind a pager.
mgene_dt_options <- function(d) {
  list(
    pageLength = MGENE_MAX, dom = "t", scrollX = TRUE, ordering = FALSE,
    columnDefs = list(
      # Hidden, NOT dropped: the click handler reads it out of the row's own data. It is also
      # searchable = FALSE so that a search box, if one is ever turned on, cannot match a
      # gene on markup the reader cannot see.
      list(targets = .mgene_col_index(d, MGENE_DRAWER_COL),
           visible = FALSE, searchable = FALSE),
      list(targets = .mgene_col_index(d, MGENE_CHEVRON_COL),
           className = "mg-control", orderable = FALSE, width = "28px")))
}

# The columns DataTables must NOT escape, as the negative index vector its `escape` argument
# takes, 1-BASED because that argument is R-side rather than JavaScript-side.
#
# STATED AS A WHITELIST OF TWO, never as escape = FALSE. The difference is the whole of this
# function: escape = FALSE would also stop escaping the `gene` column, whose contents are h5
# feature strings, and this tab's entire input is symbols the user pasted in. Two columns
# carry markup this file wrote; every other column is data and stays escaped.
mgene_dt_escape <- function(d)
  -c(.mgene_col_index(d, MGENE_CHEVRON_COL) + 1L,
     .mgene_col_index(d, MGENE_DRAWER_COL) + 1L)

# The click handler, as a STRING for app.R to wrap in DT::JS().
#
# Not wrapped here: JS() belongs to htmlwidgets, and R/plots.R ships inside pkg/OMICohort,
# whose Imports are the engine's. The same boundary .html_escape() is on the near side of.
#
# `table` is the DataTables API object DT hands a callback. The handler is delegated to the
# TABLE rather than bound to each cell, so it survives a redraw -- a per-cell binding is lost
# the moment DataTables re-renders the body, which is exactly what pressing Run again does.
mgene_dt_callback <- function(d) sprintf(
  "table.on('click', 'td.mg-control', function() {
     var tr = $(this).closest('tr');
     var row = table.row(tr);
     if (row.child.isShown()) { row.child.hide(); tr.removeClass('mg-open'); }
     else { row.child(row.data()[%d]).show(); tr.addClass('mg-open'); }
   });", .mgene_col_index(d, MGENE_DRAWER_COL))

# ---- tumor vs normal-adjacent panel ---------------------------------------------------

# Same shape as KM_SIZE_IN/forest_height_in: ONE sizing rule, used by both the screen and
# the PDF export. Two copies of a height formula drifting apart is a defect this project
# has already had (the forest's 380px literal), so the app asks this function rather than
# repeating the arithmetic. Width grows with the number of cohorts drawn; on Single query
# that is almost always one, so the default is deliberately narrow -- a lone pair of boxes
# stretched across 1400px reads as a bar chart of two bars.
# The first slot carries the y axis and the title, so it is wider than the ones after it;
# a flat per-slot width would leave the single-cohort case (the Single query tab's normal
# state) too narrow for its own title.
TN_SIZE_IN <- c(width = 4.4, extra = 2.6, height = 4.4)
tn_size_in <- function(k) c(
  width  = unname(TN_SIZE_IN[["width"]]) + unname(TN_SIZE_IN[["extra"]]) * (max(1L, k) - 1L),
  height = unname(TN_SIZE_IN[["height"]]))

TN_COL <- c(tumor = "#b2182b", normal = "#4393c3")

# Half the box width (boxwex 0.42), so the point cloud fills the box and does not spill
# past its edges at any device size.
TN_JITTER <- 0.21

# Axis unit. Every HR in this tool is per SD (get_survival standardizes within cohort), and
# tumor_normal() standardizes on the cohort's tumors for exactly that reason, so one unit
# here is one unit there. Saying so on the axis is the difference between a reader
# comparing this panel to the forest correctly and comparing it by eye to nothing.
TN_YLAB <- "SD (tumor reference)"

# Figure-title size, deliberately under forest_plot()'s 1.1. See the draw site for why.
TN_TITLE_CEX <- 0.9

# The fold change that sits UNDER the SD, for kind = "expr" only (see .tn_fc()).
#
# The axis stays per-SD; this is the companion number, because an SD alone is not
# interpretable -- it is divided by how variable that cohort's tumors are, so the same SD is
# a different amount of biology gene to gene. Breast makes the case on its own data:
# ESR1 +0.44 SD is 2.3x, ERBB2 +0.46 SD is 1.5x. Same SD, and nobody would call those the
# same result.
#
# Printed as log2FC AND the ratio, with the ratio taken as 2^log2FC rather than flipped for
# downs. "0.19x" and "5.3x" then carry their own direction, so the label needs no "higher"/
# "lower" word that could be read against the wrong box -- and no threshold has to be
# invented to decide when a change is small enough to describe as none.
.tn_fc_ratio <- function(log2fc) {
  r <- 2^log2fc
  if (r >= 1) sprintf("%.1fx", r) else sprintf("%.2fx", r)
}
.tn_fc_label <- function(x) {
  if (!identical(x$fc_status, "ok") || !length(x$log2fc) || is.na(x$log2fc)) return("")
  sprintf("   log2FC %+.2f (%s)", x$log2fc, .tn_fc_ratio(x$log2fc))
}

# Tumor vs matched normal-adjacent tissue, one slot per ESTIMABLE cohort.
#
# Cohorts that cannot be drawn are not drawn -- they are reported as prose by tn_notes(),
# beside the figure. An empty-but-present box would claim the comparison was made and
# found nothing, which is the opposite of "this cohort carries no label for it".
tn_plot <- function(tn, file = NULL) {
  stopifnot(inherits(tn, "tumor_normal"))
  est <- tn_estimable(tn)
  if (!length(est))
    stop("tn_plot(): no cohort in this result has a paired tumor/normal comparison. ",
         "Call tn_estimable() first -- the panel should not be drawn at all.")
  s <- tn_size_in(length(est))
  opened <- .open_dev(file, width = s[["width"]], height = s[["height"]])
  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    suppressWarnings(graphics::par(op))
    if (opened) grDevices::dev.off()
  })

  # One common y range across slots so two cohorts are read against each other, not each
  # against its own invisible scale.
  vals <- unlist(lapply(est, function(nm) c(tn$per_cohort[[nm]]$tumor, tn$per_cohort[[nm]]$normal)))
  ylim <- range(vals, na.rm = TRUE)
  # oma[3] reserves TWO title lines whether or not both are used. The count depends on the
  # device, and the app sizes its plotOutput from tn_size_in() before any device exists, so
  # a margin sized to the actual line count would be right in the export and wrong on
  # screen. Reserving the maximum costs ~0.15in and cannot clip.
  graphics::par(mfrow = c(1L, length(est)), mar = c(4.2, 0.6, 3.4, 0.6),
                oma = c(0, 3.6, 3.4, 0.4))

  for (i in seq_along(est)) {
    x <- tn$per_cohort[[est[i]]]
    graphics::boxplot(list(x$tumor, x$normal), ylim = ylim, xlim = c(0.4, 2.6),
                      outline = FALSE, col = unname(TN_COL[c("tumor", "normal")]),
                      border = "grey25", boxwex = 0.42, names = c("", ""), frame = FALSE,
                      las = 1L, yaxt = if (i == 1L) "s" else "n",
                      # xaxt="n": boxplot draws its own axis LINE under the boxes, which at
                      # this geometry reads as a bracket joining them. The tick-less axis()
                      # below is the one that carries the labels.
                      xaxt = "n")
    # Points over the boxes: at n=41 normals a violin is mostly smoothing artifact, and the
    # reader should be able to see how few points the narrow box is made of. Jitter is
    # seeded so redraws of the same query are identical -- an unseeded cloud that shifts
    # every reactive beat reads as the data changing.
    set.seed(1L)
    # amount=, not factor=: jitter()'s factor is scaled off the data's own range, which is
    # zero for a constant x, so the spread would come from an undocumented fallback rather
    # than from the box it has to sit inside. TN_JITTER is stated against boxwex above.
    graphics::points(jitter(rep(1, length(x$tumor)), amount = TN_JITTER), x$tumor,
                     pch = 16, cex = 0.2, col = "#00000018")
    graphics::points(jitter(rep(2, length(x$normal)), amount = TN_JITTER), x$normal,
                     pch = 16, cex = 0.35, col = "#00000055")
    graphics::axis(1, at = 1:2, labels = c("tumor", "normal"), tick = FALSE,
                   line = -0.6, cex.axis = 0.95)
    graphics::mtext(sprintf("n=%d", x$n_tumor), 1, line = 1.4, at = 1, cex = 0.62, col = "grey40")
    graphics::mtext(sprintf("n=%d", x$n_normal), 1, line = 1.4, at = 2, cex = 0.62, col = "grey40")
    # Cohort name and the two stat lines all live in the TOP MARGIN (positive mtext lines),
    # not inside the plot region. The single subtitle this replaced sat at line = -0.1, i.e.
    # just inside the region, which is invisible on a gene whose boxes sit low and a
    # collision on one whose boxes sit high: SFTPC in LUAD is 0.01x, so its normal box is
    # pinned to the top of the shared y range and the evidence line was drawn straight
    # through it. mar[3] = 3.4 lines holds all three with room to spare.
    graphics::mtext(x$cohort, 3, line = 1.9, cex = 1.05, col = "grey20")
    # The paired difference, not the gap between the two medians: see tumor_normal()'s
    # header for why only the paired quantity is invariant to cohort composition.
    #
    # Two lines, not one. The effect and the evidence for it are different questions, and
    # with the fold change appended a single line runs past the 4.4in slot and is clipped at
    # the p-value -- i.e. the number that says whether to believe it is the one that would
    # have been lost. The n here is the PAIR count, not the box counts underneath.
    graphics::mtext(sprintf("paired diff %+.2f SD%s", x$median_diff, .tn_fc_label(x)),
                    3, line = 0.95, cex = 0.62, col = "grey30")
    graphics::mtext(sprintf("%d pairs, p=%s", x$n_pair, fmt_p(x$p)),
                    3, line = 0.2, cex = 0.62, col = "grey45")
  }
  graphics::mtext(TN_YLAB, 2, outer = TRUE, line = 1.6, cex = 0.72, col = "grey30")
  # Measured against the OUTER width, because that is where the title is drawn. Using
  # .title_width() here would measure one mfrow SLOT while drawing across the whole device
  # -- and the first draft then rejoined the wrapped lines with a space, which put the
  # clipped title back exactly as .title_lines() exists to prevent. Each line is drawn.
  w_outer <- graphics::par("din")[1L] - sum(graphics::par("omi")[c(2L, 4L)])
  # TN_TITLE_CEX, not a literal, because the WRAP measurement and the DRAW must agree: sizing
  # the lines at one cex and drawing them at another silently reintroduces the clipped title
  # .title_lines() exists to prevent, and the two calls are 4 lines apart where a lone edit
  # looks complete. Below the forest's 1.1 on purpose (2026-08-13): the two figures now sit
  # side by side with their titles on the same line, and at equal weight a reader takes them
  # for two results of one query -- which is exactly the reading TN_CAPTION has to spend a
  # sentence undoing. Layout outranks caption, so the panel is made to look like the
  # annotation it is.
  tl <- .title_lines(tn$feature, "  |  tumor vs matched normal-adjacent",
                     w = w_outer, cex = TN_TITLE_CEX)
  for (i in seq_along(tl))
    graphics::mtext(tl[i], 3, outer = TRUE, line = 1.9 - 1.05 * (i - 1L),
                    cex = TN_TITLE_CEX, font = 2L, adj = 0)
  invisible(TRUE)
}

# The panel's export filename. A FUNCTION rather than a sprintf inside app.R's handler, for
# the reason mgene_export_name() records: app.R has no test harness, so a name built inline
# there can only be checked by grepping source, while this one can be called with a real
# tumor_normal object and its output read.
#
# Everything comes from the OBJECT. Nothing is passed in from the app, so there is no way for
# this name to describe a query the file does not contain -- the failure mode mq_res() carries
# feature_id to avoid, since a selector can be changed without pressing Run and a filename
# outlives the session that made it.
#
# WHAT IS DELIBERATELY ABSENT: the endpoint and the follow-up horizon, which every other
# export in this app carries (see .tau_tag). This figure has no follow-up in it at all -- it
# is expression in a tumor beside the same patient's own normal tissue -- so "_OS_tau60"
# would promise a file that changes with controls it does not depend on, which is the same
# lie .tau_tag exists to prevent, told in the other direction.
#
# WHAT IS PRESENT AND MIGHT LOOK REDUNDANT: `kind`. The same symbol drawn as expression and
# as VIPER activity are two genuinely different figures, and without it the second silently
# overwrites the first. And the COHORTS, which stand in for the cancer type: the panel is
# per-cohort, only the estimable ones are drawn, and the ids are self-describing enough
# (TCGA_BRCA) that a tissue segment would say it twice. Joined with "-" the way
# mgene_export_name() joins genes; k > 1 does not occur in the app today -- no tissue has two
# estimable tumor/normal cohorts -- but the panel is written for k slots and so is this.
#
# .safe_name() on the cohorts as well as .feature_tag() on the symbol. Cohort ids come from
# validate_cohorts() and cannot currently carry a separator, but "every filename in the app
# goes through .safe_name()" is a rule that is worth more than the one call it saves.
tn_export_name <- function(tn) {
  stopifnot(inherits(tn, "tumor_normal"))
  est <- tn_estimable(tn)
  if (!length(est))
    stop("tn_export_name(): no cohort in this result has a paired tumor/normal comparison, ",
         "so there is no figure to name. The download should not be offered at all.")
  sprintf("tumor_normal_%s_%s_%s.pdf", paste(.safe_name(est), collapse = "-"),
          .feature_tag(tn$feature), .safe_name(tn$kind))
}
# --- RPPA antibody panel ----------------------------------------------------
# NOT a forest, and deliberately not drawn like one -- but only where the difference SAYS
# something. Revised 2026-08-20 (step 51) after the first version diverged on cues that
# carried no meaning at all.
#
# In this tool a plotted row is ALWAYS a cohort -- separate patients, which is what makes
# combining them mean something. An RPPA panel's rows are antibodies on ONE patient set.
# Same grammar, opposite meaning: the single most available misreading in the whole app.
# These divergences are what hold that line, and each one is load-bearing:
#
#   no pooled diamond, ever                a diamond asserts the rows were combined
#   marker size FIXED across rows          a forest sizes by inverse-variance WEIGHT, i.e.
#                                          by contribution to a pooled estimate that does
#                                          not exist here. Equal markers cannot be read as
#                                          weights: there is no varying signal to misread.
#   n and events stated ONCE, in the       per-row n invites reading rows as independent
#     subtitle, never per row                samples
#   left header is "Antibody"              not "Cohort"
#   caption states the non-independence    the grammar cannot say it, so words must
#
# Everything else now MATCHES forest_plot(), because the first version's other differences
# were inconsistency wearing the costume of a signal. A reader who sees circles here and
# squares there, "1.28 (1.13-1.45)" here and "1.28 [1.13, 1.45]" there, a dashed grey
# reference line here and a dotted black one there, learns nothing about independence --
# they learn that the tool is untidy, and an untidy tool is one whose deliberate
# differences stop being read as deliberate. So the marker SHAPE, the CI punctuation, the
# reference line, the header rule, the stroke colours and the title placement are the
# forest's. The list above is what remains, and it is short enough to mean something.
#
# metafor::forest() is still not used here at all -- not with the diamond suppressed, not
# with arguments tuned. Passing these rows through the forest machinery is exactly the
# mistake the panel exists to prevent, and a suppressed diamond is one argument away from
# coming back. Matching its LOOK is a handful of graphics parameters; matching its CODE
# path would hand a future edit the diamond back for free.
# The caption is the only place the non-independence can actually be STATED, so it must
# never be the thing that runs off the page. Split on sentences first (deterministic, and
# where a reader would break it anyway), then wrap each sentence to the device -- the
# first draft drew it as one mtext() line and lost "...are never combined" off the right
# edge on a 7.5in render, which is the caption failing in exactly the case it exists for.
.rppa_caption_lines <- function(cex = 0.78) {
  sent <- trimws(unlist(strsplit(RPPA_PANEL_CAPTION, "(?<=\\.)\\s+", perl = TRUE)))
  sent <- sent[nzchar(sent)]
  unlist(lapply(sent, .chunk_to_width, w = .title_width(), cex = cex, font = 1L),
         use.names = FALSE)
}

RPPA_MARKER_CEX <- 1.15          # fixed for every row: see above
# The forest's square, not a shape of our own. pch 22 rather than metafor's solid 15
# because the fill is carrying the phospho/total distinction -- see .rppa_marker_bg().
RPPA_MARKER_PCH <- 22
# One constant for the right-hand header, because the margin is MEASURED from it and the
# figure is DRAWN from it. Two literals is how the panel came to reserve room for
# "HR (95% CI)" while printing "HR [95% CI]" -- a margin measured against a string the
# figure does not contain is the same defect tests/test_forest_label_fit.R exists for.
#
# It reads from the forest's constant rather than repeating its text (step 55). Step 51
# matched this header to the forest's DELIBERATELY, as part of adopting the forest's visual
# grammar; two literals that happen to be equal record that decision as a coincidence, and
# leave the panels free to drift apart the next time either is edited.
RPPA_HR_HEADER <- FOREST_HR_HEADER
RPPA_PANEL_WIDTH_IN <- 7.5

# Height must cover the caption, and the caption is a CONSTANT, so its worst-case line
# count can be bounded here even though the actual wrap is a property of the device (the
# title cannot do this -- a feature name is not a constant, which is why .title_lines()
# refuses to grow the device). Three lines is what RPPA_PANEL_CAPTION wraps to at the 7.5in
# export width; a wider screen device wraps to fewer and simply leaves a little more room.
# Was 2L until step 51 added the marker-fill key -- the constant and the caption are two
# halves of one fact, and tests/test_rppa_panel.R measures the wrap rather than trusting
# either of them.
RPPA_CAPTION_LINES_MAX <- 3L

rppa_panel_height_in <- function(k)
  2.55 + 0.42 * max(k, 1) + 0.26 * RPPA_CAPTION_LINES_MAX

rppa_panel_size_in <- function(k)
  c(width = RPPA_PANEL_WIDTH_IN, height = rppa_panel_height_in(k))

.rppa_hr_text <- function(r, show_n = FALSE) {
  if (isTRUE(r$skipped)) return("--")
  sprintf("%.2f [%.2f, %.2f]%s", r$HR, r$lo, r$hi,
          if (show_n) sprintf("  n=%d", r$n) else "")
}

# Row label. The phospho marker is a printed suffix, not a colour, so it survives
# greyscale printing and a colour-blind reader.
.rppa_marker_bg <- function(r) if (isTRUE(r$phospho)) "black" else "white"

.rppa_row_label <- function(r)
  paste0(r$antibody, if (isTRUE(r$phospho)) "  (P)" else "")

rppa_panel_plot <- function(panel, file = NULL) {
  if (!is.list(panel) || is.null(panel$rows))
    stop("rppa_panel_plot(): expected a list from rppa_panel(); got ",
         paste(class(panel), collapse = "/"))
  rows <- panel$rows
  if (!length(rows))
    stop("rppa_panel_plot(): this panel has no antibody rows to draw -- callers must ",
         "check status == 'absent' and show the message instead")
  est <- Filter(function(r) isFALSE(r$skipped), rows)
  if (!length(est))
    stop("rppa_panel_plot(): no estimable antibody in this panel -- callers must check ",
         "status == 'empty' and show the message instead")

  k <- length(rows)
  s <- rppa_panel_size_in(k)
  opened <- .open_dev(file, width = s[["width"]], height = s[["height"]])

  labs <- vapply(rows, .rppa_row_label, character(1))
  hrs  <- vapply(rows, .rppa_hr_text, character(1), show_n = isTRUE(panel$n_varies))
  # Margins measured from the strings that will actually be drawn, on the open device.
  # Guessing them in character units is what put the forest title off the page once
  # already; strwidth(units = "inches") needs font metrics, which exist now.
  lw <- max(graphics::strwidth(labs, units = "inches", cex = 0.85))
  rw <- max(graphics::strwidth(c(hrs, RPPA_HR_HEADER), units = "inches", cex = 0.85))
  # Bottom margin is set from the caption's ACTUAL line count on this device, not from a
  # constant: the same caption wraps to one line on the wide screen device and two on the
  # 7.5in export, and a fixed margin clips whichever it was not measured on.
  cap <- .rppa_caption_lines()
  # +0.45 / +0.40 rather than the label width alone: 0.28in of it is the gutter between the
  # text and the plot region (x_l / x_hr below sit at exactly that offset), and the
  # remainder is the page inset that keeps the header rule off the device edge. Guessing
  # these once already put a forest title off the page; they are measured, then spent.
  op <- graphics::par(mai = c(0.85 + 0.22 * length(cap), lw + 0.45, 0.75, rw + 0.40))
  on.exit({ graphics::par(op); if (opened) grDevices::dev.off() })

  lo <- vapply(est, `[[`, numeric(1), "lo")
  hi <- vapply(est, `[[`, numeric(1), "hi")
  xr <- range(c(lo, hi, 1))
  xr <- exp(log(xr) + c(-1, 1) * 0.08 * diff(log(xr)))
  y  <- rev(seq_len(k))

  # Headroom for the header row and the rule beneath it, the way metafor lays a forest out.
  # k + 0.85 rather than k + 0.6: the header used to sit exactly ON the top of ylim, which
  # left nowhere to put a rule between it and the first row.
  plot(NA, xlim = xr, ylim = c(0.4, k + 0.85), log = "x", axes = FALSE,
       xlab = "", ylab = "", xaxs = "i")
  graphics::abline(v = 1, lty = 3)                    # the forest's dotted refline

  # The label and CI columns live in the MARGINS, so the header rule that spans them has to
  # be drawn in user coordinates converted from inches -- the margin widths (lw, rw) are the
  # only thing that knows where those columns start and end. Derived from the same two
  # numbers the mai above was built from, so the rule cannot drift from the text it
  # underlines.
  xin  <- graphics::grconvertX(c(0, 1), "npc", "inches")
  x_l  <- graphics::grconvertX(xin[1] - lw - 0.28, "inches", "user")
  x_r  <- graphics::grconvertX(xin[2] + rw + 0.23, "inches", "user")
  x_hr <- graphics::grconvertX(xin[2] + 0.10, "inches", "user")
  old_xpd <- graphics::par(xpd = NA)

  graphics::text(x_l,  k + 0.75, "Antibody",    adj = 0, cex = 0.85, font = 2)
  graphics::text(x_r,  k + 0.75, RPPA_HR_HEADER, adj = 1, cex = 0.85, font = 2)
  graphics::segments(x_l, k + 0.45, x_r, k + 0.45)

  for (i in seq_len(k)) {
    r <- rows[[i]]
    dim_col <- if (isTRUE(r$skipped)) "grey45" else "black"
    graphics::text(x_l,  y[i], labs[i], adj = 0, cex = 0.85, col = dim_col)
    graphics::text(x_hr, y[i], hrs[i],  adj = 0, cex = 0.85, col = dim_col)
    if (isTRUE(r$skipped)) {
      graphics::text(sqrt(xr[1] * xr[2]), y[i], r$reason %||% "not estimable",
                     cex = 0.75, col = "grey45")
      next
    }
    graphics::segments(r$lo, y[i], r$hi, y[i], lwd = 1.6)
    graphics::segments(c(r$lo, r$hi), y[i] - 0.13, c(r$lo, r$hi), y[i] + 0.13, lwd = 1.6)
    # One size for every marker. Not a style choice -- see the header.
    graphics::points(r$HR, y[i], pch = RPPA_MARKER_PCH, cex = RPPA_MARKER_CEX, lwd = 1.4,
                     bg = .rppa_marker_bg(r))
  }
  graphics::par(old_xpd)

  at <- pretty(xr)
  at <- at[at > 0 & at >= xr[1] & at <= xr[2]]
  graphics::axis(1, at = at, labels = format(at, trim = TRUE), cex.axis = 0.85)
  graphics::mtext(sprintf("HR per SD of protein level (%s)", panel$endpoint),
                  side = 1, line = 2.4, cex = 0.9)

  graphics::title(main = sprintf("%s: RPPA protein in %s", panel$gene, panel$cohort),
                  line = 1.7, cex.main = 1.1)
  # ONE n for the whole panel where there IS one -- the usual case, and stating it once is
  # what stops the rows reading as independent samples.
  #
  # Where it varies (an antibody not run on part of the cohort) the RANGE is stated and
  # every row gains its own n on the right. Deliberately not a warning: this is ordinary
  # upstream missingness -- 34 of ovarian's 37 multi-antibody genes share n exactly, 3
  # differ by 10-20 patients -- and crying defect at that would teach a reader to skip the
  # line that also has to carry the real "these are the same people" claim.
  sub <- if (isTRUE(panel$n_varies))
    sprintf("n = %d-%d patients, %s %s events; rows differ where an antibody was not run",
            panel$n_range[1], panel$n_range[2],
            if (isTRUE(panel$events_varies))
              sprintf("%d-%d", panel$events_range[1], panel$events_range[2])
            else as.character(panel$events_range[1]), panel$endpoint)
  else sprintf("n = %d patients, %d %s events; the same patients in every row",
               panel$n, panel$events, panel$endpoint)
  graphics::mtext(sub, side = 3, line = 0.6, cex = 0.85, col = "grey30")
  for (i in seq_along(cap))
    graphics::mtext(cap[i], side = 1, line = 3.3 + 0.85 * i, adj = 0, cex = 0.78,
                    col = "grey35")
  invisible(panel)
}

# Filename for the panel export. Mirrors .safe_name()/.feature_tag() use elsewhere so a
# gene with a slash or a space cannot produce an unwritable path.
rppa_export_name <- function(panel)
  sprintf("rppa_panel_%s_%s_%s.pdf", .safe_name(panel$cohort),
          .safe_name(panel$gene), .safe_name(panel$endpoint))
