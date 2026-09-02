# depmap.R -- is the queried gene one that every cancer cell line needs?
#
# A DepMap "common essential" gene is one whose CRISPR knockout stops growth in nearly
# every cell line screened. That is a caution flag on a survival hit, and a specific one:
# it does not say the association is wrong, it says the gene is core machinery, so a
# hazard ratio on it is unlikely to point at a selective vulnerability. Measured on the
# live scans, 60 of 1598 TFs carry it -- CTCF, YY1, MAX, NRF1, DNMT1, CENPA, TFDP1.
# (Step 49 listed MYBL2 here. It was never on the list, in 24Q4 or 26Q1 -- an exemplar
# nobody re-derived until step 122 did. Every count in this file is re-measured, not
# carried, for exactly that reason.)
#
# It is NOT redundant with the scan's own prolif_cor, which was the obvious worry. On the
# breast OS scan the median |prolif_cor| is 0.351 for the common-essential TFs and 0.360
# for the other 1452 -- indistinguishable, so the flag is not a proliferation proxy in
# disguise. Different axis: prolif_cor asks whether the TF's activity tracks a
# proliferation signature in tumors, this asks whether cells die without the gene in a
# dish. Under 24Q4 these were 0.313 and 0.360, which read as the essential TFs being the
# LESS proliferation-correlated half; 26Q1 closes that gap, and the weaker claim is the
# one the numbers now support.
#
# This is a cell-line property and not a patient result, in the same family as R/cptac.R
# and R/tumor_normal.R. No survival path reads it.
#
# ---- WHY TWO FILES ------------------------------------------------------------------
# The published essentials list alone cannot separate "screened and not essential" from
# "never screened", and the difference is real: DepMap 26Q1 covers 18531 genes -- 641 more
# than 24Q4 carried, and 26 fewer of the ones it did. Of the 1598 TFs the scans rank, 29
# are not in 26Q1 at all -- readthrough fusions and paralogue families (BORCS8-MEF2B,
# FOXD4L3/4/6, DUX4). One silence covering both would be the kind of ambiguity this
# project refuses elsewhere, so scripts/fetch_depmap.py stages the screened set as well
# (the column header of the 440 MB matrix; 264 KB on disk).
#
# The release is NAMED in every message for the same reason: "not screened" is a fact
# about 26Q1, not about the gene -- and the 24Q4 -> 26Q1 move is what naming it buys.
# 34 of the 63 TFs 24Q4 called unscreened are screened in 26Q1, and the TFs flagged
# common essential went 48 -> 60 (13 gained, ZNF143 lost). A badge that said only
# "DepMap" would have reported all of that as though nothing had changed.
#
# ---- WHY THE KIND GATE IS NOT OPTIONAL ----------------------------------------------
# Measured over every feature the three selectors offer:
#
#   kind    total  unresolvable  essential  screened-not-ess  not screened
#   viper    1598             0         60              1509            29
#   expr    58393          5861       1827             16681         34024
#   immune      77            36          0                 0            41
#
# The immune row is measured with the kind gate BYPASSED -- it is what the badge would say
# without the gate, which is the entire argument for having one. Asked honestly, all 77
# immune features return "kind" and carry no line at all.
#
# Immune features are cell-type scores (MCP_T cells, CIBERSORT_*), and 41 of the 77 have
# no whitespace, so a shape test alone would let them through and print "not screened by
# DepMap" against a T-cell score -- true, and misleading, because it implies DepMap might
# have screened one. The gate is on the KIND, and the shape test only runs inside it.
#
# Unresolvable means the feature does not name exactly one gene: 5828 of the expression
# features are '///'-joined probe sets (ABCB6///ATG9A) and 33 are free-text microarray
# annotation carried through from the platform ("aromatase cytochrome P-450 (P-450AROM)",
# "TCR BV7"). A per-gene property cannot be attributed to either, so they get no line at
# all -- the three states below are about a GENE, and those are not genes.

DEPMAP_RAW_DIR <- function() file.path(ROOT, "data", "raw", "depmap")

# The release these files came from. Staged by scripts/fetch_depmap.py, which pins the
# md5 of both published files -- see data/raw/depmap/SOURCE.txt for the source URL and
# the CC BY 4.0 attribution. 26Q1 has no figshare DOI; it is portal-only.
DEPMAP_RELEASE <- "DepMap 26Q1"

# What DepMap asks to be cited, and the licence that makes it an obligation rather than a
# courtesy. Both are rendered on the About tab (R/about.R, "Reference data"), which is the
# only place a visitor could ever see them.
#
# DERIVED from DEPMAP_RELEASE, deliberately. A release update that changed the badge but left
# a hand-typed citation behind would credit the wrong dataset -- the same class of drift that
# step 122 found between DEPMAP_RELEASE and the staged files, and closed with a test. The
# quarter is DepMap's own ordering ("DepMap Public 26Q1"), which is not the order its release
# archive uses ("DepMap 26Q1 Public"); the citation follows what the portal asks for.
DEPMAP_RELEASE_YEAR <- "2026"   # 26Q1 was released 2026-04-01
DEPMAP_LICENCE <- "CC BY 4.0"
DEPMAP_CITATION <- sprintf("DepMap, Broad (%s). DepMap Public %s. Dataset. depmap.org",
                           DEPMAP_RELEASE_YEAR, sub("^DepMap ", "", DEPMAP_RELEASE))

# Kinds whose features are gene symbols. See the table above for why immune is absent.
#
# "rppa" is here from step 51, and it is the one entry where the kind and the feature are
# NOT the same object. An RPPA row's feature is an ANTIBODY ("mTOR_pS2448"), which
# depmap_resolves() would reject and should reject. What the RPPA tab hands in is the
# PANEL'S GENE -- the thing the user queried, gene-first -- so the symbol contract this
# list asserts still holds. R/rppa.R owns that call (rppa_depmap_note()) and owns the
# caveat that comes with it: a knockout removes the protein, so it says nothing about
# whether the phosphorylated FRACTION matters.
DEPMAP_KINDS <- c("viper", "expr", "rppa")

# Both staged files are in the published `SYMBOL (entrez)` form, so one parser serves
# both and the Entrez id survives in the file as the disambiguator.
#
# The uniqueness assertion lives HERE, where the vector is built, and is a callable
# function rather than an inline stopifnot: step 46's CPTAC defect was a by-name lookup
# silently resolving to whichever duplicate came first, and a guard the correct path
# never trips is a guard no mutation test can kill. A test can hand this one a duplicate.
.depmap_symbols <- function(lines, what) {
  x <- trimws(lines)
  x <- x[nzchar(x)]
  bad <- !grepl("^\\S.*\\s\\(\\d+\\)$", x)
  if (any(bad))
    stop(sprintf(paste("depmap: %d of %d %s entries are not `SYMBOL (entrez)` --",
                       "first %s. A line that does not parse becomes a symbol nothing",
                       "can ever match, which reads on screen as 'not screened'."),
                 sum(bad), length(x), what, paste(head(x[bad], 3), collapse = ", ")))
  sym <- sub("\\s*\\(\\d+\\)$", "", x)
  if (anyDuplicated(sym))
    stop(sprintf(paste("depmap: %d of %d %s symbols repeat (first %s). Two Entrez ids",
                       "under one symbol makes membership ambiguous, and %%in%% would",
                       "answer for whichever the file happened to list."),
                 sum(duplicated(sym)), length(sym), what,
                 paste(head(sym[duplicated(sym)], 3), collapse = ", ")))
  sym
}

# list(essential, screened, release, n_essential, n_screened).
#
# An absent directory is a legitimate state -- data/raw/ is not required to run the app --
# and returns empty vectors so depmap_note() drops the badge. Every OTHER failure is loud,
# including a present-but-malformed file: a half-read essentials list would silently
# demote common-essential genes to "screened, not essential", which is a wrong answer
# wearing the muted styling of a routine one.
depmap_lists <- function(dir = DEPMAP_RAW_DIR()) {
  ess_f <- file.path(dir, "CRISPRInferredCommonEssentials.csv")
  scr_f <- file.path(dir, "depmap_screened_genes.txt")
  empty <- list(essential = character(0), screened = character(0),
                release = DEPMAP_RELEASE, n_essential = 0L, n_screened = 0L)
  if (!file.exists(ess_f) && !file.exists(scr_f)) return(empty)
  # One present and one missing is NOT the absent state. Continuing on the essentials
  # file alone would make every gene "not screened" except the 1827, which is the one
  # failure mode this pair of files exists to prevent.
  if (!file.exists(ess_f) || !file.exists(scr_f))
    stop(sprintf(paste("depmap_lists(): %s is present but %s is missing. The screened set",
                       "is what separates 'not common essential' from 'never screened' --",
                       "run scripts/fetch_depmap.py rather than staging one of the two."),
                 basename(if (file.exists(ess_f)) ess_f else scr_f),
                 basename(if (file.exists(ess_f)) scr_f else ess_f)))

  ess <- .depmap_symbols(readLines(ess_f, warn = FALSE)[-1], "essentials")  # drop header
  scr <- .depmap_symbols(readLines(scr_f, warn = FALSE), "screened")
  if (!length(ess) || !length(scr))
    stop(sprintf("depmap_lists(): staged files hold %d essential and %d screened genes",
                 length(ess), length(scr)))
  # Essential must be a SUBSET of screened, or the three states below overlap and a gene
  # could be both "common essential" and "never screened". This is also the check that
  # catches two files taken from different releases.
  miss <- setdiff(ess, scr)
  if (length(miss))
    stop(sprintf(paste("depmap_lists(): %d common-essential genes are absent from the",
                       "screened set (first %s) -- the two files are not from one release."),
                 length(miss), paste(head(miss, 3), collapse = ", ")))
  list(essential = ess, screened = scr, release = DEPMAP_RELEASE,
       n_essential = length(ess), n_screened = length(scr))
}

# A feature names exactly one gene. HGNC symbols never contain whitespace, and '///' is
# this project's multi-gene probe-set separator, so both are shape tests rather than a
# lookup against a vocabulary -- a feature that is a gene symbol DepMap has never heard
# of must still reach the "not screened" state, not be dismissed as unresolvable.
depmap_resolves <- function(feature)
  length(feature) == 1L && !is.na(feature) && nzchar(feature) &&
    !grepl("///", feature, fixed = TRUE) && !grepl("[[:space:]]", feature)

# Which of the five states this query is in, as its own function so the wording below and
# the tests can be checked against each other. `kind` and `unresolved` carry no line.
#
#   kind        the score type is not gene-symbol-valued (immune)
#   no_data     nothing staged in data/raw/depmap
#   unresolved  the feature does not name exactly one gene
#   essential   in the release's common-essential list
#   screened    covered by the release, not common essential
#   unscreened  not covered by the release at all
depmap_state <- function(feature, kind, lists) {
  if (!length(lists$screened)) return("no_data")
  if (!isTRUE(kind %in% DEPMAP_KINDS)) return("kind")
  if (!depmap_resolves(feature)) return("unresolved")
  if (feature %in% lists$essential) return("essential")
  if (feature %in% lists$screened) return("screened")
  "unscreened"
}

# What a reader has to know before using the flag: it is about cell lines in a dish, and
# the caution it raises is about selectivity, not about whether the hazard ratio is real.
# The length cap (2 sentences / 220 chars) is enforced by info_note_html().
DEPMAP_BADGE_INFO <- paste(
  "Common essential means knocking the gene out stops growth in nearly every cancer cell",
  "line, so it is a poor selective target. It is a cell-line property, not a patient result.")

# list(text, tone, info) for the badge, or NULL when there is nothing worth saying.
# tone is one of the note-* families app.R already styles: "warn" for the flag itself,
# "muted" for the two states that merely report coverage.
depmap_note <- function(feature, kind, lists) {
  st <- depmap_state(feature, kind, lists)
  if (st %in% c("no_data", "kind", "unresolved")) return(NULL)
  rel <- lists$release
  txt <- switch(st,
    # Deliberately no cell-line count. It would be the ROW count of CRISPRGeneEffect.csv,
    # and staging reads only that file's header (see scripts/fetch_depmap.py) -- quoting
    # a number this tree cannot verify is worse than not quoting one.
    essential  = sprintf("%s: COMMON ESSENTIAL -- knockout stops growth in nearly all screened cell lines.",
                         rel),
    screened   = sprintf("%s: screened, not a common-essential gene.", rel),
    unscreened = sprintf("%s: this gene was not screened in this release.", rel))
  list(text = txt, tone = if (st == "essential") "warn" else "muted",
       info = DEPMAP_BADGE_INFO, state = st)
}
