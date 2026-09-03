# R/about.R -- the About tab: where the data came from and on whose terms, how to cite the
# tool, what licence it carries, and who made it. Added 2026-08-30 (step 110).
#
# It is a SEPARATE tab from the Guide, deliberately. The Guide explains how to read a result
# and changes when the analysis changes; this page records provenance and people and changes
# when the project's circumstances do. Merging them would tie two different clocks together.
#
# That line is what step 118 moved content ALONG. "How it was built" -- ARACNe, VIPER, Cox,
# pooling -- was on this page and is now the Guide's closing section, because it answers "how
# are these numbers produced", which is the Guide's clock, not this one's. What stayed is
# everything that is a CONDITION or an IDENTITY: the citation, the people, the licence and the
# sources' terms. The test asserts the section is on exactly one of the two pages, since the
# failure of a move is not a missing section but a duplicated one that then diverges.
#
# The same rule as R/guide.R applies and for the same reason: the figures are DERIVED. The
# cohort lists under each data source are read out of config/cohorts.tsv, so a cohort added
# to the registry appears in its source's attribution automatically -- an attribution list
# that has to be updated by hand is one that will eventually credit the wrong set.
#
# WHAT IS DELIBERATELY NOT DERIVED, AND MUST NOT BE INVENTED. Accession numbers beyond those
# the registry already carries, formal citations, contact addresses and people's names exist
# nowhere in this repository. A page that made them up would be worse than one that says
# nothing: it would be wrong in the one section a reader is most likely to act on. So each
# is either derived, or declared once in ABOUT_PEOPLE / ABOUT_CITATION below and rendered
# only when non-empty. The GEO series are the happy case -- their accession IS their cohort
# id in the registry, so the whole list is derivable and exact.


# ---- 1. where each cohort came from --------------------------------------------------

# Source families, resolved from the cohort id and the tissue. Prefix rules rather than a
# hand-kept table, so a new cohort lands in a family by construction; the test asserts the
# families PARTITION the registry, which is the invariant that matters -- no cohort may ship
# unattributed, and none may be attributed twice.
#
# The ovarian rule is the one that needs stating: those GSE/EMTAB/PMID series reached this
# project through the curatedOvarianData Bioconductor package rather than one at a time from
# GEO, so they are attributed to it. Every other GSE was taken as an individual series.
.about_family <- function(cohort, cancer_type) {
  if (grepl("^TCGA_", cohort))            return("TCGA")
  if (grepl("^CGGA_", cohort))            return("CGGA")
  if (identical(cohort, "METABRIC"))      return("METABRIC")
  if (identical(cohort, "SCANB"))         return("SCAN-B")
  if (identical(cancer_type, "ovarian"))  return("curatedOvarianData")
  "GEO series"
}

# Order is presentation, not data: the large programmatic sources first, the individually
# collected series last, because that is the order a reader needs them in.
ABOUT_FAMILY_ORDER <- c("TCGA", "METABRIC", "SCAN-B", "CGGA", "curatedOvarianData", "GEO series")

# What each source requires of anyone using it. AUTHORED, and CONSTRAINED.
#
# THE AUDIENCE IS THE VISITOR, and that is the rule these were rewritten against on
# 2026-09-01 (step 121). Everything on this page is read by a researcher who has opened the
# site, not by whoever maintains it. They cannot see this repository: a sentence naming
# `scripts/download_cohorts.sh` or `config/cohorts.tsv` tells them nothing and quietly says
# the page was not written for them. Nor do they need the REASONING behind a decision -- why
# GPL-3 was chosen belongs in this file and in SUMMARY.md, not in front of someone who wants
# to know what they may do with a hazard ratio. Each entry below should answer one question
# the visitor actually has: where did this cohort come from, and what do I owe the people who
# produced it. Provenance that only the maintainer can verify goes in a comment like this one.
#
# THE CONSTRAINT, added 2026-08-31 (step 120), because this block was wrong for weeks in the
# section a reader is most likely to act on. METABRIC's entry asserted controlled access, a
# data access committee and an approved copy -- while scripts/download_cohorts.sh fetched it
# with a bare unauthenticated curl. A terms string is not free prose: the half of it that
# describes HOW THE DATA WAS OBTAINED is a claim about this repository, and this repository
# holds the evidence. tests/test_about_derived.R section 7 now refuses any claim that a
# reader cannot obtain a source that scripts/ downloads without credentials. If a source has
# a genuinely restricted tier, say which tier and say that it is unused -- do not open the
# entry with the restriction, because the entry is read as being about the data served here.
# NOTE: these strings are rendered through .g_esc(), which escapes `&`. An HTML entity written
# here would reach the page as the literal text "&mdash;", so anything typographic in this list
# has to be a real character. Since step 133 the question is settled a third way: em dashes, en
# dashes and double hyphens are not used in anything a reader sees, in either half of this file,
# so use a colon, a comma or a new sentence. tests/test_about_derived.R checks the RENDERED page
# rather than this source, which catches the entity form and the literal form at once.
ABOUT_TERMS <- list(
  "TCGA" = paste(
    "The Cancer Genome Atlas, distributed by the NCI Genomic Data Commons. The expression and",
    "clinical data shown here are open-access. If a TCGA cohort carries a result into work you",
    "publish, the programme asks that you acknowledge the TCGA Research Network."),
  "METABRIC" = paste(
    "Molecular Taxonomy of Breast Cancer International Consortium. The expression and clinical",
    "data shown here are the public release distributed by cBioPortal; METABRIC's raw",
    "sequencing sits separately at the EGA and is not used. cBioPortal releases this dataset",
    "under the Open Database License (ODbL), so you may reuse and redistribute it provided you",
    "credit the source and place any dataset you derive from it under the same licence. Cite",
    "Curtis et al., Nature 2012 and Pereira et al., Nature Communications 2016."),
  "SCAN-B" = paste(
    "Sweden Cancerome Analysis Network - Breast, a population-based cohort released through",
    "GEO. If a SCAN-B result matters to your work, cite the SCAN-B data release itself rather",
    "than this tool."),
  "CGGA" = paste(
    "Chinese Glioma Genome Atlas, released through the CGGA portal. Cite CGGA when a glioma",
    "result rests on it. The three CGGA batches appear as separate cohorts here, and are never",
    "merged, because they were run on different platforms."),
  "curatedOvarianData" = paste(
    "The Bioconductor curatedOvarianData package, which brings published ovarian series",
    "together under one clinical schema. Each series keeps its own terms and its own citation:",
    "cite the series a result rests on, not the package alone."),
  "GEO series" = paste(
    "Individual Gene Expression Omnibus series, listed by accession below. Each was assembled",
    "and deposited by its own submitters, and the accession is how they are credited: cite it",
    "when a result rests on that series."))


# ---- 2. the people, the citation, the licence ----------------------------------------

# DECLARED, never derived -- there is nothing in this repository to derive them from, and a
# fabricated name on a published tool is a fault of a different kind from a stale number.
# The lab name is the exception and is not invented: it ships in the brand image already
# (www/logo.png reads "OMICohort / KORKMAZ LAB"; see scripts/make_logo.py).
#
# TO ADD PEOPLE: append to `members`, one entry per person, each a list(name=, role=) with
# an optional `link=`. Anything left empty simply does not render -- the page never shows a
# blank slot or a placeholder, because a visible placeholder on a live site is worse than an
# absent section.
#
# SPELLING IS AS SUPPLIED for the PEOPLE, and that is a rule and not an oversight. These names
# were given by the user on 2026-08-30 in plain ASCII. Turkish orthography would write
# Gozde/Karagoz with diacritics, and this file does NOT apply them: guessing at the spelling of
# a person's own name is inventing, in the one section where inventing is worst. Changing them
# is a one-line edit here and the page follows.
#
# The INSTITUTION is no longer in that category, and stopped being on 2026-08-30 (step 118).
# Its own wordmark now ships in this repository -- www/logo-renkli.png reads "KOC UNIVERSITESI"
# with the diacritics -- so "Koc University" here would have put text on the page contradicting
# the image beside it about the name of a real institution. This is DERIVED from the supplied
# mark, not a guess applied by pattern to every Turkish name on the page; the two personal
# names above have no such second source and are untouched.
ABOUT_PEOPLE <- list(
  lab         = "Korkmaz Lab",
  # THE RESEARCH CENTRE, added 2026-09-03 at the user's request. Its own field rather than
  # being folded into `institution`, for the reason the mark below gives: each of the three
  # identities on this page has its own mark, and each is checked against its own mark. Fold
  # two names into one string and the check on the pair becomes a check on a sentence.
  #
  # NOT EXPANDED IN THE PROSE, and that is a decision. KUTTAM stands for the Koç University
  # Research Center for Translational Medicine, and www/kuttam-1.png prints exactly that under
  # the wordmark, on the same row of the same section. Writing it out again in the sentence
  # above would be the page saying the same thing twice within one screen. The expansion is in
  # the mark's `alt` instead, which is where a reader who cannot see the image will find it.
  center      = "KUTTAM",
  institution = "Koç University",
  # Supplied by the user on 2026-08-30. NAMED vector: the name is the visible link text and
  # the value is the href, and tests/test_about_derived.R asserts the page renders exactly
  # these and no other -- an href on this page is the one thing a reader will act on without
  # checking, so a URL nobody declared must not be able to appear.
  links       = c(
    "Lab website"  = "https://research.ku.edu.tr/korkmazlab/functional-genomics-laboratory/",
    "Lab on GitHub" = "https://github.com/Korkmaz-Lab"),
  contact     = "",                # e.g. "someone@university.edu"
  members     = list(
    list(name = "Assist. Prof. Gözde Korkmaz", role = "Principal Investigator"),
    list(name = "Arda Burak Karagöz",             role = "Development & Analysis")
  )
)

# The one place the three identities are joined into a single line. CITATION.cff quotes it as
# each author's affiliation, the Zenodo creator records take it as theirs, and the People
# paragraph renders it: three consumers, one derivation, so a fourth name can never be added to
# the page and forgotten in the deposit. tests/test_release_files.R asserts the CFF still
# carries what this returns.
about_affiliation <- function(p = ABOUT_PEOPLE) {
  parts <- c(p$lab, p$center, p$institution)
  paste(parts[!vapply(parts, is.null, logical(1)) & nzchar(unlist(parts))], collapse = ", ")
}

# The institutional marks, supplied by the user on 2026-08-30 and shown in the People section
# beside the lab and university they identify. UNLIKE the three brand images, these are not
# generated: scripts/make_logo.py builds www/logo.png, favicon.png and logo_full.png and the
# brand test can rebuild and re-check them, whereas these two arrived as files and there is
# nothing to regenerate them from. That is exactly why `file` is checked against the disk in
# tests/test_about_derived.R rather than trusted: a src pointing at a name that does not exist
# renders a broken-image icon and NOTHING ERRORS -- the same silence as the favicon that was
# served but never requested (steps 96-111). It is also why the filenames here were read off
# www/ instead of typed: the lab file is `korkmazlab.jpg`, and it was described in the message
# that supplied it as `korkmaz_lab.jpg`. Only one of those two spellings is on the disk.
#
# `height` is per-image on purpose. The three marks have very different aspect ratios: the lab
# mark is square, the university's is a wide lockup roughly 4.7:1, KUTTAM's is 1.75:1. A single
# shared height would render one of them at several times the visual weight of another.
#
# KUTTAM's 104 is the one value that was MEASURED rather than eyeballed, and it had to be. Its
# file is 350x200 but the artwork occupies only 55% of that height, the rest is white margin, so
# the obvious "about half the lab mark" reading of the number renders a mark half the size it
# looks like it should be: at 46 the artwork came out 25px tall against the university's 52 and
# the lab's 92, and the centre's name was unreadable. 104 puts its ink at 57px, between the
# other two, which is what the row was checked against on screen after the number changed.
#
# NOT SETTLED HERE: the Koc University mark is an institutional trademark, and universities
# normally publish rules for it (clear space, which colour variant, whether an external site
# may host it at all). Using it on a tool built in one of its own labs is the ordinary case,
# but that is a question for the institution and not one this file can answer.
#
# `role` IS NOT DECORATION. tests/test_about_derived.R checks that the institution's name in
# the prose and the institution's own mark agree on the spelling, and until step 134 it found
# that mark by looking for "University" in the alt text. Adding a third mark whose alt reads
# "KUTTAM, Koç University Research Center for Translational Medicine" would have matched that
# search twice, and the check is written to run only on a unique match: it would have gone
# SILENT, on the commit that added the logo, with nothing failing. So the marks now say what
# they are and the test selects on that, asserting exactly one of each declared role.
ABOUT_LOGOS <- list(
  list(file = "korkmazlab.jpg",  alt = "Korkmaz Lab", role = "lab",         height = 120L),
  list(file = "kuttam-1.png",
       alt = "KUTTAM, Koç University Research Center for Translational Medicine",
       role = "center", height = 104L),
  list(file = "logo-renkli.png", alt = "Koç Üniversitesi", role = "institution", height = 52L))

# The order the page's own <h3> sections appear in, DECLARED so that it can be asserted. The
# arrangement is a decision -- citation first, then the people, then the terms, with the two
# kinds of obligation under one heading (step 118, at the user's call) -- and a decision about
# ordering is invisible in a diff of the prose it orders. Without this, a later edit that moved
# People back under the licence would read as a formatting change.
ABOUT_SECTIONS <- c("How to cite", "People", "Licence and data terms")

# THE RELEASE VERSION, DECLARED ONCE (step 135). It was written out in five places before
# this: CITATION.cff, pkg/OMICohort/DESCRIPTION, the `docker run` line on the public README,
# and twice in the Dockerfile header. Nothing held them equal, so the first release after a
# text change would have shipped an image tagged one thing while the citation metadata said
# another, and the only way to find out is for somebody to pull the tag the README names and
# get a 404. Section 12 of tests/test_release_files.R now holds every one of them to this
# string.
#
# WHAT IT IS AND IS NOT. This is the version of the SOFTWARE and of the image built from it.
# It is deliberately NOT the version of the data deposit, which has its own line in
# scripts/make_zenodo_bundle.R and its own DOI: that is the whole point of depositing the two
# separately. A wording fix on the About page is a new image and is not a new dataset.
#
# BUMPING IT: patch for text, packaging and fixes that leave every number the app reports
# unchanged; minor when the app gains or loses a capability; and if a released estimate ever
# MOVES, that is not a version bump on its own, it is a BUILD_LOG entry saying which anchors
# moved and why, with the version following from that.
ABOUT_VERSION <- "0.1.1"

# No paper yet (confirmed 2026-08-30). While `doi` is empty the page renders the "not yet
# published" form; filling `doi` in switches it to a formal citation block and nothing else
# has to change.
ABOUT_CITATION <- list(doi = "", note = "")

# GPL-3, chosen 2026-08-30 for a specific reason rather than by habit: `metafor` is
# GPL (>= 2) and it is not an optional extra -- it IS the pooling engine, so every result
# this tool reports passes through it. Licensing the whole under GPL-3 removes the question
# of whether a GPL dependency that central makes the work derivative, at no cost in an
# academic setting where the tool is to be published and university-hosted. The prose and
# figures are CC BY 4.0, which is the norm for content and lets the guide be quoted.
#
# NOT a licence for the DATA, and the page says so: none of the sources above transferred
# any right to relicense their cohorts, and this tool cannot grant what it was not given.
ABOUT_LICENSE <- list(
  code = "GPL-3.0-or-later",
  docs = "CC BY 4.0")


# ---- 3. the page ----------------------------------------------------------------------
# Uses the Guide's classes throughout rather than declaring its own: two stylesheets for two
# pages of the same document is how they start to look like different sites.

about_facts <- function(registry = COHORTS) {
  fam <- vapply(seq_len(nrow(registry)), function(i)
    .about_family(registry$cohort[[i]], registry$cancer_type[[i]]), character(1))
  by <- split(registry$cohort, factor(fam, levels = ABOUT_FAMILY_ORDER))
  by <- by[vapply(by, length, integer(1)) > 0L]
  list(n_cohorts = nrow(registry), n_types = length(unique(registry$cancer_type)),
       families = lapply(by, sort))
}

.about_sources <- function(f) {
  blocks <- vapply(names(f$families), function(nm) {
    co <- f$families[[nm]]
    paste0('<h5 class="guide-h4">', .g_esc(nm), ' <span class="guide-fine">(',
           .g_int(length(co)), if (length(co) == 1L) ' cohort)' else ' cohorts)', '</span></h5>',
           '<p>', .g_esc(ABOUT_TERMS[[nm]]), '</p>',
           '<p class="guide-fine">', paste(paste0("<code>", .g_esc(co), "</code>"),
                                           collapse = " "), '</p>')
  }, character(1))
  paste(blocks, collapse = "")
}

.about_people <- function(p) {
  if (!length(p$members) && !nzchar(p$contact) && !length(p$links))
    return(paste0('<p>Built in the <b>', .g_esc(p$lab), '</b>. ',
                  'Contributors, affiliations and a contact address are not listed here yet: ',
                  'they are declared in <code>R/about.R</code> and this section fills ',
                  'itself in from that one place.</p>'))
  mem <- if (length(p$members))
    paste0('<ul class="guide-list">', paste(vapply(p$members, function(m)
      paste0('<li><b>', .g_esc(m$name), '</b>',
             if (!is.null(m$role) && nzchar(m$role)) paste0('<b>, </b>', .g_esc(m$role)) else '',
             '</li>'), character(1)), collapse = ""), '</ul>') else ""
  lnk <- if (length(p$links))
    paste0('<p>', paste(vapply(seq_along(p$links), function(i)
      paste0('<a href="', .g_esc(p$links[[i]]), '">', .g_esc(names(p$links)[[i]]), '</a>'),
      character(1)), collapse = " &middot; "), '</p>') else ""
  con <- if (nzchar(p$contact))
    paste0('<p>Contact: <code>', .g_esc(p$contact), '</code></p>') else ""
  # The institution is rendered only if declared -- same rule as every other field here. It is
  # the one identity with no second source in the repository to check it against: `lab` is
  # cross-checked against the wordmark in scripts/make_logo.py, this came from the user.
  # The centre and the institution each render only if declared, same rule as every other field
  # here, and in that order: lab inside centre inside university, which is also the order of the
  # marks below them.
  at <- paste0(vapply(c("center", "institution"), function(k)
    if (!is.null(p[[k]]) && nzchar(p[[k]])) paste0(', ', .g_esc(p[[k]])) else '',
    character(1)), collapse = "")
  paste0('<p>Built in the <b>', .g_esc(p$lab), '</b>', at, '.</p>', mem, lnk, con)
}

# The marks, rendered from ABOUT_LOGOS and from nothing else. Deliberately NOT wrapped in
# links: the href set on this page is asserted to be exactly ABOUT_PEOPLE$links, and quietly
# turning an image into a second route to a URL would make that assertion weaker while still
# passing it. A mark identifies; a link is a claim about where to go, and those are declared
# in one place.
.about_logos <- function(logos = ABOUT_LOGOS) {
  if (!length(logos)) return("")
  paste0('<div class="guide-marks">',
         paste(vapply(logos, function(g)
           sprintf('<img src="%s" class="guide-mark" alt="%s" style="height:%dpx">',
                   .g_esc(g$file), .g_esc(g$alt), as.integer(g$height)), character(1)),
           collapse = ""),
         '</div>')
}

about_html <- function(f, people = ABOUT_PEOPLE, citation = ABOUT_CITATION,
                       license = ABOUT_LICENSE, logos = ABOUT_LOGOS) {
  cite <- if (nzchar(citation$doi))
    paste0('<p>If OMICohort contributed to your work, please cite it: <code>',
           .g_esc(citation$doi), '</code>.</p>')
  else
    paste0('<p><b>There is no paper yet.</b> Until there is, cite OMICohort by name along ',
           'with the lab and the year you used it. Please cite the <b>data sources</b> below ',
           'separately: they carry their own requirements, and this tool does not stand in ',
           'for them. When a paper exists this section becomes a formal citation and nothing ',
           'else on the page changes.</p>')

  paste0(
'<div class="guide">',
'<h2 class="guide-h1">About OMICohort</h2>',
# The lead says what is NOT here as well as what is. Since step 118 the method left this page
# for the Guide, and a reader who arrives on About looking for it needs to be sent one click
# away rather than concluding the tool does not document it.
'<p class="guide-lead">A survival-analysis tool for molecular scores across <b>',
.g_int(f$n_cohorts), ' public cohorts</b> in <b>', .g_int(f$n_types),
' cancer types</b>. This page records who built it, how to cite it, and the terms the code ',
'and the data carry. <b>How the numbers are produced</b> is on the Guide tab.</p>',

'<h3 class="guide-h2">How to cite</h3>', cite,

'<h3 class="guide-h2">People</h3>',
.about_people(people),
.about_logos(logos),

# ONE section, TWO obligations, and the sub-headings are what keep them from reading as one
# (step 118, at the user's call). Merging them entirely was the alternative and was not taken:
# the licence points outward -- what this tool grants a reader -- and the source terms point
# inward -- what a reader owes the people whose patients these were. The third licence bullet
# is exactly the seam, and it only makes sense if the two halves are still distinguishable.
'<h3 class="guide-h2">Licence and data terms</h3>',
#'<p>Two different obligations, kept under one heading because a reader looking for either ',
#'should find both. The <b>licence</b> is what this tool grants you. The <b>source terms</b> ',
#'are what you owe the people whose data it was; nothing in the licence below alters them.</p>',

# The REASONING for GPL-3 (metafor is GPL >= 2 and is the pooling engine, so every reported
# result passes through it) deliberately does NOT appear on the page any more -- step 121. A
# visitor needs to know what the licence lets them do; why it was picked is a maintainer's
# question and lives in SUMMARY.md and in tests/test_about_derived.R section 4, which still
# holds metafor to being a hard Import.
'<h4 class="guide-h3">Licence</h4>',
'<ul class="guide-list">',
'<li><b>The code:</b> ', .g_esc(license$code), '. You may use, modify and redistribute it, ',
'including for commercial purposes, provided anything you distribute stays under the same ',
'licence.</li>',
'<li><b>This documentation, and the figures on these pages:</b> ', .g_esc(license$docs),
'. Reuse them freely with credit.</li>',
'<li><b>The data is covered by neither.</b> The cohorts below were produced by other groups, ',
'who set their own terms; nothing in the licences above changes what you owe them. Those ',
'terms are listed for each source.</li>',
'</ul>',

'<h4 class="guide-h3">Data sources and their terms</h4>',
# Stated POSITIVELY and only positively, which is a constraint and not a style note: this
# sentence is held by tests/test_about_derived.R section 7b, which reads the rendered page for
# claims that a source needed credentials. Phrase the correction as a denial -- "none of these
# required a data access committee" -- and it contains the very words the guard looks for, so
# the guard either fires on correct text or has to learn to parse negation. Say what is true
# instead. (The draft written in step 120 made exactly that mistake, and named `scripts/` at
# the reader as well, which is the step-121 fault: see the audience rule above ABOUT_TERMS.)
'<p>Every cohort below is <b>publicly available</b>, and each one was produced by a study ',
'that should be cited in its own right. OMICohort shows <b>summary statistics only</b> ',
'(hazard ratios, confidence intervals and survival curves) and never displays ',
'or distributes individual patient records. The list is exactly the set of cohorts this ',
'tool analyses.</p>',
.about_sources(f),

# ---- reference data --------------------------------------------------------------------
# NOT part of .about_sources(f), and that separation is load-bearing. The paragraph above
# ends "the list is exactly the set of cohorts this tool analyses", and ABOUT_TERMS is held
# to the registry in both directions by tests/test_about_derived.R -- a source described
# there that no cohort uses is an error. DepMap and the Firehose RPPA release are neither
# cohorts nor optional: DepMap ships under CC BY 4.0, which makes crediting it an OBLIGATION
# the page was silently failing to meet, and the Firehose release states its own citation.
# So they get their own block, after the cohorts, with the reason they are separate said out
# loud rather than left for the reader to infer from an odd-looking entry in a cohort list.
#
# The DepMap strings are the CONSTANTS from R/depmap.R, never retyped. The badge beside a
# gene and the credit on this page name one release or they are a bug; step 122 found that
# exact drift between the release constant and the staged files and closed it with a test.
# Section 9 there now holds this block to every staged input that declares an obligation.
'<h4 class="guide-h3">Reference data</h4>',
'<p>Two of the inputs behind these pages are <b>not patient cohorts</b>. They describe genes ',
'and proteins rather than people, they carry no survival information, and no hazard ratio is ',
'computed from them: they only annotate a gene you have already asked about. Each is ',
'listed here because each asks to be cited in its own right.</p>',
'<ul class="guide-list">',
'<li><b>', .g_esc(DEPMAP_RELEASE), '</b>: the source of the <i>common essential</i> ',
'note that can appear beside a gene, which means knocking that gene out stops growth in ',
'nearly every cancer cell line screened. It is a property of cell lines in a dish, not of ',
'the patients above. Released under ', .g_esc(DEPMAP_LICENCE), ', which asks that you credit ',
'it wherever the note carries into your own work. Cite: ', .g_esc(DEPMAP_CITATION), '</li>',
'<li><b>TCGA RPPA, Broad GDAC Firehose run stddata__2016_01_28</b>: the protein and ',
'phospho-protein measurements behind the Protein (RPPA) tab. The patients are the TCGA ',
'cohorts credited above; this is the analysis-ready release those measurements were taken ',
'from, and it asks to be cited alongside them. Cite: Broad Institute TCGA Genome Data ',
'Analysis Center (2016). Analysis-ready standardized TCGA data from Broad GDAC Firehose ',
'stddata__2016_01_28. Broad Institute of MIT and Harvard. doi:10.7908/C11G0KM9</li>',
'</ul>',

# "How it was built" stood here from step 110 until step 118, when it moved to the end of the
# Guide. It had ended by pointing AT the Guide for the same material, which is the tell that it
# was on the wrong page: a section whose closing move is to send the reader elsewhere for the
# rest of itself belongs where the rest is. What it must not do is come back and sit in both
# places -- tests/test_about_derived.R holds it to exactly one.

'</div>')
}
