<!--
  This is the repository's landing page. GitHub resolves README in the order
  .github/README.md -> README.md -> docs/README.md, so this file wins.

  The README.md at the root is a different document: the internal development
  manual, ~1,100 lines, written for whoever maintains the engine. It stays where
  it is because BUILD_LOG.md refers to it by name 51 times and three test files
  assert numbers against it by path; renaming it to make room here would either
  falsify an append-only log or start a rename that touches ~80 prose sites for
  no gain. Two files, two audiences, one repository.

  HOUSE RULE FOR THIS PAGE (2026-09-03): no em dashes, no en dashes, no double
  hyphens. Section 11 of tests/test_release_files.R enforces it on prose lines.
  Use commas, colons, semicolons or a new sentence instead.
-->

# OMICohort

**Does this molecular score predict patient outcome?** OMICohort asks that question across
**53 patient cohorts** in **7 cancer types**, covering **10,430 patients**
and **4,227 recorded events**.

It is a survival analysis tool for cancer genomics, built in the Korkmaz Lab at Koç
University.

## The problem it is built around

A molecular marker that separates survival curves in one published cohort very often fails to
do so in the next. Part of that is biology: cohorts differ in stage, treatment era, platform
and subtype composition. Part of it is method, and two habits in particular make a result look
stronger than it is.

**Splitting patients at the median.** Cutting a continuous score into "high" and "low" throws
away most of the information in it and makes the answer depend on where the cut happened to
fall. A marker can look decisive in one cohort and disappear in another simply because the
median moved.

**Merging cohorts into one pile.** Pooling patients from several studies and fitting a single
model gives a number that is not comparable to the per-study numbers. A hazard ratio estimated
from pooled patients and a pooled hazard ratio are different quantities, and only the second
one carries the between-study variation that decides whether a finding travels.

OMICohort avoids both. Every cohort is fitted on its own, on the continuous score, and the
separate fits are then combined.

## How a cohort is fitted

Within each cohort the score is standardised, so the effect reported is the change in hazard
per one standard deviation of that score in that cohort. The model is Cox proportional
hazards, fitted on the cohort alone. Where a tissue has a subtype that shifts the baseline
hazard, such as IDH status in glioma or PAM50 class in breast, the fit is stratified on it
rather than left to average over it.

Overall survival and disease-free survival are treated as separate questions and are never
mixed together.

## How cohorts are combined

The per-cohort log hazard ratios are combined by random-effects meta-analysis (REML, using
`metafor`), which estimates how much the true effect varies between cohorts instead of
assuming that it does not vary at all. Alongside the pooled estimate the app reports I², the
share of the observed variation that is between-cohort rather than sampling noise. A pooled
estimate with high heterogeneity is not the same claim as a tight one, so the interface shows
both rather than only the headline number.

The forest plot is the primary output. It shows every cohort, its own confidence interval and
its weight, so a pooled result that rests almost entirely on one large study is visible as
exactly that.

## Kaplan-Meier curves are for looking, not for deciding

The app does draw a median-split Kaplan-Meier curve with a log-rank test, because a curve is
how most readers see a survival difference. That split is **for visualisation only**. The
number OMICohort reports, and the number the meta-analysis is built from, is always the
continuous per-standard-deviation estimate.

## Regulon activity, not expression

Unlike the public Kaplan-Meier plotters, the score being queried does not have to be a gene's
expression.

A transcription factor does its work as a protein, and how much of that protein is active is
set by things mRNA cannot see: phosphorylation, cofactor availability, whether it is in the
nucleus at all. A TF can be flat at the mRNA level and still be driving a tumour.

So OMICohort can ask a different question: is this TF's regulatory programme switched on? That
is read off the behaviour of the genes the TF controls, in two steps. **ARACNe** infers a
tissue-specific regulatory network from tumour expression, and **VIPER** scores each patient
for how consistently that TF's targets move together. Eight tissue-specific networks are in
use, covering **1,598 transcription factors**.

Expression, VIPER activity, multi-gene signatures and RPPA protein and phospho measurements
all go through the same engine, because every score is standardised within its cohort before
fitting.

## Beyond one query at a time

Every transcription factor in a tissue's network has already been scanned against survival,
for each endpoint, so a candidate can be placed among all of them rather than judged only on
whether it is significant on its own. Rankings are corrected for multiple testing with
Benjamini-Hochberg, and hits are counted at q < 0.05.

Where a tissue has enough cohorts, the scan is also run as a discovery and validation split:
one set of cohorts proposes and a held-out set tests. A regulator that clears FDR in discovery
and then holds up in validation is a stronger claim than one that clears FDR once.

## What a result does not mean

An association with survival is not a mechanism, a driver, or a target. These cohorts are
observational, treatment is uncontrolled, and a marker can track outcome because it tracks
stage or subtype rather than anything causal. The tool is built to make an association honest
and reproducible, not to tell you what it means.

## Run it

The image carries the code **and** the roughly 1.76 GB of data it serves from, so there is
nothing to download separately and nothing to configure.

```bash
docker run --rm -p 7654:7654 ghcr.io/korkmaz-lab/omicohort:0.1.0
```

Then open **http://localhost:7654**. The same command works on Intel, AMD and Apple Silicon
machines: the tag covers both architectures, and the two builds were checked against each
other down to the rendered figures.

The container runs unprivileged, needs no network access once it has started, and writes
nothing outside its own temporary directory.

## The data

OMICohort serves **derived** matrices: expression standardised within cohort, inferred VIPER
activity, and the clinical fields needed to fit survival. It redistributes no
controlled-access data.

The cohorts come from TCGA, METABRIC, SCAN-B, CGGA, curatedOvarianData and individual GEO
series. Each one's own terms are stated beside it on the **About** tab. Two reference inputs
make attribution a licence condition rather than a courtesy: **DepMap Public 26Q1** (CC BY
4.0) and the **Broad GDAC Firehose** RPPA release (`doi:10.7908/C11G0KM9`). Both are credited
there and in `CITATION.cff`.

This repository holds the code. The data is deposited separately so that the two can be cited
and versioned on their own terms.

## Citing

There is no paper yet. Until there is, cite the software: GitHub's *Cite this repository*
button reads `CITATION.cff`.

## Licence

Code is **GPL-3.0-or-later** (`LICENSE`); the prose and figures are **CC BY 4.0**. Neither
covers the patient data, and this project cannot grant rights it was not given.

## People

Built in the [Korkmaz Lab](https://research.ku.edu.tr/korkmazlab/functional-genomics-laboratory/)
at Koç University: Assist. Prof. Gözde Korkmaz (principal investigator) and Arda Burak Karagöz
(development and analysis).
