<!--
  This is the repository's landing page. GitHub resolves README in the order
  .github/README.md -> README.md -> docs/README.md, so this file wins.

  The README.md at the root is a different document: the internal development
  manual, ~1,100 lines, written for whoever maintains the engine. It stays where
  it is because BUILD_LOG.md refers to it by name 51 times and three test files
  assert numbers against it by path; renaming it to make room here would either
  falsify an append-only log or start a rename that touches ~80 prose sites for
  no gain. Two files, two audiences, one repository.
-->

# OMICohort

**Does this molecular score predict outcome?** — asked across **53 patient cohorts** in
**7 cancer types**, covering **10,430 patients** and **4,227 recorded events**.

Every cohort is fitted on its own and the fits are combined by meta-analysis. Cohorts are
never merged into a single pile, because a hazard ratio from pooled patients is a different
quantity from a pooled hazard ratio, and only the second one is comparable across studies.

Unlike the public Kaplan–Meier plotters, **the score you query does not have to be a gene's
expression.**

## Regulon activity, not expression

A transcription factor does its work as a protein, and how much of that protein is active is
set by things mRNA does not see — phosphorylation, cofactor availability, whether it is in the
nucleus at all. A TF can be flat at the mRNA level and still be driving a tumour.

So OMICohort asks a different question: *is this TF's regulatory programme switched on?* That
is read off the behaviour of the genes the TF controls, in two steps — **ARACNe** infers a
tissue-specific network from tumour expression, and **VIPER** scores each patient for how
consistently that TF's targets move together. Eight tissue-specific networks are in use,
covering **1,598 transcription factors**.

Expression, VIPER activity, multi-gene signatures and RPPA protein/phospho measurements all
go through the same engine, because every score is z-scored within its cohort before fitting.

## Run it

The image carries the code **and** the ~1.76 GB of data it serves from, so there is nothing to
download separately and nothing to configure.

```bash
docker run --rm -p 7654:7654 ghcr.io/korkmaz-lab/omicohort:0.1.0
```

Then open **http://localhost:7654**.

That tag is a multi-architecture manifest covering **`linux/amd64` and `linux/arm64`**, so the
same command pulls the right image on an Intel or AMD server and on an Apple Silicon Mac —
there is no `--platform` flag to get right. The two are not merely both built: the survival
fits, the pooled estimates and the rendered figures were compared across them and are
identical, byte for byte.

To serve on a different port, set `PORT` and map it to match:

```bash
docker run --rm -e PORT=8080 -p 8080:8080 ghcr.io/korkmaz-lab/omicohort:0.1.0
```

The container runs unprivileged, needs no network access once started, and writes nothing
outside its own temporary directory — so `--rm` loses nothing.

### Building the image yourself

This repository holds the code, not the data. A clone is therefore not enough to build on
its own: unpack the data bundle over it first, so that the following exist under the
repository root —

| from the bundle | what it is |
|---|---|
| `data/processed/` | the per-cohort expression and VIPER matrices, and `clinical.db` |
| `data/raw/depmap/`, `data/raw/cptac/` | the two annotation inputs read at runtime |
| `networks/*/regulon_*.rds`, `networks/*/tf_list_*.txt` | the 8 ARACNe networks |
| `results/*/*.csv` | the per-tissue scan tables the Multiple-query tab browses |

— and then:

```bash
docker build -t omicohort:0.1.0 .
```

The build runs `scripts/docker_selfcheck.R` inside the image and **fails by name** if anything
the app opens is missing, so an incomplete bundle stops the build rather than producing an
image that starts and then serves a tissue with a cohort quietly absent. You can run that same
check against a plain directory before building:

```bash
Rscript scripts/docker_selfcheck.R
```

## Running without Docker

R 4.5.2 with nine packages — `shiny`, `DT`, `shinythemes`, `RSQLite`, `DBI`, `data.table`,
`metafor`, `survival` (CRAN) and `rhdf5` (Bioconductor 3.22). R must be **started** in the
project directory, because `.Rprofile` is read at startup and turns off Shiny's automatic
sourcing of `R/`:

```bash
Rscript -e "shiny::runApp('.', port = 7654)"
```

## What is in this repository

Deliberately small: this is the tree the container image is built from, and nothing else.

| | |
|---|---|
| `app.R` | The Shiny application |
| `R/` | The analysis engine — the 12 files the app loads |
| `config/cohorts.tsv` | The cohort registry: the single source of truth for what exists |
| `www/` | Branding assets the page renders |
| `Dockerfile`, `.dockerignore` | How the image is built, and what goes into it |
| `scripts/docker_selfcheck.R` | Run inside the build; fails it by name if anything the app opens is absent |
| `CITATION.cff`, `LICENSE` | How to cite this, and the terms |

What is **not** here is as deliberate. The pipeline that staged the cohorts and ran the
genome-wide scans, the test suite, the packaged R library, and the development record are
kept but not published: they need the >10 GB raw tree to do anything, and several of them
hard-code the absolute paths of the machine they were written on. The scan tables and the
matrices travel with the data bundle, where output data belongs.

## Data

The tool ships **derived** matrices: expression z-scored within cohort, inferred VIPER
activity, and the clinical fields needed to fit survival. It redistributes no controlled-access
data. Each cohort's own terms are stated beside it on the **About** tab, and the two reference
inputs that make attribution a licence condition rather than a courtesy — **DepMap Public 26Q1**
(CC BY 4.0) and the **Broad GDAC Firehose** RPPA release (`doi:10.7908/C11G0KM9`) — are credited
there and in `CITATION.cff`.

## Citing

There is no paper yet. Until there is, cite the software — GitHub's *Cite this repository*
button reads `CITATION.cff`.

## Licence

Code is **GPL-3.0-or-later** (`LICENSE`); the prose and figures are **CC BY 4.0**. Neither
covers the patient data, and this project cannot grant rights it was not given.

## People

Built in the [Korkmaz Lab](https://research.ku.edu.tr/korkmazlab/functional-genomics-laboratory/)
at Koç University — Assist. Prof. Gözde Korkmaz (principal investigator) and
Arda Burak Karagöz (development and analysis).
