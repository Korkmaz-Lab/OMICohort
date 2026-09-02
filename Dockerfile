# OMICohort -- multi-cohort survival analysis, as a self-contained image.
# (2026-09-01, step 124.)
#
# WHAT IS IN HERE. Code plus the ~1.76 GB of data the running app opens, and nothing else:
# the working tree is ~19 GB, of which 14 GB is pipeline input (raw TCGA downloads), 1.7 GB
# is a GEO download cache, and the rest is prior-state archives. .dockerignore is a
# WHITELIST for that reason -- see the header there.
#
# WHY THE DATA IS BAKED IN RATHER THAN MOUNTED. Mounting means the image and the data can
# disagree about which release they are, and the app has no way to notice: it would fit
# whatever matrices it finds. Baking makes `docker run` reproduce one published state.
#
# BUILD (from the project root, with the data tree present):
#   docker build -t omicohort:0.1.0 .
# RUN:
#   docker run --rm -p 7654:7654 omicohort:0.1.0     ->  http://localhost:7654

# R is pinned to the version this project is developed and tested on. rocker/r-ver
# additionally pins a CRAN snapshot, which is what makes a rebuild in six months install
# the same package versions rather than whatever is current then.
FROM rocker/r-ver:4.5.2

LABEL org.opencontainers.image.title="OMICohort" \
      org.opencontainers.image.description="Multi-cohort Cox meta-analysis of molecular scores across 53 cancer cohorts" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.source="https://github.com/Korkmaz-Lab/OMICohort"

# System libraries. cairo/freetype/fontconfig are NOT optional decoration: without them
# png() has no device in a headless container and every Kaplan-Meier curve renders as an
# error box, while the app itself starts and looks fine. fonts-dejavu gives the devices a
# real font to draw with -- without one, axis labels come out as boxes.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcairo2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff5-dev \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        zlib1g-dev \
        fonts-dejavu-core \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Pin CRAN to a dated snapshot, deriving the Ubuntu codename from the base image rather
# than writing it down -- a hardcoded "noble" silently falls back to source builds the day
# rocker moves to the next LTS, which turns a 3-minute install into a 40-minute one.
# DO NOT ESCAPE THE DOLLARS HERE. Measured on this machine (docker 29.7.2) rather than
# recalled: in a RUN, `$FOO` set earlier in the same command expands normally, and `\$FOO`
# is what produces a LITERAL "$FOO". A first version of this file escaped them on the
# belief that Docker blanks undeclared names -- the shell then wrote a repo URL containing
# the characters "$REPO" and curl failed with "Bad hostname". Docker substitutes the names
# it knows (ARG/ENV, so CRAN_SNAPSHOT below) and leaves the rest for the shell, which is
# where `. /etc/os-release` puts UBUNTU_CODENAME.
#
# The curl is the point of the extra line: a wrong codename yields a URL that 404s into a
# silent source-build fallback -- same image, forty minutes instead of three, no error
# anywhere. Checking the index turns that into a build failure.
ARG CRAN_SNAPSHOT=2026-09-01
RUN . /etc/os-release && \
    REPO="https://packagemanager.posit.co/cran/__linux__/$UBUNTU_CODENAME/${CRAN_SNAPSHOT}" && \
    echo "using CRAN repo: $REPO" && \
    curl -fsS -o /dev/null "$REPO/src/contrib/PACKAGES.gz" && \
    echo "options(repos = c(CRAN = '$REPO'))" >> "$R_HOME/etc/Rprofile.site"

# The nine packages the running app loads. `viper` is deliberately absent: VIPER activity
# is precomputed into viper_*.h5 by the pipeline, so the served app never needs it, and
# leaving it out keeps the whole Bioconductor single-cell dependency chain out of the image.
RUN install2.r --error --skipinstalled --ncpus -1 \
        shiny \
        DT \
        shinythemes \
        RSQLite \
        DBI \
        data.table \
        metafor \
        survival

# rhdf5 is the one Bioconductor dependency. Rhdf5lib vendors HDF5 itself, so no system
# libhdf5 is required and the version is fixed by the Bioconductor release, not by apt.
ARG BIOC_VERSION=3.22
RUN install2.r --error BiocManager && \
    R -q -e "BiocManager::install('rhdf5', version = '${BIOC_VERSION}', ask = FALSE, update = FALSE)"

# R.utils, and it gets its own line because the REASON is the interesting part. Nothing in
# this codebase names it: `data.table::fread()` needs it to read a .gz, and it asks for it at
# CALL TIME. R/cptac.R reads cptac_proteomics.tsv.gz during app startup, so without this the
# image builds cleanly, passes the self-check, and then dies on boot with an error about a
# package no static scan of the source could ever have found. It was found by starting the
# app, which is why the build now does that too (see the smoke test below).
RUN install2.r --error --skipinstalled --ncpus -1 R.utils

WORKDIR /srv/omicohort

# Data first, code second: the matrices are ~1.76 GB and change once a release, the code
# changes constantly. In this order an edit to app.R rebuilds a 200 KB layer instead of
# re-copying the whole data tree.
COPY data/      /srv/omicohort/data/
COPY networks/  /srv/omicohort/networks/
COPY results/   /srv/omicohort/results/
COPY config/    /srv/omicohort/config/
COPY www/       /srv/omicohort/www/
COPY R/         /srv/omicohort/R/
COPY scripts/docker_selfcheck.R /srv/omicohort/scripts/
COPY .Rprofile app.R /srv/omicohort/

# FAIL THE BUILD, not the visitor. .dockerignore is a whitelist and Docker's re-include
# rules are subtle; every way of getting one wrong produces an image that builds, starts,
# and serves a page with one cohort quietly missing. This derives the required set from
# app.R and config/cohorts.tsv and names anything absent. It also proves .Rprofile is read
# from this WORKDIR and that png() has a device -- both invisible until a user hits them.
RUN Rscript scripts/docker_selfcheck.R

# Run unprivileged. The app writes only to tempdir(), so the tree stays root-owned and
# read-only to the process -- which is the intent, not an oversight: nothing in
# /srv/omicohort should be mutable at runtime.
RUN useradd --create-home --shell /usr/sbin/nologin --uid 10001 omicohort
USER omicohort

# START THE APP AND FETCH A PAGE. The self-check above proves the FILES are present and the
# nine declared packages are installed; it cannot prove the app boots, because the thing that
# broke it -- data.table reaching for R.utils at call time -- is named nowhere in the source.
# A list of dependencies can only ever be as complete as the reading that produced it. Serving
# one real request is the check that does not depend on my reading being right.
RUN set -eu; \
    Rscript -e "shiny::runApp('.', port = 7654, host = '127.0.0.1', launch.browser = FALSE)" \
      > /tmp/smoke.log 2>&1 & \
    pid=$!; \
    ok=0; \
    for i in $(seq 1 60); do \
      if curl -fsS -o /tmp/smoke.html http://127.0.0.1:7654/ 2>/dev/null; then ok=1; break; fi; \
      kill -0 "$pid" 2>/dev/null || break; \
      sleep 2; \
    done; \
    kill "$pid" 2>/dev/null || true; \
    if [ "$ok" -ne 1 ]; then echo "SMOKE TEST FAILED -- the app never served a page:"; cat /tmp/smoke.log; exit 1; fi; \
    grep -q "OMICohort" /tmp/smoke.html || { echo "served a page that is not OMICohort:"; head -40 /tmp/smoke.html; exit 1; }; \
    echo "smoke test OK: served $(wc -c < /tmp/smoke.html) bytes containing the app title"

EXPOSE 7654

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
    CMD curl -fsS "http://127.0.0.1:${PORT:-7654}/" > /dev/null || exit 1

# R must START here for .Rprofile (shiny.autoload.r = FALSE) to be read and for the app's
# ROOT resolution -- getOption("omicohort.root") then $OMICOHORT_ROOT then getwd() -- to
# land on this directory. Do not add --vanilla.
CMD ["Rscript", "-e", "shiny::runApp('.', port = as.integer(Sys.getenv('PORT', '7654')), host = '0.0.0.0', launch.browser = FALSE)"]
