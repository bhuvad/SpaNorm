#!/bin/sh
# Build the generated reference pages (git-ignored): R pages from man/*.Rd via altdoc into
# reference/r/, Python pages via quartodoc into reference/python/, after building the shared fitted
# objects the examples load (reference-examples.yml decides which examples run; they execute during
# `quarto render`). Run from anywhere:
#   bash site/scripts/build_reference.sh && quarto render site
# Local settings (RETICULATE_PYTHON) come from _environment.local; PYTHON overrides them.
set -e
cd "$(dirname "$0")/.."

if [ -f _environment.local ]; then
  # exported, so that Rscript/python children see RETICULATE_PYTHON and LD_PRELOAD too
  set -a
  # shellcheck disable=SC1091
  . ./_environment.local 2>/dev/null || true
  set +a
fi
PYTHON="${PYTHON:-${RETICULATE_PYTHON:-python}}"

# Python tools: pip install -r requirements.txt (quartodoc pinned; needs griffe < 2)
# altdoc goes into the project-local library, never the user's
mkdir -p _rlib
# A FIXED altdoc version: the Rd conversion uses its private altdoc:::.rd2qmd.
ALTDOC_VERSION=0.7.3
export ALTDOC_VERSION
Rscript -e 'v <- Sys.getenv("ALTDOC_VERSION"); ok <- requireNamespace("altdoc", lib.loc = "_rlib", quietly = TRUE) && as.character(packageVersion("altdoc", lib.loc = "_rlib")) == v; if (!ok) { if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org"); remotes::install_version("altdoc", version = v, lib = "_rlib", repos = "https://cloud.r-project.org", upgrade = "never") }'

# phase timer: prints the seconds each phase took
t_last=$(date +%s)
phase() { t=$(date +%s); echo "== phase '$1': $((t - t_last)) s"; t_last=$t; }

# The shared fitted objects, built ONCE per site build and loaded by every executed example
# (R: _cache/HumanDLPFC-fit.rds; Python: _cache/HumanDLPFC-fit.h5ad, same data and genes).
# Files in _cache/ are reused within a build and across local builds; delete them to refit.
mkdir -p _cache
[ -f _cache/HumanDLPFC-fit.rds ] || Rscript scripts/build-shared-fit-r.R
phase "shared R fit"
[ -f _cache/HumanDLPFC.h5ad ] || Rscript scripts/export-h5ad.R _cache/HumanDLPFC.h5ad
[ -f _cache/HumanDLPFC-fit.h5ad ] || "$PYTHON" scripts/build-shared-fit-py.py
phase "shared Python fit"

Rscript scripts/build-reference-r.R
phase "R reference pages (altdoc)"

# quartodoc needs griffe < 2 (quartodoc 0.11 does not support griffe 2): see requirements.txt
"$PYTHON" -m quartodoc build --config _quarto.yml
# turn the docstring Examples into executable cells (reference-examples.yml decides which run)
"$PYTHON" scripts/build-reference-py.py
phase "Python reference pages (quartodoc)"

# Quarto resolves {{< include >}} before its pre-render step runs, so copy the included
# files (spanormpy docs, NEWS.md) here too. SPANORMPY_DIR: see copy-assets.sh.
sh scripts/copy-assets.sh

# every R export and every spanormpy.__all__ name must have a generated page
"$PYTHON" scripts/check_reference_coverage.py
