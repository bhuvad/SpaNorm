#!/bin/sh
# Build the generated reference pages (git-ignored): R pages from man/*.Rd via altdoc into
# reference/r/, Python pages via quartodoc into reference/python/. Run from anywhere:
#   bash site/scripts/build_reference.sh && quarto render site
# Local settings (RETICULATE_PYTHON) come from _environment.local; PYTHON overrides them.
set -e
cd "$(dirname "$0")/.."

if [ -f _environment.local ]; then
  # shellcheck disable=SC1091
  . ./_environment.local 2>/dev/null || true
fi
PYTHON="${PYTHON:-${RETICULATE_PYTHON:-python}}"

# Python tools: pip install -r requirements.txt (quartodoc pinned; needs griffe < 2)
# altdoc goes into the project-local library, never the user's
mkdir -p _rlib
# A FIXED altdoc version: the Rd conversion uses its private altdoc:::.rd2qmd.
ALTDOC_VERSION=0.7.3
export ALTDOC_VERSION
Rscript -e 'v <- Sys.getenv("ALTDOC_VERSION"); ok <- requireNamespace("altdoc", lib.loc = "_rlib", quietly = TRUE) && as.character(packageVersion("altdoc", lib.loc = "_rlib")) == v; if (!ok) { if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org"); remotes::install_version("altdoc", version = v, lib = "_rlib", repos = "https://cloud.r-project.org", upgrade = "never") }'
Rscript scripts/build-reference-r.R

# quartodoc needs griffe < 2 (quartodoc 0.11 does not support griffe 2): see requirements.txt
"$PYTHON" -m quartodoc build --config _quarto.yml

# Quarto resolves {{< include >}} before its pre-render step runs, so copy the included
# files (spanormpy docs, NEWS.md) here too. SPANORMPY_DIR: see copy-assets.sh.
sh scripts/copy-assets.sh

# every R export and every spanormpy.__all__ name must have a generated page
"$PYTHON" scripts/check_reference_coverage.py
