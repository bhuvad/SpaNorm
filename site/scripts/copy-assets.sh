#!/bin/sh
# Pre-render step (run by quarto from the site/ directory). Copies files that live outside the
# project into it, since Quarto only publishes files inside the project. All copies are
# git-ignored.
#
# SPANORMPY_DIR: checkout of https://github.com/BhuvaLab/spanormpy holding docs/r-to-python.md and
# docs/r-discrepancies.md. Default: a sibling of this repository (../spanormpy relative to the
# repository root), which is where it lives locally. CI checks it out elsewhere and sets the variable.
set -e
SPANORMPY_DIR="${SPANORMPY_DIR:-../../spanormpy}"   # relative to site/

cp ../vignettes/SpaNormWorkflow.png SpaNormWorkflow.png
cp ../man/figures/logo.png logo.png
# total Bioconductor downloads, computed at build time (replaces the rank badge)
python3 scripts/download-stats.py

mkdir -p _includes
cp ../NEWS.md _includes/NEWS.md
for f in r-to-python r-discrepancies; do
  if [ ! -f "$SPANORMPY_DIR/docs/$f.md" ]; then
    echo "copy-assets.sh: $SPANORMPY_DIR/docs/$f.md not found; set SPANORMPY_DIR to a spanormpy checkout" >&2
    exit 1
  fi
done
cp "$SPANORMPY_DIR/docs/r-to-python.md" _includes/r-to-python.md
cp "$SPANORMPY_DIR/docs/r-discrepancies.md" _includes/r-discrepancies.md
# the two pages are shown on one site page, so the cross-link becomes an anchor
sed -i 's|(r-discrepancies.md)|(#discrepancies-from-r-spanorm)|' _includes/r-to-python.md
