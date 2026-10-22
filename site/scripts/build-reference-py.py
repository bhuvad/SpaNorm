"""Make the Python reference pages show the output of their examples.

Run after `quartodoc build` (build_reference.sh does). quartodoc renders a docstring's
`Examples` section as a plain ```python block with the `>>>` prompts. This turns each example
(blocks of `>>>` lines separated by blank lines) into executable ```{python} cells, adds a
hidden setup cell that loads the shared fitted AnnData (built once per site build by
build-shared-fit-py.py) as `adata`, sets the figure size, and prints a note about it. Which pages
run is decided by reference-examples.yml (python:); a page that is listed but has no Examples
section, or a page that is not listed, fails the build. A failing example fails the render
(Quarto stops and names the page).

Usage: python scripts/build-reference-py.py     (any working directory)
"""
import re
import sys
from pathlib import Path

import yaml

site = Path(__file__).resolve().parent.parent
pages_dir = site / "reference" / "python"
spec = yaml.safe_load((site / "reference-examples.yml").read_text())["python"]

pages = sorted(p.stem for p in pages_dir.glob("*.qmd") if p.stem != "index")
problems = [f"not listed in reference-examples.yml (python): {p}" for p in pages if p not in spec]
problems += [f"listed in reference-examples.yml but no such page: {p}" for p in spec if p not in pages]

FRONT = """---
execute:
  eval: true
  error: false
  fig-width: 7
  fig-height: 5
  fig-dpi: 144
---

"""
SETUP = '''```{{python}}
#| include: false
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
plt.rcParams.update({{"figure.dpi": 144, "savefig.bbox": "tight"}})
{shared}```

'''
SHARED = '''import anndata as ad
# the shared fitted AnnData, built once per site build (scripts/build-shared-fit-py.py)
adata = ad.read_h5ad("../../_cache/HumanDLPFC-fit.h5ad")
'''
NOTE_SHARED = (
    "::: {.callout-note}\nThe code below runs when the site is built, on a small pre-fitted subset of "
    "`HumanDLPFC` (301 genes x 4015 spots; `spanorm(sample_p=0.05, df_tps=2)`, `spanorm_svg()` and "
    "`spanorm_pca()` already applied) that is fitted once per build and loaded as `adata`. The outputs "
    "are those of the code shown.\n:::\n\n"
)


def to_cells(block: str) -> str:
    """`>>>`-prompted lines -> executable cells, one per blank-line separated group."""
    cells, cur = [], []
    for line in block.split("\n"):
        if not line.strip():
            if cur:
                cells.append(cur)
                cur = []
            continue
        m = re.match(r"^(>>> ?|\.\.\. ?)(.*)$", line)
        cur.append(m.group(2) if m else line)
    if cur:
        cells.append(cur)
    return "\n\n".join("```{python}\n" + "\n".join(c) + "\n```" for c in cells)


for page in pages:
    if page not in spec:
        continue
    path = pages_dir / f"{page}.qmd"
    text = path.read_text(encoding="utf-8")
    if "```{python}\n#| include: false\nimport warnings" in text:
        continue  # already processed
    sp = spec[page]
    m = re.search(r"^## Examples[^\n]*\n", text, re.M)
    if m is None:
        if sp.get("run"):
            problems.append(f"{page}: run: true but the docstring has no Examples section")
        else:
            path.write_text(text, encoding="utf-8")
        continue
    head, tail = text[: m.end()], text[m.end():]
    nxt = re.search(r"^## ", tail, re.M)
    body, rest = (tail[: nxt.start()], tail[nxt.start():]) if nxt else (tail, "")
    if sp.get("run"):
        body = re.sub(r"```python\n(.*?)\n```", lambda mm: to_cells(mm.group(1)), body, flags=re.S)
        note = NOTE_SHARED if sp.get("shared") else ""
        text = FRONT + SETUP.format(shared=SHARED if sp.get("shared") else "") + head + "\n" + note + body + rest
    else:
        reason = sp.get("reason")
        if not reason:
            problems.append(f"{page}: run: false needs a reason")
            continue
        note = f"::: {{.callout-note}}\nThis example is shown but not run when the site is built: {reason}.\n:::\n\n"
        text = "---\nexecute:\n  eval: false\n---\n\n" + head + "\n" + note + body + rest
    path.write_text(text, encoding="utf-8")

if problems:
    print("\n".join(problems), file=sys.stderr)
    sys.exit(1)
print(f"python reference examples prepared: {sum(bool(spec[p].get('run')) for p in pages)} run, "
      f"{sum(not spec[p].get('run') for p in pages)} not run")
