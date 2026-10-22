"""Check the RENDERED reference pages against reference-examples.yml (run after `quarto render`).

* every page listed `run: true` must show at least one output block (unless `output: optional`)
  and no error output, so a page can never silently render with an empty or failed example;
* every page listed `run: false` must carry the visible "shown but not run" note;
* every listed page must exist in the rendered site.

Usage: python scripts/check_reference_outputs.py [site_dir]   (default: _site next to scripts/)
Exits 1 and lists the problems if any.
"""
import sys
from pathlib import Path

import yaml

site = Path(__file__).resolve().parent.parent
root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else site / "_site"
spec = yaml.safe_load((site / "reference-examples.yml").read_text())
problems, n_out = [], 0

for lang, folder in (("r", "r"), ("python", "python")):
    for page, sp in spec[lang].items():
        f = root / "reference" / folder / f"{page}.html"
        if not f.exists():
            problems.append(f"{lang} {page}: rendered page {f.relative_to(root)} not found")
            continue
        html = f.read_text(encoding="utf-8", errors="replace")
        if "cell-output-error" in html or "Traceback (most recent call last)" in html:
            problems.append(f"{lang} {page}: error output on the page")
        if sp.get("run"):
            has_out = "cell-output" in html
            n_out += has_out
            if not has_out and sp.get("output") != "optional":
                problems.append(f"{lang} {page}: run: true but the page shows no output")
        elif "shown but not run when the site is built" not in html and "has no example" not in sp.get("reason", ""):
            problems.append(f"{lang} {page}: excluded but the page lacks the 'not run' note")

if problems:
    print("\n".join(problems), file=sys.stderr)
    sys.exit(1)
print(f"reference outputs OK: {n_out} pages show example output, "
      f"{sum(not sp.get('run') for l in ('r', 'python') for sp in spec[l].values())} excluded")
