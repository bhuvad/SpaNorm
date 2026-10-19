"""Check that every R export (NAMESPACE) and every spanormpy.__all__ name has a generated
reference page. Run after the reference pages are built (build_reference.sh calls it).

Usage: python scripts/check_reference_coverage.py     (any working directory)
Exits 1 and lists the missing pages if any.
"""
import importlib
import re
import sys
from pathlib import Path

site = Path(__file__).resolve().parent.parent
repo = site.parent
missing = []

# R: every export() must be an \alias{} of an Rd file that has a page in reference/r/
alias_page = {}
for rd in (repo / "man").glob("*.Rd"):
    for alias in re.findall(r"^\\alias\{([^}]+)\}", rd.read_text(encoding="utf-8"), re.M):
        alias_page[alias] = rd.stem
exports = re.findall(r"^export\(([^)]+)\)", (repo / "NAMESPACE").read_text(), re.M)
for name in exports:
    name = name.strip('"')
    page = alias_page.get(name)
    if page is None or not (site / "reference" / "r" / f"{page}.qmd").exists():
        missing.append(f"R export {name}: no page reference/r/{page or '<no Rd alias>'}.qmd")

# Python: every name in __all__ (and in the submodule `pl`'s __all__) has a page
pkg = importlib.import_module("spanormpy")


def need(prefix, names):
    for name in names:
        obj = getattr(importlib.import_module(f"spanormpy.{prefix}".rstrip(".")) if prefix else pkg, name)
        if hasattr(obj, "__path__") or type(obj).__name__ == "module":
            need(name, getattr(obj, "__all__", []))
        else:
            page = f"{prefix + '.' if prefix else ''}{name}.qmd"
            if not (site / "reference" / "python" / page).exists():
                missing.append(f"Python {prefix + '.' if prefix else ''}{name}: no page reference/python/{page}")


need("", pkg.__all__)

if missing:
    print("\n".join(missing), file=sys.stderr)
    sys.exit(1)
print(f"reference coverage OK: {len(exports)} R exports, spanormpy.__all__ = {len(pkg.__all__)} names")
