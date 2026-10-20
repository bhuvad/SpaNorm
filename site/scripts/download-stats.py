#!/usr/bin/env python3
"""Write site/downloads.svg: a badge with SpaNorm's total Bioconductor downloads.

The count is computed at build time by summing the per-year totals in Bioconductor's published
statistics (https://bioconductor.org/packages/stats/bioc/SpaNorm/SpaNorm_stats.tab), instead of
showing Bioconductor's rank badge. If the statistics cannot be fetched the badge says
"unavailable" and the build carries on. Standard library only. Run from the site/ directory
(scripts/copy-assets.sh does this before every render).
"""
import csv
import datetime
import io
import sys
import urllib.request

PKG = "SpaNorm"
URL = f"https://bioconductor.org/packages/stats/bioc/{PKG}/{PKG}_stats.tab"
OUT = "downloads.svg"


def total_downloads():
    # bioconductor.org answers 403 to urllib's default User-Agent
    req = urllib.request.Request(URL, headers={"User-Agent": "SpaNorm-site-build (https://github.com/BhuvaLab/SpaNorm)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        text = r.read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    # one summary row per year (Month == "all"); the monthly rows add up to the same total
    yearly = [int(x["Nb_of_downloads"]) for x in rows if x["Month"] == "all"]
    total = sum(yearly) if yearly else sum(int(x["Nb_of_downloads"]) for x in rows)
    years = sorted({x["Year"] for x in rows})
    return total, years


def badge(label, value, colour, title):
    w_label = int(len(label) * 6.3 + 12)
    w_value = int(len(value) * 6.6 + 12)
    w = w_label + w_value
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="20" role="img" aria-label="{label}: {value}">
<title>{title}</title>
<linearGradient id="s" x2="0" y2="100%"><stop offset="0" stop-color="#bbb" stop-opacity=".1"/><stop offset="1" stop-opacity=".1"/></linearGradient>
<clipPath id="r"><rect width="{w}" height="20" rx="3" fill="#fff"/></clipPath>
<g clip-path="url(#r)"><rect width="{w_label}" height="20" fill="#555"/><rect x="{w_label}" width="{w_value}" height="20" fill="{colour}"/><rect width="{w}" height="20" fill="url(#s)"/></g>
<g fill="#fff" text-anchor="middle" font-family="Verdana,Geneva,DejaVu Sans,sans-serif" font-size="11">
<text x="{w_label / 2:.1f}" y="14">{label}</text><text x="{w_label + w_value / 2:.1f}" y="14">{value}</text></g>
</svg>
"""


def main():
    label = "Bioconductor downloads"
    try:
        total, years = total_downloads()
        value = f"{total:,} total"
        title = (f"{PKG}: {total:,} downloads from Bioconductor ({years[0]}-{years[-1]}), "
                 f"computed {datetime.date.today().isoformat()}")
        svg = badge(label, value, "#0c6fb4", title)
        print(f"download-stats: {PKG} {total:,} Bioconductor downloads ({years[0]}-{years[-1]})")
    except Exception as e:  # network or format problems must not break the site build
        print(f"download-stats: could not compute downloads ({e}); writing placeholder", file=sys.stderr)
        svg = badge(label, "unavailable", "#9f9f9f", f"{PKG} Bioconductor download statistics unavailable")
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(svg)


if __name__ == "__main__":
    main()
