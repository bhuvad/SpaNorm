"""Check that every internal link and image in the rendered site (_site/*.html) resolves.

Usage: python site/scripts/check_links.py [site_dir]   (default: _site next to this script's parent)
Exits 1 and lists the broken links if any.
"""
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlparse

root = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "_site").resolve()


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        for key in ("href", "src"):
            if a.get(key):
                self.urls.append(a[key])


broken = []
n = 0
for page in sorted(root.rglob("*.html")):
    parser = Links()
    parser.feed(page.read_text(encoding="utf-8", errors="replace"))
    for url in parser.urls:
        u = urlparse(url)
        if u.scheme or u.netloc or not u.path or url.startswith(("#", "mailto:", "data:")):
            continue
        n += 1
        target = (root / u.path.lstrip("/")) if u.path.startswith("/") else (page.parent / unquote(u.path))
        target = target.resolve()
        if target.is_dir():
            target = target / "index.html"
        if not target.exists():
            broken.append(f"{page.relative_to(root)}: {url}")

print(f"checked {n} internal links in {len(list(root.rglob('*.html')))} pages")
if broken:
    print(f"{len(broken)} broken:")
    print("\n".join(sorted(set(broken))))
    sys.exit(1)
