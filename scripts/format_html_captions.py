#!/usr/bin/env python3
"""Normalize rendered Quarto figure/table captions for the HTML site.

Quarto controls the caption label and numbering, but the rendered HTML exposes
the label as plain text. This post-render step keeps the generated pages static
while enforcing the house style:

  <strong>Figura X.</strong> caption text
  <strong>Tabla X.</strong> caption text
"""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"

LABEL = r"(?:Figura|Tabla)&nbsp;[A-Za-z0-9]+(?:\.[0-9]+)*"
CAPTION_LABEL_RE = re.compile(
    rf"(<figcaption\b[^>]*>\s*)(?!<strong>)(?P<label>{LABEL})(?::|\.)\s*(?:\.\s*)?",
    flags=re.MULTILINE,
)
STRONG_DUPLICATE_DOT_RE = re.compile(
    rf"(<strong>{LABEL}\.</strong>)\s+\.\s*"
)
TITLE_LABEL_RE = re.compile(rf'(title="(?P<label>{LABEL})(?::|\.)\s*(?:\.\s*)?)')
SEARCH_LABEL_RE = re.compile(
    r"(?P<label>(?:Figura|Tabla)\u00a0[A-Za-z0-9]+(?:\.[0-9]+)*)(?::|\.)\s*(?:\.\s*)?"
)


def normalize_html(text: str) -> str:
    text = CAPTION_LABEL_RE.sub(
        lambda match: f"{match.group(1)}<strong>{match.group('label')}.</strong> ",
        text,
    )
    text = STRONG_DUPLICATE_DOT_RE.sub(lambda match: f"{match.group(1)} ", text)
    text = TITLE_LABEL_RE.sub(
        lambda match: f'title="{match.group("label")}. ',
        text,
    )
    return text


def normalize_search_json(text: str) -> str:
    return SEARCH_LABEL_RE.sub(lambda match: f"{match.group('label')}. ", text)


def rewrite(path: Path, transform) -> bool:
    original = path.read_text(encoding="utf-8")
    updated = transform(original)
    if updated == original:
        return False
    path.write_text(updated, encoding="utf-8")
    return True


def main() -> None:
    if not DOCS.exists():
        raise SystemExit(f"Rendered output directory not found: {DOCS}")

    changed = []
    for path in sorted(DOCS.rglob("*.html")):
        if rewrite(path, normalize_html):
            changed.append(path)

    search_json = DOCS / "search.json"
    if search_json.exists() and rewrite(search_json, normalize_search_json):
        changed.append(search_json)

    print(f"[caption-format] normalized captions in {len(changed)} rendered file(s)")


if __name__ == "__main__":
    main()
