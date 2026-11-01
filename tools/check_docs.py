#!/usr/bin/env python3
"""Find missing local Markdown link targets in the tracked repository docs.

This check is offline. It ignores external URLs, file fragments, code fences
and display-math blocks delimited by \\[ ... \\] or double dollar signs.
"""

import argparse
import posixpath
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

INLINE_LINK = re.compile(r"!?\[[^\]\n]+\]\(\s*(<[^>]+>|[^\s)]+)")
REFERENCE_LINK = re.compile(r"^\s{0,3}\[[^\]]+\]:\s*(<[^>]+>|\S+)")
FENCE = re.compile(r"^\s*(" + chr(96) + r"{3,}|~{3,})")


def tracked_files(root: Path) -> set[str]:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"],
        check=True,
        stdout=subprocess.PIPE,
    )
    return {
        p.decode("utf-8", errors="surrogateescape")
        for p in result.stdout.split(b"\0") if p
    }


def link_targets(document: str, text: str, files: set[str]) -> list[str]:
    """Return path and line diagnostics for missing local links."""
    directories = {"."}
    for name in files:
        parent = posixpath.dirname(name)
        while parent and parent != ".":
            directories.add(parent)
            parent = posixpath.dirname(parent)
    missing = []
    fence_char = None
    fence_width = 0
    math_end = None
    for line_number, line in enumerate(text.splitlines(), 1):
        if math_end is not None:
            stripped = line.lstrip()
            if not stripped.startswith(math_end):
                continue
            line = stripped[len(math_end):]
            math_end = None
        marker = FENCE.match(line)
        if marker:
            fence = marker.group(1)
            if fence_char is None:
                fence_char, fence_width = fence[0], len(fence)
            elif fence[0] == fence_char and len(fence) >= fence_width:
                fence_char = None
            continue
        if fence_char is not None:
            continue
        if line.strip() in (r"\[", "$$"):
            math_end = r"\]" if line.strip() == r"\[" else "$$"
            continue
        urls = [m.group(1) for m in INLINE_LINK.finditer(line)]
        ref = REFERENCE_LINK.match(line)
        if ref:
            urls.append(ref.group(1))
        for raw in urls:
            url = raw.strip("<>")
            parsed = urlsplit(url)
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            path = unquote(parsed.path)
            if path.startswith("/"):
                continue  # absolute site paths are not repository file links
            target = posixpath.normpath(posixpath.join(posixpath.dirname(document), path))
            if target.startswith("../") or target == "..":
                missing.append(f"{document}:{line_number}: outside repository: {url}")
            elif target not in files and target not in directories:
                missing.append(f"{document}:{line_number}: missing {target} (from {url})")
    return missing


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        files = tracked_files(root)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"Cannot enumerate tracked files: {exc}", file=sys.stderr)
        return 2
    docs = sorted(name for name in files if name.endswith(".md"))
    errors = []
    for name in docs:
        text = (root / name).read_text(encoding="utf-8")
        errors.extend(link_targets(name, text, files))
    if errors:
        for err in errors:
            print(err, file=sys.stderr)
        print(f"Found {len(errors)} missing links in {len(docs)} Markdown files.", file=sys.stderr)
        return 1
    print(f"Checked local links in {len(docs)} Markdown files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
