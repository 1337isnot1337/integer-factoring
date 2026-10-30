"""Find possible overlap in current GitHub work and this checkout's research notes.

This read-only search returns leads, not a novelty certificate or reservation.
It requires the GitHub CLI and repository access for the remote sources.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPO = "1337isnot1337/integer-factoring"
DISPLAY_LIMIT = 30


def github_pages(endpoint: str) -> list[dict]:
    result = subprocess.run(
        ["gh", "api", "--paginate", "--slurp", endpoint],
        check=True, capture_output=True, text=True, timeout=45,
    )
    pages = json.loads(result.stdout)
    if not isinstance(pages, list) or any(not isinstance(page, list) for page in pages):
        raise ValueError("unexpected GitHub pagination response")
    return [item for page in pages for item in page]


def matches(text: str, terms: list[str]) -> bool:
    text = text.casefold()
    return any(term.casefold() in text for term in terms)


def snippet(text: str, terms: list[str]) -> str:
    for line in text.splitlines():
        if matches(line, terms):
            return " ".join(line.split())[:240]
    return ""


def search(terms: list[str], root: Path = ROOT) -> list[tuple[str, list[tuple[str, str]]]]:
    pulls = github_pages(f"repos/{REPO}/pulls?state=open&per_page=100")
    issues = github_pages(f"repos/{REPO}/issues?state=open&per_page=100")
    comments = github_pages(f"repos/{REPO}/issues/8/comments?per_page=100")
    issues = [issue for issue in issues if "pull_request" not in issue]
    notes = sorted((root / "research").rglob("*.md"))

    def records(items: list[dict]) -> list[tuple[str, str]]:
        result = []
        for item in items:
            text = (item.get("title") or "") + "\n" + (item.get("body") or "")
            if matches(text, terms):
                label = f"{item['html_url']} (updated {item['updated_at']})"
                result.append((label, snippet(text, terms)))
        return result

    local = []
    for path in notes:
        content = path.read_text(encoding="utf-8")
        label = str(path.relative_to(root))
        if matches(label + "\n" + content, terms):
            local.append((label, snippet(content, terms)))
    return [
        (f"Open pull requests ({len(pulls)} scanned)", records(pulls)),
        (f"Open issues ({len(issues)} scanned)", records(issues)),
        (f"Coordination comments ({len(comments)} scanned)", records(comments)),
        (f"Local research Markdown files ({len(notes)} scanned)", local),
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("terms", nargs="+", help="case-insensitive phrases; any match is shown")
    args = parser.parse_args()
    if any(not term.strip() for term in args.terms):
        parser.error("search terms must not be empty")
    try:
        groups = search(args.terms)
    except (OSError, subprocess.SubprocessError, ValueError, KeyError) as exc:
        print(f"Overlap search incomplete: {exc}", file=sys.stderr)
        return 1
    print(f"Snapshot: {datetime.now(timezone.utc).isoformat()}")
    print("Matches are leads, not active claims or proof of novelty; check linked updates.")
    print("PR comments, closed work, other worktrees and non-Markdown files are not scanned.")
    for title, found in groups:
        print(f"\n{title}: {len(found)} matches")
        for label, excerpt in found[:DISPLAY_LIMIT]:
            print(f"- {label}\n  {excerpt}")
        if len(found) > DISPLAY_LIMIT:
            print(f"  ... {len(found) - DISPLAY_LIMIT} further matches omitted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
