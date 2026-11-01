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

# Keep the three connections separate: nested pagination would otherwise miss
# later comment/reply pages, and large nested queries exceed GitHub's node limit.
DISCUSSIONS_QUERY = """
query($owner: String!, $name: String!, $cursor: String) {
  repository(owner: $owner, name: $name) {
    discussions(first: 100, after: $cursor) {
      nodes { id title body html_url: url updated_at: updatedAt }
      pageInfo { hasNextPage endCursor }
    }
  }
}
"""
COMMENTS_QUERY = """
query($id: ID!, $cursor: String) {
  node(id: $id) {
    ... on Discussion {
      comments(first: 100, after: $cursor) {
        nodes { id body html_url: url updated_at: updatedAt replies { totalCount } }
        pageInfo { hasNextPage endCursor }
      }
    }
  }
}
"""
REPLIES_QUERY = """
query($id: ID!, $cursor: String) {
  node(id: $id) {
    ... on DiscussionComment {
      replies(first: 100, after: $cursor) {
        nodes { id body html_url: url updated_at: updatedAt }
        pageInfo { hasNextPage endCursor }
      }
    }
  }
}
"""


def github_graphql_pages(query: str, variables: dict, path: tuple[str, ...]) -> list[dict]:
    """Read one connection completely, rejecting errors and stalled cursors."""
    records = []
    cursor = None
    seen_cursors = set()
    while True:
        result = subprocess.run(
            ["gh", "api", "graphql", "--input", "-"],
            input=json.dumps({"query": query, "variables": {**variables, "cursor": cursor}}),
            check=True, capture_output=True, text=True, timeout=45,
        )
        response = json.loads(result.stdout)
        if not isinstance(response, dict):
            raise ValueError("unexpected GitHub GraphQL response")
        if response.get("errors"):
            raise ValueError(f"GitHub GraphQL errors: {response['errors']}")
        try:
            connection = response["data"]
            for key in path:
                connection = connection[key]
            nodes = connection["nodes"]
            info = connection["pageInfo"]
            has_next = info["hasNextPage"]
            next_cursor = info["endCursor"]
        except (KeyError, TypeError) as exc:
            raise ValueError("missing GitHub GraphQL connection or page info") from exc
        if (not isinstance(nodes, list) or any(not isinstance(node, dict) for node in nodes)
                or not isinstance(has_next, bool)):
            raise ValueError("unexpected GitHub GraphQL connection shape")
        records.extend(nodes)
        if not has_next:
            return records
        if not isinstance(next_cursor, str) or not next_cursor or next_cursor in seen_cursors:
            raise ValueError("GitHub GraphQL pagination did not advance")
        seen_cursors.add(next_cursor)
        cursor = next_cursor


def github_discussions() -> tuple[list[dict], list[dict]]:
    """Return discussion posts and all comments/replies, even for nonmatching parents."""
    owner, name = REPO.split("/", 1)
    discussions = github_graphql_pages(
        DISCUSSIONS_QUERY, {"owner": owner, "name": name}, ("repository", "discussions"))
    comments = []
    for discussion in discussions:
        parents = github_graphql_pages(
            COMMENTS_QUERY, {"id": discussion["id"]}, ("node", "comments"))
        comments.extend(parents)
        for comment in parents:
            if comment["replies"]["totalCount"]:
                comments.extend(github_graphql_pages(
                    REPLIES_QUERY, {"id": comment["id"]}, ("node", "replies")))
    return discussions, comments


def github_pages(endpoint: str) -> list[dict]:
    result = subprocess.run(
        ["gh", "api", "--paginate", endpoint],
        check=True, capture_output=True, text=True, timeout=45,
    )
    # --paginate emits successive JSON values. Parse that stream directly so
    # this also works with CLI releases predating the optional --slurp flag.
    decoder = json.JSONDecoder()
    remaining = result.stdout.lstrip()
    records = []
    if not remaining:
        raise ValueError("empty GitHub pagination response")
    while remaining:
        page, end = decoder.raw_decode(remaining)
        if not isinstance(page, list) or any(not isinstance(item, dict) for item in page):
            raise ValueError("unexpected GitHub pagination response")
        records.extend(page)
        remaining = remaining[end:].lstrip()
    return records


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
    discussions, discussion_comments = github_discussions()
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
        (f"Discussion posts ({len(discussions)} scanned)", records(discussions)),
        (f"Discussion comments and replies ({len(discussion_comments)} scanned)", records(discussion_comments)),
        (f"Historical issue #8 comments ({len(comments)} scanned)", records(comments)),
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
    print("All accessible discussion posts, comments and replies are scanned, including closed discussions.")
    print("PR comments, other issue comments, closed issues/PRs, other worktrees and non-Markdown files are not scanned.")
    print("This is a live paginated scan, not an atomic snapshot of concurrent edits.")
    for title, found in groups:
        print(f"\n{title}: {len(found)} matches")
        for label, excerpt in found[:DISPLAY_LIMIT]:
            print(f"- {label}\n  {excerpt}")
        if len(found) > DISPLAY_LIMIT:
            print(f"  ... {len(found) - DISPLAY_LIMIT} further matches omitted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
