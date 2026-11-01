"""Offline regressions for the live coordination-board search."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from tools import find_overlap


def page(path, nodes, cursor=None):
    data = {"nodes": nodes, "pageInfo": {
        "hasNextPage": cursor is not None, "endCursor": cursor}}
    for key in reversed(path):
        data = {key: data}
    return subprocess.CompletedProcess([], 0, json.dumps({"data": data}), "")


def record(identifier, body, **extra):
    return dict(id=identifier, body=body, html_url=f"https://example/{identifier}",
                updated_at="2026-11-01", **extra)


class DiscussionOverlapTests(unittest.TestCase):
    def test_search_finds_later_discussion_comment_and_reply_pages(self):
        # Neither parent matches: a reply-only claim must still be discovered.
        responses = [
            page(("repository", "discussions"), [record("d0", "unrelated")], "d-next"),
            page(("repository", "discussions"), [record("d1", "board")]),
            page(("node", "comments"), []),
            page(("node", "comments"), [record("c0", "unrelated", replies={"totalCount": 0})], "c-next"),
            page(("node", "comments"), [record("c1", "board", replies={"totalCount": 2})]),
            page(("node", "replies"), [record("r0", "unrelated")], "r-next"),
            page(("node", "replies"), [record("r1", "GROUPING claim")]),
        ]
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(find_overlap, "github_pages", return_value=[]), \
             patch.object(find_overlap.subprocess, "run", side_effect=responses) as run:
            groups = find_overlap.search(["grouping"], Path(directory))
        found = [label for _, items in groups for label, _ in items]
        self.assertEqual(found, ["https://example/r1 (updated 2026-11-01)"])
        self.assertEqual(run.call_count, 7)
        payloads = [json.loads(call.kwargs["input"]) for call in run.call_args_list]
        self.assertEqual([p["variables"].get("cursor") for p in payloads],
                         [None, "d-next", None, None, "c-next", None, "r-next"])
        self.assertEqual([p["variables"].get("id") for p in payloads],
                         [None, None, "d0", "d1", "d1", "c1", "c1"])
        self.assertTrue(all(call.args[0] == ["gh", "api", "graphql", "--input", "-"]
                            for call in run.call_args_list))

    def test_discussion_title_and_body_are_searched(self):
        discussion = record("d", "a body-only phrase", title="A title-only phrase")
        for term in ("TITLE-ONLY", "body-only"):
            with self.subTest(term=term), tempfile.TemporaryDirectory() as directory, \
                 patch.object(find_overlap, "github_pages", return_value=[]), \
                 patch.object(find_overlap.subprocess, "run", side_effect=[
                     page(("repository", "discussions"), [discussion]),
                     page(("node", "comments"), [])]):
                groups = find_overlap.search([term], Path(directory))
            self.assertEqual(sum(len(items) for _, items in groups), 1)

    def test_partial_graphql_error_is_not_a_successful_empty_search(self):
        response = page(("repository", "discussions"), [])
        data = json.loads(response.stdout)
        data["errors"] = [{"message": "Discussion access denied"}]
        response.stdout = json.dumps(data)
        with patch.object(find_overlap, "github_pages", return_value=[]), \
             patch.object(find_overlap.subprocess, "run", return_value=response), \
             patch.object(sys, "argv", ["find_overlap.py", "grouping"]), \
             contextlib.redirect_stderr(io.StringIO()) as err, \
             contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(find_overlap.main(), 1)
            self.assertIn("incomplete", err.getvalue())
            self.assertIn("Discussion access denied", err.getvalue())
            self.assertEqual(out.getvalue(), "")

    def test_reply_timeout_is_not_silently_skipped(self):
        responses = [
            page(("repository", "discussions"), [record("d", "board")]),
            page(("node", "comments"), [record("c", "claim", replies={"totalCount": 1})]),
            subprocess.TimeoutExpired(["gh", "api"], 45),
        ]
        with patch.object(find_overlap.subprocess, "run", side_effect=responses):
            with self.assertRaises(subprocess.TimeoutExpired):
                find_overlap.github_discussions()

    def test_invalid_connections_and_cursors_fail_closed(self):
        cases = [
            {"data": {"repository": None}},
            {"data": {"repository": {"discussions": {"nodes": [None], "pageInfo": {}}}}},
            {"data": {"repository": {"discussions": {"nodes": [], "pageInfo": {
                "hasNextPage": True, "endCursor": None}}}}},
        ]
        for payload in cases:
            with self.subTest(payload=payload), patch.object(find_overlap.subprocess, "run",
                    return_value=subprocess.CompletedProcess([], 0, json.dumps(payload), "")):
                with self.assertRaises(ValueError):
                    find_overlap.github_discussions()
        repeated = page(("repository", "discussions"), [], "same")
        with patch.object(find_overlap.subprocess, "run", return_value=repeated) as run:
            with self.assertRaises(ValueError):
                find_overlap.github_discussions()
        self.assertEqual(run.call_count, 2)


if __name__ == "__main__":
    unittest.main()
