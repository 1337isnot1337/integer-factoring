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


class OverlapSearchTests(unittest.TestCase):
    def test_pagination_is_flattened_and_bad_shape_is_rejected(self):
        result = subprocess.CompletedProcess([], 0, json.dumps([[{"id": 1}], [{"id": 2}]]), "")
        with patch.object(find_overlap.subprocess, "run", return_value=result) as run:
            self.assertEqual(find_overlap.github_pages("example"), [{"id": 1}, {"id": 2}])
            self.assertIn("--paginate", run.call_args.args[0])
            self.assertIn("--slurp", run.call_args.args[0])
        result.stdout = json.dumps({"message": "unexpected response"})
        with patch.object(find_overlap.subprocess, "run", return_value=result):
            with self.assertRaises(ValueError):
                find_overlap.github_pages("example")

    def test_search_reports_matches_from_each_source_and_excludes_prs_from_issues(self):
        def record(url, title="", body=""):
            return {"html_url": url, "title": title, "body": body, "updated_at": "2026-10-30"}

        pr = record("https://example/pr/1", "Grouping proof")
        issue = record("https://example/issue/2", body="Grouping question")
        pr_as_issue = {**pr, "pull_request": {}}
        comment = record("https://example/comment/3", body="GROUPING active")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "research").mkdir()
            (root / "research" / "note.md").write_text("A grouping obstacle\n")
            with patch.object(find_overlap, "github_pages", side_effect=[[pr], [issue, pr_as_issue], [comment]]):
                groups = find_overlap.search(["grouping"], root)
        self.assertEqual([len(items) for _, items in groups], [1, 1, 1, 1])
        self.assertEqual(groups[1][1][0][0].split()[0], "https://example/issue/2")
        self.assertEqual(groups[3][1][0][0], "research/note.md")

    def test_failure_is_not_reported_as_no_overlap(self):
        error = subprocess.CalledProcessError(1, ["gh", "api"])
        with patch.object(find_overlap, "search", side_effect=error), \
             patch.object(sys, "argv", ["find_overlap.py", "grouping"]), \
             contextlib.redirect_stderr(io.StringIO()) as err, \
             contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(find_overlap.main(), 1)
            self.assertIn("incomplete", err.getvalue())
            self.assertEqual(out.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
