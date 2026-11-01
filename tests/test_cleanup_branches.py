"""Exercise cleanup against a real local remote, without GitHub writes."""

import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

PATH = Path(__file__).resolve().parents[1] / "tools" / "cleanup_branches.py"
SPEC = importlib.util.spec_from_file_location("cleanup_branches", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.repo = self.root / "work"
        self.remote = self.root / "remote.git"
        self.git("init", "-b", "main", str(self.repo))
        self.git("-C", str(self.repo), "config", "user.name", "Test")
        self.git("-C", str(self.repo), "config", "user.email", "test@example.com")
        self.git("init", "--bare", str(self.remote))
        self.git("-C", str(self.repo), "remote", "add", "origin", str(self.remote))
        old_env = dict(os.environ, GIT_AUTHOR_DATE="2020-01-01T00:00:00Z",
                       GIT_COMMITTER_DATE="2020-01-01T00:00:00Z")
        self.git("-C", str(self.repo), "commit", "--allow-empty", "-m", "base", env=old_env)
        self.base = self.git("-C", str(self.repo), "rev-parse", "HEAD")
        for name in ("stale", "stale-two", "protected", "open-head", "open-base"):
            self.git("-C", str(self.repo), "branch", name)
        self.git("-C", str(self.repo), "commit", "--allow-empty", "-m", "current")
        self.git("-C", str(self.repo), "branch", "recent")
        self.git("-C", str(self.repo), "checkout", "-b", "unique", self.base)
        self.git("-C", str(self.repo), "commit", "--allow-empty", "-m", "unique", env=old_env)
        self.unique = self.git("-C", str(self.repo), "rev-parse", "HEAD")
        self.git("-C", str(self.repo), "checkout", "main")
        self.git("-C", str(self.repo), "push", "origin", "--all")
        self.addCleanup(patch.stopall)
        patch.object(MODULE, "ROOT", self.repo).start()
        patch.object(MODULE, "api", side_effect=self.api).start()

    def git(self, *args, **kwargs):
        return subprocess.check_output(
            ["git", *args], text=True, stderr=subprocess.STDOUT, **kwargs,
        ).strip()

    def refs(self):
        return dict(line.split() for line in self.git(
            "--git-dir", str(self.remote), "for-each-ref",
            "--format=%(refname:short) %(objectname)", "refs/heads",
        ).splitlines())

    def api(self, repository, endpoint, paginate=False):
        self.assertEqual(repository, "owner/repo")
        if not endpoint:
            return {"default_branch": "main"}
        if endpoint.startswith("/pulls"):
            return [{"head": {"ref": "open-head"}, "base": {"ref": "open-base"}}]
        self.assertTrue(endpoint.startswith("/branches"))
        return [{"name": name, "commit": {"sha": sha}, "protected": name == "protected"}
                for name, sha in self.refs().items()]

    def test_preview_never_deletes(self):
        before = self.refs()
        MODULE.cleanup("owner/repo")
        self.assertEqual(self.refs(), before)

    def test_apply_preserves_open_protected_recent_and_unique_work(self):
        MODULE.cleanup("owner/repo", apply=True)
        self.assertEqual(set(self.refs()), {
            "main", "protected", "open-head", "open-base", "recent", "unique",
        })
        self.assertEqual(self.git("-C", str(self.repo), "rev-parse", "stale"), self.base)

    def test_new_open_pr_aborts_batch(self):
        before = self.refs()
        original = MODULE.open_pr_branches
        with patch.object(MODULE, "open_pr_branches", side_effect=[
            original("owner/repo"), {"stale"},
        ]):
            with self.assertRaisesRegex(RuntimeError, "in use"):
                MODULE.cleanup("owner/repo", apply=True)
        self.assertEqual(self.refs(), before)

    def test_concurrent_push_rejects_entire_atomic_deletion(self):
        original = MODULE.run

        def race(*command):
            if command[:2] == ("git", "push"):
                self.git("--git-dir", str(self.remote), "update-ref", "refs/heads/stale", self.unique)
            return original(*command)

        with patch.object(MODULE, "run", side_effect=race):
            with self.assertRaises(subprocess.CalledProcessError):
                MODULE.cleanup("owner/repo", apply=True)
        self.assertEqual(self.refs()["stale"], self.unique)
        self.assertEqual(self.refs()["stale-two"], self.base)


if __name__ == "__main__":
    unittest.main()
