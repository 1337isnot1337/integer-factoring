"""Regression tests for evidence binding and failed-build invalidation."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tools.check_k6_receipt import (
    ENDPOINTS, LOG_NAMES, PROTECTED, REFERENCE_COMMIT,
    VerificationError, verify_receipt,
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class K6ReceiptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "root"
        self.build = Path(self.temporary.name) / "build.fresh"
        self.logs = Path(self.temporary.name) / "logs"
        for path in (self.root, self.build, self.logs):
            path.mkdir()
        self.proof_name = "variants/k6/FactoringK6/Helper.lean"
        sources = {
            "variants/k6/FactoringK6/Certified.lean":
                "module\n" + "\n".join(f"#print axioms {name}" for name in sorted(ENDPOINTS)) + "\n",
            self.proof_name: "theorem helper : True := trivial\n",
        }
        hashes = {}
        for name, text in sources.items():
            self.write(self.root / name, text)
            self.write(self.build / name.removeprefix("variants/k6/"), text)
            hashes[name] = digest(text.encode())
        patch_text = (
            "diff --git a/FactoringAux/Schedule.lean b/FactoringAux/Schedule.lean\n"
            "--- a/FactoringAux/Schedule.lean\n+++ b/FactoringAux/Schedule.lean\n"
            "@@ -1 +1 @@\n-def degree := 16\n+def degree := 12\n"
        )
        self.write(self.root / "variants/k6/reference.patch", patch_text)
        self.write(self.build / "FactoringAux/Schedule.lean", "def degree := 12\n")
        protected = {}
        for name in PROTECTED:
            text = f"fixed reference input: {name}\n"
            self.write(self.build / name, text)
            protected[name] = digest(text.encode())
        manifest = {
            "reference_commit": REFERENCE_COMMIT,
            "patch_sha256": digest(patch_text.encode()),
            "patched_reference_sha256": {
                "FactoringAux/Schedule.lean": digest(b"def degree := 12\n"),
            },
            "protected_reference_sha256": protected,
            "source_sha256": hashes,
        }
        self.write(self.logs / "sources.json", json.dumps(manifest))
        self.log_contents = {
            "patch-check.log": "Checking patch FactoringAux/Schedule.lean...\n",
            "patch-apply.log": "Applied patch FactoringAux/Schedule.lean cleanly.\n",
            "toolchain.log": "Lean (version 4.35.0, x86_64-unknown-linux-gnu, commit test, Release)\n",
            "cache.log": "Completed successfully!\n",
            "build.log": "Built FactoringK6.Certified\nBuild completed successfully (42 jobs).\n",
            "Certified.log": "\n".join(
                f"info: FactoringK6/Certified.lean:10:0: '{name}' depends on axioms: "
                "[propext, Classical.choice, Quot.sound]"
                for name in sorted(ENDPOINTS)
            ) + "\n",
        }
        self.log_contents["build.log"] = (
            "Built FactoringK6.Certified\n" + self.log_contents["Certified.log"]
            + "Build completed successfully (42 jobs).\n"
        )
        for name, text in self.log_contents.items():
            self.write(self.logs / name, text)

    @staticmethod
    def write(path: Path, text: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def verify(self) -> dict:
        return verify_receipt(self.root, self.build, self.logs)

    def assert_invalidated(self) -> None:
        self.assertFalse((self.logs / "axioms.json").exists())
        self.assertFalse((self.logs / "verification.json").exists())

    def test_success_binds_sources_logs_and_axiom_receipt(self) -> None:
        result = self.verify()
        saved = json.loads((self.logs / "verification.json").read_text())
        self.assertEqual(saved, result)
        self.assertEqual(result["status"], "success")
        self.assertEqual(result["source_manifest_sha256"], digest((self.logs / "sources.json").read_bytes()))
        for name in LOG_NAMES:
            self.assertEqual(result["log_sha256"][name], digest((self.logs / name).read_bytes()))
        for name, expected in result["build_source_sha256"].items():
            self.assertEqual(expected, digest((self.build / name).read_bytes()))
        for name, expected in result["root_proof_sha256"].items():
            self.assertEqual(expected, digest((self.root / name).read_bytes()))
        self.assertEqual(result["axioms_receipt_sha256"], digest((self.logs / "axioms.json").read_bytes()))
        self.assertEqual({entry["declaration"] for entry in result["checked_declarations"]}, ENDPOINTS)

    def test_changed_copied_proof_invalidates_previous_success(self) -> None:
        self.verify()
        self.write(self.build / "FactoringK6/Helper.lean", "different copied proof\n")
        with self.assertRaisesRegex(VerificationError, "Source hash mismatch"):
            self.verify()
        self.assert_invalidated()

    def test_changed_root_proof_invalidates_previous_success(self) -> None:
        self.verify()
        self.write(self.root / self.proof_name, "new root proof after copying\n")
        with self.assertRaisesRegex(VerificationError, "Source hash mismatch"):
            self.verify()
        self.assert_invalidated()

    def test_changed_reference_inputs_cannot_reuse_a_success(self) -> None:
        for name in ("FactoringAux/Schedule.lean", "FactoringAux/Definitions.lean"):
            with self.subTest(name=name):
                original = (self.build / name).read_text()
                self.verify()
                self.write(self.build / name, "changed reference input\n")
                with self.assertRaisesRegex(VerificationError, "Source hash mismatch"):
                    self.verify()
                self.assert_invalidated()
                self.write(self.build / name, original)

    def test_missing_endpoint_cannot_reuse_a_success(self) -> None:
        self.verify()
        missing = "Factoring.K6.certified_improvement"
        text = "\n".join(line for line in self.log_contents["Certified.log"].splitlines() if missing not in line)
        self.write(self.logs / "Certified.log", text)
        with self.assertRaisesRegex(VerificationError, "Endpoint axiom reports"):
            self.verify()
        self.assert_invalidated()

    def test_extra_endpoint_cannot_expand_the_verified_contract(self) -> None:
        self.verify()
        self.write(self.logs / "Certified.log", self.log_contents["Certified.log"]
                   + "'Factoring.K6.unrequested_endpoint' depends on axioms: [propext]\n")
        with self.assertRaisesRegex(VerificationError, "Endpoint axiom reports"):
            self.verify()
        self.assert_invalidated()

    def test_extra_axiom_rejects_an_otherwise_successful_build(self) -> None:
        self.verify()
        self.write(self.logs / "Certified.log", self.log_contents["Certified.log"].replace(
            "Quot.sound]", "Quot.sound, unprovedBudget]", 1
        ))
        with self.assertRaisesRegex(VerificationError, "Unexpected axioms"):
            self.verify()
        self.assert_invalidated()

    def test_stale_endpoint_report_cannot_be_mixed_with_new_build(self) -> None:
        self.verify()
        self.write(self.logs / "build.log", self.log_contents["build.log"].replace(
            "[propext, Classical.choice, Quot.sound]", "[propext, Quot.sound]", 1
        ))
        with self.assertRaisesRegex(VerificationError, "axiom reports disagree"):
            self.verify()
        self.assert_invalidated()

    def test_failed_or_incomplete_build_invalidates_success(self) -> None:
        failed_logs = (
            "Built FactoringK6.MachineBridge\n",
            "error: theorem failed\nBuild completed successfully (42 jobs).\n",
            "warning: declaration uses 'sorry'\nBuild completed successfully (42 jobs).\n",
        )
        for text in failed_logs:
            with self.subTest(text=text):
                self.write(self.logs / "build.log", self.log_contents["build.log"])
                self.verify()
                self.write(self.logs / "build.log", text)
                with self.assertRaises(VerificationError):
                    self.verify()
                self.assert_invalidated()

    def test_partial_receipt_write_does_not_leave_a_success(self) -> None:
        self.verify()
        import tools.check_k6_receipt as checker
        original = checker.atomic_write

        def fail_second(path: Path, data: bytes) -> None:
            if path.name == "verification.json":
                raise OSError("simulated second-write failure")
            original(path, data)

        with patch.object(checker, "atomic_write", side_effect=fail_second):
            with self.assertRaisesRegex(OSError, "second-write"):
                self.verify()
        self.assert_invalidated()


if __name__ == "__main__":
    unittest.main()
