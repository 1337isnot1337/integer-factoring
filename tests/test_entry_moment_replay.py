"""Tests for cross-platform entry-moment replay and python -O tamper checks."""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import json
import lzma
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "tools") not in sys.path:
    sys.path.insert(0, str(ROOT / "tools"))

import check_entry_moment_cycles as entry_cycles  # noqa: E402


class EntryMomentReplayTests(unittest.TestCase):
    def test_saved_flags_replay_passes_and_preserves_frozen_receipt(self) -> None:
        receipt = ROOT / "evidence" / "entry-moment-cycles.json"
        before = receipt.read_bytes()
        entry_cycles.replay(receipt, saved_flags=True)
        self.assertEqual(receipt.read_bytes(), before)

    def test_required_constant_matches_high_precision_decimal_root(self) -> None:
        sources = json.loads(
            (ROOT / "evidence" / "clustered-power-relations" / "receipt.json").read_text()
        )["certificates"]
        source = next(item for item in sources if item["N"] == 1000036000099)
        flags = lzma.decompress((ROOT / source["flags_file"]).read_bytes())
        got = entry_cycles.required_constant(flags, window=90480, order=36)
        period, good = len(flags), sum(flags)
        counts = entry_cycles.cyclic_window_counts(flags, 90480)
        lower = sum((-(period * c - good * 90480)) ** 72 for c in counts if period * c < good * 90480)
        with localcontext() as ctx:
            ctx.prec = 120
            unit_rhs = Decimal(36 * 90480 * good) ** 36 * Decimal(period) ** 37
            ref = float(((Decimal(lower) / unit_rhs).ln() / Decimal(36)).exp())
        self.assertEqual(got, ref)

    def test_replay_rejects_exact_and_decimal_corruptions_under_python_optimize(self) -> None:
        receipt_data = json.loads((ROOT / "evidence" / "entry-moment-cycles.json").read_text())

        bad_int = copy.deepcopy(receipt_data)
        bad_int["cycles"][0]["sufficient_windows"][0]["lower_tail"]["lower_tail_numerator"] += 1

        bad_dec = copy.deepcopy(receipt_data)
        bad_dec["cycles"][1]["matched_mean_windows"][1]["required_constant_decimal"] += 1e-8

        with tempfile.TemporaryDirectory() as tmp:
            for name, payload in (("bad_int.json", bad_int), ("bad_dec.json", bad_dec)):
                target = Path(tmp) / name
                target.write_text(json.dumps(payload))
                proc = subprocess.run(
                    [
                        sys.executable,
                        "-O",
                        "-B",
                        str(ROOT / "tools" / "check_entry_moment_cycles.py"),
                        "--replay",
                        str(target),
                        "--saved-flags",
                    ],
                    cwd=str(ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertNotEqual(proc.returncode, 0, msg=name)
                self.assertIn("ValueError", proc.stderr)


if __name__ == "__main__":
    unittest.main()
