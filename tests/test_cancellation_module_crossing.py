"""Exact certificate check for a complete coupled-auxiliary crossing."""

from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import check_cancellation_module_crossing as example


class CoupledCrossingTests(unittest.TestCase):
    def test_strict_certificate_gap_and_literal_coefficients(self):
        report = example.check(maximum_word=2)
        self.assertEqual(report["B_module"], "1625/108")
        self.assertEqual(report["B_min_plus_lower"], "425/27")
        self.assertLess(Fraction(report["B_module"]), Fraction(report["B_min_plus_lower"]))
        self.assertGreater(report["literal_coefficient_rank_inequalities"], 0)


if __name__ == "__main__":
    unittest.main()
