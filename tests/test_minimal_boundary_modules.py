"""Exact controls for minimal-representation boundary rate certificates."""

from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import check_minimal_boundary_modules as minimal


class MinimalBoundaryModuleTests(unittest.TestCase):
    def test_observable_quotient_removes_fixed_precision_rate_failure(self):
        result = minimal.hidden_component_control()
        self.assertEqual(result["minimal_dimension"], 1)
        self.assertTrue(all(row["minimal_module_rate"] == "1/2" for row in result["records"]))

    def test_minimal_boundary_rates_need_a_limsup(self):
        result = minimal.parity_control()
        self.assertEqual(result["boundary_partition"], [str(int(length % 2 == 0)) for length in range(9)])

    def test_physical_choices_are_grouped_not_counted(self):
        self.assertEqual(minimal.grouped_physical_control()["rate"], "1")

    def test_rational_reduction_and_finite_window(self):
        result = minimal.embedded_controls(2)
        self.assertEqual(result["finite_window_inequalities"], 12)
        self.assertEqual(result["literal_scalar_word_comparisons"], 170)


if __name__ == "__main__":
    unittest.main()
