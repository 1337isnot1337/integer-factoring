"""Finite boundary checks for the exact schedule.

The infinite growth and contraction estimates are proved in recurrence.md;
these tests catch integer-rounding and cutoff mistakes in the utility.
"""

from fractions import Fraction
import unittest

from tools.parameters import admissible_depth, schedule


class ParameterTests(unittest.TestCase):
    def test_threshold_boundaries(self) -> None:
        self.assertEqual(admissible_depth(0), 0)
        levels = schedule(4)
        for depth in range(1, len(levels)):
            threshold = levels[depth].k ** 40
            self.assertEqual(admissible_depth(threshold - 1), depth - 1)
            self.assertEqual(admissible_depth(threshold), depth)
            self.assertEqual(admissible_depth(threshold + 1), depth)

    def test_coarse_bounds_and_summable_error_certificate(self) -> None:
        rho = Fraction(281, 312)
        product_bound = Fraction(1)
        levels = schedule(6)
        for level in levels:
            k, g, d = level.k, level.grid_side, level.d
            self.assertGreaterEqual(k, 48 * 2 ** d)
            self.assertGreaterEqual(k, d + 48)
            self.assertLessEqual(level.slack_denominator, k ** 2)
            self.assertGreaterEqual(g, k ** 2)
            self.assertLessEqual(g, k ** 8)
            self.assertLessEqual(2 ** (level.logarithm_bound - 1), k)
            self.assertLess(k, 2 ** level.logarithm_bound)
            product_bound *= (rho + Fraction(1, g) + Fraction(1, level.slack_denominator)) / rho
            self.assertLess(product_bound, Fraction(11, 10))
        total_error_bound = Fraction(1, 47) + Fraction(1, 1728)
        self.assertLess(total_error_bound / rho, Fraction(1, 11))

    def test_input_validation(self) -> None:
        with self.assertRaises(ValueError):
            admissible_depth(-1)
        with self.assertRaises(ValueError):
            schedule(-1)


if __name__ == "__main__":
    unittest.main()
