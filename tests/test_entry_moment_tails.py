"""Asymmetric regressions distinguish the lower tail from the upper tail."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from check_cofactor_orbit_moments import moment_record
from check_entry_moment_cycles import lower_tail_record, required_constant


class EntryMomentTailTests(unittest.TestCase):
    def test_asymmetric_three_cycle(self):
        # Mean 1/3: the two misses contribute 1 each, while the hit's
        # upper-tail contribution would be 4. Denominator is 3^3.
        flags = bytes([1, 0, 0])
        record = lower_tail_record(flags, 1, 1, 1)
        self.assertEqual(record["lower_tail_numerator"], 2)
        self.assertAlmostEqual(required_constant(flags, 1, 1), 2 / 9)

    def test_independent_lower_tail_and_trivial_bound(self):
        flags = bytes([1, 0, 0, 1, 0])
        for order, window in ((1, 1), (3, 4), (3, 7)):
            with self.subTest(order=order, window=window):
                record = lower_tail_record(flags, window, order, 1)
                other = moment_record(flags, order, window)
                self.assertEqual(record["lower_tail_numerator"], other["lower_tail_scaled_sum"])
                self.assertLessEqual(record["lower_tail_numerator"],
                                     len(flags) * (sum(flags) * window) ** (2 * order))


if __name__ == "__main__":
    unittest.main()
