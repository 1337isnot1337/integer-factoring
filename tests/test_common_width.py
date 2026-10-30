"""Exact rounding and composed-volume certificates.

The corresponding universal arithmetic statements are checked in
formal/CommonWidth.lean; these tests also check the executable schedule.
"""

from fractions import Fraction
import unittest

from tools.parameters import Level, schedule, width_schedule, width_step


class CommonWidthTests(unittest.TestCase):
    def test_common_width_is_exact_maximum_allowed_by_pair_bound(self):
        # Synthetic small levels exercise both rounding operations, without
        # claiming the asymptotic grouping threshold holds at these sizes.
        remainders = set()
        for k in range(1, 8):
            for d in range(1, 5):
                for g in range(1, 8):
                    level = Level(0, k, d, 1, g)
                    for width in range(7):
                        step = width_step(level, width)
                        total = Fraction(281, 312) * g ** 4 * width + g ** 3 * width + 3 * g ** 4
                        self.assertEqual(step.ordered_total_bound, total.numerator // total.denominator)
                        numerator = d * step.ordered_total_bound + g ** 4 * width
                        denominator = d * k ** 2
                        remainders.add(numerator % denominator != 0)
                        self.assertLessEqual(denominator * step.output_width, numerator)
                        self.assertGreater(denominator * (step.output_width + 1), numerator)
                        bound = (Fraction(281, 312) + Fraction(1, g) + Fraction(1, d)) * step.input_volume + 3 * g ** 4
                        self.assertLessEqual(step.output_volume, bound)
        self.assertEqual(remainders, {False, True})

    def test_schedule_composes_in_evaluation_order(self):
        levels = schedule(4)
        # Thousands of bits make accidental floating-point conversion fail.
        top_width = 2 ** 4096 + 117
        steps = width_schedule(levels, top_width)
        self.assertEqual([step.level.d for step in steps], [3, 2, 1, 0])
        bound = Fraction(steps[0].input_volume)
        for index, step in enumerate(steps):
            if index:
                previous = steps[index - 1]
                self.assertEqual(step.input_width, previous.output_width)
                self.assertEqual(step.input_volume, previous.output_volume)
            g, d = step.level.grid_side, step.level.slack_denominator
            bound = (Fraction(281, 312) + Fraction(1, g) + Fraction(1, d)) * bound + 3 * g ** 4
            self.assertLessEqual(step.output_volume, bound)
        self.assertEqual(width_schedule([], 0), [])
        self.assertEqual(width_schedule(levels, 0)[0].input_width, 0)

    def test_invalid_inputs(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                width_step(schedule(1)[0], value)
            with self.assertRaises(ValueError):
                width_schedule([], value)
        for level in (Level(0, 0, 1, 1, 1), Level(0, 1, 0, 1, 1), Level(0, 1, 1, 1, 0)):
            with self.assertRaises(ValueError):
                width_step(level, 1)


if __name__ == "__main__":
    unittest.main()
