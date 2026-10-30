"""Finite checks of channel capacity and the common integer bound."""

import unittest

from tools.check_crossing import exact_rank_formula
from tools.exact_widths import common_total_bound, piece_widths, seven_width


class ExactWidthsTests(unittest.TestCase):
    def test_integer_log_boundaries(self):
        for exponent in range(201):
            width = seven_width(exponent)
            self.assertLessEqual(7 ** exponent, 2 ** width)
            if width:
                self.assertGreater(7 ** exponent, 2 ** (width - 1))

    def test_piece_widths_cover_channels_and_common_bound(self):
        for label_width in range(61):
            blocks = (label_width + 2) // 3
            total_bound = common_total_bound(5, label_width)
            for i in range(blocks + 1):
                row, col = piece_widths(label_width, i)
                for width, j in ((row, i), (col, blocks - i)):
                    self.assertTrue(width == 2 * label_width or
                                    2 ** width >= exact_rank_formula(blocks, j))
                self.assertLessEqual(4 * (total_bound - 5 ** 3 * label_width),
                                     5 ** 4 * (8 * blocks + seven_width(blocks) + 1))
                self.assertLessEqual(5 ** 4 * (row + col),
                                     4 * (total_bound - 5 ** 3 * label_width) + 3)

    def test_invalid_inputs(self):
        for call in (lambda: seven_width(-1), lambda: piece_widths(-1, 0),
                     lambda: piece_widths(3, 2), lambda: common_total_bound(0, 1),
                     lambda: common_total_bound(1, -1)):
            with self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main()
