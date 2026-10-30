"""Modular decoder comparisons with the independent polynomial coefficients."""

from itertools import product
import unittest

from tools.channel_decoder import channel_factor, piece_factor
from tools.check_crossing import original_coefficients
from tools.check_grid_level import _expanded_channel


class ChannelDecoderTests(unittest.TestCase):
    def test_selected_width_factors_in_both_orientations(self):
        # Small widths select the fallback; width six includes both choices.
        cases = ((0, (0, 0, 0, 0)), (1, (0, 1, 1, 0)),
                 (2, (1, 2, 3, 0)), (3, (1, 2, 3, 0)),
                 (6, (9, 14, 23, 16)))
        for width, values in cases:
            t = (width + 2) // 3
            blocks = tuple(tuple((v >> (3 * s)) & 7 for v in values) for s in range(t))
            coefficients = original_coefficients(blocks)
            for i, orientation, modulus in product(range(t + 1), ("row", "column"), (6, 101)):
                x, y, u, v = values
                p, q, r, s = (x, y, u, v) if orientation == "row" else (x, u, y, v)
                j = i if orientation == "row" else t - i
                bits = min(2 * width, 4 * t + (73 * j + 25) // 26)
                total = sum(piece_factor(width, j, label, p, q, modulus, left=True)
                            * piece_factor(width, j, label, r, s, modulus, left=False)
                            for label in range(1 << bits))
                self.assertEqual(total % modulus, 4 ** t * coefficients[i] % modulus)
                self.assertEqual(piece_factor(width, j, 1 << bits, p, q, modulus, left=True), 0)

    def test_fallback_large_width_matches_polynomial(self):
        # At i=t the trivial factorization is selected, with ten blocks.
        width, index, modulus = 30, 10, 45
        p, q, r, s = 0x1234567, 0x2345678, 0x3456789, 0x0123456
        label = (p << width) | q
        blocks = tuple(tuple((v >> (3 * b)) & 7 for v in (p, q, r, s)) for b in range(index))
        self.assertEqual(piece_factor(width, index, label, p, q, modulus, left=True), 1)
        self.assertEqual(piece_factor(width, index, label, r, s, modulus, left=False),
                         4 ** index * original_coefficients(blocks)[index] % modulus)

    def test_individual_channels_against_integer_reference(self):
        for width, index in ((0, 0), (3, 0), (3, 1), (6, 1), (7, 2)):
            t = (width + 2) // 3
            limit = 16 ** t * 7 ** index
            labels = (0, 1, limit // 3, limit - 1, limit, limit + 1)
            for label, first, second, left, modulus in product(labels, (0, (1 << width) - 1),
                                                                (0, (1 << width) - 1),
                                                                (False, True), (4, 6, 101)):
                expected = _expanded_channel(width, index, label, first, second, left) % modulus
                self.assertEqual(channel_factor(width, index, label, first, second, modulus, left=left), expected)

    def test_channel_sum_against_polynomial_and_transpose(self):
        for width, index, values in ((3, 1, (1, 2, 3, 0)), (6, 1, (9, 14, 23, 16))):
            t = (width + 2) // 3
            blocks = tuple(tuple((v >> (3 * s)) & 7 for v in values) for s in range(t))
            expected = 4 ** t * original_coefficients(blocks)[index]
            for orientation, modulus in product(("row", "column"), (4, 9, 101)):
                x, y, u, v = values
                p, q, r, s = (x, y, u, v) if orientation == "row" else (x, u, y, v)
                j = index if orientation == "row" else t - index
                total = sum(channel_factor(width, j, label, p, q, modulus, left=True)
                            * channel_factor(width, j, label, r, s, modulus, left=False)
                            for label in range(16 ** t * 7 ** j))
                self.assertEqual(total % modulus, expected % modulus)

    def test_validation(self):
        for args in ((-1, 0, 0, 0, 0, 2), (3, 2, 0, 0, 0, 2),
                     (3, 0, -1, 0, 0, 2), (3, 0, 0, 8, 0, 2), (3, 0, 0, 0, 0, 1)):
            with self.assertRaises(ValueError):
                channel_factor(*args, left=True)


if __name__ == "__main__":
    unittest.main()
