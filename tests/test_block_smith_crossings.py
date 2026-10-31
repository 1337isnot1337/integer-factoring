"""Independent finite controls for the block-profile implementation."""

from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import check_block_smith_crossings as blocks


class BlockSmithTests(unittest.TestCase):
    def test_visible_pivots_do_not_certify_missing_tail_as_zero(self):
        matrix = [[1, 0], [0, 1 << 7]]
        self.assertIsNone(blocks.smith_at_precision(matrix, 7))
        orders, precision = blocks.smith_profile(matrix, 7)
        self.assertEqual(orders, (0, 7))
        self.assertGreater(precision, 7)

    def test_elementary_blocks_keep_the_exact_product_certificate(self):
        report = blocks.check_family("auxiliary elementary", *blocks.local_tools.elementary(2), maximum_block=2)
        self.assertEqual([Fraction(row["B"]) for row in report["records"]],
                         [Fraction(256, 27), Fraction(256, 27) ** 2])
        self.assertGreater(report["independent_toeplitz_ranks"], 0)


if __name__ == "__main__":
    unittest.main()
