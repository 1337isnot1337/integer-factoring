"""Finite tests joining crossing factors to a complete four-node grid."""

from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.check_grid_level import (
    channel_value,
    compressed_channel_value,
    make_channels,
    run_suite,
    scaled_piece,
    signed_tables,
    verify_grid_level,
)


class GridLevelTests(unittest.TestCase):
    def test_complete_signed_and_modular_suite(self):
        report = run_suite()
        self.assertEqual(report["case_count"], 8)
        for case in report["cases"]:
            certificate = case["certificate"]
            self.assertTrue(certificate["valid"])
            self.assertEqual(certificate["multiplier"], 8)
            self.assertEqual(sum(certificate["direct_pieces"]), 8 * certificate["original_value"])
            self.assertTrue(all(certificate["fallback_row"]))
            self.assertTrue(all(certificate["fallback_column"]))

    def test_active_compressed_channels_without_large_grid(self):
        # At lambda=1,2 the cap 2*lambda selects the fallback. Lambda=3
        # exercises an actual compressed rank-8/rank-56-derived channel.
        for orientation, index in (("row", 0), ("column", 1)):
            channels = make_channels(3, index, orientation)
            self.assertFalse(channels.fallback)
            self.assertEqual(channels.channel_bits, 4)
            for x, y, u, v in product(range(8), repeat=4):
                self.assertEqual(channel_value(channels, x, y, u, v), scaled_piece(3, index, x, y, u, v))

    def test_zero_width(self):
        certificate = verify_grid_level(((2,), (-3,), (4,), (5,)), 0)
        self.assertTrue(certificate.valid)
        self.assertEqual(certificate.original_value, -120)
        self.assertEqual(certificate.multiplier, 1)

    def test_active_rank56_mixed_piece_without_grid_materialization(self):
        # lambda=6,t=2,i=1 has compressed width11<12 and uses both
        # cheap rank8 and residual rank56 channels. No 64^8 grid is built.
        for x, y, u in ((0, 0, 0), (1, 2, 3), (9, 14, 23), (63, 17, 4), (21, 42, 0)):
            v = x ^ y ^ u
            expected = scaled_piece(6, 1, x, y, u, v)
            for orientation in ("row", "column"):
                self.assertEqual(compressed_channel_value(6, 1, orientation, x, y, u, v), expected)

    def test_state_cap_and_invalid_inputs(self):
        with self.assertRaisesRegex(ValueError, "cap"):
            verify_grid_level((), 3)
        with self.assertRaisesRegex(ValueError, "cap"):
            verify_grid_level((), 10 ** 100)
        with self.assertRaises(ValueError):
            verify_grid_level((), 1)
        with self.assertRaises(ValueError):
            verify_grid_level(signed_tables(1, 0), 1, moduli=(1,))
        with self.assertRaises(ValueError):
            make_channels(1, 2, "row")
        with self.assertRaises(ValueError):
            make_channels(1, 0, "diagonal")


if __name__ == "__main__":
    unittest.main()
