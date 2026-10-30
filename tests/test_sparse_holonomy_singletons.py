"""Weighted singleton and sparse Gaussian regression controls."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import check_sparse_holonomy_singletons as holonomy


class SparseHolonomySingletonTests(unittest.TestCase):
    def test_binary_pairing_profile(self):
        weights = tuple(1 << j for j in range(8))
        self.assertEqual(holonomy.pairing_profile(weights, 255),
                         (63, 63, 207, 207, 243, 243, 252, 252, 255))

    def test_nonbinary_weights_and_exact_unit_recurrence(self):
        weights = (2, 5, 9, 14, 20)
        matrix = holonomy.singleton_matrix(weights, 50)
        orders, _ = holonomy.smith.smith_profile(matrix, 50)
        self.assertEqual(orders, holonomy.pairing_profile(weights, 50))
        self.assertEqual(holonomy.unit_recurrence(weights), 2)

    def test_filtered_gaussian_encoding(self):
        self.assertEqual(holonomy.gaussian_check(4, 3), 225)


if __name__ == "__main__":
    unittest.main()
