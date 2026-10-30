import random
import unittest
from itertools import product

from tools.dense_grouping import (
    additive_balance_certificate,
    bounded_edge_group_size,
    bounded_edge_union_budget,
    dense_group_size,
    dense_schedule,
    first_good_grouping,
    good_grouping,
    grouping_statistics,
    has_dense_widths,
    strong_balance_certificate,
    union_budget_certificate,
)


class DenseGroupingTests(unittest.TestCase):
    def test_exact_union_budgets_at_small_and_large_parameters(self):
        values = (2, 3, 7, 48, 127, 1024, 2**127 + 1, 2**1024 - 1)
        for groups, accuracy in product(values, repeat=2):
            certificate = union_budget_certificate(groups, accuracy)
            self.assertTrue(certificate["tail_exponent_ge_4L"])
            self.assertTrue(certificate["events_lt_two_pow_3L_minus_5"])
        self.assertEqual(dense_group_size(2, 2), 55296)

    def test_bounded_edge_union_budget(self):
        values = (2, 3, 7, 48, 127, 1024, 2**127 + 1, 2**1024 - 1)
        for groups, accuracy in product(values, repeat=2):
            certificate = bounded_edge_union_budget(groups, accuracy)
            self.assertTrue(certificate["tail_exponent_ge_4L"])
            self.assertTrue(certificate["events_lt_two_pow_3L_minus_7"])
            self.assertEqual(bounded_edge_group_size(groups, accuracy) * 4, dense_group_size(groups, accuracy))

    def test_degree_balance_bridge_at_density_boundary(self):
        # Deliberately nonconstant degrees, with all weights in {4,5,6}.
        # Every group contains four vertices of each degree type.
        types = [2, 3] * 12
        weights = [
            [0 if u == v else types[u] + types[v] for v in range(24)]
            for u in range(24)
        ]
        coloring = [u // 8 for u in range(24)]
        self.assertTrue(has_dense_widths(weights, 6))
        self.assertTrue(strong_balance_certificate(weights, coloring, 3, 2))
        self.assertTrue(additive_balance_certificate(weights, coloring, 3, 2, 6))
        self.assertTrue(good_grouping(weights, coloring, 3, 2))
        # Independent unordered-edge calculation catches ordered-total and
        # double-counting errors in the checker.
        sizes, pairs, vertex_sums, degree_sums, total = grouping_statistics(
            weights, coloring, 3
        )
        self.assertEqual(sizes, [8, 8, 8])
        self.assertEqual(total, 2 * sum(weights[u][v] for u in range(24) for v in range(u + 1, 24)))
        for a, b in ((0, 1), (0, 2), (1, 2)):
            direct = sum(
                weights[u][v]
                for u in range(24)
                for v in range(u + 1, 24)
                if {coloring[u], coloring[v]} == {a, b}
            )
            self.assertEqual(pairs[a][b], direct)
            self.assertEqual(pairs[a][b], sum(vertex_sums[u][b] for u in range(24) if coloring[u] == a))
        self.assertEqual(degree_sums, [920] * 3)
        self.assertFalse(strong_balance_certificate(weights, [0] * 24, 3, 2))
        self.assertFalse(good_grouping(weights, [0] * 24, 3, 2))

    def test_finite_search_checked_by_unordered_edge_verifier(self):
        # All dense boundary graphs on four vertices; this finite test is
        # not an existence proof at the theorem's parameter threshold.
        edges = [(u, v) for u in range(4) for v in range(u + 1, 4)]
        for values in product((2, 3), repeat=len(edges)):
            weights = [[0] * 4 for _ in range(4)]
            for (u, v), value in zip(edges, values):
                weights[u][v] = weights[v][u] = value
            coloring = first_good_grouping(weights, 2, 2)
            self.assertIsNotNone(coloring)
            degrees = [sum(row) for row in weights]
            total = 2 * sum(values)
            pair_weight = sum(value for (u, v), value in zip(edges, values) if coloring[u] != coloring[v])
            self.assertLessEqual(8 * pair_weight, 2 * total + 4 * max(degrees))

    def test_density_is_required_and_zero_width_is_separate(self):
        path = [[int(abs(u - v) == 1) for v in range(6)] for u in range(6)]
        self.assertFalse(has_dense_widths(path, 1))
        zero = [[0] * 12 for _ in range(12)]
        with self.assertRaises(ValueError):
            has_dense_widths(zero, 0)
        self.assertTrue(good_grouping(zero, [u % 3 for u in range(12)], 3, 2))

    def test_bounded_edge_form_does_not_refute_path_obstruction(self):
        path = [[int(abs(u - v) == 1) for v in range(5)] for u in range(5)]
        coloring = (0, 0, 1, 1, 2)
        self.assertFalse(has_dense_widths(path, 1))
        self.assertIsNone(first_good_grouping(path, 4, 2))
        self.assertFalse(good_grouping(path, coloring, 4, 2, degree_bound=2))
        # The bounded-edge theorem uses the larger budget P=g*B=5,
        # so it has a different conclusion from the path's P=2 problem.
        self.assertTrue(good_grouping(path, coloring, 4, 2, degree_bound=5))

    def test_additive_balance_on_sparse_graphs(self):
        # Four disjoint copies with equal color counts. This case satisfies
        # the additive certificates while violating every positive density
        # hypothesis, so the two theorem interfaces are distinguished.
        size = 24
        weights = [[int(u != v and u // 6 == v // 6) for v in range(size)] for u in range(size)]
        coloring = [u % 3 for u in range(size)]
        self.assertFalse(has_dense_widths(weights, 1))
        self.assertTrue(additive_balance_certificate(weights, coloring, 3, 2, 1))
        self.assertTrue(good_grouping(weights, coloring, 3, 2, degree_bound=size))

    def test_actual_piece_widths_force_dense_grid_weights(self):
        rng = random.Random(92813)
        for width in (1, 2, 3, 4, 8, 31, 512):
            blocks = (width + 2) // 3
            for mode in ("original", "exact", "rational"):
                rows = 5
                weights = [[0] * rows for _ in range(rows)]
                for u in range(rows):
                    for v in range(u + 1, rows):
                        value = rows * width
                        for _ in range(rows * (rows - 1) // 2):
                            index = rng.randrange(blocks + 1)
                            if mode == "original":
                                extra = 3 * index
                            elif mode == "exact":
                                extra = (7**index - 1).bit_length()
                            else:
                                extra = (73 * index + 25) // 26
                            value += min(2 * width, 4 * blocks + extra)
                        weights[u][v] = weights[v][u] = value
                self.assertTrue(has_dense_widths(weights, rows**2 * width))

    def test_schedule_coarse_bounds_and_growth(self):
        schedule = dense_schedule(12)
        for d, (groups, accuracy, size) in enumerate(schedule):
            length = groups.bit_length() + accuracy.bit_length() + 8
            self.assertGreaterEqual(groups, 48 * 2**d)
            self.assertLessEqual(accuracy, groups**2)
            self.assertLessEqual(length, groups)
            self.assertLessEqual(groups, size)
            self.assertLessEqual(size, groups**8)
            if d + 1 < len(schedule):
                self.assertEqual(schedule[d + 1][0], size**2)


if __name__ == "__main__":
    unittest.main()
