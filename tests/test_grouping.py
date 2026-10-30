"""Independent finite checks of the exact grouping certificates."""

from itertools import combinations, product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.grouping import (  # noqa: E402
    NoGroupingFound,
    SearchLimitExceeded,
    WeightedGraph,
    _bounded_power,
    demonstration,
    find_grouping,
    sufficient_node_count,
    verify_grouping,
)


def direct_spec(n, edges, k, d, colors, p):
    """Separate specification, directly using the ordered-edge formula."""

    u = 2 * sum(weight for _, _, weight in edges)
    for a in range(k):
        if k * sum(color == a for color in colors) > 2 * n:
            return False
        for b in range(a + 1, k):
            y = sum(
                weight
                for v, w, weight in edges
                if (colors[v] == a and colors[w] == b) or (colors[v] == b and colors[w] == a)
            )
            if d * k * k * y > d * u + n * p:
                return False
    return True


class GroupingTests(unittest.TestCase):
    def test_normalization_and_ordered_total(self):
        graph = WeightedGraph(4, ((0, 1, 2), (1, 0, 3), (2, 3, 4), (1, 2, 0)))
        self.assertEqual(graph.edges, ((0, 1, 5), (2, 3, 4)))
        self.assertEqual(graph.degrees, (5, 5, 4, 4))
        self.assertEqual(graph.ordered_total, 18)
        self.assertEqual(graph.max_degree, 5)

    def test_all_small_weighted_graphs_against_independent_spec(self):
        # Every loopless labeled graph on <=4 vertices with weights in {0,1,2}.
        # The checker and the returned first solution are compared to direct
        # exhaustive inequalities, not merely to the implementation's result.
        graph_count = 0
        coloring_count = 0
        for n in range(5):
            endpoints = tuple(combinations(range(n), 2))
            for weights in product(range(3), repeat=len(endpoints)):
                graph = WeightedGraph(n, tuple((u, v, w) for (u, v), w in zip(endpoints, weights)))
                graph_count += 1
                for k, d in ((2, 2), (3, 2), (3, 4)):
                    first = None
                    for colors in product(range(k), repeat=n):
                        expected = direct_spec(n, graph.edges, k, d, colors, graph.max_degree)
                        actual = verify_grouping(graph, k, d, colors)
                        self.assertEqual(actual.valid, expected, (graph, k, d, colors))
                        if expected and first is None:
                            first = colors
                        coloring_count += 1
                    if first is None:
                        with self.assertRaises(NoGroupingFound):
                            find_grouping(graph, k, d)
                    else:
                        found = find_grouping(graph, k, d)
                        self.assertEqual(found.certificate.colors, first)
                        self.assertTrue(found.certificate.valid)
        self.assertEqual(graph_count, 761)
        self.assertEqual(coloring_count, 131513)

    def test_matching_star_and_zero_examples(self):
        for graph in (
            WeightedGraph(6, ((0, 1, 9), (2, 3, 9), (4, 5, 9))),
            WeightedGraph(5, tuple((0, v, 7) for v in range(1, 5))),
            WeightedGraph(5, ()),
        ):
            result = find_grouping(graph, 3, 2)
            self.assertTrue(direct_spec(graph.vertex_count, graph.edges, 3, 2, result.certificate.colors, graph.max_degree))
            self.assertFalse(result.certificate.sufficient_threshold_met)

    def test_nonexistence_is_distinct_from_a_search_cap(self):
        with self.assertRaisesRegex(NoGroupingFound, "9 colorings"):
            find_grouping(WeightedGraph(2, ((0, 1, 1),)), 3, 2)
        with self.assertRaisesRegex(SearchLimitExceeded, "no colorings were enumerated"):
            find_grouping(WeightedGraph(30, ()), 3, 2, max_candidates=100)

    def test_empty_graph(self):
        result = find_grouping(WeightedGraph(0, ()), 3, 2)
        self.assertEqual(result.certificate.colors, ())
        self.assertEqual(result.total_colorings, 1)
        self.assertEqual(result.candidates_examined, 1)
        huge_groups = find_grouping(WeightedGraph(0, ()), 10 ** 100, 2)
        self.assertEqual(huge_groups.certificate.colors, ())
        self.assertTrue(huge_groups.certificate.valid)

    def test_invalid_graphs(self):
        bad_cases = (
            (-1, ()), (True, ()), (2, ((0, 0, 1),)),
            (2, ((0, 0, 0),)), (2, ((0, 2, 1),)),
            (2, ((-1, 1, 1),)), (2, ((0, 1, -1),)),
            (2, ((0, 1, 1.5),)), (2, ((0, 1, True),)),
            (2, ((0, 1),)), (2, ((False, 1, 1),)),
        )
        for n, edges in bad_cases:
            with self.subTest(n=n, edges=edges), self.assertRaises(ValueError):
                WeightedGraph(n, edges)

    def test_invalid_parameters_and_colorings(self):
        graph = WeightedGraph(2, ((0, 1, 3),))
        for k, d in ((1, 2), (2, 1), (True, 2), (2, 2.5)):
            with self.subTest(k=k, d=d), self.assertRaises(ValueError):
                find_grouping(graph, k, d)
        for colors in ((), (0,), (0, 0, 0), (0, 2), (False, 0), (0, -1)):
            with self.subTest(colors=colors), self.assertRaises(ValueError):
                verify_grouping(graph, 2, 2, colors)
        with self.assertRaises(ValueError):
            verify_grouping(graph, 2, 2, (0, 1), degree_bound=2)
        with self.assertRaises(ValueError):
            find_grouping(graph, 2, 2, max_candidates=0)

    def test_certificate_detects_both_failure_types(self):
        graph = WeightedGraph(2, ((0, 1, 1),))
        same = verify_grouping(graph, 3, 2, (0, 0))
        split = verify_grouping(graph, 3, 2, (0, 1))
        self.assertFalse(same.valid)
        self.assertFalse(split.valid)
        self.assertIn("group 0", same.violations[0])
        self.assertIn("pair (0, 1)", split.violations[0])

    def test_safe_power_and_threshold(self):
        for base in range(2, 8):
            for exponent in range(10):
                for cap in (1, 10, 100, 1000):
                    self.assertEqual(_bounded_power(base, exponent, cap), min(base ** exponent, cap + 1))
        self.assertEqual(_bounded_power(2, 10 ** 100, 100), 101)
        self.assertEqual(sufficient_node_count(3, 2), 864)
        self.assertEqual(sufficient_node_count(2, 2), 384)

    def test_demo_contains_recomputed_certificates(self):
        report = demonstration()
        for name, example in report["examples"].items():
            if name != "nonexistence":
                self.assertTrue(example["result"]["certificate"]["valid"])
        self.assertIn("none of the 9", report["examples"]["nonexistence"]["result"])


if __name__ == "__main__":
    unittest.main()
