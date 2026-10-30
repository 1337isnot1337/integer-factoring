"""Check the obstruction against exhaustive search and actual group weights."""

import unittest

from tools.check_grouping_obstruction import path_graph, path_threshold, run_suite
from tools.grouping import verify_grouping


class GroupingObstructionTests(unittest.TestCase):
    def test_exhaustive_nonexistence_and_boundary_witnesses(self):
        report = run_suite()
        self.assertEqual([case["exhausted_colorings"] for case in report["cases"]],
                         [27, 81, 1024, 65536])
        for case in report["cases"]:
            self.assertGreater(case["strict_obstruction_margin"], 0)
            self.assertTrue(case["boundary_witness"]["valid"])

    def test_exact_rounding_and_contiguous_construction(self):
        # Include integral equality (K=4,D=2,g=6): the obstruction is strict,
        # so it must not exclude this actual valid grouping.
        self.assertEqual(path_threshold(4, 2), 6)
        for k in range(3, 20):
            for d in (2, 3, 10, 100):
                first = path_threshold(k, d)
                graph = path_graph(first - 1)
                self.assertGreater(d * k * k,
                                   d * graph.ordered_total + graph.vertex_count * graph.max_degree)
                for g in (first, first + 1, first + k):
                    block = (g + k - 1) // k
                    certificate = verify_grouping(path_graph(g), k, d,
                                                  tuple(u // block for u in range(g)))
                    self.assertTrue(certificate.valid, (k, d, g, certificate.violations))

    def test_invalid_inputs(self):
        for value in (True, 2, -1, 3.0):
            with self.assertRaises(ValueError):
                path_graph(value)
        for k, d in ((2, 2), (3, 1), (True, 2), (3, 2.0)):
            with self.assertRaises(ValueError):
                path_threshold(k, d)


if __name__ == "__main__":
    unittest.main()
