"""Check grouped, padded clique tables against direct network contractions."""

from dataclasses import replace
from itertools import product
import random
import unittest

from tools.check_grid_level import make_channels, network_contraction, original_contraction, signed_tables
from tools.coarsened_network import Edge, Network, clique_value, network_value
from tools.grouping import WeightedGraph, verify_grouping


def signed_network():
    rng = random.Random(837)
    states = ((0, 1),) * 4
    weights = tuple((lambda z, a=a: a[z]) for a in
                    [tuple(rng.randint(-3, 3) for _ in range(2)) for _ in states])
    edges = []
    for first, second, bits in ((0, 1, 1), (0, 2, 1), (0, 3, 0),
                                (1, 2, 1), (1, 3, 1), (2, 3, 1)):
        left = tuple(tuple(rng.randint(-2, 2) for _ in range(1 << bits)) for _ in states[first])
        right = tuple(tuple(rng.randint(-2, 2) for _ in range(1 << bits)) for _ in states[second])
        edges.append(Edge(first, second, bits,
                          lambda z, x, a=left: a[z][x],
                          lambda z, x, a=right: a[z][x]))
    return Network(states, weights, tuple(edges))


def grid_row_network(tables, piece):
    """Reconstruct the width-one 2-by-2 row network from its node tables."""
    channels = make_channels(1, piece, "row")
    states = tuple(product(range(2), repeat=5))

    def entry(table, a, b, c):
        return table[(a * 2 + b) * 2 + c]

    def top(z):
        a, c, x, y, e = z
        return entry(tables[0], a, c, x) * entry(tables[1], a, y, e)

    def bottom(z):
        f, c, u, v, e = z
        return entry(tables[2], c, u, f) * entry(tables[3], v, e, f)

    def left(z, label):
        _, c, x, y, e = z
        return int(c == (label & 1) and e == ((label >> 1) & 1)) * channels.left[label >> 2][2 * x + y]

    def right(z, label):
        _, c, u, v, e = z
        return int(c == (label & 1) and e == ((label >> 1) & 1)) * channels.right[label >> 2][2 * u + v]

    edge = Edge(0, 1, 2 + channels.channel_bits, left, right)
    return Network((states, states), (top, bottom), (edge,))


class CoarsenedNetworkTests(unittest.TestCase):
    def test_signed_network_internal_crossing_and_empty_group(self):
        network = signed_network()
        expected = network_value(network)
        self.assertEqual(expected, -720)
        first = network.edges[0]
        swapped = replace(first, left=first.right, right=first.left)
        self.assertEqual(network_value(replace(network, edges=(swapped,) + network.edges[1:])), 1440)
        # Includes an internal edge, reversed edge orientations, a zero-bit
        # factor, multiple factors per group pair, and an empty group.
        grouping = (1, 0, 1, 2)
        self.assertEqual(clique_value(network, grouping, 4, 2), expected)
        for modulus in (4, 6, 9, 15):
            self.assertEqual(clique_value(network, grouping, 4, 2, modulus), expected % modulus)
        with self.assertRaises(ValueError):
            clique_value(network, grouping, 4, 1)

    def test_actual_grid_pieces_through_coarsened_tables(self):
        tables = signed_tables(1, 100)
        coarsened_pieces = []
        for piece in range(2):
            network = grid_row_network(tables, piece)
            direct = network_contraction(tables, make_channels(1, piece, "row"))
            self.assertEqual(network_value(network), direct)
            coarsened_pieces.append(clique_value(network, (1, 0), 2, 6))
            # Split groups in reverse order, with two padding bits; then put
            # both nodes together, producing an internal edge and empty group.
            for grouping, width in (((1, 0), 6), ((0, 0), 0)):
                self.assertEqual(clique_value(network, grouping, 2, width), direct)
                for modulus in (4, 9):
                    self.assertEqual(clique_value(network, grouping, 2, width, modulus), direct % modulus)
        self.assertEqual(sum(coarsened_pieces), 8 * original_contraction(tables, 1))

    def test_bad_grouping_and_modulus(self):
        network = signed_network()
        for call in (lambda: clique_value(network, (0,), 2, 2),
                     lambda: clique_value(network, (0, 0, 0, 2), 2, 2),
                     lambda: clique_value(network, (0, 0, 0, 0), 2, 0, 1),
                     lambda: network_value(network, 1)):
            with self.assertRaises(ValueError):
                call()

    def test_verified_grouping_and_derived_common_width(self):
        network = signed_network()
        graph = WeightedGraph(len(network.states),
                              tuple((e.first, e.second, e.bits) for e in network.edges))
        grouping, groups, discrepancy = (0, 1, 1, 0), 2, 2
        cert = verify_grouping(graph, groups, discrepancy, grouping)
        self.assertTrue(cert.valid)
        self.assertFalse(cert.sufficient_threshold_met)
        width = ((discrepancy * graph.ordered_total + graph.vertex_count * graph.max_degree)
                 // (discrepancy * groups ** 2))
        self.assertEqual(width, 4)
        self.assertEqual(cert.pair_weights, ((0, 1, width),))
        self.assertEqual(clique_value(network, grouping, groups, width), -720)
        self.assertEqual(clique_value(network, grouping, groups, width, 9), -720 % 9)


if __name__ == "__main__":
    unittest.main()
