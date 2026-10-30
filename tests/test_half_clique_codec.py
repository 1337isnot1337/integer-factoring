"""Independent finite tests of the ordinary half-rate table adapter."""

from itertools import product
import random
import unittest

from tools.half_clique_codec import (
    EqualityNetwork, MessageEdge, of_word, pinst_value, source_clique_value,
    source_tables_to_pinst, unflat,
)


class HalfCliqueCodecTests(unittest.TestCase):
    def test_off_diagonal_flatten_round_trip(self):
        for k in range(2, 5):
            for width in (1, 2):
                for vertex in range(k):
                    for index in range(1 << ((k - 1) * width)):
                        blocks = unflat(k, width, vertex, index)
                        self.assertEqual(blocks[vertex], ())
                        self.assertTrue(all(len(block) == width for j, block in enumerate(blocks)
                                            if j != vertex))
                        self.assertEqual(of_word(tuple(bit for block in blocks for bit in block)), index)

    def test_zero_diagonal_codec_on_arbitrary_signed_tables(self):
        rng = random.Random(53117)
        for k, width in ((2, 1), (2, 2), (3, 1)):
            for _ in range(8):
                source = tuple(tuple(rng.randrange(-4, 5) for _ in range(1 << (k * width)))
                               for _ in range(k))
                target = source_tables_to_pinst(source, width)
                self.assertEqual(source_clique_value(source, width), pinst_value(target, width))

    @staticmethod
    def oriented_network():
        # Group 1 contains gates 0 and 2, so edge 0--1 is reversed in group order.
        # Internal equality enforces z2=1-z0; gate 1 fixes the accepted boundary to integer 2.
        return EqualityNetwork(((1, 1), (1, 0), (1, 1)), (
            MessageEdge(0, 1, 1, (0, 1), (0, 1)),
            MessageEdge(0, 2, 1, (1, 0), (0, 1)),
            MessageEdge(1, 2, 1, (1, 0), (0, 1)),
        ))

    def test_grouping_internal_equalities_reversed_orientation_and_empty_group(self):
        network = self.oriented_network()
        for k, width in ((2, 2), (2, 3), (3, 2)):
            tables = network.group_tables((1, 0, 1), k, width)
            self.assertEqual(network.value(), 1)
            self.assertEqual(pinst_value(tables, width), 1)
            self.assertEqual(max(tables[0]), 1)
            if k == 3:
                self.assertEqual(tables[2][0], 1)
                self.assertEqual(sum(tables[2]), 1)

    def test_padding_bits_are_forced_to_zero(self):
        tables = self.oriented_network().group_tables((1, 0, 1), 2, 3)
        self.assertEqual(tables[0], (0, 0, 1, 0, 0, 0, 0, 0))
        self.assertEqual(tables[1], (0, 1, 1, 0, 0, 0, 0, 0))
        self.assertEqual(pinst_value(tables, 3), 1)

    def test_one_endpoint_chunk_reversal_is_detected(self):
        network = EqualityNetwork(((1, 0), (1, 0), (1, 1)), self.oriented_network().edges)
        tables = network.group_tables((1, 0, 1), 2, 2)
        corrupted = list(tables[0])
        corrupted[1], corrupted[2] = corrupted[2], corrupted[1]
        # Only the fixed gate's endpoint is reordered; the neighboring row keeps its ordering.
        self.assertEqual(network.value(), 1)
        self.assertEqual(pinst_value(tables, 2), 1)
        self.assertEqual(pinst_value((tuple(corrupted), tables[1]), 2), 0)

    def test_all_binary_message_maps_and_composite_moduli(self):
        # Independent state enumeration versus edge-label contraction. The two-edge path has
        # 2^8 endpoint-message maps and two opposite grouping orientations.
        for bits in product((0, 1), repeat=8):
            network = EqualityNetwork(((1, 1), (1, 1), (1, 1)), (
                MessageEdge(0, 1, 1, bits[:2], bits[2:4]),
                MessageEdge(1, 2, 1, bits[4:6], bits[6:]),
            ))
            exact = network.value()
            for groups in ((0, 1, 0), (1, 0, 1)):
                tables = network.group_tables(groups, 2, 2)
                self.assertEqual(pinst_value(tables, 2), exact)
                for modulus in (4, 6, 8, 9):
                    self.assertEqual(pinst_value(tables, 2, modulus), exact % modulus)

    def test_insufficient_width_is_rejected(self):
        with self.assertRaises(ValueError):
            self.oriented_network().group_tables((1, 0, 1), 2, 1)


if __name__ == "__main__":
    unittest.main()
