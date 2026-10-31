"""Behavioral controls for the cancellation-module order certificate."""

from fractions import Fraction
from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import check_cancellation_module_automata as module


class CancellationModuleTests(unittest.TestCase):
    def test_dense_nilpotent_cancellation_is_not_a_path_cost(self):
        one = module.ONE
        matrix = ((one, one), (one, one))
        for precision in (1, 2, 4):
            automaton = module.ModuleAutomaton(((matrix,),), precision)
            self.assertEqual(module.minplus_bound(automaton, (0,) * 7), 0)
            self.assertEqual(automaton.word_bound((0,) * 7), 6 * precision)
            self.assertIsNone(module.literal_word_order(automaton, (0, 0)))
            transfer, terminal = automaton.weighted_matrix((0,), Fraction(2))
            self.assertEqual(module.weighted_value(transfer, terminal, 7), Fraction(1, 2 ** (6 * precision)))

    def test_full_tail_prevents_false_exact_zero(self):
        one, zero = module.ONE, module.ZERO
        actual_row = (one, module.rat(1 << 3))
        matrix = ((zero, zero), (one, zero))
        self.assertEqual(module.rational_vector_product((one, zero), matrix), (zero, zero))
        automaton = module.ModuleAutomaton(((matrix,),), 3, initial=actual_row)
        self.assertEqual(automaton.word_bound((0,)), 3)
        self.assertEqual(module.literal_word_order(automaton, (0,)), 3)

    def test_terminal_cancellation_and_precision_hierarchy(self):
        one, zero = module.ONE, module.ZERO
        channels = (
            (((one, one), (zero, one)),),
            (((one, zero), (one, module.rat(5, 3))),),
        )
        automata = [module.ModuleAutomaton(channels, precision, terminal=(one, one))
                    for precision in (1, 2, 3)]
        for length in range(1, 5):
            for word in product(range(2), repeat=length):
                bounds = [automaton.word_bound(word) for automaton in automata]
                self.assertEqual(bounds, sorted(bounds))
                exact = module.literal_word_order(automata[0], word)
                for precision, bound in enumerate(bounds, 1):
                    self.assertEqual(min(bound, precision), precision if exact is None else min(exact, precision))

    def test_free_boundary_word_counts_are_submultiplicative(self):
        one, zero = module.ONE, module.ZERO
        channels = (
            (((one, one), (one, one)),),
            (((one, zero), (zero, module.rat(2))),),
        )
        for precision in (1, 2, 3):
            automaton = module.ModuleAutomaton(channels, precision, free_boundaries=True)
            transfer, terminal = automaton.weighted_matrix((1, -1), Fraction(3, 2))
            counts = {length: module.weighted_value(transfer, terminal, length) for length in range(1, 5)}
            for a, b in ((1, 1), (1, 2), (2, 1), (2, 2)):
                self.assertLessEqual(counts[a + b], counts[a] * counts[b])


if __name__ == "__main__":
    unittest.main()
