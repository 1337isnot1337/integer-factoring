"""Tests for exact DVR row lattices and completed finite module orbits."""

import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import check_cancellation_module_automata as module
import check_exact_module_orbits as exact
import check_minimal_boundary_modules as minimal


class ExactModuleOrbitTests(unittest.TestCase):
    def test_rational_rank_deficient_lattice(self):
        row = ((1, 3), (2, 3))
        basis = exact.hermite_lattice((row, tuple(module.rshift(value, 2) for value in row)), 2)
        self.assertEqual(basis, (((1, 1), (2, 1)),))
        exact.membership(row, basis)
        with self.assertRaises(AssertionError):
            exact.membership((module.ZERO, module.ONE), basis)

    def test_nonunit_pivots_and_regular_remainders(self):
        basis = exact.hermite_lattice((((4, 1), (3, 1)), (module.ZERO, (2, 1))), 2)
        self.assertEqual(basis, (((4, 1), (1, 1)), (module.ZERO, (2, 1))))
        with self.assertRaises(AssertionError):
            exact.membership((module.ONE, module.ZERO), basis)

    def test_companion_closure_includes_rank_one_state(self):
        factors = module.local_factors(exact.companion_two(), 1, 0, (0, 0), False)
        reduced = minimal.regularize(
            minimal.minimize(factors["channels"], (module.ONE, module.ZERO),
                             (module.ONE, module.ZERO)), factors["exponents"])
        orbit = exact.ExactModuleOrbit(reduced["channels"])
        self.assertEqual(len(orbit.states), 7)
        self.assertEqual(sorted(set(map(len, orbit.states))), [1, 2])
        self.assertEqual(orbit.word_order((2,)), 0)
        with self.assertRaisesRegex(AssertionError, "state budget"):
            exact.ExactModuleOrbit(reduced["channels"], maximum_states=1)

    def test_zero_channel_is_an_exact_dead_transition(self):
        orbit = exact.ExactModuleOrbit(((((module.ZERO,),),),))
        self.assertEqual(orbit.word_order(()), 0)
        self.assertIsNone(orbit.word_order((0,)))
        self.assertEqual(len(orbit.states), 1)


if __name__ == "__main__":
    unittest.main()
