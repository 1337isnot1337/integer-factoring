#!/usr/bin/env python3
"""Literal polynomial controls for floor-preserving row transvections."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_f4_recentered_lattices import random_crossing
from f4_recentered_lattices import RecenteredCrossingLattice, RecenteredFiber, assignment, dot, multiply_matrix
from f4_smith_flag_fibers import SmithFlagChannels
from f4_polynomial_smith_flags import PolynomialSmithFlagFiber, random_transvections


def run():
    rng = random.Random(2026103197)
    counts = {'polynomial_flags': 0, 'transvections': 0, 'literal_coefficients': 0,
              'feasible_witnesses': 0, 'identity_flags': 0, 'negative_controls': 0}
    records = []
    for q, degree in ((2, 2), (2, 3), (3, 2)):
        values = random_crossing(q, degree, rng)
        heights = [tuple(rng.randrange(-1, 3) for _ in range(q)) for _ in range(8)]
        lattice = RecenteredCrossingLattice(values, q, degree, heights, raise_cap=4)
        floors = tuple(max(0, order-2) for order in lattice.orders)
        fiber = RecenteredFiber(lattice, floors, compression='tables')
        channels = SmithFlagChannels(lattice, fiber)
        identity = PolynomialSmithFlagFiber(channels, ())
        assert identity.equations == fiber.equations
        counts['identity_flags'] += 1
        operations_record = []
        for flag in range(5):
            operations = random_transvections(channels, rng, sweeps=2, extra_degree=3)
            mixed = PolynomialSmithFlagFiber(channels, operations)
            counts['polynomial_flags'] += 1
            counts['transvections'] += len(operations)
            operations_record.append(operations)
            assert mixed.check_witness(values, lattice.orders) == lattice.profiles(values)
            free_low, free_high = rng.getrandbits(fiber.variables), rng.getrandbits(fiber.variables)
            original_low = sum(((free_low >> i) & 1) << j for i, j in enumerate(fiber.free))
            original_high = sum(((free_high >> i) & 1) << j for i, j in enumerate(fiber.free))
            sample = lattice.seed(*assignment(fiber.base_pivots, original_low, original_high))
            for side, (left, right) in enumerate(mixed.frames):
                matrix = lattice.local_matrix(sample, side)
                transformed = multiply_matrix(multiply_matrix(left, matrix), right)
                for i in range(lattice.size):
                    row = side*lattice.size+i
                    for power in range(floors[row], lattice.caps[row]):
                        for j in range(lattice.size):
                            a, b, rhs = mixed.coefficient_channels[row][power-floors[row]][j]
                            assert dot(a, b, free_low, free_high) ^ rhs == (transformed[i][j] >> (2*power)) & 3
                            counts['literal_coefficients'] += 1
            for trial in range(8):
                targets = lattice.orders if trial == 0 else tuple(
                    rng.randrange(floor, min(cap, order+1)+1)
                    for floor, cap, order in zip(floors, lattice.caps, lattice.orders))
                pivots = mixed.solve(targets)
                if pivots is not None:
                    mixed.check_witness(mixed.sample(pivots, rng), targets)
                    counts['feasible_witnesses'] += 1
        bad_operations = (((0, 0, 0, 1, 1),), ((0, 0, 1, -1, 1),), ((0, 0, 1, 1, 0),))
        for operations in bad_operations:
            try:
                PolynomialSmithFlagFiber(channels, operations)
            except ValueError:
                counts['negative_controls'] += 1
            else:
                raise AssertionError('An invalid transvection was accepted.')
        records.append({'alphabet': q, 'degree': degree, 'heights': heights,
                        'orders': lattice.orders, 'floor_variables': fiber.variables,
                        'operations': operations_record})
    return {'scope': 'Exact finite F4 polynomial row-flag controls; no exponent improvement.',
            'checks': counts, 'cases': records,
            'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in ('f4_polynomial_smith_flags.py', 'f4_smith_flag_fibers.py',
                                           'f4_recentered_lattices.py', Path(__file__).name)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['checks']))


if __name__ == '__main__':
    main()
