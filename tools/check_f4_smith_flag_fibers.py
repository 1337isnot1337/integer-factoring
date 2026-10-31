#!/usr/bin/env python3
"""Literal controls for common-floor-preserving F4 Smith row flags."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_f4_recentered_lattices import random_crossing
from f4_recentered_lattices import RecenteredCrossingLattice, RecenteredFiber, assignment, dot, multiply_matrix
from f4_smith_flag_fibers import SmithFlagChannels, add_scaled


def run():
    rng = random.Random(2026103195)
    counts = {'literal_mixed_coefficients': 0, 'complete_coefficient_channel_comparisons': 0,
              'feasible_witnesses': 0, 'nonbinary_row_flags': 0, 'negative_controls': 0}
    records = []
    for q, degree in ((2, 2), (2, 3), (3, 2)):
        values = random_crossing(q, degree, rng)
        heights = [(0,)*q]*8
        lattice = RecenteredCrossingLattice(values, q, degree, heights, raise_cap=3)
        floors = tuple(max(0, order-1) for order in lattice.orders)
        fiber = RecenteredFiber(lattice, floors, compression='tables')
        channels = SmithFlagChannels(lattice, fiber)
        counts['complete_coefficient_channel_comparisons'] += 1
        identity = channels.mix(channels.identity_mixes())
        assert identity.equations == fiber.equations
        mixtures = []
        for flag in range(4):
            matrices = channels.random_mixes(rng)
            mixed = channels.mix(matrices)
            assert any(c in (2, 3) for matrix in matrices for row in matrix for c in row)
            counts['nonbinary_row_flags'] += 1
            free_low, free_high = rng.getrandbits(fiber.variables), rng.getrandbits(fiber.variables)
            original_low = sum(((free_low >> i) & 1) << j for i, j in enumerate(fiber.free))
            original_high = sum(((free_high >> i) & 1) << j for i, j in enumerate(fiber.free))
            sample = lattice.seed(*assignment(fiber.base_pivots, original_low, original_high))
            for side, (left, right) in enumerate(mixed.frames):
                matrix = lattice.local_matrix(sample, side)
                transformed = multiply_matrix(multiply_matrix(left, matrix), right)
                for row in range(lattice.size):
                    first = side*lattice.size+row
                    for power in range(floors[first], lattice.caps[first]):
                        for column in range(lattice.size):
                            equation = (0, 0, 0)
                            for j, coefficient in enumerate(matrices[side][row]):
                                if coefficient:
                                    equation = add_scaled(equation,
                                        channels.channels[side*lattice.size+j][power-floors[first]][column], coefficient)
                            a, b, rhs = equation
                            assert dot(a, b, free_low, free_high) ^ rhs == (transformed[row][column] >> (2*power)) & 3
                            counts['literal_mixed_coefficients'] += 1
            found = 0
            for trial in range(8):
                targets = lattice.orders if trial == 0 else tuple(
                    rng.randrange(floor, min(cap, order+1)+1)
                    for floor, cap, order in zip(floors, lattice.caps, lattice.orders))
                pivots = mixed.solve(targets)
                if pivots is None:
                    continue
                sample = mixed.sample(pivots, rng)
                mixed.check_witness(sample, targets)
                found += 1
                counts['feasible_witnesses'] += 1
            mixtures.append({'matrices': matrices, 'feasible_witnesses': found})
        bad = [[list(row) for row in matrix] for matrix in channels.identity_mixes()]
        bad[0][0] = [0]*lattice.size
        try:
            channels.mix(bad)
        except ValueError:
            counts['negative_controls'] += 1
        else:
            raise AssertionError('A singular frame was accepted.')
        for side in range(2):
            orders = lattice.orders[side*lattice.size:(side+1)*lattice.size]
            pair = next(((i, j) for i in range(lattice.size) for j in range(lattice.size) if orders[i] != orders[j]), None)
            if pair:
                bad = [[list(row) for row in matrix] for matrix in channels.identity_mixes()]
                bad[side][pair[0]][pair[1]] = 1
                try:
                    channels.mix(bad)
                except ValueError:
                    counts['negative_controls'] += 1
                else:
                    raise AssertionError('A mixture across unequal orders was accepted.')
                break
        records.append({'alphabet': q, 'degree': degree, 'floor_variables': fiber.variables,
                        'orders': lattice.orders, 'blocks': channels.blocks, 'flags': mixtures})
    return {'scope': 'Exact finite common-floor coordinate and F4 row-frame controls; no exponent improvement.',
            'checks': counts, 'cases': records,
            'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in ('f4_smith_flag_fibers.py', 'f4_recentered_lattices.py', Path(__file__).name)}}


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
