#!/usr/bin/env python3
"""Independent endpoint, Toeplitz-rank, and recentering controls over F4."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_block_smith_crossings import binary_rank
from f4_recentered_lattices import (
    RecenteredCrossingLattice, RecenteredFiber, add_equation, assignment,
    dot, f4_rank, pack, polynomial_multiply, regular_frame, scale_row, unpack,
)
from check_polynomial_crossing_interpolation import (
    column_elimination, polynomial_matrix_product, smith_jet_changes, solve_in_basis,
)


ROOT = Path(__file__).resolve().parent.parent


def scalar_product(a, b):
    # Literal multiplication followed by omega**2 = omega+1.
    a0, a1 = a & 1, a >> 1
    b0, b1 = b & 1, b >> 1
    return (a0*b0 ^ a1*b1) | ((a0*b1 ^ a1*b0 ^ a1*b1) << 1)


def toeplitz_rank(matrix, precision):
    """Expand multiplication by every F4 coefficient into its binary 2x2 map."""
    n = len(matrix)
    rows = []
    for power in range(precision):
        for row in matrix:
            low = high = 0
            for source_power in range(power+1):
                for j, polynomial in enumerate(row):
                    c = (polynomial >> (2*(power-source_power))) & 3
                    for basis in range(2):
                        out = scalar_product(c, 1 << basis)
                        position = 2*(source_power*n+j)+basis
                        low |= (out & 1) << position
                        high |= ((out >> 1) & 1) << position
            rows.extend((low, high))
    return binary_rank(rows)


def random_crossing(q, degree, rng, binary=False):
    values = []
    for entry in range(q**4):
        x, y, u, v = ((entry//q**i) % q for i in range(4))
        a = rng.getrandbits(degree)
        b = 0 if binary else rng.getrandbits(degree)
        values.append((int(x == v and y == u) ^ a ^ (a << 1), b ^ (b << 1)))
    return tuple(values)


def direct_binary_fiber(lattice, targets):
    """Expand every basis vector literally in independent binary Smith frames."""
    frames = []
    base_blocks = []
    for side in range(2):
        source = [[unpack(value)[0] for value in row]
                  for row in lattice.local_matrix(lattice.base, side)]
        precision = lattice.precisions[side]
        orders, left, right = smith_jet_changes(source, precision)
        frames.append((left, right))
        base_blocks.append(polynomial_matrix_product(
            polynomial_matrix_product(left, source, precision), right, precision))

    def encode(blocks):
        value = offset = 0
        for side, matrix in enumerate(blocks):
            for i, row in enumerate(matrix):
                target = targets[side*lattice.size+i]
                for polynomial in row:
                    value ^= (polynomial & ((1 << target)-1)) << offset
                    offset += target
        return value

    columns = []
    for entry, bits in enumerate(lattice.labels):
        for power in range(lattice.degree):
            blocks = []
            for side, (left, right) in enumerate(frames):
                matrix = [[0]*lattice.size for _ in range(lattice.size)]
                i, j = lattice.location(bits, side)
                h = lattice.heights[4*side:4*side+4]
                shift = lattice.scalar_shifts[side] + sum(row[label] for row, label in zip(h, bits))
                for exponent in (power, power+1):
                    at = shift + (-exponent if side == 0 else exponent)
                    assert at >= 0
                    matrix[i][j] ^= 1 << at
                precision = lattice.precisions[side]
                blocks.append(polynomial_matrix_product(
                    polynomial_matrix_product(left, matrix, precision), right, precision))
            columns.append(encode(blocks))
    basis, kernel = column_elimination(columns)
    solution = solve_in_basis(basis, encode(base_blocks))
    return None if solution is None else lattice.variables-len(basis)


def run():
    rng = random.Random(2026103191)
    checks = {'scalar_products': 0, 'packed_roundtrips': 0,
              'affine_equations': 0, 'literal_endpoints': 0,
              'independent_toeplitz_ranks': 0, 'full_tensor_replays': 0,
              'recentered_nonbinary_references': 0, 'negative_controls': 0}
    for a, b in product(range(4), repeat=2):
        low, high = scale_row(a & 1, a >> 1, b)
        assert low | (high << 1) == scalar_product(a, b)
        assert polynomial_multiply(a, b) == scalar_product(a, b)
        checks['scalar_products'] += 1
    for _ in range(100):
        pair = rng.getrandbits(30), rng.getrandbits(30)
        assert unpack(pack(pair)) == pair
        checks['packed_roundtrips'] += 1
    for n in range(1, 9):
        low, high = rng.getrandbits(n), rng.getrandbits(n)
        pivots = {}
        rows = []
        for _ in range(2*n):
            a, b = rng.getrandbits(n), rng.getrandbits(n)
            rhs = dot(a, b, low, high)
            assert add_equation(pivots, a, b, rhs)
            rows.append((a, b, rhs))
        solution = assignment(pivots)
        assert all(dot(a, b, *solution) == rhs for a, b, rhs in rows)
        checks['affine_equations'] += len(rows)

    records = []
    cases = ((2, 2), (2, 3), (3, 2), (3, 3))
    for q, degree in cases:
        heights = [tuple(rng.randrange(-1, 2) for _ in range(q)) for _ in range(8)]
        values = random_crossing(q, degree, rng)
        assert any(b for a, b in values)
        lattice = RecenteredCrossingLattice(values, q, degree, heights, raise_cap=3)
        reference = lattice.check_witness(values, lattice.orders)
        checks['full_tensor_replays'] += 1
        for side in range(2):
            matrix = lattice.local_matrix(values, side)
            orders = [v+lattice.scalar_shifts[side] for v in reference[side]]
            for precision in sorted({1, 2, 3, max(orders)+1}):
                assert toeplitz_rank(matrix, precision) == 2*sum(max(0, precision-v) for v in orders)
                checks['independent_toeplitz_ranks'] += 1
        floors = tuple(max(0, value-2) for value in lattice.orders)
        fiber = RecenteredFiber(lattice, floors)
        target_records = []
        retained = None
        for trial in range(8):
            targets = tuple(rng.randrange(floor, min(cap, order+1)+1)
                            for floor, cap, order in zip(floors, lattice.caps, lattice.orders))
            direct, compressed = lattice.solve(targets), fiber.solve(targets)
            assert (direct is None) == (compressed is None)
            if direct is None:
                target_records.append({'feasible': False})
                continue
            assert lattice.variables-len(direct) == fiber.variables-len(compressed)
            sample = fiber.sample(compressed, rng)
            profiles = lattice.check_witness(sample, targets)
            checks['literal_endpoints'] += q**4
            checks['full_tensor_replays'] += 1
            target_records.append({'feasible': True, 'dimension': fiber.variables-len(compressed),
                                   'profiles': profiles})
            if any(b for a, b in sample):
                retained = sample
        if retained is None:
            retained = values
        second = RecenteredCrossingLattice(retained, q, degree, heights, raise_cap=2)
        assert second.check_witness(retained, second.orders) == lattice.profiles(retained)
        checks['recentered_nonbinary_references'] += 1
        checks['full_tensor_replays'] += 1
        checks['literal_endpoints'] += q**4
        records.append({'alphabet': q, 'degree': degree, 'variables': lattice.variables,
                        'floor_variables': fiber.variables, 'reference_profiles': reference,
                        'targets': target_records, 'heights': heights,
                        'nonbinary_frame_polynomial_entries': sum(bool(unpack(v)[1]) for frames in second.frames
                                                                  for matrix in frames for row in matrix
                                                                  for v in row)})
        try:
            lattice.solve(tuple(cap+1 for cap in lattice.caps))
        except ValueError:
            checks['negative_controls'] += 1
        else:
            raise AssertionError('Out-of-window targets were silently accepted.')

    # Independent full basis expansion uses the older binary jet routine and
    # avoids the optional Weyl-constructor import chain altogether.
    values = random_crossing(2, 3, rng, binary=True)
    heights = [(0, 1), (0, -1), (1, 0), (0, 1),
               (0, 0), (1, 0), (0, -1), (0, 1)]
    new = RecenteredCrossingLattice(values, 2, 3, heights, raise_cap=3)
    for trial in range(20):
        targets = tuple(rng.randrange(max(0, order-2), cap+1)
                        for order, cap in zip(new.orders, new.caps))
        first, second = direct_binary_fiber(new, targets), new.solve(targets)
        assert (first is None) == (second is None)
        if first is not None:
            assert first == new.variables-len(second)
    comparisons = {'status': 'passed', 'complete_fiber_queries': 20,
                   'method': 'Literal basis expansion in independent binary Smith jets.'}
    return {'scope': 'Exact finite F4 reference-frame and full-entry affine controls; no new factoring bound.',
            'checks': checks, 'cases': records, 'binary_predecessor_comparison': comparisons,
            'source_sha256': {name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('f4_recentered_lattices.py', Path(__file__).name)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'comparison': result['binary_predecessor_comparison']}))


if __name__ == '__main__':
    main()
