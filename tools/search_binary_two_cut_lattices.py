#!/usr/bin/env python3
"""Synthesize unrestricted binary crossings in two fixed Smith row lattices.

Every entry is BASE + (1+w) times a free F4 polynomial. The literal
crossing endpoint is therefore fixed. Both complete weighted cut matrices
are transformed by exact regular F2 polynomial bases, and row-order
targets become affine binary equations. Returned witnesses use complete
coefficient arrays; floating pressure optimization is exploratory.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize_scalar

from native_rectangular_smith_f4 import RectangularOracle
from search_binary_affine_crossings import critical, rate, physical_energies
from search_binary_polynomial_neighborhood import BASE, HEIGHTS
from screen_full_alphabet_monomial_pressure import optimize
from weyl_lattice_filter_synthesis import regular_smith_basis, polynomial_multiply, polynomial_degrees


def labels(entry):
    return tuple((entry >> i) & 1 for i in range(4))


def flatten(entry, side):
    x, y, u, v = labels(entry)
    return (2*x+u, 2*y+v) if side == 0 else (2*x+y, 2*u+v)


def weighted_monomial(entry, power, side, degree, heights):
    weight = sum(a*b for a, b in zip(labels(entry), heights[4*side:4*side+4]))
    shift = -sum(min(0, h) for h in heights[4*side:4*side+4])
    exponent = weight + shift + (degree-power if side == 0 else power)
    assert exponent >= 0
    return 1 << exponent


def xor_all(values):
    value = 0
    for item in values:
        value ^= item
    return value


class TwoCutLattice:
    def __init__(self, degree=5, heights=HEIGHTS, baseline_shift=0):
        if (int(baseline_shift) != baseline_shift or baseline_shift < 0
                or degree < 5 + baseline_shift or len(heights) != 8
                or any(int(h) != h for h in heights)):
            raise ValueError('Use a nonnegative integer baseline shift, degree at least five plus that shift, and eight integer heights.')
        self.degree, self.heights, self.variables = degree, tuple(heights), 16*degree
        self.baseline_shift = int(baseline_shift)
        self.base = tuple(value << self.baseline_shift for value in BASE)
        self.equations, self.orders, self.scalar_shifts = [], [], []
        self.bases = []
        for side in range(2):
            source = [[0]*4 for _ in range(4)]
            for entry, value in enumerate(self.base):
                i, j = flatten(entry, side)
                source[i][j] = xor_all(weighted_monomial(entry, power, side, degree, heights)
                                       for power in polynomial_degrees(value))
            orders, left, right, diagonal = regular_smith_basis(source)
            shift = -sum(min(0, h) for h in heights[4*side:4*side+4]) + (degree if side == 0 else 0)
            self.scalar_shifts.append(shift)
            self.orders.extend(orders)
            table = defaultdict(int)
            for entry in range(16):
                i, j = flatten(entry, side)
                outer = [[polynomial_multiply(left[a][i], right[j][b]) for b in range(4)] for a in range(4)]
                for power in range(degree):
                    coefficient = weighted_monomial(entry, power, side, degree, heights) ^ weighted_monomial(entry, power+1, side, degree, heights)
                    bit = 1 << (entry*degree+power)
                    for a in range(4):
                        for b in range(4):
                            transformed = polynomial_multiply(outer[a][b], coefficient)
                            for exponent in polynomial_degrees(transformed):
                                table[a, b, exponent] ^= bit
            for a in range(4):
                row = []
                keys = {(b, exponent) for aa, b, exponent in table if aa == a}
                keys |= {(b, exponent) for b in range(4) for exponent in polynomial_degrees(diagonal[a][b])}
                for b, exponent in sorted(keys, key=lambda key: (key[1], key[0])):
                    row.append((exponent, table[a, b, exponent], (diagonal[a][b] >> exponent) & 1))
                self.equations.append(row)
            self.bases.append((left, right))
        self.orders = tuple(self.orders)

    def solve(self, targets):
        pivots = {}
        for target, equations in zip(targets, self.equations):
            for exponent, row, rhs in equations:
                if exponent >= target:
                    break
                while row:
                    pivot = row.bit_length()-1
                    if pivot not in pivots:
                        pivots[pivot] = (row, rhs)
                        break
                    old, value = pivots[pivot]
                    row ^= old; rhs ^= value
                if not row and rhs:
                    return None
        return pivots

    def vector(self, pivots, free_bits=0, homogeneous=False):
        value = free_bits
        for pivot, (row, rhs) in sorted(pivots.items()):
            assert not (value >> pivot) & 1
            bit = ((row & value).bit_count() & 1) ^ (0 if homogeneous else rhs)
            value |= bit << pivot
        return value

    def sample(self, pivots, rng, colored=True):
        free = [i for i in range(self.variables) if i not in pivots]
        first = sum(rng.randrange(2) << i for i in free)
        second = sum(rng.randrange(2) << i for i in free) if colored else 0
        return self.seed(self.vector(pivots, first), self.vector(pivots, second, True))

    def seed(self, low, high=0):
        values = []
        mask = (1 << self.degree)-1
        for entry, base in enumerate(self.base):
            a, b = (low >> (entry*self.degree)) & mask, (high >> (entry*self.degree)) & mask
            values.append((base ^ a ^ (a << 1), b ^ (b << 1)))
        return {'values': values, 'q_x': 2, 'q_u': 2, 'dimensions': [2]*4,
                'denominator_degree': 0, 'strict_endpoint': True}

    def profiles(self, seed):
        oracle = RectangularOracle(seed['values'], 2, 2, endpoint_check=True)
        _, _, orders = oracle.evaluate(list(self.heights)+[1.], (.5, .5), (.5, .5), 1, True)
        assert all(value == int(value) for value in orders)
        return tuple(map(int, orders[:4])), tuple(map(int, orders[4:]))

    def check_witness(self, seed, targets):
        from weyl_lattice_filter_synthesis import matrix_multiply
        for side, (left, right) in enumerate(self.bases):
            for plane in range(2):
                matrix = [[0]*4 for _ in range(4)]
                for entry, pair in enumerate(seed['values']):
                    i, j = flatten(entry, side)
                    matrix[i][j] = xor_all(weighted_monomial(entry, power, side, self.degree, self.heights)
                                           for power in polynomial_degrees(pair[plane]))
                transformed = matrix_multiply(matrix_multiply(left, matrix), right)
                for i, row in enumerate(transformed):
                    for value in row:
                        assert value & ((1 << max(0, targets[4*side+i]))-1) == 0
        profiles = self.profiles(seed)
        for side in range(2):
            bounds = sorted(value-self.scalar_shifts[side] for value in targets[4*side:4*side+4])
            for k in range(1, 5):
                assert sum(profiles[side][:k]) >= sum(bounds[:k])
        return profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=5)
    parser.add_argument('--baseline-shift', type=int, default=0,
                        help='Clear this many negative Laurent powers before the affine synthesis.')
    parser.add_argument('--targets', type=int, default=2000)
    parser.add_argument('--samples', type=int, default=40)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--top', type=int, default=30)
    parser.add_argument('--seed', type=int, default=310779)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    started = time.monotonic(); deadline = started+args.seconds
    lattice = TwoCutLattice(args.degree, baseline_shift=args.baseline_shift)
    energies = physical_energies(lattice.heights)
    baseline = lattice.profiles(lattice.seed(0))
    assert tuple(value-lattice.scalar_shifts[i//4] for i, value in enumerate(lattice.orders)) == baseline[0]+baseline[1]
    base_solution = lattice.solve(lattice.orders); assert base_solution is not None
    result = {'status': 'Exact two-cut lattice feasibility with exploratory pressure optimization; no global exclusion.',
              'degree': args.degree, 'baseline_laurent_clearing': args.baseline_shift,
              'base_coefficients': lattice.base, 'height_ray': lattice.heights,
              'base_profiles': baseline, 'row_orders': lattice.orders, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'variables': lattice.variables, 'base_lattice_dimension': lattice.variables-len(base_solution),
              'targets_attempted': 0, 'feasible_targets': 0, 'profile_queries': 0, 'independent_witness_replays': 0,
              'classes': [], 'optimized': []}
    classes, seen_targets = {}, set()
    scales = (.1, .2, .29277969091981376, .4, .6, 1.)
    def save():
        result['classes'] = sorted(classes.values(), key=lambda row: row['critical_grid'])[:500]
        result['profile_class_count'] = len(classes)
        result['elapsed_seconds'] = time.monotonic()-started
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    for trial in range(args.targets):
        if time.monotonic() > deadline:
            break
        targets = list(lattice.orders)
        if trial:
            for _ in range(rng.randint(1, 4)):
                index = rng.randrange(8); targets[index] -= rng.randint(1, 3)
            order = list(range(8)); rng.shuffle(order)
            for index in order:
                for _ in range(5):
                    proposed = targets.copy(); proposed[index] += 1
                    if lattice.solve(proposed) is None:
                        break
                    targets = proposed
        targets = tuple(targets)
        if targets in seen_targets:
            continue
        seen_targets.add(targets); result['targets_attempted'] += 1
        solution = lattice.solve(targets)
        if solution is None:
            continue
        result['feasible_targets'] += 1
        for sample in range(args.samples):
            seed = lattice.sample(solution, rng, colored=sample%3 != 0)
            profiles = lattice.check_witness(seed, targets) if sample == 0 else lattice.profiles(seed)
            result['profile_queries'] += 1
            result['independent_witness_replays'] += int(sample == 0)
            if profiles not in classes:
                pressures = [critical(profiles, energies, scale) for scale in scales]
                classes[profiles] = {'values': seed['values'], 'profiles': profiles, 'targets': targets,
                                     'dimension': lattice.variables-len(solution), 'count': 0,
                                     'critical_grid': min(pressures)}
            classes[profiles]['count'] += 1
        if trial%25 == 0:
            save();print(json.dumps({k: result[k] for k in ('targets_attempted','feasible_targets','profile_queries','profile_class_count','elapsed_seconds')}), flush=True)
    save()
    nrng = np.random.default_rng(args.seed)
    for row in result['classes'][:args.top]:
        if time.monotonic() > deadline:
            break
        fixed = minimize_scalar(lambda x: rate(row['profiles'], energies, math.exp(x)), bounds=(-6., 2.), method='bounded')
        scale = math.exp(fixed.x)
        seed = lattice.seed(0);seed['values'] = row['values']
        outcome = optimize(seed, cut=1, rng=nrng, restarts=5, maxfun=800, height_bound=48,
                           deadline=min(deadline,time.monotonic()+40),
                           starts=[[height*scale for height in lattice.heights]+[scale]])
        result['optimized'].append({'values': row['values'], 'fixed_ray_rho': float(fixed.fun), 'optimization': outcome})
        save()
    result['complete_requested_target_trials'] = trial+1 == args.targets
    save()


if __name__ == '__main__':
    main()
