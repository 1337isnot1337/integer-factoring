#!/usr/bin/env python3
"""Exact finite binary row-flag search in both unrestricted crossing cuts.

The 315 complete F2 flags in dimension four are enumerated at each cut.
Their polynomial frames are the fixed regular frames of the degree-five
compressed seed. Each flag pair gives affine coefficient equations, not
a floating feasibility decision. This is a scoped finite search.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import random
import time

from search_binary_two_cut_lattices import (TwoCutLattice, BASE, weighted_monomial,
                                           flatten, xor_all)
from weyl_lattice_filter_synthesis import matrix_multiply, polynomial_degrees, polynomial_multiply
from search_binary_affine_crossings import critical, physical_energies


def span(rows):
    values = {0}
    for row in rows:
        values |= {v ^ row for v in tuple(values)}
    return values


@lru_cache(maxsize=1)
def flags():
    representatives = {}
    for rows in product(range(1, 16), repeat=4):
        if len(span(rows)) != 16:
            continue
        key = tuple(sum(1 << value for value in span(rows[k:])) for k in (1, 2, 3))
        representatives.setdefault(key, rows)
    assert len(representatives) == 315
    return tuple(representatives.values())


def insert(pivots, equation, variables):
    row, rhs = equation & ((1 << variables)-1), equation >> variables
    while row:
        pivot = row.bit_length()-1
        if pivot not in pivots:
            pivots[pivot] = (row, rhs)
            return True
        old, value = pivots[pivot]
        row ^= old; rhs ^= value
    return rhs == 0


def equation_rows(lattice, side, raw=False):
    degree, heights, variables = lattice.degree, lattice.heights, lattice.variables
    left, right = lattice.bases[side]
    source = [[0]*4 for _ in range(4)]
    for entry, value in enumerate(lattice.base):
        i, j = flatten(entry, side)
        source[i][j] = xor_all(weighted_monomial(entry, p, side, degree, heights)
                               for p in polynomial_degrees(value))
    base = matrix_multiply(matrix_multiply(left, source), right)
    table = {}
    for i in range(4):
        for j in range(4):
            for power in polynomial_degrees(base[i][j]):
                table[i, j, power] = 1 << variables
    for entry in range(16):
        i, j = flatten(entry, side)
        for power in range(degree):
            coefficient = weighted_monomial(entry, power, side, degree, heights) ^ weighted_monomial(entry, power+1, side, degree, heights)
            bit = 1 << (entry*degree+power)
            for a in range(4):
                for b in range(4):
                    value = polynomial_multiply(polynomial_multiply(left[a][i], right[j][b]), coefficient)
                    for exponent in polynomial_degrees(value):
                        key = a, b, exponent
                        table[key] = table.get(key, 0) ^ bit
    if raw:
        return table
    result = {}
    for mask in range(1, 16):
        combined = {}
        for (i, j, power), value in table.items():
            if (mask >> i) & 1:
                key = j, power
                combined[key] = combined.get(key, 0) ^ value
        result[mask] = tuple(sorted(((power, value) for (j, power), value in combined.items() if value)))
    return result


def flag_systems(lattice, equations, targets):
    systems = []
    for frame in flags():
        pivots = {}
        valid = True
        for row, target in zip(frame, targets):
            for exponent, equation in equations[row]:
                if exponent >= target:
                    break
                if not insert(pivots, equation, lattice.variables):
                    valid = False;break
            if not valid:
                break
        systems.append(pivots if valid else None)
    return systems


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=5)
    parser.add_argument('--baseline-shift', type=int, default=0)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic(); rng = random.Random(310780)
    lattice = TwoCutLattice(args.degree, baseline_shift=args.baseline_shift)
    equations = [equation_rows(lattice, side) for side in range(2)]
    base_systems = [flag_systems(lattice, equations[side], lattice.orders[4*side:4*side+4]) for side in range(2)]
    result = {'status': 'Complete stated F2 flag pairs in two fixed polynomial frames; no global exclusion.',
              'degree': args.degree, 'baseline_laurent_clearing': args.baseline_shift,
              'flags_per_cut': len(flags()), 'base_row_orders': lattice.orders,
              'height_ray': lattice.heights, 'variants': [], 'profile_representatives': []}
    classes = {}
    for improved in range(-1, 8):
        systems = list(base_systems)
        if improved >= 0:
            side = improved//4; targets = list(lattice.orders[4*side:4*side+4]); targets[improved%4] += 1
            systems[side] = flag_systems(lattice, equations[side], targets)
        attempted = feasible = queries = 0
        for i, first in enumerate(systems[0]):
            for j, second in enumerate(systems[1]):
                attempted += 1
                if first is None or second is None:
                    continue
                pivots = dict(first)
                if not all(insert(pivots, row | (rhs << lattice.variables), lattice.variables)
                           for row, rhs in second.values()):
                    continue
                feasible += 1
                for _ in range(args.samples):
                    seed = lattice.sample(pivots, rng)
                    profiles = lattice.profiles(seed);queries += 1
                    if profiles not in classes:
                        energies = physical_energies(lattice.heights)
                        classes[profiles] = {'values': seed['values'], 'profiles': profiles,
                                             'flags': [flags()[i], flags()[j]], 'improved_row': improved,
                                             'dimension': lattice.variables-len(pivots),
                                             'minimum_sampled_critical_pressure': min(critical(profiles, energies, s) for s in (.1,.2,.3,.4,.6,1.))}
        row = {'improved_row': improved, 'flag_pairs_decided': attempted,
               'feasible_pairs': feasible, 'exact_profile_queries': queries}
        result['variants'].append(row)
        result['profile_representatives'] = sorted(classes.values(), key=lambda row: row['minimum_sampled_critical_pressure'])
        result['seconds'] = time.monotonic()-start
        args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({**row,'profile_classes':len(classes),'seconds':result['seconds']}),flush=True)


if __name__ == '__main__':
    main()
