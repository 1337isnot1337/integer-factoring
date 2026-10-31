#!/usr/bin/env python3
"""Exact affine synthesis with nonconstant regular row flags over F4.

The catalogues are finite. A successful affine solve is checked against the
complete source tensor and both actual Smith profiles. No floating comparison
decides feasibility. This module extends, but does not overwrite, the fixed
constant-flag interface.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import hashlib
import json
import math
from pathlib import Path
import random
import time


from f4_lattice_linear import MUL, coefficient, dot, encode, eliminate
from search_binary_f4_flags import FlagSearch, insert
from search_binary_affine_crossings import critical, rate, physical_energies


def constant_row(row):
    return tuple(((value & 1), (value >> 1)) for value in row)


def polynomial_row(constant, higher=()):
    values = [list(pair) for pair in constant_row(constant)]
    for degree, row in higher:
        if degree < 1:
            raise ValueError('Higher row coefficients need a positive degree.')
        for i, value in enumerate(row):
            values[i][0] ^= (value & 1) << degree
            values[i][1] ^= (value >> 1) << degree
    return tuple(map(tuple, values))


def field_terms(pair):
    union = pair[0] | pair[1]
    while union:
        bit = union & -union
        degree = bit.bit_length()-1
        yield degree, ((pair[0] >> degree) & 1) | (((pair[1] >> degree) & 1) << 1)
        union ^= bit


def regular_frame(frame):
    rows = [encode([a & 1 | ((b & 1) << 1) for a, b in row]) for row in frame]
    return len(eliminate(rows)) == 4


class PolynomialFlagSearch(FlagSearch):
    @lru_cache(maxsize=None)
    def polynomial_equations(self, side, vector):
        combined = {}
        for (i, j, power), value in self.tables[side].items():
            for shift, scalar in field_terms(vector[i]):
                key = j, power+shift
                first, second = combined.get(key, (0, 0))
                combined[key] = (first ^ (value if scalar & 1 else 0),
                                 second ^ (value if scalar & 2 else 0))
        return tuple(sorted((power, equation) for (_, power), equation in combined.items()
                            if equation[0] | equation[1]))

    def polynomial_system(self, side, frame, targets, initial=None):
        if not regular_frame(frame):
            raise ValueError('The constant row frame is singular.')
        pivots = {} if initial is None else dict(initial)
        for vector, target in zip(frame, targets):
            for power, equation in self.polynomial_equations(side, vector):
                if power >= target:
                    break
                if not insert(pivots, equation, self.lattice.variables):
                    return None
        return pivots

    def assignment(self, seed):
        planes = []
        for plane in range(2):
            packed = 0
            for entry, pair in enumerate(seed['values']):
                remainder = pair[plane] ^ (self.lattice.base[entry] if plane == 0 else 0)
                quotient = 0
                for power in range(self.lattice.degree):
                    if (remainder >> power) & 1:
                        quotient |= 1 << power
                        remainder ^= 3 << power
                if remainder:
                    raise ValueError('The witness is not in the declared endpoint fiber.')
                packed |= quotient << (entry*self.lattice.degree)
            planes.append(packed)
        return planes[0] | (1 << self.lattice.variables), planes[1]

    def polynomial_verify(self, seed, frames, targets):
        augmented = self.assignment(seed)
        for side, frame in enumerate(frames):
            assert regular_frame(frame)
            for i, vector in enumerate(frame):
                for power, equation in self.polynomial_equations(side, vector):
                    if power >= targets[4*side+i]:
                        break
                    assert dot(equation, augmented) == 0
        profiles = self.lattice.profiles(seed)
        for side in range(2):
            bounds = sorted(targets[4*side+i]-self.lattice.scalar_shifts[side] for i in range(4))
            for k in range(1, 5):
                assert sum(profiles[side][:k]) >= sum(bounds[:k])
        return profiles


def lower_jet_frames(power=1, field_order=4):
    """All I+z^power N with N strictly lower triangular over F2 or F4."""
    if power < 1 or field_order not in (2, 4):
        raise ValueError('Use positive jet degree and field order two or four.')
    positions = [(i, j) for i in range(4) for j in range(i)]
    for coefficients in product(range(field_order), repeat=len(positions)):
        rows = [[0]*4 for _ in range(4)]
        for (i, j), value in zip(positions, coefficients):
            rows[i][j] = value
        yield tuple(polynomial_row(tuple(int(i == j) for j in range(4)),
                                   ((power, rows[i]),)) for i in range(4))


def frame_digest(frame):
    return hashlib.sha256(json.dumps(frame, separators=(',', ':')).encode()).hexdigest()


def intersect(first, second, variables):
    pivots = dict(first)
    for a, b, rhs in second.values():
        if not insert(pivots, (a | ((rhs & 1) << variables),
                               b | ((rhs >> 1) << variables)), variables):
            return None
    return pivots


def target_catalogue(lattice, transfers):
    targets = {lattice.orders}
    for i in range(8):
        row = list(lattice.orders); row[i] += 1; targets.add(tuple(row))
    if transfers:
        for i, j in product(range(8), repeat=2):
            if i == j:
                continue
            for step in range(1, transfers+1):
                row = list(lattice.orders); row[i] += step; row[j] -= step
                if all(row[4*s:4*s+4] == sorted(row[4*s:4*s+4]) for s in range(2)):
                    targets.add(tuple(row))
    return sorted(targets)


def main():
    # SciPy is optional for the exact checkers that import this module.
    from scipy.optimize import minimize_scalar
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=12)
    parser.add_argument('--baseline-shift', type=int, default=4)
    parser.add_argument('--jet-power', type=int, default=1)
    parser.add_argument('--field-order', type=int, choices=(2, 4), default=4)
    parser.add_argument('--transfers', type=int, default=3)
    parser.add_argument('--samples', type=int, default=4)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic(); deadline = start+args.seconds
    engine = PolynomialFlagSearch(args.degree, args.baseline_shift)
    lattice = engine.lattice
    frames = tuple(lower_jet_frames(args.jet_power, args.field_order))
    targets = target_catalogue(lattice, args.transfers)
    energies = physical_energies(lattice.heights)
    # This numerical filter defines the schedule only; it proves no exclusion
    # of omitted targets and has no role in the exact feasibility decisions.
    def pressure(t):
        profiles = [tuple(t[4*s+i]-lattice.scalar_shifts[s] for i in range(4)) for s in range(2)]
        return critical(profiles, energies, .29277969091981376)
    base_pressure = pressure(lattice.orders)
    targets = [lattice.orders]+sorted(t for t in targets if t != lattice.orders and pressure(t) <= base_pressure+1e-12)
    rng = random.Random(311631)
    result = {'status': 'Finite exact nonconstant regular-flag synthesis; no factoring improvement claimed.',
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'degree': args.degree, 'baseline_shift': args.baseline_shift,
              'height_ray': lattice.heights, 'jet_power': args.jet_power,
              'field_order': args.field_order, 'frames_per_cut': len(frames),
              'catalogue': 'I + z^jet_power N, all strictly lower triangular N',
              'target_count': len(targets), 'targets': targets,
              'variants': [], 'witnesses': [], 'completed': False}
    systems = {}; classes = {}
    def save():
        result['elapsed_seconds'] = time.monotonic()-start
        result['witnesses'] = list(classes.values())
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    for target in targets:
        if time.monotonic() >= deadline:
            break
        common = engine.common_system(target)
        common_key = None if common is None else tuple(sorted(common.items()))
        valid = []
        interrupted = False
        for side in range(2):
            key = side, target[4*side:4*side+4], common_key
            if key not in systems:
                survivors = []
                if common is not None:
                    for i, frame in enumerate(frames):
                        if time.monotonic() >= deadline:
                            interrupted = True; break
                        system = engine.polynomial_system(side, frame, key[1], common)
                        if system is not None:
                            survivors.append((i, system))
                if interrupted:
                    break
                systems[key] = survivors
            valid.append(systems[key])
        if interrupted:
            break
        feasible = comparisons = samples = 0
        for i, first in valid[0]:
            for j, second in valid[1]:
                if time.monotonic() >= deadline:
                    interrupted = True; break
                comparisons += 1
                pivots = intersect(first, second, lattice.variables)
                if pivots is None:
                    continue
                feasible += 1
                for _ in range(args.samples):
                    seed = engine.sample(pivots, rng)
                    profiles = engine.polynomial_verify(seed, (frames[i], frames[j]), target)
                    samples += 1
                    if profiles not in classes:
                        fit = minimize_scalar(lambda x: rate(profiles, energies, math.exp(x)),
                                              bounds=(-6., 2.), method='bounded')
                        classes[profiles] = {'values': seed['values'], 'profiles': profiles,
                                             'target': target, 'frames': [frames[i], frames[j]],
                                             'frame_indices': [i, j],
                                             'affine_dimension': lattice.variables-len(pivots),
                                             'fixed_ray_rho': float(fit.fun)}
            if interrupted:
                break
        if interrupted:
            break
        result['variants'].append({'target': target,
                                    'common_rank': None if common is None else len(common),
                                    'single_cut_survivors': [len(v) for v in valid],
                                    'pairs_decided': len(frames)**2,
                                    'pair_intersections': comparisons,
                                    'feasible_pairs': feasible, 'samples_checked': samples})
        save()
        print(json.dumps({'targets_done': len(result['variants']), 'survivors': [len(v) for v in valid],
                          'feasible_pairs': feasible, 'profile_classes': len(classes),
                          'elapsed_seconds': result['elapsed_seconds']}), flush=True)
    result['completed'] = len(result['variants']) == len(targets)
    save()


if __name__ == '__main__':
    main()
