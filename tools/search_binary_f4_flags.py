#!/usr/bin/env python3
"""Search nonbinary constant flags in the two unrestricted crossing frames.

Coefficient equations are exact over F4; every polynomial crossing has
the literal permutation endpoint. Targets, sampled flags, and numerical
pressure comparisons are finite and explicitly recorded.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import random
import time

from scipy.optimize import minimize_scalar

from f4_lattice_linear import MUL, INV, add, coefficient, dot, eliminate, encode, scale
from search_binary_two_cut_lattices import TwoCutLattice, BASE
from search_binary_flag_lattices import equation_rows
from search_binary_affine_crossings import critical, rate, physical_energies


def projective(dimension):
    return tuple(values for values in product(range(4), repeat=dimension)
                 if next((value for value in values if value), 0) == 1)


@lru_cache(maxsize=1)
def flags():
    result = []
    def visit(rows, pivots):
        if len(rows) == 4:
            result.append(tuple(reversed(rows)));return
        free = [i for i in range(4) if i not in pivots]
        for values in projective(len(free)):
            vector = [0]*4
            for i, value in zip(free, values):
                vector[i] = value
            visit(rows+[tuple(vector)], eliminate(list(pivots.values())+[encode(vector)]))
    visit([], {})
    assert len(result) == 8925
    return tuple(result)


def insert(pivots, equation, variables):
    mask = (1 << variables)-1
    rhs = coefficient(equation, variables)
    row = equation[0] & mask, equation[1] & mask
    while row[0] | row[1]:
        pivot = (row[0] | row[1]).bit_length()-1
        value = coefficient(row, pivot)
        if pivot not in pivots:
            inv = INV[value]
            pivots[pivot] = (*scale(row, inv), MUL[rhs][inv])
            return True
        a, b, old_rhs = pivots[pivot]
        row = add(row, scale((a, b), value));rhs ^= MUL[value][old_rhs]
    return rhs == 0


class FlagSearch:
    def __init__(self, degree, baseline_shift=0):
        self.lattice = TwoCutLattice(degree, baseline_shift=baseline_shift)
        self.tables = [equation_rows(self.lattice, side, True) for side in range(2)]

    @lru_cache(maxsize=None)
    def row(self, side, vector):
        combined = {}
        for (i, j, power), value in self.tables[side].items():
            scalar = vector[i]
            if not scalar:
                continue
            key = j, power
            first, second = combined.get(key, (0, 0))
            combined[key] = (first ^ (value if scalar & 1 else 0),
                             second ^ (value if scalar & 2 else 0))
        return tuple(sorted((power, equation) for (j, power), equation in combined.items() if equation[0] | equation[1]))

    def common_system(self, targets):
        """Necessary two-cut minimum constraints, independent of both flags."""
        pivots = {}
        for side in range(2):
            minimum = min(targets[4*side:4*side+4])
            for i in range(4):
                vector = tuple(int(i == j) for j in range(4))
                for power, equation in self.row(side, vector):
                    if power >= minimum:
                        break
                    if not insert(pivots, equation, self.lattice.variables):
                        return None
        return pivots

    def system(self, side, frame, targets, initial=None):
        pivots = {} if initial is None else dict(initial)
        for vector, target in zip(frame, targets):
            for power, equation in self.row(side, vector):
                if power >= target:
                    break
                if not insert(pivots, equation, self.lattice.variables):
                    return None
        return pivots

    def sample(self, pivots, rng):
        assignment = [0]*self.lattice.variables
        for i in range(self.lattice.variables):
            if i not in pivots:
                assignment[i] = rng.randrange(4)
        planes = encode(assignment)
        for pivot, (a, b, rhs) in sorted(pivots.items()):
            value = dot((a, b), planes) ^ rhs
            assignment[pivot] = value
            planes = planes[0] ^ ((value & 1) << pivot), planes[1] ^ ((value >> 1) << pivot)
        return self.lattice.seed(*planes)

    def verify(self, seed, frame_pair, targets):
        # Recover the finite coefficient assignments from BASE+(1+w)Q.
        plane_assignments = []
        for plane in range(2):
            packed = 0
            for entry, pair in enumerate(seed['values']):
                remainder = pair[plane] ^ (self.lattice.base[entry] if plane == 0 else 0)
                quotient = 0
                for power in range(self.lattice.degree):
                    if (remainder >> power) & 1:
                        quotient |= 1 << power; remainder ^= 3 << power
                assert remainder == 0
                packed |= quotient << (entry*self.lattice.degree)
            plane_assignments.append(packed)
        augmented = plane_assignments[0] | (1 << self.lattice.variables), plane_assignments[1]
        for side, frame in enumerate(frame_pair):
            for i, vector in enumerate(frame):
                for power, equation in self.row(side, vector):
                    if power >= targets[4*side+i]:
                        break
                    assert dot(equation, augmented) == 0
        profiles = self.lattice.profiles(seed)
        for side in range(2):
            ordered = sorted(targets[4*side+i]-self.lattice.scalar_shifts[side] for i in range(4))
            for k in range(1, 5):
                assert sum(profiles[side][:k]) >= sum(ordered[:k])
        return profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=5)
    parser.add_argument('--baseline-shift', type=int, default=0)
    parser.add_argument('--flags', type=int, default=384)
    parser.add_argument('--samples', type=int, default=8)
    parser.add_argument('--max-transfer', type=int, default=3)
    parser.add_argument('--target-box', type=int, default=0,
                        help='Also test all zero-sum eight-row changes in this integer box.')
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=310781)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.max_transfer < 1:
        raise ValueError('The maximum row-target transfer must be positive.')
    if args.target_box < 0:
        raise ValueError('The row-target box radius must be nonnegative.')
    rng = random.Random(args.seed);start = time.monotonic();deadline = start+args.seconds
    engine = FlagSearch(args.degree,args.baseline_shift);lattice = engine.lattice
    all_flags = flags(); selected = rng.sample(range(len(all_flags)), min(args.flags,len(all_flags)))
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    selected = sorted(set(selected+[all_flags.index(identity)]))
    frames = [all_flags[i] for i in selected]
    energies = physical_energies(lattice.heights)
    targets = {lattice.orders}
    for first, second in product(range(8),repeat=2):
        if first == second:
            continue
        for step in range(1,args.max_transfer+1):
            target = list(lattice.orders);target[first] += step;target[second] -= step
            if any(sorted(target[4*s:4*s+4]) != target[4*s:4*s+4] for s in range(2)):
                continue
            targets.add(tuple(target))
    if args.target_box:
        radius = args.target_box
        for prefix in product(range(-radius,radius+1),repeat=7):
            last = -sum(prefix)
            if abs(last) > radius:
                continue
            target = tuple(a+b for a,b in zip(lattice.orders,prefix+(last,)))
            if all(tuple(sorted(target[4*s:4*s+4])) == target[4*s:4*s+4] for s in range(2)):
                targets.add(target)
    def target_pressure(target):
        profiles = [tuple(target[4*s+i]-lattice.scalar_shifts[s] for i in range(4)) for s in range(2)]
        return critical(profiles, energies, .29277969091981376)
    base_pressure = target_pressure(lattice.orders)
    raw_target_count = len(targets)
    targets = sorted((t for t in targets if target_pressure(t) <= base_pressure+1e-12), key=target_pressure)
    targets = [lattice.orders]+[t for t in targets if t != lattice.orders]
    result = {'status':'Finite exact F4 flags and polynomial lattice feasibility; exploratory pressures.',
              'degree':args.degree,'baseline_laurent_clearing':args.baseline_shift,
              'flag_indices':selected,'available_flags':len(all_flags),
              'target_count':len(targets),'random_seed':args.seed,'variants':[],'profiles':[],
              'maximum_row_target_transfer':args.max_transfer,
              'zero_sum_target_box_radius':args.target_box,
              'raw_sorted_target_count':raw_target_count,
              'feasibility_pruning':'Necessary two-cut minimum constraints are imposed before either flag.'}
    systems = {}
    classes = {}
    for target in targets:
        if time.monotonic() >= deadline:
            break
        common = engine.common_system(target)
        common_key = None if common is None else tuple(sorted(common.items()))
        pair = []
        for side in range(2):
            key = side, target[4*side:4*side+4], common_key
            if key not in systems:
                systems[key] = ([None]*len(frames) if common is None else
                                [engine.system(side, frame, key[1], common) for frame in frames])
            pair.append(systems[key])
        # A pair containing an inconsistent single-cut system is already
        # decided. Only intersect the surviving affine systems explicitly.
        attempted = len(frames)**2
        feasible = queries = 0
        valid = [[(i, value) for i, value in enumerate(side) if value is not None] for side in pair]
        for i, first in valid[0]:
            for j, second in valid[1]:
                pivots = dict(first)
                if not all(insert(pivots, (a | ((rhs & 1) << lattice.variables), b | ((rhs >> 1) << lattice.variables)), lattice.variables)
                           for a,b,rhs in second.values()):
                    continue
                feasible += 1
                for sample in range(args.samples):
                    seed = engine.sample(pivots,rng)
                    profiles = engine.verify(seed,(frames[i],frames[j]),target) if sample == 0 else lattice.profiles(seed)
                    queries += 1
                    if profiles not in classes:
                        fixed = minimize_scalar(lambda x:rate(profiles,energies,math.exp(x)),bounds=(-6.,2.),method='bounded')
                        classes[profiles]={'values':seed['values'],'profiles':profiles,'target':target,
                                          'flag_indices':[selected[i],selected[j]],'fixed_ray_rho':float(fixed.fun),
                                          'log_s':math.exp(fixed.x),'dimension':lattice.variables-len(pivots)}
        result['variants'].append({'targets':target,'flag_pairs_decided':attempted,
                                   'common_constraint_rank':None if common is None else len(common),
                                   'separately_feasible_cut_flags':[len(side) for side in valid],
                                   'feasible_pairs':feasible,'profile_queries':queries})
        result['profiles']=sorted(classes.values(),key=lambda row:row['fixed_ray_rho'])
        result['seconds']=time.monotonic()-start
        result['complete_requested_targets']=len(result['variants'])==len(targets)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'targets':len(result['variants']),'pairs':attempted,'feasible':feasible,
                          'classes':len(classes),'best':result['profiles'][0]['fixed_ray_rho'] if classes else None,
                          'seconds':result['seconds']}),flush=True)


if __name__=='__main__':
    main()
