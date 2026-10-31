#!/usr/bin/env python3
"""Finite five-symbol crossing searches with exact F4 frame recentering.

The coefficient systems and complete profiles are exact.  Pressure/root
optimization is exploratory floating point, and no minimum is certified.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

from f4_recentered_lattices import RecenteredCrossingLattice, RecenteredFiber, add_equation


def root_for_profiles(profiles, energy_x, energy_y):
    """Numerical objective only; copied from the complete-alphabet predecessor.

    This grid and local minimization are not a certificate of the global minimum.
    A negative result here must pass a separate rigorous pressure check before use.
    """
    import numpy as np
    from scipy.optimize import brentq, minimize_scalar
    from scipy.special import logsumexp
    first, second = (np.asarray(values, dtype=float) for values in profiles)
    ex = np.asarray(energy_x, dtype=float)
    ey = np.asarray(energy_y, dtype=float)
    ex_min, ey_min = min(ex), min(ey)
    ex, ey = ex-ex_min, ey-ey_min

    def root(log_scale):
        scale = math.exp(log_scale)
        constant = logsumexp(-scale*first)+logsumexp(-scale*second)+scale*(ex_min+ey_min)
        if constant < 1e-12:
            return 1.

        def pressure(r):
            return constant-r*(logsumexp(-scale*ex/r)+logsumexp(-scale*ey/r))

        upper = 2.
        while pressure(upper) > 0 and upper < 64:
            upper *= 2
        if pressure(upper) > 0:
            return 32.
        return brentq(pressure, 1e-8, upper, xtol=2e-11)/2

    grid = np.linspace(-8., 1., 25)
    values = [root(point) for point in grid]
    candidates = [(value, float(point)) for point, value in zip(grid, values)]
    for index in sorted(range(len(grid)), key=lambda i: values[i])[:3]:
        lo, hi = grid[max(0, index-1)], grid[min(len(grid)-1, index+1)]
        if lo == hi:
            continue
        result = minimize_scalar(root, bounds=(lo, hi), method='bounded', options={'xatol': 1e-9})
        candidates.append((float(result.fun), float(result.x)))
    value, log_scale = min(candidates)
    scale = math.exp(log_scale)
    critical = logsumexp(-scale*first)+logsumexp(-scale*second)
    critical -= math.sqrt(2)*(logsumexp(-scale*np.asarray(energy_x)/math.sqrt(2))+
                             logsumexp(-scale*np.asarray(energy_y)/math.sqrt(2)))
    return {'rho': value, 'scale': scale, 'critical_pressure': float(critical)}


def five_letter_source():
    """The retained reflection code, computed directly in its compressed variable.

    This is the M=8 formula in reflection_dual_rail_gaussian.py, restricted to
    five letters, divided by its common w^16 factor, and written in z=w^2.
    A small Laurent-Pfaffian recurrence avoids an unrelated source import chain.
    """
    letters = (0, 2, 4, 8, 12)
    kernel = [[frozenset() for _ in range(16)] for _ in range(16)]

    def put(leg1, i, leg2, j, exponents):
        a, b = 4*leg1+i, 4*leg2+j
        assert a != b and not kernel[a][b]
        terms = set()
        for exponent in exponents:
            terms.symmetric_difference_update((exponent,))
        kernel[a][b] = kernel[b][a] = frozenset(terms)

    for leg in range(4):
        for i in range(2):
            put(leg, i, leg, 3-i, (0, 8 if leg in (2, 3) else -8))
    for i in range(4):
        for j in range(4):
            power = 2*(j & ~i)
            put(3, i, 1, j, (power, power+1))
            put(2, i, 0, j, (power, power+1))
            power = 2*(3-(i | j))
            put(3, i, 2, j, (power, power+1))
            power = -8+2*(i & j)
            put(1, i, 0, j, (power, power+1))
            if i == j:
                put(3, i, 0, j, (0,))
                put(1, i, 2, j, (0,))

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return frozenset((0,))
        if mask.bit_count() % 2:
            return frozenset()
        first = mask & -mask
        i = first.bit_length()-1
        rest = remaining = mask ^ first
        value = set()
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length()-1
            for power in kernel[i][j]:
                value.symmetric_difference_update(power+k for k in pfaffian(rest ^ chosen))
        return frozenset(value)

    values = []
    for index in range(5**4):
        mask = sum(letters[(index//5**i) % 5] << (4*i) for i in range(4))
        powers = tuple(power+8 for power in pfaffian(mask))
        assert all(0 <= power <= 18 for power in powers)
        values.append((sum(1 << power for power in powers), 0))
    assert max(a.bit_length()-1 for a, _ in values) == 18
    energy = (0, 10, 6, 2, 8)
    half = tuple(value//2 for value in energy)
    zero = (0,)*5
    return tuple(values), (zero, zero, half, half, half, half, zero, zero), energy, energy


def brauer_source(q=5):
    """Strict degree-one crossing from the two Brauer pairings R and E.

    R_xyuv=[x=y][u=v] and E_xyuv=[x=v][y=u].  The source is R+z(E+R).
    Coefficient searches may allow a larger degree than this reference.
    """
    values = []
    for index in range(q**4):
        x, y, u, v = ((index//q**leg) % q for leg in range(4))
        r, e = int(x == y and u == v), int(x == v and y == u)
        values.append((r | ((e ^ r) << 1), 0))
    energy = (0, 12, 6, 0, 8)
    half = tuple(value//2 for value in energy)
    zero = (0,)*q
    return tuple(values), (zero, zero, half, half, half, half, zero, zero), energy, energy


def increment(fiber, pivots, row, target):
    result = dict(pivots)
    for power, low, high, rhs in fiber.equations[row]:
        if power < target:
            continue
        if power > target:
            break
        if not add_equation(result, low, high, rhs):
            return None
    return result


def digest(values):
    return sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rounds', type=int, default=4)
    parser.add_argument('--trials-per-round', type=int, default=5000)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--lower-window', type=int, default=4)
    parser.add_argument('--upper-window', type=int, default=5)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--source', choices=('reflection-five', 'brauer-five'), default='reflection-five')
    parser.add_argument('--degree', type=int, help='Allowed polynomial degree; defaults to 18 for reflection and 6 for Brauer.')
    parser.add_argument('--paired-energies', help='Five even integer paired energies; default is the natural source gauge.')
    parser.add_argument('--seed', type=int, default=2026103192)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if min(args.rounds, args.trials_per_round, args.samples, args.lower_window, args.upper_window) < 1:
        parser.error('All counts and windows must be positive.')
    started = time.monotonic()
    deadline = started+args.seconds
    rng = random.Random(args.seed)
    values, heights, energy_x, energy_y = (five_letter_source() if args.source == 'reflection-five' else brauer_source())
    if args.paired_energies:
        energy_x = energy_y = tuple(map(int, args.paired_energies.split(',')))
        if len(energy_x) != 5 or any(value % 2 for value in energy_x):
            parser.error('The selected gauge requires five even integer paired energies.')
        half = tuple(value//2 for value in energy_x)
        zero = (0,)*5
        heights = (zero, zero, half, half, half, half, zero, zero)
    initial_digest = digest(values)
    q = 5
    degree = args.degree if args.degree is not None else (18 if args.source == 'reflection-five' else 6)
    if degree < max(a.bit_length()-1 for a, b in values):
        parser.error('The coefficient degree cannot be smaller than the source degree.')
    retained = []
    rounds = []
    used_references = set()
    total_queries = 0
    best = None
    script_root = Path(__file__).resolve().parent
    source_names = (Path(__file__).name, 'f4_recentered_lattices.py',
                    'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
                    'weighted_smith_oracle_f4.py')
    report = {'scope': 'Exact F4 coefficient fibers and literal endpoint/profile replays; finite stochastic searches and numerical roots only.',
              'arguments': vars(args) | {'output': str(args.output)}, 'alphabet': q, 'degree': degree,
              'reference_family': args.source, 'initial_source_sha256': initial_digest,
              'heights_infinity_then_zero': heights, 'paired_energies': (energy_x, energy_y),
              'source_sha256': {name: sha256((script_root/name).read_bytes()).hexdigest() for name in source_names}}
    if args.source == 'reflection-five':
        report.update(source_letters=(0, 2, 4, 8, 12), removed_common_monomial=16, removed_exponent_gcd=2)

    def save(completed=False):
        report.update(rounds=rounds, retained_witnesses=sorted(retained, key=lambda r: r['rho'])[:24],
                      complete_profile_queries=total_queries, best_recorded_rho=best,
                      completed_requested_schedule=completed, elapsed_seconds=time.monotonic()-started,
                      strict_pressure_candidate_found=any(r['critical_pressure'] < -1e-9 for r in retained))
        args.output.write_text(json.dumps(report, indent=2)+'\n')

    for round_index in range(args.rounds):
        if time.monotonic() >= deadline:
            break
        round_started = time.monotonic()
        reference_digest = digest(values)
        used_references.add(reference_digest)
        lattice = RecenteredCrossingLattice(values, q, degree, heights, args.upper_window)
        floors = tuple(max(0, order-args.lower_window) for order in lattice.orders)
        fiber = RecenteredFiber(lattice, floors)
        reference_profiles = lattice.check_witness(values, lattice.orders)
        reference_cost = root_for_profiles(reference_profiles, energy_x, energy_y)
        if best is None:
            best = reference_cost['rho']
        record = {'round': round_index, 'reference_values': values, 'reference_digest': reference_digest,
                  'nonbinary_reference': any(b for a, b in values), 'reference_profiles': reference_profiles,
                  'reference_cost': reference_cost, 'fiber_variables': fiber.variables,
                  'floor_targets': floors, 'ceiling_targets': lattice.caps, 'reference_targets': lattice.orders,
                  'construction_seconds': time.monotonic()-round_started, 'completed_trials': 0,
                  'profile_classes': [], 'literal_replays': 1}
        rounds.append(record)
        print(json.dumps({'event': 'constructed', 'round': round_index,
                          'nonbinary_reference': record['nonbinary_reference'],
                          'fiber_variables': fiber.variables, 'rho': reference_cost['rho'],
                          'seconds': time.monotonic()-started}), flush=True)
        parents = [lattice.orders]
        seen = set()
        classes = {}
        candidates = []
        for trial in range(args.trials_per_round):
            if time.monotonic() >= deadline:
                break
            record['completed_trials'] = trial+1
            if trial and trial % 100 == 0:
                record.update(distinct_target_vectors=len(seen),
                              profile_classes=[{'profiles': p, 'samples': count} for p, count in sorted(classes.items())],
                              elapsed_seconds=time.monotonic()-round_started)
                save()
                print(json.dumps({'event': 'iteration_checkpoint', 'round': round_index, 'trial': trial,
                                  'targets': len(seen), 'profile_classes': len(classes),
                                  'best_rho': best, 'seconds': time.monotonic()-started}), flush=True)
            targets = list(lattice.orders if trial == 0 else rng.choice(parents))
            if trial:
                for row in rng.sample(range(2*q*q), rng.randrange(2, q*q+1)):
                    targets[row] = max(floors[row], targets[row]-rng.randrange(1, args.lower_window+1))
            pivots = fiber.solve(targets)
            assert pivots is not None
            ordering = list(range(2*q*q))
            rng.shuffle(ordering)
            if trial % 4 == 1:
                ordering.sort(key=lambda row: targets[row]+rng.random()*args.lower_window)
            for row in ordering:
                while targets[row] < lattice.caps[row]:
                    proposal = increment(fiber, pivots, row, targets[row])
                    if proposal is None:
                        break
                    pivots = proposal
                    targets[row] += 1
            targets = tuple(targets)
            if targets in seen:
                continue
            seen.add(targets)
            if len(parents) < 128:
                parents.append(targets)
            elif rng.randrange(4) == 0:
                parents[rng.randrange(1, len(parents))] = targets
            for sample_index in range(args.samples):
                sample = fiber.sample(pivots, rng)
                profiles = lattice.profiles(sample)
                total_queries += 1
                if profiles in classes:
                    classes[profiles] += 1
                    # Same-profile but different F4 frames are a central purpose of this search.
                    if sample_index == 0 and len(candidates) < 64 and any(b for a, b in sample):
                        cost = root_for_profiles(profiles, energy_x, energy_y)
                        candidates.append({'values': sample, 'profiles': profiles, 'targets': targets,
                                           'round': round_index, 'trial': trial, 'digest': digest(sample), **cost})
                    continue
                classes[profiles] = 1
                cost = root_for_profiles(profiles, energy_x, energy_y)
                witness = {'values': sample, 'profiles': profiles, 'targets': targets,
                           'round': round_index, 'trial': trial, 'digest': digest(sample),
                           'fiber_dimension': fiber.variables-len(pivots), **cost}
                if any(b for a, b in sample):
                    candidates.append(witness)
                if len(retained) < 24 or cost['rho'] < max(r['rho'] for r in retained):
                    assert lattice.check_witness(sample, targets) == profiles
                    record['literal_replays'] += 1
                    retained.append(witness)
                    retained[:] = sorted(retained, key=lambda r: r['rho'])[:24]
                if cost['rho'] < best-1e-11:
                    best = cost['rho']
                    print(json.dumps({'event': 'new_best', 'round': round_index, 'trial': trial, **cost}), flush=True)
            if trial % 100 == 0:
                record.update(distinct_target_vectors=len(seen),
                              profile_classes=[{'profiles': p, 'samples': count} for p, count in sorted(classes.items())],
                              elapsed_seconds=time.monotonic()-round_started)
                save()
                print(json.dumps({'event': 'progress', 'round': round_index, 'trial': trial,
                                  'targets': len(seen), 'profile_classes': len(classes),
                                  'best_rho': best, 'seconds': time.monotonic()-started}), flush=True)
        record.update(distinct_target_vectors=len(seen),
                      profile_classes=[{'profiles': p, 'samples': count} for p, count in sorted(classes.items())],
                      completed_requested_trials=record['completed_trials'] == args.trials_per_round,
                      elapsed_seconds=time.monotonic()-round_started)
        candidates.sort(key=lambda r: (r['rho'], r['digest']))
        candidates = [r for r in candidates if r['digest'] not in used_references]
        if candidates:
            next_source = candidates[min(round_index, min(3, len(candidates)-1))]
            assert lattice.check_witness(next_source['values'], next_source['targets']) == next_source['profiles']
            record['literal_replays'] += 1
            values = next_source['values']
            record['next_reference_digest'] = next_source['digest']
        else:
            record['next_reference_digest'] = reference_digest
        save()
    completed = len(rounds) == args.rounds and all(r['completed_requested_trials'] for r in rounds)
    save(completed)
    print(json.dumps({'event': 'finished', 'completed': completed, 'rounds': len(rounds),
                      'queries': total_queries, 'best_rho': best,
                      'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
