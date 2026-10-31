#!/usr/bin/env python3
"""Growing code-coset filters with the complete multiplier envelope."""

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import random
from time import monotonic

import numpy as np

from boolean_multiplier_pressure import MultiplierEnvelope
from boolean_multiplier_profiles import code_group_signatures_for_types, signature_patterns
from search_symmetric_tiled_filters import TiledPressure


CODE_MATRICES = (
    (2, (1, 2, 2, 1, 3, 1, 3, 3)),
    (3, (1, 2, 4, 3, 7, 1, 1, 2)),
    (3, (1, 2, 4, 1, 2, 2, 3, 4)),
    (3, (1, 2, 4, 6, 2, 5, 1, 6)),
    (4, (1, 2, 4, 8, 9, 2, 7, 5)),
    (4, (1, 2, 4, 8, 4, 1, 15, 1)),
    (4, (1, 2, 4, 8, 8, 7, 4, 6)),
    (4, (1, 2, 4, 8, 13, 6, 2, 7)),
    (5, (1, 2, 4, 8, 16, 3, 23, 26)),
    (5, (1, 2, 4, 8, 16, 17, 18, 1)),
    (6, (1, 2, 4, 8, 16, 32, 47, 26)),
    (3, (1, 2, 3, 4, 5, 6, 7, 1)),
    (4, (8, 9, 10, 11, 12, 13, 14, 15)),
    (4, (3, 6, 12, 8, 1, 5, 10, 4)),
    (4, (1, 2, 4, 8, 3, 6, 12, 9)),
)


def digest(value):
    return sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def build_families(n, rank, seed):
    energies = tuple(range(1, 2 * n, 2))
    result = []
    for matrix_index, (matrix_rank, columns) in enumerate(CODE_MATRICES):
        if matrix_rank != rank:
            continue
        types = tuple(sorted(set(columns)))
        for count in sorted({8, 16, n // 2, n}):
            if count > n or count % 8:
                continue
            selected = energies[-count:]
            for name in ('contiguous', 'interleaved', 'reflected', 'shuffled'):
                ordered = list(selected)
                if name == 'reflected':
                    ordered = [value for pair in zip(selected[:count // 2], selected[count // 2:][::-1])
                               for value in pair]
                elif name == 'shuffled':
                    random.Random(seed + 104729 * matrix_index + 1009 * n + count).shuffle(ordered)
                if name == 'contiguous':
                    slots = [ordered[(index * count) // 8:((index + 1) * count) // 8]
                             for index in range(8)]
                else:
                    slots = [ordered[index::8] for index in range(8)]
                groups = [sorted(value for index, column in enumerate(columns) if column == group
                                 for value in slots[index]) for group in types]
                assert sorted(value for group in groups for value in group) == list(selected)
                result.append({
                    'matrix_index': matrix_index, 'rank': rank, 'source_columns': columns,
                    'active_modes': count, 'partition': name, 'column_types': types,
                    'actual_mode_energy_groups': groups,
                    'spectator_mode_energies': energies[:-count],
                })
    return tuple(result)


@lru_cache(maxsize=None)
def group_data(rank, types):
    signatures = code_group_signatures_for_types(rank, types)
    states = np.array([[(state >> bit) & 1 for bit in range(len(types))]
                       for state in range(1 << len(types))], dtype=bool)
    patterns = np.array([signature_patterns(signatures, rank, syndrome)
                         for syndrome in range(1 << rank)], dtype=np.int16)
    syndromes = []
    for state in range(len(states)):
        syndrome = 0
        for bit, column in enumerate(types):
            if state & (1 << bit):
                syndrome ^= column
        syndromes.append(syndrome)
    assert all(((int(patterns[syndrome, state]).bit_count() & 1) == (syndromes[state] == syndrome))
               for syndrome in range(1 << rank) for state in range(len(states)))
    return states, patterns, np.array(syndromes, dtype=np.int16)


def parity_probabilities(groups, n, tau, charge):
    """Both parity probabilities, retaining tiny tails without 1-exp loss."""
    even, odd = [], []
    for group in groups:
        arguments = tau * (np.asarray(group, dtype=float) / n + charge)
        with np.errstate(divide='ignore'):
            exponential = np.exp(-np.abs(arguments))
            log_magnitude = float(np.sum(np.log1p(-exponential) - np.log1p(exponential)))
        small = float(-np.expm1(log_magnitude) / 2)
        large = 1 - small
        negative = int((arguments < 0).sum()) & 1
        even.append(small if negative else large)
        odd.append(large if negative else small)
    return np.array(even), np.array(odd)


class GlobalCodePressure(TiledPressure):
    def __init__(self, n, rank, families):
        self.n, self.rank, self.families = n, rank, families
        self.modes, self.blocks, self.calls = n, [list(range(1, 2 * n, 2))], 0
        self.envelope = MultiplierEnvelope(rank)
        self.chosen_detail = {'kind': 'baseline'}

    @lru_cache(maxsize=2)
    def physical_supports(self, tau):
        result = []
        for family in self.families:
            states, _, syndromes = group_data(self.rank, tuple(family['column_types']))
            even, odd = parity_probabilities(family['actual_mode_energy_groups'], self.n, tau, 0)
            weights = np.prod(np.where(states, odd[None, :], even[None, :]), axis=1)
            result.extend(np.bincount(syndromes, weights=weights, minlength=1 << self.rank))
        return np.array(result)

    def distributions(self, tau, charge):
        result = []
        for family in self.families:
            states, patterns, _ = group_data(self.rank, tuple(family['column_types']))
            even, odd = parity_probabilities(family['actual_mode_energy_groups'], self.n, tau, charge)
            weights = np.prod(np.where(states, odd[None, :], even[None, :]), axis=1)
            result.extend(np.bincount(patterns[syndrome], weights=weights,
                                      minlength=1 << (self.rank + 1))
                          for syndrome in range(1 << self.rank))
        return np.array(result), self.physical_supports(float(tau))

    def relative_cost(self, tau, charge):
        layers, supports = self.distributions(tau, charge)
        costs, profiles, divergences = self.envelope.evaluate(layers, supports, tau)
        row = int(np.argmin(costs))
        cost = float(costs[row])
        if profiles[row] < 0:
            self.chosen_detail = {'kind': 'baseline'}
            return 0.0, np.array([-1]), np.array([0.0])
        family_index, syndrome = divmod(row, 1 << self.rank)
        self.chosen_detail = {
            'family_index': family_index, 'syndrome': syndrome,
            'support_probability': float(supports[row]),
            'layer_probabilities': layers[row].tolist(),
            **self.envelope.describe(int(profiles[row]), float(divergences[row])),
        }
        index = row * (len(self.envelope.profiles) + len(self.envelope.fixed)) + int(profiles[row])
        return cost, np.array([index]), np.array([cost])

    def release(self):
        self.physical_supports.cache_clear()
        self.envelope.exponential_arrays.cache_clear()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lengths', default='64,128,256')
    parser.add_argument('--ranks', default='2,3,4,5,6')
    parser.add_argument('--iterations', type=int, default=35)
    parser.add_argument('--seed', type=int, default=720491)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {
        'status': 'Running the listed outer searches.',
        'scope': 'Complete multiplier envelope for the explicit code supports; outer scale and charge optimization is exploratory.',
        'source_code_matrices': CODE_MATRICES, 'seed': args.seed,
        'optimizer_iterations': args.iterations,
        'outer_bounds': {'log_tau': [-4, float(np.log(15))], 'scaled_charge': [-6, 6]},
        'automata': {}, 'optimized_runs': [], 'strict_crossing_candidate_found': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for rank in map(int, args.ranks.split(',')):
        for n in map(int, args.lengths.split(',')):
            if n < 16 or n & (n - 1):
                raise ValueError('Use power-of-two physical lengths at least sixteen.')
            families = build_families(n, rank, args.seed)
            if not families:
                raise ValueError('The requested rank has no source matrix.')
            run = GlobalCodePressure(n, rank, families)
            automaton = run.envelope.automaton
            result['automata'][str(rank)] = {
                'sha256': digest(automaton), 'pattern_count': automaton['pattern_count'],
                'fixed_profiles': len(automaton['fixed_profiles']),
                'continuous_profiles': len(automaton['continuous_profiles']),
                'negative_gap_startup_rows': sum(len(tail['startup']) for tail in automaton['singleton_tails']),
                'negative_gap_periodic_rows': sum(len(tail['cycle']) for tail in automaton['singleton_tails']),
                'maximum_period': max(tail['period'] for tail in automaton['singleton_tails']),
            }
            started = monotonic()
            row = run.optimize(args.seed, args.iterations)
            row.update(n=n, rank=rank, elapsed_seconds=monotonic() - started,
                       support_families=families, support_families_sha256=digest(families),
                       selected_multiplier=run.chosen_detail)
            result['optimized_runs'].append(row)
            result['strict_crossing_candidate_found'] |= row['full_critical_pressure'] < -1e-9
            args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
            print('GLOBAL', rank, n, len(families) * (1 << rank), row['rho'],
                  row['relative_log_cost'], row['elapsed_seconds'], flush=True)
            run.release()
    result['status'] = 'All requested genuine-rho searches completed; no all-code pressure exclusion is inferred.'
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    main()
