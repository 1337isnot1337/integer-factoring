#!/usr/bin/env python3
"""Certified finite search of nine-term matchings from the first low-rank flag."""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from low_rank_two_copy_updates import LowRankTwoCopyFrames
from native_matrix_smith_f4 import MatrixOracle
from nine_term_flag_templates import ANTITONE, component_data, updates
from search_entangled_two_copy_updates import compact, convolution, prefix_dominated
from search_f4_recentered_crossings import five_letter_source, root_for_profiles
from two_copy_crossing_couplings import CoupledCopies


def recipes(count, seed):
    patterns = [(1,)*9, (2,)*9, tuple(1+i % 3 for i in range(9)),
                tuple(1+(i//3) % 3 for i in range(9)),
                tuple(1+(i+i//3) % 3 for i in range(9)),
                (1, 2, 1, 3, 1, 2, 1, 3, 1)]
    generated = 0
    for shape in ('binomial', 'endpoint_power'):
        for zero in range(16, 22):
            for coefficients in patterns:
                if generated == count:
                    return
                yield {'kind': 'antitone_grid', 'pairing': list(ANTITONE),
                       'zero_orders': [zero]*9, 'coefficients': list(coefficients), 'shape': shape}
                generated += 1
    rng = random.Random(seed)
    data = component_data()
    while generated < count:
        pairing = list(ANTITONE)
        if generated % 3:
            rng.shuffle(pairing)
        zero_orders = []
        for i, j in enumerate(pairing):
            a, b = data[i], data[j]
            maximum = 9+a['balanced_height_sum']+b['balanced_height_sum']
            lower = min(16, maximum)
            zero_orders.append(rng.randrange(lower, maximum+1))
        yield {'kind': 'seeded_variable_low_orders', 'pairing': pairing,
               'zero_orders': zero_orders, 'coefficients': [rng.randrange(1, 4) for _ in range(9)],
               'shape': rng.choice(('binomial', 'endpoint_power'))}
        generated += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', type=int, default=4000)
    parser.add_argument('--seconds', type=float, default=5400.)
    parser.add_argument('--seed', type=int, default=311026220)
    parser.add_argument('--replay-every', type=int, default=200)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.cases < 1 or args.replay_every < 1 or args.seconds <= 0:
        parser.error('Use positive counts and a positive time budget.')
    started = time.monotonic(); deadline = started+args.seconds
    source, heights, ex, ey = five_letter_source()
    compiler = CoupledCopies(5, source, 18, heights, precision=128)
    frames = LowRankTwoCopyFrames(compiler)
    # Exact pure-function cache; it changes no polynomial arithmetic or precision.
    original_transform = frames.transformed_vectors
    cached = lru_cache(maxsize=512)(original_transform)
    frames.transformed_vectors = lambda vectors, side: cached(tuple(map(tuple, vectors)), side)
    baseline_profiles = compiler.baseline_control()
    energy = ([a+b for a in ex for b in ex], [a+b for a in ey for b in ey])
    baseline = root_for_profiles(baseline_profiles, *energy)
    reference = convolution(baseline_profiles)
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'nine_term_flag_templates.py', 'low_rank_two_copy_updates.py',
             'native_certified_low_rank_smith_f4.py', 'native_certified_low_rank_smith_f4.cpp',
             'native_certified_rank_one_smith_f4.cpp', 'native_rectangular_smith_f4.cpp',
             'entangled_two_copy_frames.py', 'pure_two_copy_frames.py', 'native_tensor_frame_f4.py',
             'native_tensor_frame_f4.cpp', 'two_copy_crossing_couplings.py',
             'search_f4_recentered_crossings.py', 'search_entangled_two_copy_updates.py')
    record = {'scope': 'Exact complete profiles through tracked-error elimination and finite-jet transfer; finite template search only. Numerical rates remain exploratory.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names},
              'physical_alphabet': 25, 'complete_degree_bound': 36, 'frame_precision': 128,
              'baseline': baseline, 'baseline_profiles': compact(baseline_profiles),
              'paired_energies': energy, 'cases': [], 'profile_classes': [], 'unresolved': [],
              'original_matrix_replays': [], 'completed_requested_schedule': False}
    signatures = {}; counts = Counter(); best = {'kind': 'unchanged reference', **baseline}
    smallest_margin = None

    def save():
        record.update(elapsed_seconds=time.monotonic()-started, counters=dict(counts), best=best,
                      minimum_certified_pivot_error_margin=smallest_margin,
                      strict_pressure_candidate_found=best['critical_pressure'] < -1e-9)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(record, indent=2)+'\n'); temporary.replace(args.output)

    save()
    for index, recipe in enumerate(recipes(args.cases, args.seed)):
        if time.monotonic() >= deadline:
            record['time_budget_reached'] = True; break
        terms = updates(**{k: v for k, v in recipe.items() if k != 'kind'})
        witness_digest = sha256(json.dumps(terms, separators=(',', ':')).encode()).hexdigest()
        failure = None
        for precision in (256, 512, 1024):
            try:
                profiles, certificates = frames.profiles(terms, work_precision=precision)
                break
            except ArithmeticError as error:
                counts['unresolved_precision_attempts'] += 1
                failure = str(error)
        else:
            record['unresolved'].append({'index': index, **recipe, 'reason': failure})
            save(); continue
        counts['certified_complete_profiles'] += 2
        margins = [c['minimum_certified_pivot_error_margin'] for c in certificates]
        smallest_margin = min(margins+[smallest_margin] if smallest_margin is not None else margins)
        frame_margin = min(128-max(side)-2*compiler.scalars[cut] for cut, side in enumerate(profiles))
        assert frame_margin > 0
        signature = tuple(tuple(side) for side in profiles)
        if signature not in signatures:
            dominated = prefix_dominated(convolution(profiles), reference)
            item = {'first_case': index, 'profiles': compact(profiles),
                    'all_scale_prefix_dominated_by_reference': dominated}
            if not dominated:
                item.update(root_for_profiles(profiles, *energy)); counts['numerical_pressure_searches'] += 1
            signatures[signature] = len(record['profile_classes']); record['profile_classes'].append(item)
        profile_class = signatures[signature]; item = record['profile_classes'][profile_class]
        dominated = item['all_scale_prefix_dominated_by_reference']
        counts['dominated_queries' if dominated else 'nondominated_queries'] += 1
        improving = not dominated and item['rho'] < best['rho']-1e-11
        if index % args.replay_every == 0 or improving:
            original = tuple(compiler.integer_orders(MatrixOracle(matrix, 625, precision=128).profile()['orders'],
                                                     2*compiler.scalars[side])
                             for side, matrix in enumerate(frames.original_sum_matrices(terms)))
            assert original == profiles, (index, 'Complete original matrix disagreement')
            counts['original_matrix_profile_replays'] += 2
            record['original_matrix_replays'].append({'case': index, 'profiles': compact(original),
                                                       'literal_updates_sha256': witness_digest})
        case = {'index': index, **recipe, 'literal_updates_sha256': witness_digest,
                'profile_class': profile_class, 'work_precision': precision,
                'minimum_pivot_error_margin': min(margins), 'frame_tail_margin': frame_margin}
        record['cases'].append(case)
        if improving:
            best = {'kind': 'certified nine-term update', **case, **item, 'updates': terms}
            print(json.dumps({'improvement': best}), flush=True)
        if index % 10 == 0:
            save()
            print(json.dumps({'case': index, 'classes': len(signatures), 'best_rho': best['rho'],
                              'counters': dict(counts), 'elapsed_seconds': time.monotonic()-started}), flush=True)
    else:
        record['completed_requested_schedule'] = True
    save()


if __name__ == '__main__':
    main()
