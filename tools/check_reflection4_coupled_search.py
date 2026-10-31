#!/usr/bin/env python3
"""Replay all completed coordinate frames and independent winner profiles."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import flatten
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from search_reflection_physical_bases import packed_polynomial
from square_zero_coupled_reflection_frames import construct
from weighted_smith_oracle_f4 import WeightedSmithOracle


def rebuild(case, base, bell_data):
    if case['kind'] == 'identity':
        return base if case['coordinates'] == 'original' else bell_data[0]
    generator = [[0] * 4 for _ in range(4)]
    generator[case['target']][case['source']] = case['scalar']
    return construct(base, generator, pole=case['pole_order'],
                     restore=case['coordinates'] == 'original', bell_data=bell_data)


def rational_profiles(seed, variables):
    scale = float(variables[24])
    heights = [Fraction(str(float(value) / scale)).limit_denominator(1000)
               for value in variables[:24]]
    records = []
    for side in (0, 1):
        cut, infinity = 1 - side, side == 0
        h = [[Fraction(0)] + heights[side * 12 + leg * 3:side * 12 + (leg + 1) * 3]
             for leg in range(4)]
        rows = [h[0][a] + h[2 if cut else 1][b] for a in range(4) for b in range(4)]
        columns = [h[1 if cut else 2][a] + h[3][b] for a in range(4) for b in range(4)]
        matrix, _, _ = flatten(seed, cut)
        exact = WeightedSmithOracle([[packed_polynomial(value) for value in row] for row in matrix],
                                    infinity=infinity)
        expected = exact.profile(rows, columns, Fraction(1))[0]
        native = MatrixOracle(matrix, 16)
        observed = sorted(native.profile(rows, columns, scale=1, infinity=infinity)['orders'])
        assert max(abs(x - float(y)) for x, y in zip(observed, sorted(expected))) < 1e-8
        records.append({'cut': cut, 'infinity': infinity,
                        'row_weights': [[x.numerator, x.denominator] for x in rows],
                        'column_weights': [[x.numerator, x.denominator] for x in columns],
                        'complete_orders': [[x.numerator, x.denominator] for x in expected]})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/reflection4-coupled-frame-search.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    report = json.loads(raw)
    assert report['completed'] and report['complete_catalogue_optimization']
    assert report['completed_cases'] == len(report['cases']) == report['planned_coordinate_cases'] == 206
    assert report['admissible_frames'] == 102
    source_dir = Path(__file__).parent
    for name, digest in report['source_sha256'].items():
        assert sha256((source_dir / name).read_bytes()).hexdigest() == digest, name
    base = tensor(4)
    bell_data = bell_normalize(base)
    seeds = []
    counts = {'reconstructed_candidates': 0, 'endpoint_entries': 0,
              'completed_optimizer_loops': 0, 'optimizer_successes': 0}
    expected_keys = {('identity', coordinate) for coordinate in ('original', 'bell')}
    for entry in report['admissibility_catalogue']:
        for pole in entry['allowed_k']:
            for coordinate in ('original', 'bell'):
                expected_keys.add((entry['target'], entry['source'], entry['scalar'], pole, coordinate))
    observed_keys = set()
    for case in report['cases']:
        key = (('identity', case['coordinates']) if case['kind'] == 'identity' else
               (case['target'], case['source'], case['scalar'], case['pole_order'], case['coordinates']))
        assert key not in observed_keys
        observed_keys.add(key)
        seed = rebuild(case, base, bell_data)
        digest = sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest()
        assert digest == case['seed_sha256']
        assert seed['degree'] == case['degree'] == 24
        for index, (a, b) in enumerate(seed['values']):
            x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
            assert (a.bit_count() % 2, b.bit_count() % 2) == (int(x == v and y == u), 0)
            counts['endpoint_entries'] += 1
        counts['reconstructed_candidates'] += 1
        counts['completed_optimizer_loops'] += len(case['optimizations'])
        counts['optimizer_successes'] += sum(bool(x['success']) for x in case['optimizations'])
        seeds.append(seed)
    assert observed_keys == expected_keys
    selections = {}
    winners = []
    for mode in ('balanced', 'asymmetric'):
        for nontrivial in (False, True):
            eligible = [(case[mode]['value'], i) for i, case in enumerate(report['cases'])
                        if mode in case and (not nontrivial or case['kind'] != 'identity')]
            value, index = min(eligible)
            winners.append({'mode': mode, 'nonidentity_only': nontrivial,
                            'case_index': index, 'value': value})
            selections[(mode, index)] = value
    assert min(w['value'] for w in winners if w['mode'] == 'balanced') == report['best_rho']
    assert min(w['value'] for w in winners if w['mode'] == 'asymmetric') == report['best_support_value']
    replays = []
    for (mode, index), expected in selections.items():
        stored = report['cases'][index][mode]
        params = np.asarray(stored['params'], float)
        variables = np.r_[params / math.sqrt(2), math.pi / 4] if mode == 'balanced' else params
        oracle = AsymmetricPressure(seeds[index], 1)
        value, details = oracle.evaluate(variables, gradient=False)
        if mode == 'balanced':
            value /= math.sqrt(2)
        assert abs(value - expected) < 2e-8
        replays.append({'case_index': index, 'mode': mode, 'value': value,
                        'absolute_replay_error': abs(value - expected),
                        'rational_profiles': rational_profiles(seeds[index], variables)})
    counts['numerical_winner_replays'] = len(replays)
    counts['independent_complete_rational_profiles'] = 2 * len(replays)
    result = {'status': 'passed', 'scope': 'complete candidate accounting and independent winner profiles; no global optimization claim',
              'input_sha256': sha256(raw).hexdigest(), 'elapsed_search_seconds': report['elapsed_seconds'],
              'counts': counts, 'case_statuses': dict(Counter(x['status'] for x in report['cases'])),
              'winners': winners, 'replays': replays,
              'source_sha256': {name: sha256((source_dir / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'square_zero_coupled_reflection_frames.py',
                                             'weighted_smith_oracle_f4.py', 'native_matrix_smith_f4.py',
                                             'native_matrix_smith_f4.cpp')}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], **counts, 'winners': winners}))


if __name__ == '__main__':
    main()
