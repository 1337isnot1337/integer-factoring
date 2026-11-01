#!/usr/bin/env python3
"""Compact frozen root-worker logs and replay retained mathematical witnesses.

The physical-family archive is a summary, not a byte-for-byte copy of each
optimizer log. It retains every completed candidate's digest and scores,
the full winning records, and original-byte/source/continuation provenance.
The smaller arithmetic-target logs are retained in full. Numerical pressure
replays do not certify a rate, global optimum, or factoring improvement.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import gzip
from hashlib import sha256
import importlib.util
import inspect
import json
import math
from pathlib import Path

import numpy as np

FAMILIES = {
    'basis4': 'reflection-physical-basis-search',
    'asymmetric4': 'reflection-asymmetric-basis-search',
    'basis16': 'reflection16-local-basis-search',
    'polynomial4': 'physical-polynomial-frame-search',
    'squarezero4': 'reflection4-square-zero-frame-search',
    'affine4': 'affine-coupled-frame-search',
}
EXTRAS = (
    'reflection-independent-charge-search',
    'dense-affine-rank-target-search',
    'affine-hasse-rank-target-search',
    'dense-affine-pairing-targets',
    'reflection16-filtered-laurent-series-search',
    'truncated-degree16-affine-lattice-search',
    'truncated-degree22-affine-lattice-search',
)
ROOT = Path(__file__).resolve().parent.parent


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def digest(value):
    return sha256(canonical(value)).hexdigest()


def tensor_digest(values):
    return sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()


def score_records(row):
    """Return the named objective, stored value, and exact saved call record."""
    result = []
    for key in ('best', 'balanced', 'asymmetric'):
        value = row.get(key)
        if not isinstance(value, dict):
            continue
        if key == 'balanced' or (key == 'best' and row.get('phase') == 'balanced') or 'rho' in value:
            metric, number = 'balanced_rho', value.get('rho', value.get('value'))
        else:
            metric, number = 'asymmetric_support', value.get('support_value', value.get('value'))
        if number is not None:
            result.append((metric, float(number), key, value))
    return result


def compact_run(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    omit = {'cases', 'screen_records', 'admissibility'}
    header = {k: v for k, v in data.items() if k not in omit}
    ledger, winners = [], {}
    for index, row in enumerate(data['cases']):
        scores = score_records(row)
        ledger.append([index, digest(row), row.get('status', 'saved completed candidate'),
                       {metric: value for metric, value, _, _ in scores}])
        for metric, value, field, _ in scores:
            if metric not in winners or value < winners[metric]['score']:
                winners[metric] = {'index': index, 'score': value, 'field': field, 'record': row}
    omitted = {}
    for key in sorted(omit - {'cases'}):
        if key in data:
            omitted[key] = {'records': len(data[key]), 'canonical_sha256': digest(data[key])}
    return {'path': str(path.relative_to(ROOT)), 'original_sha256': sha256(raw).hexdigest(),
            'original_bytes': len(raw), 'header': header,
            'ledger_fields': ['row_index', 'canonical_record_sha256', 'terminal_status', 'scores'],
            'candidate_ledger': ledger, 'winning_records': winners, 'omitted_arrays': omitted}


def create_archive(path):
    physical = {kind: [compact_run(ROOT / 'evidence' / (name + suffix + '.json'))
                       for suffix in ('', '-continuation')]
                for kind, name in FAMILIES.items()}
    extras = {}
    for name in EXTRAS:
        source = ROOT / 'evidence' / (name + '.json')
        raw = source.read_bytes()
        extras[name] = {'path': str(source.relative_to(ROOT)), 'original_sha256': sha256(raw).hexdigest(),
                        'original_bytes': len(raw), 'data': json.loads(raw)}
    # The polynomial-span input already has a separately owned complete archive.
    span_hash = extras['affine-hasse-rank-target-search']['data']['span_sha256']
    matches = [p for p in (ROOT / 'evidence').glob('*.json.gz')
               if sha256(p.read_bytes()).hexdigest() == span_hash]
    assert len(matches) == 1, 'The exact Hasse-span input must be available.'
    inputs = [matches[0], ROOT / 'evidence/dense-affine-scalar-content.json',
              ROOT / 'evidence/dense-rank-one-coupled-frames.json']
    data = {'schema_version': 1,
            'scope': 'Compact physical-search summary plus complete smaller root-worker logs. '
                     'Original optimizer traces are not embedded; byte hashes bind their local originals. '
                     'All completed physical candidate score records and full family winners are retained.',
            'physical_families': physical, 'other_runs': extras,
            'input_bindings': [{'path': str(p.relative_to(ROOT)), 'sha256': sha256(p.read_bytes()).hexdigest(),
                                'bytes': p.stat().st_size} for p in inputs]}
    path.write_bytes(gzip.compress(canonical(data) + b'\n', mtime=0))


def source_checks(archive):
    sources = {}
    for runs in archive['physical_families'].values():
        for run in runs:
            for name, expected in run['header']['source_sha256'].items():
                assert sources.setdefault(name, expected) == expected
    for run in archive['other_runs'].values():
        for name, expected in run['data']['source_sha256'].items():
            assert sources.setdefault(name, expected) == expected
    sources['continue_physical_frame_search.py'] = next(iter(archive['physical_families'].values()))[1]['header']['continuation']['harness_sha256']
    for name, expected in sources.items():
        assert sha256((ROOT / 'tools' / name).read_bytes()).hexdigest() == expected, name
    for item in archive['input_bindings']:
        assert sha256((ROOT / item['path']).read_bytes()).hexdigest() == item['sha256'], item['path']
    return sources


def check_continuations(archive):
    reports = []
    for kind, (prefix, continuation) in archive['physical_families'].items():
        old, new = prefix['header'], continuation['header']
        meta = new['continuation']
        assert sha256((ROOT/'tools/continue_physical_frame_search.py').read_bytes()).hexdigest() == meta['harness_sha256']
        assert old['completed'] is False and new['completed'] is True
        assert meta['prefix_sha256'] == prefix['original_sha256']
        assert meta['previous_saved_cases'] == len(prefix['candidate_ledger'])
        assert meta['previous_checkpoint_seconds'] == old['elapsed_seconds']
        assert meta['requested_total_seconds'] == old['requested_seconds']
        assert abs(meta['requested_continuation_seconds'] + old['elapsed_seconds'] - old['requested_seconds']) < 1e-7
        assert abs(new['combined_checkpoint_seconds'] - old['elapsed_seconds'] - new['elapsed_seconds']) < 1e-7
        assert new['combined_saved_cases'] == sum(len(run['candidate_ledger']) for run in (prefix, continuation))
        assert old['source_sha256'] == new['source_sha256']
        worker = ROOT / 'tools' / meta['original_worker']
        assert sha256(worker.read_bytes()).hexdigest() == meta['original_worker_sha256']
        spec = importlib.util.spec_from_file_location('_provenance_' + kind, worker)
        module = importlib.util.module_from_spec(spec)
        assert spec and spec.loader
        spec.loader.exec_module(module)
        instrumented = inspect.getsource(module.main)
        for change in meta['initialization_changes']:
            assert instrumented.count(change['old']) == 1
            instrumented = instrumented.replace(change['old'], change['new'])
        assert sha256(instrumented.encode()).hexdigest() == meta['instrumented_main_sha256']
        winners = {}
        for run in (prefix, continuation):
            for metric, candidate in run['winning_records'].items():
                index = candidate['index']
                assert digest(candidate['record']) == run['candidate_ledger'][index][1]
                assert candidate['score'] == min(row[3][metric] for row in run['candidate_ledger'] if metric in row[3])
                if metric not in winners or candidate['score'] < winners[metric]['score']:
                    winners[metric] = dict(candidate, source_path=run['path'])
        reports.append({'kind': kind, 'prefix_cases': len(prefix['candidate_ledger']),
                        'continuation_cases': len(continuation['candidate_ledger']),
                        'combined_saved_cases': new['combined_saved_cases'],
                        'requested_seconds': old['requested_seconds'],
                        'saved_elapsed_seconds': new['combined_checkpoint_seconds'],
                        'stopping_reason': new['stopping_reason'],
                        'metadata_best_rho': new.get('best_rho'),
                        'metadata_best_support': new.get('best_support_value'),
                        'winning_scores': {k: v['score'] for k, v in winners.items()},
                        'winners': winners,
                        'prefix_sha256': prefix['original_sha256'],
                        'continuation_sha256': continuation['original_sha256'],
                        'random_stream': meta['random_stream'],
                        'counters_prefix': old.get('counters', {}),
                        'counters_continuation': new.get('counters', {})})
    return reports


@lru_cache(maxsize=None)
def base_tensor(interval):
    from reflection_dual_rail_gaussian import tensor
    return tensor(interval)


@lru_cache(maxsize=1)
def affine_catalogue():
    from search_affine_coupled_reflection_frames import catalogue
    return catalogue(base_tensor(4))


def construct_physical(kind, row):
    if kind in ('basis4', 'asymmetric4', 'basis16'):
        from search_reflection_physical_bases import change_basis
        return change_basis(base_tensor(8 if kind == 'basis16' else 4), row['basis_x'], row['basis_y'])
    if kind == 'polynomial4':
        from physical_polynomial_frames import apply_word
        return apply_word(base_tensor(4), row['word'])
    if kind == 'squarezero4':
        from square_zero_coupled_reflection_frames import construct
        return construct(base_tensor(4), row['generator'], pole=row['pole'])
    if kind == 'affine4':
        from search_affine_coupled_reflection_frames import affine_seed
        return affine_seed(base_tensor(4), affine_catalogue(), row['word'])
    raise ValueError(kind)


def check_endpoint(seed):
    a, b = seed['q_x'], seed['q_u']
    assert len(seed['values']) == a * a * b * b
    for index, (lo, hi) in enumerate(seed['values']):
        x, rest = index % a, index // a
        y, rest = rest % b, rest // b
        u, v = rest % b, rest // b
        assert (lo.bit_count() % 2, hi.bit_count() % 2) == (int(x == v and y == u), 0)
    return len(seed['values'])


def replay_call(seed, record, field, metric, kind=''):
    from asymmetric_crossing_pressure import AsymmetricPressure
    from native_rectangular_smith_f4 import RectangularOracle
    value = record[field]
    cut = record.get('cut_at_infinity', record.get('cut', value.get('cut', 1)))
    if kind == 'basis4':
        oracle = RectangularOracle(seed['values'], seed['q_x'], seed['q_u'])
        rho, _, _, orders, _, _ = oracle.types(value['params'], cut)
        score, expected, stored = rho, value['rho'], value['scaled_orders']
    else:
        params = np.asarray(value.get('params', value.get('variables')), dtype=float)
        oracle = AsymmetricPressure(seed, cut)
        if metric == 'balanced_rho':
            score, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4], gradient=False)
            score /= math.sqrt(2)
        else:
            score, details = oracle.evaluate(params, gradient=False)
        expected = value.get('rho', value.get('support_value', value.get('value')))
        orders = details['scaled_orders']
        stored = value.get('details', value)['scaled_orders']
    observed_array, stored_array = np.asarray(orders), np.asarray(stored)
    assert observed_array.shape == stored_array.shape
    assert observed_array.size == 2 * seed['q_x'] * seed['q_u']
    error = max(abs(float(a) - float(b)) for a, b in zip(observed_array.flat, stored_array.flat))
    assert error < 2e-7 and abs(score - expected) < 2e-8, (kind, score, expected, error)
    return {'metric': metric, 'recorded': expected, 'recomputed': score,
            'complete_profile_orders': np.asarray(orders).size,
            'maximum_saved_order_difference': error, 'orders': orders,
            'tensor_sha256': tensor_digest(seed['values']), 'endpoint_entries': check_endpoint(seed)}


def replay_physical(summaries):
    controls = []
    for summary in summaries:
        kind = summary['kind']
        for metric, candidate in summary['winners'].items():
            row = candidate['record']
            seed = construct_physical(kind, row)
            if 'seed_sha256' in row:
                assert tensor_digest(seed['values']) == row['seed_sha256']
            replay = replay_call(seed, row, candidate['field'], metric, kind)
            replay.update(family=kind, source_path=candidate['source_path'], row_index=candidate['index'])
            controls.append(replay)
    return controls


def replay_extras(archive):
    from asymmetric_crossing_pressure import AsymmetricPressure
    from native_rectangular_smith_f4 import RectangularOracle
    from search_physical_code_gauges import start_params
    from screen_affine_finite_places import bitplane_rank, rows
    runs = {name: item['data'] for name, item in archive['other_runs'].items()}
    result = {}
    # The independent character masks are finite and exact; rerun all 3^6.
    from itertools import product
    from search_reflection_independent_charges import masked
    data = runs['reflection-independent-charge-search']
    degrees, seen = Counter(), set()
    hashes = {row['seed_sha256'] for row in data['screen_records']}
    for tail in product(range(3), repeat=6):
        seed = masked(base_tensor(4), (0,) + tail[:3], (0,) + tail[3:])
        degrees[seed['degree']] += 1
        seen.add(tensor_digest(seed['values']))
    assert seen == hashes and len(seen) == data['unique_tensors']
    assert {str(k): v for k, v in degrees.items()} == data['degree_histogram']
    calls = []
    for field, metric in (('balanced', 'balanced_rho'), ('asymmetric', 'asymmetric_support')):
        candidate = min((r for r in data['cases'] if field in r), key=lambda r: r[field]['value'])
        seed = masked(base_tensor(4), candidate['first_charge'], candidate['second_charge'])
        calls.append(replay_call(seed, candidate, field, metric))
    result['independent_charges'] = {'exact_masks': 729, 'unique_tensors': len(seen),
                                     'complete_character_entries': 729 * 256, 'winners': calls}
    # Retained dense witnesses are checked against every imposed rectangle.
    data = runs['dense-affine-rank-target-search']
    assert len(data['cases']) == data['planned_targets'] == data['completed_targets']
    assert sum(r['feasible'] for r in data['cases']) == data['feasible_targets']
    for seed in data['retained_seeds']:
        assert tensor_digest(seed['values']) == seed['seed_sha256']
        check_endpoint(seed)
        finite = [[bitplane_rank(rows(seed['values'], cut, point)) for cut in (0, 1)] for point in (2, 3)]
        assert finite == seed['finite_ranks']
        record = data['cases'][seed['target_index']]
        assert record['seed_sha256'] == seed['seed_sha256'] and record['feasible']
        from screen_affine_finite_places import evaluate
        if record['preserve_zero']:
            assert [evaluate(v, 0) for v in seed['values']] == [evaluate(v, 0) for v in base_tensor(4)['values']]
        for point, cut, first, second in record['constraints']:
            keep_rows, keep_columns = {5*x for x in first}, {5*x for x in second}
            for index, value in enumerate(seed['values']):
                x, y, u, v = ((index // 4**j) % 4 for j in range(4))
                i, j = (4*x+y, 4*u+v) if cut == 0 else (4*x+u, 4*y+v)
                if i not in keep_rows and j not in keep_columns:
                    assert evaluate(value, point) == 0
    best = min(data['retained_seeds'], key=lambda r: r['natural_rho'] if r['natural_rho'] is not None else math.inf)
    natural = start_params(np.asarray([0, 12, 4, 16]) * .2326459924883792, .2326459924883792)
    oracle = AsymmetricPressure(best, 1)
    scores = [oracle.evaluate(np.r_[natural*factor/math.sqrt(2), math.pi/4], gradient=False)[0]/math.sqrt(2)
              for factor in (1, .5, .2)]
    assert abs(min(scores)-best['natural_rho']) < 2e-8
    result['dense_rank_targets'] = {'planned_and_completed': len(data['cases']),
        'feasible_targets': data['feasible_targets'], 'unique_seeds': data['unique_seeds'],
        'retained_witnesses_replayed': len(data['retained_seeds']),
        'finite_ranks_replayed': 4*len(data['retained_seeds']),
        'winner_target': best['target_index'], 'winner_rho': min(scores)}
    # The Hasse search exhausted its time, not its finite target list.
    data = runs['affine-hasse-rank-target-search']
    assert len(data['cases']) == data['completed_targets'] < data['planned_targets']
    assert not any(r['feasible'] for r in data['cases']) and not data['retained_seeds']
    from search_affine_hasse_rank_targets import expanded, packed, protection, rectangle, add_column
    from solve_dense_affine_scalar_content import Span
    span_path = next(ROOT/item['path'] for item in archive['input_bindings'] if item['sha256'] == data['span_sha256'])
    source = json.loads(gzip.decompress(span_path.read_bytes()))
    basis = [tuple(int(v, 16) for v in row['bitplanes']) for row in source['basis']]
    transformed = [expanded(v) for v in basis]
    target = expanded(packed(base_tensor(4)['values']))
    for row in (data['cases'][0], data['cases'][-1]):
        fixed, finite = protection(*row['protected_cutoffs']), 0
        for condition in row['conditions']:
            finite |= rectangle(*condition)
        mask = fixed | finite
        system = object.__new__(Span)
        system.pivots, system.increments = {}, []
        for j, column in enumerate(transformed):
            add_column(system, (column[0]&mask, column[1]&mask), j)
        residual, _ = system.reduce((target[0]&finite, target[1]&finite))
        assert [hex(v) for v in residual] == row['residual_bitplanes']
        assert len(system.pivots) == row['constraint_rank'] and any(residual)
    result['hasse_targets'] = {'completed_targets': len(data['cases']), 'planned_targets': data['planned_targets'],
        'feasible_recorded': 0, 'fresh_first_and_last_negative_systems': 2,
        'queue_exhausted': False, 'span_dimension': len(basis)}
    # Full natural-profile replay for the 64-element filtered Laurent family.
    from search_filtered_laurent_series16 import filtered_basis, seed_from_parameters
    from reflection_laurent_physical_frames import bell_normalize
    from reflection_dual_rail_gaussian import physical_bit_energies
    data = runs['reflection16-filtered-laurent-series-search']
    base = base_tensor(8)
    bell = bell_normalize(base)
    weights = physical_bit_energies(8)
    energy = np.asarray([sum(v for j, v in enumerate(weights) if i>>j&1) for i in range(16)])
    spaces = [filtered_basis(bell[1]['first_factor'], bell[1]['inverse_factor'], list(map(int, energy/2)), p)
              for p in range(1, 5)]
    assert spaces == data['linear_spaces']
    assert len(data['screen']) == 4**sum(r['dimension'] for r in spaces) == 64
    baseline_orders = RectangularOracle(base['values'], 16, 16).prefixes(start_params(energy, 1), 1)['orders']
    assert baseline_orders == data['original_natural_profiles']
    for row in data['screen']:
        seed = seed_from_parameters(base, bell, spaces, row['coefficients'])
        assert tensor_digest(seed['values']) == row['seed_sha256']
        check_endpoint(seed)
        orders = RectangularOracle(seed['values'], 16, 16).prefixes(start_params(energy, 1), 1)['orders']
        assert orders == row['natural_profiles']
        assert orders[:256] == baseline_orders[:256]
    best = min((r for r in data['refinements'] if 'best' in r), key=lambda r: r['best']['rho'])
    seed = seed_from_parameters(base, bell, spaces, best['coefficients'])
    result['filtered_laurent16'] = {'frames': 64, 'natural_orders': 64*512,
        'unchanged_infinity_orders': 64*256, 'field_dimensions': [r['dimension'] for r in spaces],
        'refinements': len(data['refinements']), 'queue_exhausted': data['refinement_queue_exhausted'],
        'winner': replay_call(seed, best, 'best', 'balanced_rho')}
    # Rebuild the two regular Smith frames and replay the best full tensor.
    from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice
    from truncate_reflection_crossings import truncate_even_crossing
    from search_full_alphabet_two_cut_lattices import root_for_profiles
    for degree in (16, 22):
        data = runs[f'truncated-degree{degree}-affine-lattice-search']
        base = truncate_even_crossing(base_tensor(4), degree)
        assert tensor_digest(base['values']) == data['source_tensor_sha256']
        lattice = FullAlphabetTwoCutLattice(base['values'], 4, degree, data['heights_infinity_then_zero'],
                                          raise_cap=data['arguments']['upper_window'])
        best = min(data['retained_witnesses'], key=lambda r: r['rho'])
        assert tensor_digest(best['values']) == best['tensor_sha256']
        observed = lattice.check_witness(best['values'], best['targets'])
        assert [list(v) for v in observed] == best['profiles']
        cost = root_for_profiles(observed, [0,12,4,16], [0,12,4,16])
        assert abs(cost['rho']-best['rho']) < 2e-10
        assert sum(r['samples'] for r in data['classes']) == data['complete_profile_queries']
        result[f'truncated_degree{degree}'] = {'completed_trials': data['completed_trials'],
            'planned_trials': data['arguments']['trials'], 'completed_requested_trials': data['completed_requested_trials'],
            'distinct_targets': data['distinct_target_vectors'], 'profile_classes': len(data['classes']),
            'complete_profile_queries': data['complete_profile_queries'], 'winner_trial': best['trial'],
            'winner_rho': cost['rho'], 'exact_integer_orders': best['profiles'],
            'tensor_sha256': best['tensor_sha256'], 'endpoint_entries': 256}
    pairing = runs['dense-affine-pairing-targets']
    assert len(pairing['records']) == 18 and not any(r['consistent'] for r in pairing['records'])
    result['dense_pairing_targets'] = {'archived_exact_negative_records': 18,
        'scope': 'Source-bound original records retained; these 18 systems are not independently rerun by this checker.'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=ROOT/'evidence/root-worker-search-archive.json.gz')
    parser.add_argument('--output', type=Path, default=ROOT/'evidence/root-worker-search-replay.json')
    parser.add_argument('--build-archive', action='store_true', help='Read the preserved local raw logs and replace the compact archive.')
    args = parser.parse_args()
    if args.build_archive:
        create_archive(args.archive)
    raw = args.archive.read_bytes()
    archive = json.loads(gzip.decompress(raw))
    sources = source_checks(archive)
    summaries = check_continuations(archive)
    physical = replay_physical(summaries)
    extras = replay_extras(archive)
    result = {'schema_version': 1,
        'scope': 'Source-bound provenance and fresh exact coefficient/endpoint controls. Complete native weighted profiles '
                 'use exact F4 coefficients but floating valuation comparisons and numerical pressure; no threshold certificate.',
        'archive_sha256': sha256(raw).hexdigest(), 'archive_bytes': len(raw),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256': sources, 'input_bindings': archive['input_bindings'],
        'physical_families': [{k:v for k,v in row.items() if k != 'winners'} for row in summaries],
        'physical_winner_replays': physical, 'other_replays': extras,
        'score_claim': 'No new complete factoring exponent; no numerical winner crosses its claimed threshold.'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'archive_bytes': len(raw), 'source_files': len(sources),
                      'physical_families': len(summaries),
                      'saved_physical_candidates': sum(r['combined_saved_cases'] for r in summaries),
                      'fresh_physical_winners': len(physical), 'other_replays': list(extras)}))


if __name__ == '__main__':
    main()
