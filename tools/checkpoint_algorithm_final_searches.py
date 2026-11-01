#!/usr/bin/env python3
"""Portable archive and scoped fresh replay of completed algorithm searches.

The archive preserves original JSON bytes, including historical path strings.
Replays resolve inputs inside the archive and never require those paths.
Numerical optimization records are evidence of a finite search, not exclusions.
"""
from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np

from affine_place_profile_oracle import AffinePlaceOracle, bitplanes, polynomial, primitive_seed
from asymmetric_mixed_place_pressure import ScaledMixedSupportingPressure
from entrywise_reflection_constraints import solve_target
from f16_identity_operator_profiles import INV, evaluate_one, multiply
from finite_jet_smith_fibers import FiniteJetFiber
from finite_precision_smith_frames import build_lattice
from fused_six_vertex_crossings import pscale
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_entrywise_reflection_jets import basis_catalogue, tasks
from search_finite_jet_smith_fibers import templates
from search_fused_physical_bases import ScaledHeightPressure
from search_adaptive_fixed_type_gauges import LocalFixedTypeRates
from hot_code_charge_profiles import natural_profiles, validate_code
from search_hot_downsets import pressure as natural_pressure
from weighted_hot_code_oracle import HotWeightedOracle
from hot_compound_pressure import evaluate as weighted_pressure
from certify_hot_charge_pressure import exact_comparison


REPOSITORY_REPORTS = (
    'entrywise-reflection-jet-search',
    'finite-jet-smith-fiber-search',
    'finite-jet-smith-fiber-infinity4-search',
    'finite-jet-smith-fiber-infinity6-search',
    'affine-mixed-place-search',
    'dense-affine-mixed-place-search',
    'affine-quadratic-mixed-search',
    'entrywise-bell-jet-mixed-search',
    'entrywise-bell-jet-asymmetric-search',
    'finite-fiber-full-mixed-search',
    'finite-fiber-asymmetric-mixed-search',
    'truncated-reflection4-mixed-search',
    'affine-asymmetric-mixed-place-search',
    'reflection4-adaptive-gauge-pilot',
    'reflection16-adaptive-gauge-pilot',
    'truncated-reflection16-pressure-search',
    'reflection16-coupled-frame-search',
    'reflection16-coupled-frame-search-continuation',
    'reflection16-square-zero-search',
    'reflection16-square-zero-search-continuation',
)


def raw_record(path):
    raw = path.read_bytes()
    return {'original_path': str(path), 'sha256': sha256(raw).hexdigest(),
            'bytes': len(raw), 'raw_json': raw.decode()}


def build_archive(path, include_hot):
    files = [Path('evidence') / (name + '.json') for name in REPOSITORY_REPORTS]
    files.append(Path.cwd().parent / 'hot16-mixed-place-pressure-search.json')
    files.append(Path.cwd().parent / 'coherent-hot-refinement.json')
    if include_hot:
        files += [Path('evidence/hot-downset-natural-search.json'),
                  Path('evidence/hot24-weighted-cone-search.json'),
                  Path.cwd().parent / 'hot-charge-larger-pilot.json']
    archive = {'scope': __doc__, 'reports': {}, 'inputs': {}, 'input_aliases': {}}
    for filename in files:
        record = raw_record(filename)
        data = json.loads(record['raw_json'])
        archive['reports'][filename.stem] = record
        for key in ('input', 'catalogue', 'rank_screen'):
            value = data.get('arguments', {}).get(key)
            if value is not None:
                source = Path(value)
                archive['inputs'].setdefault(value, raw_record(source))
    coherent_input = Path.cwd().parent / 'coherent-hot-refinement-input.json'
    archive['inputs'][str(coherent_input)] = raw_record(coherent_input)
    archive['input_aliases']['coherent-hot-refinement'] = str(coherent_input)
    raw = (json.dumps(archive, separators=(',', ':')) + '\n').encode()
    path.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))


def endpoints(seed):
    q = seed['q_x']
    assert seed['q_u'] == q and len(seed['values']) == q ** 4
    for index, value in enumerate(seed['values']):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        assert evaluate_one(polynomial(value)) == int(x == v and y == u)
    return len(seed['values'])


def source_binding(report):
    for name, expected in report.get('source_sha256', {}).items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name


def accounting(name, report):
    keys = ('elapsed_seconds', 'attempted_trials', 'completed_trials', 'completed_targets',
            'feasible_targets', 'unique_seeds', 'distinct_seeds', 'completed_optimizations',
            'unresolved', 'unresolved_pressure_calls', 'completed', 'budget_completed',
            'completed_budget', 'trial_list_completed', 'target_queue_exhausted',
            'complete_template_list', 'completed_plan', 'completed_screens',
            'completed_screen_list', 'completed_refinement_list', 'best_grid_rho')
    result = {key: report[key] for key in keys if key in report}
    result['requested_seconds'] = report.get('arguments', {}).get('seconds')
    if 'records' in report:
        result['records'] = len(report['records'])
    if 'queries' in report and isinstance(report['queries'], list):
        result['queries'] = len(report['queries'])
    if 'iterations' in report and isinstance(report['iterations'], list):
        result['iterations'] = len(report['iterations'])
    if 'continuation' in report:
        result['continuation'] = report['continuation']
    return result


def replay_entrywise(report, counts):
    base = tensor(4)
    queue = list(tasks(basis_catalogue(base)))
    assert report['target_queue_exhausted'] and report['completed_targets'] == len(queue) == len(report['cases'])
    assert report['feasible_targets'] == sum(row['feasible'] for row in report['cases'])
    assert report['unique_seeds'] == len({row['seed_sha256'] for row in report['cases'] if row['feasible']})
    for saved in report['retained_seeds']:
        target = queue[saved['target_index']]
        restored = solve_target(base, target['cap'], target['conditions'], *target['cutoffs'],
                                trace_power=16, protection_mode=target['mode'],
                                prescribed_jets=target['prescribed_jets'])['seed']
        assert digest(restored['values']) == digest(saved['seed']['values']) == saved['seed_sha256']
        assert restored['degree'] == saved['actual_degree']
        counts['entrywise_exact_reconstructions'] += 1
        counts['endpoints'] += endpoints(restored)


def replay_fibers(reports, counts):
    base = tensor(4)
    h, zero = [0, 6, 2, 8], [0] * 4
    lattice = build_lattice(base['values'], 4, 24, [zero, zero, h, h, h, h, zero, zero], 8)
    matrices = basis_catalogue(base)
    matrices['alternating'] = [[int(j == (i ^ 1)) for j in range(4)] for i in range(4)]
    matrices['upper'] = [[int(i == j) ^ (2 if j == i + 1 else 0) for j in range(4)] for i in range(4)]
    queue = list(templates(matrices))
    for report in reports:
        assert report['complete_template_list'] and len(report['records']) == len(queue) == report['completed_targets']
        assert report['feasible_targets'] == sum(row['feasible'] for row in report['records'])
        fiber = FiniteJetFiber(lattice, tuple(report['common_floors']))
        for saved in report['retained_seeds']:
            common = {int(k): value for k, value in saved['common_factor'].items()}
            factor = pscale(common, INV[evaluate_one(common)])
            values = [bitplanes(multiply(polynomial(value), factor)) for value in saved['seed']['values']]
            assert digest(values) == saved['tensor_sha256']
            raw = dict(saved['seed'], values=values, degree=saved['degree_cap'])
            counts['endpoints'] += endpoints(raw) + endpoints(saved['seed'])
            fiber.verify_jets(values, queue[saved['target_index']][1])
            restored_profiles = lattice.check_witness(values, saved['row_targets'])
            assert [list(row) for row in restored_profiles] == saved['two_place_profiles']
            restored, content = primitive_seed(raw)
            assert content == common and restored['degree'] == saved['actual_degree']
            assert restored['values'] == [tuple(value) for value in saved['seed']['values']]
            oracle = AffinePlaceOracle(raw)
            orders = oracle.log_profile_with_gradient(0, (2, 1), [0] * 16)[3]
            assert orders.shape == (16,) and np.array_equal(orders, saved['finite_row_orders'])
            counts['fiber_exact_reconstructions'] += 1
            counts['integer_profile_channels'] += 48


def replay_mixed(report, inputs, counts):
    winners = list(report['winners'].values())
    if not winners:
        return None
    metric = 'support_value' if 'support_value' in winners[0]['best'] else 'rho'
    winner = min(winners, key=lambda row: row['best'][metric])
    competitors = [row for row in report['records'] if row['case_index'] == winner['case_index'] and row.get('best')]
    assert winner['best'][metric] == min(row['best'][metric] for row in competitors)
    assert winner['best'][metric] == min(row['best'][metric] for row in report['records'] if row.get('best'))
    source = inputs[report['arguments']['input']]
    case = source['cases'][winner['case_index']]
    seed, content = primitive_seed(case['seed'])
    assert content == {int(k): value for k, value in winner['common_factor'].items()}
    counts['endpoints'] += endpoints(seed)
    oracle = AffinePlaceOracle(seed)
    cls = ScaledMixedSupportingPressure if metric == 'support_value' else ScaledHeightPressure
    partition = winner.get('partition', winner.get('balanced_input_partition'))
    assert partition == case['partitions'][winner['partition_index']]
    pressure = cls(oracle, partition, 'full', False, True, float('inf'))
    pressure.evaluate(np.asarray(winner['best']['optimizer_parameters']))
    restored = pressure.best
    assert abs(restored[metric] - winner['best'][metric]) < 3e-9
    assert len(winner['best']['place_records']) == len(restored['place_records']) == 2
    for cut, (saved_places, new_places) in enumerate(zip(winner['best']['place_records'], restored['place_records'])):
        assert len(saved_places) == len(new_places) == len(pressure.groups[cut])
        for saved, new, group in zip(saved_places, new_places, pressure.groups[cut]):
            assert saved['point'] == new['point'] == (list(group['point']) if group['point'] is not None else None)
            assert saved['conjugate_labels'] == new['conjugate_labels'] == group['labels']
            order_key = 'complete_orders' if 'complete_orders' in saved else 'restored_orders'
            a, b = np.asarray(saved[order_key]), np.asarray(new[order_key])
            assert a.shape == b.shape == (oracle.alphabet ** 2,)
            assert np.max(np.abs(a - b)) < 3e-8
            counts['numerical_profile_channels'] += len(a)
    counts['mixed_winner_replays'] += 1
    return {'label': winner['label'], 'metric': metric, 'value': restored[metric],
            'degree': seed['degree'], 'tensor_sha256': digest(seed['values'])}


def replay_hot(reports, summaries, counts):
    if 'hot-downset-natural-search' in reports:
        report = reports['hot-downset-natural-search']
        assert report['completed'] and report['elapsed_seconds'] >= report['arguments']['seconds']
        assert report['completed_trials'] == sum(row['completed'] for row in report['records'])
        winner = report['records'][report['best_record_index']]
        assert winner['pressure']['rho'] == min(row['pressure']['rho'] for row in report['records'] if row['completed'])
        for row in report['records']:
            if row['completed']:
                validate_code(row['letters'], report['arguments']['modes'])
                assert len(row['orders_w']) == len(row['letters']) ** 2
        profile = natural_profiles(winner['letters'])
        assert profile['orders_w'] == winner['orders_w']
        restored = natural_pressure(profile)
        assert abs(restored['rho'] - winner['pressure']['rho']) < 1e-10
        summaries['hot-downset-natural-search']['best_replay'] = restored
        counts['integer_profile_channels'] += len(profile['orders_w'])
    if 'hot24-weighted-cone-search' in reports:
        report = reports['hot24-weighted-cone-search']
        assert report['completed'] and report['budget_completed']
        winner = report['best']
        assert winner['rho'] == min(row['rho'] for row in report['queries'])
        oracle = HotWeightedOracle(report['letters'])
        parameters = np.asarray(winner['rational_parameters'], dtype=np.int64)
        denominator = winner['rational_denominator']
        heights = np.zeros(sum(oracle.dimensions))
        heights[oracle.active_heights] = parameters[:-1] / denominator
        scale = parameters[-1] / denominator
        matrix, bounds = oracle.cone(1.)
        assert max(abs(parameters)) < 2 ** 50
        residue = matrix[:, oracle.active_heights].astype(np.int64) @ parameters[:-1] - bounds.astype(np.int64) * parameters[-1]
        assert int(np.min(residue)) == winner['minimum_cone_numerator'] >= 0
        assert all(row['minimum_cone_numerator'] >= 0 for row in report['queries'])
        rho, _, _, restored = weighted_pressure(oracle, heights, scale)
        assert abs(rho - winner['rho']) < 1e-10
        assert len(winner['orders']) == len(restored['orders']) == oracle.qx * oracle.qy
        assert np.max(np.abs(np.asarray(winner['orders']) - restored['orders'])) < 1e-8
        summaries['hot24-weighted-cone-search']['best_replay'] = {'rho': rho, 'scale': scale}
        counts['numerical_profile_channels'] += len(restored['orders'])
    if 'hot-charge-larger-pilot' in reports:
        report = reports['hot-charge-larger-pilot']
        rows = []
        certified = json.loads(Path('evidence/hot64-charge-rational-pressure.json').read_text())
        for saved in report['pilots']:
            profile = saved['profile']
            q = saved['size']
            validate_code(profile['x_letters'], profile['modes'])
            assert profile['x_letters'] == profile['y_letters'] and len(profile['orders_w']) == q * q
            assert sum(block['size'] for block in profile['blocks']) == q * q
            restored = natural_pressure(profile)
            assert abs(restored['rho'] - saved['rho']) < 1e-10
            if q == 64:
                assert profile['x_letters'] == certified['profile']['x_letters']
                assert profile['y_letters'] == certified['profile']['y_letters']
                assert profile['orders_w'] == certified['profile']['orders_w']
            if q == 96:
                frozen = json.loads(Path('evidence/hot96-charge-profile.json').read_text())
                assert profile == frozen
                exact = json.loads(Path('evidence/hot96-charge-rational-pressure.json').read_text())
                source_binding(exact)
                assert exact['profile'] == frozen
                for key in ('certificate', 'weaker_certificate'):
                    saved_certificate = exact[key]
                    replayed = exact_comparison(profile, *saved_certificate['rational_base'],
                                                *saved_certificate['support_parameter'])
                    assert replayed == saved_certificate and replayed['passed']
                    counts['rational_pressure_comparisons'] += 1
            rows.append({'size': q, 'rho': restored['rho'], 'scale': restored['scale'],
                         'complete_channels': q * q, 'profile_seconds': saved['seconds'],
                         'fresh_profile_recomputation': False})
        summaries['hot-charge-larger-pilot']['completed_sizes'] = rows
        summaries['hot-charge-larger-pilot']['scope'] = ('Complete integer profiles from the original audited constructor run. '
            'This checkpoint repeats their numerical pressure and size checks, not their large native eliminations; '
            'the separate q64 rational certificate has an independent full replay.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('evidence/algorithm-final-searches.json.gz'))
    parser.add_argument('--output', type=Path, default=Path('evidence/algorithm-final-searches-checkpoint.json'))
    parser.add_argument('--build', action='store_true', help='Freeze original reports and inputs from their historical paths.')
    parser.add_argument('--include-hot', action='store_true', help='Include the final hot searches when freezing.')
    args = parser.parse_args()
    if args.build:
        build_archive(args.archive, args.include_hot)
    began = time.monotonic()
    archive = json.loads(gzip.decompress(args.archive.read_bytes()))
    reports, inputs = {}, {}
    for kind, destination in (('reports', reports), ('inputs', inputs)):
        for name, record in archive[kind].items():
            raw = record['raw_json'].encode()
            assert len(raw) == record['bytes'] and sha256(raw).hexdigest() == record['sha256']
            destination[name] = json.loads(raw)
    counts = {'endpoints': 0, 'entrywise_exact_reconstructions': 0, 'fiber_exact_reconstructions': 0,
              'integer_profile_channels': 0, 'numerical_profile_channels': 0, 'mixed_winner_replays': 0,
              'adaptive_catalogue_profiles': 0, 'rational_pressure_comparisons': 0}
    summaries = {}
    for name, report in reports.items():
        source_binding(report)
        if 'input_sha256' in report:
            input_path = report['arguments'].get('input', archive.get('input_aliases', {}).get(name))
            assert archive['inputs'][input_path]['sha256'] == report['input_sha256']
        for argument, hash_key in (('catalogue', 'catalogue_sha256'), ('rank_screen', 'rank_screen_sha256')):
            if hash_key in report:
                assert archive['inputs'][report['arguments'][argument]]['sha256'] == report[hash_key]
        summaries[name] = accounting(name, report)
        if 'attempted_trials' in report and 'records' in report:
            assert report['attempted_trials'] == len(report['records'])
        if 'winners' in report and 'input_sha256' in report and name != 'coherent-hot-refinement':
            assert report['completed_trials'] == sum(row['completed'] for row in report['records'])
            assert report['budget_completed'] or report.get('trial_list_completed')
            summaries[name]['best_replay'] = replay_mixed(report, inputs, counts)
        if name == 'coherent-hot-refinement':
            assert report['completed'] and len(report['records']) == 6
            summaries[name]['completed_optimizer_sequences'] = sum(row['completed'] for row in report['records'])
            summaries[name]['scope'] = 'Provenance and counters only; no fresh optimization or winner replay.'
        if 'continuation' in report:
            lineage = report['continuation']
            assert sha256(Path(__file__).with_name(lineage['original_worker']).read_bytes()).hexdigest() == lineage['original_worker_sha256']
            assert sha256(Path(__file__).with_name('continue_interrupted_worker.py').read_bytes()).hexdigest() == lineage['continuation_harness_sha256']
            prefix_name = Path(lineage['prefix_path']).stem
            assert archive['reports'][prefix_name]['sha256'] == lineage['prefix_sha256']
            previous = reports[prefix_name]
            assert previous['elapsed_seconds'] == lineage['previous_checkpoint_seconds']
            assert report['elapsed_seconds'] >= lineage['requested_continuation_seconds']
            assert abs(lineage['requested_total_seconds'] - lineage['previous_checkpoint_seconds']
                       - lineage['requested_continuation_seconds']) < 1e-8
            summaries[name]['combined_recorded_seconds'] = previous['elapsed_seconds'] + report['elapsed_seconds']
    replay_entrywise(reports['entrywise-reflection-jet-search'], counts)
    replay_fibers([reports[name] for name in REPOSITORY_REPORTS if name.startswith('finite-jet-smith-fiber')], counts)
    for name in ('reflection4-adaptive-gauge-pilot', 'reflection16-adaptive-gauge-pilot'):
        report = reports[name]
        oracle = LocalFixedTypeRates(tensor(report['arguments']['interval']), report['type_x'], report['type_y'])
        assert report['completed_trials'] == sum(row['optimization_completed'] for row in report['records'])
        for saved in (profile for side in report['profiles'] for profile in side):
            restored = oracle.catalogue_record(saved['side'], saved['parameters'])
            a, b = np.asarray(restored['orders']), np.asarray(saved['orders'])
            assert a.shape == b.shape == (oracle.n,) and np.max(np.abs(a - b)) < 3e-8
            assert abs(restored['mean_shift'] - saved['mean_shift']) < 3e-8
            counts['adaptive_catalogue_profiles'] += 1
            counts['numerical_profile_channels'] += oracle.n
    replay_hot(reports, summaries, counts)
    result = {'scope': 'Original reports and inputs are byte-preserved. Exact retained constructions and complete integer profiles are replayed where counted; floating optimized costs are diagnostics only.',
              'archive': {'path': str(args.archive), 'sha256': sha256(args.archive.read_bytes()).hexdigest(),
                          'bytes': args.archive.stat().st_size},
              'source_sha256': {Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()},
              'summaries': summaries, 'fresh_controls': counts, 'replay_seconds': time.monotonic() - began,
              'passed': True}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': True, 'archive_bytes': args.archive.stat().st_size, 'fresh_controls': counts}))


if __name__ == '__main__':
    main()
