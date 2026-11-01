#!/usr/bin/env python3
"""Archive two completed finite-fiber workers and replay retained witnesses."""
from __future__ import annotations
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from affine_place_profile_oracle import primitive_seed
from asymmetric_crossing_pressure import AsymmetricPressure
from f4_two_cut_smith_fibers import F4FullAlphabetLattice
from reflection_laurent_physical_frames import digest
from search_affine_coupled_reflection_frames import balanced
from search_full_alphabet_two_cut_lattices import root_for_profiles
from search_reflection_physical_bases import change_basis


def read(path):
    raw = path.read_bytes()
    return gzip.decompress(raw) if path.suffix == '.gz' else raw


def bind(report):
    assert report['completed'] and report['elapsed_seconds'] >= report['arguments']['seconds']
    for name, expected in report['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name
    raw = Path(report['arguments']['input']).read_bytes()
    assert sha256(raw).hexdigest() == report['input_sha256']
    return json.loads(raw)


def endpoints(seed):
    assert seed['q_x'] == seed['q_u'] == 4
    for index, (first, second) in enumerate(seed['values']):
        x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
        assert (first.bit_count() % 2, second.bit_count() % 2) == (int(x == v and y == u), 0)
    return len(seed['values'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basis-input', type=Path, default=Path('evidence/finite-fiber-physical-basis-search.json.gz'))
    parser.add_argument('--fiber-input', type=Path, default=Path('evidence/deflated-f4-smith-fiber-search.json.gz'))
    parser.add_argument('--output', type=Path, default=Path('evidence/finite-fiber-search-checkpoint.json'))
    args = parser.parse_args()
    raw = [read(path) for path in (args.basis_input, args.fiber_input)]
    basis, fiber = map(json.loads, raw)
    basis_sources, fiber_sources = bind(basis), bind(fiber)
    assert basis['attempted_cases'] == len(basis['records']) + basis['counters']['duplicate_bases']
    assert basis['completed_trials'] == sum(row['completed'] for row in basis['records'])
    assert fiber['completed_trials'] == sum(fiber['distinct_target_vectors'])
    assert fiber['complete_profile_queries'] == fiber['arguments']['samples'] * fiber['completed_trials']
    assert len(fiber['records']) == sum(fiber['distinct_profile_classes'])
    archives = []
    for name, content in zip(('finite-fiber-physical-basis-search.json.gz', 'deflated-f4-smith-fiber-search.json.gz'), raw):
        compressed = gzip.compress(content, compresslevel=9, mtime=0)
        assert gzip.decompress(compressed) == content
        path = args.output.with_name(name); path.write_bytes(compressed)
        archives.append({'path': str(path), 'sha256': sha256(compressed).hexdigest(),
                         'uncompressed_sha256': sha256(content).hexdigest(),
                         'compressed_bytes': len(compressed), 'uncompressed_bytes': len(content)})
    endpoint_count, channel_count, basis_replays = 0, 0, []
    for key, winner in sorted(basis['winners'].items()):
        source = basis_sources['cases'][winner['source_index']]['seed']
        seed = change_basis(source, winner['basis_x'], winner['basis_y'])
        assert digest(seed['values']) == winner['tensor_sha256'] and seed['degree'] == winner['degree']
        endpoint_count += endpoints(seed)
        engine = AsymmetricPressure(seed)
        parameters = np.asarray(winner['best']['parameters'])
        value, _, details = balanced(engine, parameters) if winner['mode'] == 'balanced' else engine.evaluate(parameters)
        assert abs(value - winner['best']['value']) < 2e-10
        assert np.max(np.abs(np.asarray(details['scaled_orders']) - winner['best']['type_details']['scaled_orders'])) < 1e-9
        channel_count += sum(map(len, details['scaled_orders']))
        basis_replays.append({'key': key, 'tensor_sha256': winner['tensor_sha256'], 'mode': winner['mode'],
                              'value': value, 'degree': seed['degree'], 'scaled_orders': details['scaled_orders']})
    heights = fiber['heights_infinity_then_zero']
    lattices = [F4FullAlphabetLattice(case['seed']['values'], 4, case['seed']['degree'], heights,
                                     fiber['arguments']['upper_window']) for case in fiber_sources['cases']]
    indexed = {(row['trial'], row['source_index'], row['sample']): row for row in fiber['records']}
    h = list(map(int, fiber['arguments']['heights'].split(',')))
    witness_replays = []
    for witness in fiber['retained_witnesses']:
        key = witness['trial'], witness['source_index'], witness['sample']
        assert all(witness[name] == value for name, value in indexed[key].items())
        source = fiber_sources['cases'][witness['source_index']]['seed']
        raw_seed = source | {'values': witness['raw_values']}
        assert digest(raw_seed['values']) == witness['raw_tensor_sha256']
        endpoint_count += endpoints(raw_seed)
        lattices[key[1]].check_witness(raw_seed['values'], witness['targets'])
        seed, common = primitive_seed(raw_seed)
        assert common == {int(k): v for k, v in witness['common_factor'].items()}
        assert digest(seed['values']) == witness['tensor_sha256'] and seed['degree'] == witness['actual_degree']
        endpoint_count += endpoints(seed)
        profiles = lattices[key[1]].profiles(seed['values'])
        assert [list(row) for row in profiles] == witness['profiles']
        cost = root_for_profiles(profiles, [2 * value for value in h], [2 * value for value in h])
        assert abs(cost['rho'] - witness['cost']['rho']) < 2e-10
        channel_count += sum(map(len, profiles))
        witness_replays.append({'trial': key[0], 'source_index': key[1], 'sample': key[2],
                                'tensor_sha256': witness['tensor_sha256'], 'degree': seed['degree'],
                                'profiles': profiles, 'rho': cost['rho']})
    names = (Path(__file__).name, 'search_explicit_physical_bases.py', 'search_iterated_f4_smith_fibers.py',
             'f4_two_cut_smith_fibers.py', 'f4_regular_smith_frames.py', 'search_reflection_physical_bases.py',
             'asymmetric_crossing_pressure.py', 'search_affine_coupled_reflection_frames.py',
             'search_full_alphabet_two_cut_lattices.py')
    report = {'scope': 'Two completed thirty-minute searches; retained witnesses freshly reconstructed and completely reprofiled. Final pressures are floating diagnostics.',
              'archives': archives,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'basis': {'seconds': basis['elapsed_seconds'], 'attempts': basis['attempted_cases'],
                        'retained_trials': len(basis['records']), 'completed_sequences': basis['completed_trials'],
                        'counters': basis['counters'], 'winner_replays': basis_replays},
              'fiber': {'seconds': fiber['elapsed_seconds'], 'attempts': fiber['attempted_trials'],
                        'distinct_targets': fiber['completed_trials'], 'complete_profile_queries': fiber['complete_profile_queries'],
                        'recorded_profile_classes': len(fiber['records']), 'previous_full_replays': fiber['exact_full_replays'],
                        'best_by_source': fiber['best_by_source'], 'retained_witness_replays': witness_replays},
              'fresh_controls': {'endpoints': endpoint_count, 'complete_profile_channels': channel_count,
                                 'full_fiber_witnesses': len(witness_replays), 'basis_winners': len(basis_replays)},
              'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'fresh_controls': report['fresh_controls'],
                      'basis_attempts': basis['attempted_cases'], 'fiber_targets': fiber['completed_trials']}))


if __name__ == '__main__':
    main()
