#!/usr/bin/env python3
"""Bind the completed dense pressure archive to the published exact catalogue."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def digest_file(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('evidence/dense-rank-one-pressure-search.json'))
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    archive = json.loads(args.archive.read_text())
    catalogue = json.loads(args.catalogue.read_text())
    assert archive['completed'] and archive['complete_screen']
    assert archive['published_catalogue_sha256'] == digest_file(args.catalogue)
    records_digest = sha256(json.dumps(catalogue['records'], separators=(',', ':')).encode()).hexdigest()
    assert archive['catalogue_records_sha256'] == records_digest
    # The old complete generator/candidate records are bound independently
    # of the old header, which tested only through12 and had its own timing.
    assert archive['catalogue_sha256'] == '1ecf137146da4013ab8aedf7c7278f97e81e82d193e0596e1f63ed13499cbd58'
    assert archive['catalogue_sha256'] != archive['published_catalogue_sha256']
    source_checks = {}
    for name, recorded in archive['source_sha256'].items():
        current = digest_file(Path(__file__).with_name(name))
        assert current == recorded, (name, recorded, current)
        source_checks[name] = current
    candidates = {(row['index'], candidate['pole']): candidate
                  for row in catalogue['records'] for candidate in row['candidates']}
    seen = set()
    for row in archive['screen']:
        key = row['generator_index'], row['pole']
        assert key in candidates and key not in seen
        seen.add(key)
        assert row['seed_sha256'] == candidates[key]['seed_sha256']
        assert 0 <= row['profile_index'] < len(archive['profiles'])
    assert seen == set(candidates)
    for row in archive['refinements']:
        exact = catalogue['records'][row['generator_index']]
        assert row['u'] == exact['u'] and row['v'] == exact['v']
        assert row['seed_sha256'] == candidates[row['generator_index'], row['pole']]['seed_sha256']

    from asymmetric_crossing_pressure import AsymmetricPressure
    from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms, digest
    from native_rectangular_smith_f4 import RectangularOracle
    from reflection_dual_rail_gaussian import tensor
    from reflection_laurent_physical_frames import bell_normalize

    base = tensor(4)
    _, factors = bell_normalize(base)
    replays = []
    for mode in ('balanced', 'asymmetric'):
        best = min((row for row in archive['refinements'] if mode in row),
                   key=lambda row: row[mode]['value'])
        terms = rank_one_terms(base, best['u'], best['v'], factors)
        changed = transformed_from_terms(base, terms, best['pole'])
        assert digest(changed['values']) == best['seed_sha256']
        params = best[mode]['params']
        if mode == 'balanced':
            value = RectangularOracle(changed['values'], 4, 4).types(params)[0]
        else:
            value = AsymmetricPressure(changed, 1).evaluate(params, gradient=False)[0]
        error = abs(value - best[mode]['value'])
        assert error < 1e-10
        replays.append({'mode': mode, 'generator_index': best['generator_index'], 'pole': best['pole'],
                        'reported_value': best[mode]['value'], 'replayed_value': value,
                        'absolute_difference': error})
    result = {'status': 'passed', 'scope': 'exact archive/catalogue/source binding and fresh numerical winner replays; no rate certificate',
              'archive_sha256': digest_file(args.archive),
              'historical_input_catalogue_sha256': archive['catalogue_sha256'],
              'published_catalogue_sha256': digest_file(args.catalogue),
              'identical_generator_candidate_records_sha256': records_digest,
              'screen_records_bound': len(seen), 'refinements_bound': len(archive['refinements']),
              'all_recorded_search_source_hashes_match_current': True,
              'search_source_sha256': source_checks, 'fresh_winner_replays': replays,
              'checker_sha256': digest_file(Path(__file__))}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'search_source_sha256'}))


if __name__ == '__main__':
    main()
