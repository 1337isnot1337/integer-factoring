#!/usr/bin/env python3
"""Bind and replay the completed five-direction q16 affine search."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from asymmetric_crossing_pressure import AsymmetricPressure
from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms, digest
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from search_filtered_reflection16_affine import aligned_directions, affine_terms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path,
                        default=Path('evidence/filtered-reflection16-affine-search.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    raw = args.archive.read_bytes()
    archive = json.loads(raw)
    assert archive['completed'] and archive['complete_requested_screen']
    for name, expected in archive['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name
    base = tensor(8)
    _, factors = bell_normalize(base)
    terms = []
    for exact, stored in zip(aligned_directions(), archive['directions']):
        assert exact == (stored['u'], stored['v'])
        common = sorted(set(stored['allowed_poles_by_nonzero_scalar'][0]).intersection(
            *map(set, stored['allowed_poles_by_nonzero_scalar'][1:])))
        assert common == stored['common_allowed_poles']
        terms.append(rank_one_terms(base, *exact, factors))
    expected_proposals = 63 * sum(len(set(row['common_allowed_poles']).intersection(
        archive['requested_poles'])) for row in archive['directions'])
    assert expected_proposals == archive['counters']['proposed_sums'] == 1386
    assert len(archive['screen']) + archive['counters']['duplicate_tensors'] == expected_proposals
    known, hashes = {}, set()
    for row in archive['screen']:
        key = row['direction'], row['pole'], tuple(row['parameters'])
        assert key not in known
        assert row['pole'] in archive['directions'][row['direction']]['common_allowed_poles']
        assert any(row['parameters']) and all(value in range(4) for value in row['parameters'])
        changed = transformed_from_terms(base, affine_terms(terms[row['direction']], row['parameters']),
                                         row['pole'])
        actual = digest(changed['values'])
        assert actual == row['seed_sha256'] and actual not in hashes
        assert changed['degree'] == row['degree'] == 96
        assert changed['monomial_content'] == row['monomial_content'] == 0
        assert 0 <= row['profile_index'] < len(archive['profiles'])
        assert row['natural_rho'] == archive['profiles'][row['profile_index']]['rho']
        hashes.add(actual)
        known[key] = row
    assert len(hashes) == 211
    for row in archive['refinements']:
        original = known[row['direction'], row['pole'], tuple(row['parameters'])]
        assert row['seed_sha256'] == original['seed_sha256']
    replays = []
    for mode in ('balanced', 'asymmetric'):
        best = min((row for row in archive['refinements'] if row['mode'] == mode),
                   key=lambda row: row['best_support'])
        changed = transformed_from_terms(base, affine_terms(terms[best['direction']], best['parameters']),
                                         best['pole'])
        value = AsymmetricPressure(changed, 1).evaluate(best['best_parameters'], gradient=False)[0]
        error = abs(value - best['best_support'])
        assert error < 1e-9
        replays.append({'mode': mode, 'direction': best['direction'], 'pole': best['pole'],
                        'parameters': best['parameters'], 'seed_sha256': best['seed_sha256'],
                        'reported_support': best['best_support'], 'replayed_support': value,
                        'balanced_rho': value / math.sqrt(2) if mode == 'balanced' else None,
                        'is_original_seed': best['seed_sha256'] == digest(base['values']),
                        'absolute_difference': error})
    result = {
        'status': 'passed',
        'scope': 'complete retained polynomial-tensor and archive/source binding; fresh numerical winner replay only',
        'archive_sha256': sha256(raw).hexdigest(),
        'all_recorded_search_source_hashes_match_current': True,
        'declared_affine_parameter_choices': expected_proposals,
        'distinct_complete_tensors_reconstructed': len(hashes),
        'complete_polynomial_entries_reconstructed': len(hashes) * 16 ** 4,
        'refinement_records_bound': len(archive['refinements']),
        'winner_replays': replays,
        'elapsed_seconds': time.monotonic() - started,
        'source_sha256': {Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
