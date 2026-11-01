#!/usr/bin/env python3
"""Bind and independently replay the completed constant-companion catalogue."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from asymmetric_crossing_pressure import AsymmetricPressure
from check_dense_coupled_reflection_frames import coefficients, apply_matrix, multiply
from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms, digest
from native_rectangular_smith_f4 import RectangularOracle
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from search_dense_coupled_companions import translate_terms


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


def reference_companion(seed, factors, left, right, pole, constant):
    """Separate exponent dictionaries; four complete matrices for f=w^-k+c."""
    size = seed['q_x']
    values = [coefficients(value) for value in seed['values']]
    values = apply_matrix(values, factors['inverse_factor'], 0, size)
    values = apply_matrix(values, list(map(list, zip(*factors['first_factor']))), 3, size)
    generator = [[multiply(u, v) for v in right] for u in left]
    dual = list(map(list, zip(*generator)))
    for leg, matrix in enumerate((generator, dual, generator, dual)):
        values = apply_matrix(values, matrix, leg, size, -pole, identity=True)
        constant_matrix = [[multiply(constant, value) for value in row] for row in matrix]
        values = apply_matrix(values, constant_matrix, leg, size, identity=True)
    values = apply_matrix(values, factors['first_factor'], 0, size)
    values = apply_matrix(values, list(map(list, zip(*factors['inverse_factor']))), 3, size)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path,
                        default=Path('evidence/dense-coupled-companion-pressure.json'))
    parser.add_argument('--catalogue', type=Path,
                        default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    archive = json.loads(args.archive.read_text())
    catalogue = json.loads(args.catalogue.read_text())
    assert archive['completed'] and archive['complete_companion_screen']
    assert archive['catalogue_sha256'] == file_digest(args.catalogue)
    for name, expected in archive['source_sha256'].items():
        assert file_digest(Path(__file__).with_name(name)) == expected, name
    base = tensor(4)
    _, factors = bell_normalize(base)
    assert catalogue['source_seed_sha256'] == digest(base['values'])
    expected = {(record['index'], pole, constant)
                for record in catalogue['records']
                for pole in record['pole_orders'] for constant in (1, 2, 3)}
    indexed = {}
    for row in archive['screen']:
        key = row['generator_index'], row['pole'], row['constant']
        assert key in expected and key not in indexed
        assert 0 <= row['profile_index'] < len(archive['profiles'])
        indexed[key] = row
    assert set(indexed) == expected and len(expected) == 21528
    assert archive['counters']['duplicate_tensors'] == 0
    hashes = set()
    for generator in catalogue['records']:
        terms = rank_one_terms(base, generator['u'], generator['v'], factors)
        for constant in (1, 2, 3):
            translated = translate_terms(terms, constant)
            for pole in generator['pole_orders']:
                changed = transformed_from_terms(base, translated, pole)
                actual = digest(changed['values'])
                assert actual == indexed[generator['index'], pole, constant]['seed_sha256']
                hashes.add(actual)
    assert len(hashes) == len(expected)
    for row in archive['refinements']:
        key = row['generator_index'], row['pole'], row['constant']
        assert row['seed_sha256'] == indexed[key]['seed_sha256']
        generator = catalogue['records'][row['generator_index']]
        assert row['u'] == generator['u'] and row['v'] == generator['v']
    controls = []
    selected = [0, 100, 900, 2715, 4000, 5354]
    selected.append(next(record['index'] for record in catalogue['records']
                         if max(record['pole_orders']) == 8))
    for index in dict.fromkeys(selected):
        generator = catalogue['records'][index]
        pole = max(generator['pole_orders'])
        terms = rank_one_terms(base, generator['u'], generator['v'], factors)
        for constant in (1, 2, 3):
            reference = reference_companion(base, factors, generator['u'], generator['v'],
                                            pole, constant)
            changed = transformed_from_terms(base, translate_terms(terms, constant), pole)
            assert reference == [coefficients(value) for value in changed['values']]
            controls.append({'generator_index': index, 'pole': pole, 'constant': constant})
    replays = []
    for mode in ('balanced', 'asymmetric'):
        best = min((row for row in archive['refinements'] if mode in row),
                   key=lambda row: row[mode]['value'])
        terms = rank_one_terms(base, best['u'], best['v'], factors)
        changed = transformed_from_terms(base, translate_terms(terms, best['constant']),
                                         best['pole'])
        assert digest(changed['values']) == best['seed_sha256']
        parameters = best[mode]['params']
        if mode == 'balanced':
            value = RectangularOracle(changed['values'], 4, 4).types(parameters)[0]
        else:
            value = AsymmetricPressure(changed, 1).evaluate(parameters, gradient=False)[0]
        error = abs(value - best[mode]['value'])
        assert error < 1e-10
        replays.append({'mode': mode, 'generator_index': best['generator_index'],
                        'pole': best['pole'], 'constant': best['constant'],
                        'reported_value': best[mode]['value'], 'replayed_value': value,
                        'absolute_difference': error})
    result = {
        'status': 'passed',
        'scope': 'complete exact companion reconstruction and independent Laurent controls; numerical winner replay only',
        'archive_sha256': file_digest(args.archive),
        'catalogue_sha256': file_digest(args.catalogue),
        'all_recorded_search_source_hashes_match_current': True,
        'complete_tensor_hashes_reconstructed': len(hashes),
        'screen_records_bound': len(indexed),
        'refinement_records_bound': len(archive['refinements']),
        'independent_translation_controls': controls,
        'independent_laurent_tensor_entries': 256 * len(controls),
        'fresh_winner_replays': replays,
        'elapsed_seconds': time.monotonic() - started,
        'source_sha256': {name: file_digest(Path(__file__).with_name(name))
                          for name in (Path(__file__).name,
                                       'check_dense_coupled_reflection_frames.py')},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('source_sha256', 'independent_translation_controls')}))


if __name__ == '__main__':
    main()
