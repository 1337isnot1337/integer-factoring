#!/usr/bin/env python3
"""Exact finite F4-place rank screen of a frozen q16 affine prefix.

Ranks select candidates for later complete profiles.  No pressure or
contraction certificate follows from these point ranks alone.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms, digest
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, scale
from screen_affine_finite_places import bitplane_rank, evaluate, rank
from search_filtered_reflection16_affine import affine_terms


def point_value(value, point, masks):
    lo, hi = value
    def binary(mask):
        p0, p1, p2 = ((mask & selected).bit_count() & 1 for selected in masks)
        return ((p0 ^ p2) if point == 2 else (p0 ^ p1)) | ((p1 ^ p2) << 1)
    return binary(lo) ^ F4_PRODUCT[2][binary(hi)]


def rows(entries, point, cut, masks):
    output = [(0, 0)] * 256
    for index, polynomial in entries.items():
        value = point_value(polynomial, point, masks)
        if not value:
            continue
        x, y, u, v = index & 15, (index >> 4) & 15, (index >> 8) & 15, index >> 12
        row, column = (16 * x + y, 16 * u + v) if cut == 0 else (16 * x + u, 16 * y + v)
        a, b = output[row]
        output[row] = a ^ ((value & 1) << column), b ^ ((value >> 1) << column)
    return output


def field_power(value, exponent):
    result = 1
    for _ in range(exponent):
        result = F4_PRODUCT[result][value]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.snapshot.read_bytes()
    snapshot = json.loads(raw)
    for name, expected in snapshot['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    base = tensor(8)
    _, factors = bell_normalize(base)
    masks = [sum(1 << exponent for exponent in range(residue, base['degree'] + 1, 3))
             for residue in range(3)]
    base_entries = {i: value for i, value in enumerate(base['values']) if value != (0, 0)}
    terms, cache = [], {}
    controls = 0
    for direction, record in enumerate(snapshot['directions']):
        expanded = rank_one_terms(base, record['u'], record['v'], factors)
        terms.append(expanded)
        for point in (2, 3):
            for polynomial in list(expanded[1].values())[:16]:
                assert point_value(polynomial, point, masks) == evaluate(polynomial, point)
                controls += 1
            for cut in (0, 1):
                translated = [rows(entries, point, cut, masks) for entries in expanded]
                first = [(a ^ c, b ^ d) for (a, b), (c, d) in zip(translated[1], translated[4])]
                cache[direction, point, cut] = (translated[0], first, translated[2], translated[3])
    histograms, retained = Counter(), []
    started = time.monotonic()
    independent_ranks = 0
    for index, record in enumerate(snapshot['screen']):
        observed = []
        for point in (2, 3):
            pole_value = field_power(point, (-record['pole']) % 3)
            a, b, c = record['parameters']
            coefficients = [F4_PRODUCT[a][pole_value],
                            F4_PRODUCT[b][F4_PRODUCT[pole_value][pole_value]], c]
            pair = []
            for cut in (0, 1):
                original, *directions = cache[record['direction'], point, cut]
                matrix = list(original)
                for coefficient, direction_rows in zip(coefficients, directions):
                    if coefficient:
                        matrix = [(a ^ c, b ^ d) for (a, b), (c, d) in
                                  zip(matrix, [scale(row, coefficient) for row in direction_rows])]
                actual = bitplane_rank(matrix)
                if index == 0:
                    dense = [[((a >> column) & 1) | (((b >> column) & 1) << 1)
                              for column in range(256)] for a, b in matrix]
                    assert rank(dense) == actual
                    independent_ranks += 1
                pair.append(actual)
            observed.append(pair)
        signature = tuple(value for pair in observed for value in pair)
        histograms[signature] += 1
        if min(signature) < 256:
            retained.append({'snapshot_index': index, 'direction': record['direction'],
                             'pole': record['pole'], 'parameters': record['parameters'],
                             'seed_sha256': record['seed_sha256'],
                             'natural_rho': record.get('natural_rho'), 'ranks': observed})
    retained.sort(key=lambda record: (min(value for pair in record['ranks'] for value in pair),
                                      sum(value for pair in record['ranks'] for value in pair)))
    reconstruction_entries = 0
    for record in retained[:4]:
        changed = transformed_from_terms(base, affine_terms(terms[record['direction']], record['parameters']),
                                         record['pole'])
        assert digest(changed['values']) == record['seed_sha256']
        entries = {i: value for i, value in enumerate(changed['values']) if value != (0, 0)}
        for point, expected in zip((2, 3), record['ranks']):
            for cut in (0, 1):
                assert bitplane_rank(rows(entries, point, cut, masks)) == expected[cut]
                reconstruction_entries += 65536
    names = (Path(__file__).name, 'search_filtered_reflection16_affine.py',
             'screen_affine_finite_places.py', 'dense_coupled_reflection_frames.py',
             'reflection_dual_rail_gaussian.py', 'reflection_laurent_physical_frames.py')
    result = {'scope': 'exact point ranks of the stated frozen prefix; no complete-profile or pressure claim',
              'snapshot_sha256': sha256(raw).hexdigest(), 'snapshot_bytes': len(raw),
              'screened_tensors': len(snapshot['screen']), 'snapshot_completed': snapshot.get('completed', False),
              'rational_places': [2, 3], 'singular_tensors': len(retained),
              'independent_polynomial_evaluations': controls,
              'independent_dense_rank_controls': independent_ranks,
              'full_reconstruction_entries': reconstruction_entries,
              'histogram': [{'ranks': list(signature), 'count': count} for signature, count in sorted(histograms.items())],
              'retained': retained, 'elapsed_seconds': time.monotonic() - started,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps({'screened': len(snapshot['screen']), 'singular': len(retained),
                      'best_ranks': retained[0]['ranks'] if retained else None,
                      'elapsed_seconds': result['elapsed_seconds']}))


if __name__ == '__main__':
    main()
