#!/usr/bin/env python3
"""Intersect the full dense-frame affine span with finite low-rank targets.

Targets are scalar multiples of one of the three Bell pairings at each of
the two nontrivial F4 points. The optional zero-place constraint keeps the
original constant coefficient. This finite exact search does not assert
any pressure improvement.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import time

from affine_place_profile_oracle import primitive_seed
from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, bell_normalize, scale
from physical_split_crossings import transform_leg
from screen_affine_finite_places import evaluate
from solve_dense_affine_scalar_content import Span, pack


def include(span, vector, index):
    reduced, witness = span.reduce(vector, (1 << index, 0))
    if not any(reduced):
        return
    position = ((reduced[0] | reduced[1]) & -(reduced[0] | reduced[1])).bit_length() - 1
    c = ((reduced[0] >> position) & 1) + 2 * ((reduced[1] >> position) & 1)
    inv = (0, 1, 3, 2)[c]
    span.pivots[position] = scale(reduced, inv), scale(witness, inv)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--span', type=Path, default=Path('evidence/dense-affine-scalar-content.json'))
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw_span, raw_catalogue = args.span.read_bytes(), args.catalogue.read_bytes()
    source, catalogue = json.loads(raw_span), json.loads(raw_catalogue)
    assert source['complete_catalogue'] and source['catalogue_sha256'] == sha256(raw_catalogue).hexdigest()
    base = tensor(4)
    bell, factors = bell_normalize(base)
    all_candidates = [(r, c) for r in catalogue['records'] for c in r['candidates']]
    complete = next(r for r in source['results'] if r['points'] == [0, 2, 3])
    basis, terms_cache = [], {}
    for entry in complete['rank_increments']:
        r, c = all_candidates[entry['candidate']]
        if r['index'] not in terms_cache:
            terms_cache[r['index']] = rank_one_terms(base, r['u'], r['v'], factors)
        changed = transformed_from_terms(base, terms_cache[r['index']], c['pole'])
        assert sha256(json.dumps(changed['values'], separators=(',', ':')).encode()).hexdigest() == c['seed_sha256']
        basis.append({'generator_index': r['index'], 'pole': c['pole'],
                      'seed_sha256': c['seed_sha256'], 'values': changed['values']})
    templates = []
    for kind in ('xy_uv', 'xu_yv', 'xv_yu'):
        values = []
        for index in range(256):
            x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
            value = int((x == y and u == v) if kind == 'xy_uv' else
                        (x == u and y == v) if kind == 'xu_yv' else (x == v and y == u))
            values.append((value, 0))
        # A target value is not itself a strict path, so use the raw linear
        # maps rather than the endpoint-checking bell_restore wrapper.
        values, dimensions = transform_leg(values, [4] * 4, 0, factors['first_factor'])
        values, _ = transform_leg(values, dimensions, 3, list(map(list, zip(*factors['inverse_factor']))))
        templates.append((kind, [v[0] + 2 * v[1] for v in values]))
    began = time.monotonic()
    records, seeds = [], []
    for preserve_zero in (False, True):
        points = [0, 2, 3] if preserve_zero else [2, 3]
        system = Span(points, base)
        for index, candidate in enumerate(basis):
            column = pack(evaluate(v, p) for p in points for v in candidate['values'])
            include(system, add(column, system.target), index)
        target = pack((0 if p == 0 else evaluate(v, p)) for p in points for v in base['values'])
        for name_a, a in templates:
            for name_b, b in templates:
                trial = deepcopy(system)
                include(trial, pack(v if p == 2 else 0 for p in points for v in a), len(basis))
                include(trial, pack(v if p == 3 else 0 for p in points for v in b), len(basis) + 1)
                residual, word = trial.reduce(target)
                record = {'preserve_zero': preserve_zero, 'point2_pairing': name_a,
                          'point3_pairing': name_b, 'consistent': not any(residual),
                          'span_rank': len(system.pivots), 'augmented_span_rank': len(trial.pivots)}
                if not any(residual):
                    coefficients = [((word[0] >> i) & 1) + 2 * ((word[1] >> i) & 1)
                                    for i in range(len(basis) + 2)]
                    values = list(base['values'])
                    terms = []
                    for candidate, coefficient in zip(basis, coefficients):
                        if not coefficient:
                            continue
                        values = [add(value, scale(add(old, new), coefficient))
                                  for value, old, new in zip(values, base['values'], candidate['values'])]
                        terms.append({k: v for k, v in dict(candidate, coefficient=coefficient).items() if k != 'values'})
                    from physical_split_crossings import F4_PRODUCT
                    for point, template, scalar in ((2, a, coefficients[-2]), (3, b, coefficients[-1])):
                        assert [evaluate(v, point) for v in values] == [F4_PRODUCT[scalar][v] for v in template]
                    if preserve_zero:
                        assert [evaluate(v, 0) for v in values] == [evaluate(v, 0) for v in base['values']]
                    assert [evaluate(v, 1) for v in values] == [evaluate(v, 1) for v in base['values']]
                    changed, common = primitive_seed(dict(base, values=values))
                    digest = sha256(json.dumps(changed['values'], separators=(',', ':')).encode()).hexdigest()
                    record.update(scalars=coefficients[-2:], word=terms, common_factor=common,
                                  primitive_degree=changed['degree'], seed_sha256=digest)
                    seeds.append(dict(changed, target_record=len(records), seed_sha256=digest))
                records.append(record)
    report = {'schema_version': 1, 'scope': __doc__, 'span_sha256': sha256(raw_span).hexdigest(),
              'catalogue_sha256': sha256(raw_catalogue).hexdigest(), 'basis_candidates': len(basis),
              'records': records, 'seeds': seeds, 'elapsed_seconds': time.monotonic() - began,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in
                                (Path(__file__).name, 'solve_dense_affine_scalar_content.py', 'affine_place_profile_oracle.py')}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'basis_candidates': len(basis), 'targets_tested': len(records), 'feasible_targets': len(seeds),
                      'elapsed_seconds': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
