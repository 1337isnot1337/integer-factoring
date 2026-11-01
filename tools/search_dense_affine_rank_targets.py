#!/usr/bin/env python3
"""Exact complementary-rectangle constraints in the dense affine frame span.

At a specified finite place, requiring a cut matrix to vanish outside
selected rows and columns forces rank at most their combined cardinality.
These are linear constraints on the full affine tensor. Natural pressure
values are bounded numerical diagnostics only, not rank-rate certificates.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
import math
from pathlib import Path
import time

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, bell_normalize, scale
from screen_affine_finite_places import bitplane_rank, evaluate, rows
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import save
from solve_dense_affine_scalar_content import Span, pack


def include(span, vector, index):
    reduced, witness = span.reduce(vector, (1 << index, 0))
    if not any(reduced):
        return
    pivot = ((reduced[0] | reduced[1]) & -(reduced[0] | reduced[1])).bit_length() - 1
    c = ((reduced[0] >> pivot) & 1) + 2 * ((reduced[1] >> pivot) & 1)
    inv = (0, 1, 3, 2)[c]
    span.pivots[pivot] = scale(reduced, inv), scale(witness, inv)


def rectangle_mask(point, cut, row_labels, column_labels):
    selected_rows = {5 * i for i in row_labels}
    selected_columns = {5 * i for i in column_labels}
    shift = {0: 0, 2: 256, 3: 512}[point]
    mask = 0
    for index in range(256):
        x, y, u, v = ((index // 4 ** j) % 4 for j in range(4))
        row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
        if row not in selected_rows and column not in selected_columns:
            mask |= 1 << (shift + index)
    return mask


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source_path = Path('evidence/dense-affine-scalar-content.json')
    catalogue_path = Path('evidence/dense-rank-one-coupled-frames.json')
    raw_source, raw_catalogue = source_path.read_bytes(), catalogue_path.read_bytes()
    source, catalogue = json.loads(raw_source), json.loads(raw_catalogue)
    assert source['complete_catalogue'] and source['catalogue_sha256'] == sha256(raw_catalogue).hexdigest()
    base = tensor(4)
    _, factors = bell_normalize(base)
    points = (0, 2, 3)
    target = pack(evaluate(v, p) for p in points for v in base['values'])
    all_candidates = [(r, c) for r in catalogue['records'] for c in r['candidates']]
    complete = next(r for r in source['results'] if r['points'] == list(points))
    basis, cache = [], {}
    for entry in complete['rank_increments']:
        r, c = all_candidates[entry['candidate']]
        if r['index'] not in cache:
            cache[r['index']] = rank_one_terms(base, r['u'], r['v'], factors)
        changed = transformed_from_terms(base, cache[r['index']], c['pole'])
        assert sha256(json.dumps(changed['values'], separators=(',', ':')).encode()).hexdigest() == c['seed_sha256']
        column = pack(evaluate(v, p) for p in points for v in changed['values'])
        basis.append({'generator_index': r['index'], 'pole': c['pole'], 'seed_sha256': c['seed_sha256'],
                      'values': changed['values'], 'column': add(column, target)})
    supports = []
    for a, b in ((2, 2), (2, 3), (3, 2), (3, 3), (3, 4), (4, 3), (4, 4)):
        for first in combinations(range(4), a):
            for second in combinations(range(4), b):
                if set(first) | set(second) == set(range(4)):
                    supports.append((first, second))
    tasks = []
    for preserve_zero in (False, True):
        for point, cut, support in product((2, 3), (0, 1), supports):
            tasks.append((preserve_zero, [(point, cut, *support)]))
        for first, second in product(supports, repeat=2):
            tasks.append((preserve_zero, [(2, 0, *first), (3, 1, *second)]))
    names = (Path(__file__).name, 'solve_dense_affine_scalar_content.py', 'dense_coupled_reflection_frames.py',
             'asymmetric_crossing_pressure.py', 'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'schema_version': 1, 'scope': __doc__, 'requested_seconds': args.seconds,
              'source_span_sha256': sha256(raw_source).hexdigest(), 'catalogue_sha256': sha256(raw_catalogue).hexdigest(),
              'basis': [{k: v for k, v in r.items() if k not in ('values', 'column')} for r in basis],
              'supports': supports, 'planned_targets': len(tasks), 'cases': [], 'retained_seeds': [],
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'feasible_targets': 0, 'unique_seeds': 0, 'unresolved_pressure_calls': 0}
    seen = set()
    began = time.monotonic()
    natural = start_params(np.asarray([0, 12, 4, 16]) * .2326459924883792, .2326459924883792)
    for index, (preserve_zero, constraints) in enumerate(tasks):
        if time.monotonic() - began >= args.seconds:
            break
        mask = (1 << 256) - 1 if preserve_zero else 0
        for point, cut, first, second in constraints:
            mask |= rectangle_mask(point, cut, first, second)
        rhs = (target[0] & mask, target[1] & mask)
        if preserve_zero:
            high = ~((1 << 256) - 1)
            rhs = rhs[0] & high, rhs[1] & high
        system = object.__new__(Span)
        system.pivots, system.increments = {}, []
        for j, candidate in enumerate(basis):
            column = candidate['column']
            include(system, (column[0] & mask, column[1] & mask), j)
        remainder, word = system.reduce(rhs)
        record = {'target_index': index, 'preserve_zero': preserve_zero, 'constraints': constraints,
                  'constraint_entries': mask.bit_count(), 'constraint_rank': len(system.pivots),
                  'feasible': not any(remainder), 'remainder_bitplanes': [hex(v) for v in remainder]}
        if not any(remainder):
            report['feasible_targets'] += 1
            record['coefficient_bitplanes'] = [hex(v) for v in word]
            values = list(base['values'])
            for j, candidate in enumerate(basis):
                c = ((word[0] >> j) & 1) + 2 * ((word[1] >> j) & 1)
                if c:
                    values = [add(v, scale(add(old, new), c))
                              for v, old, new in zip(values, base['values'], candidate['values'])]
            observed = pack(evaluate(v, p) for p in points for v in values)
            expected = (target[0] & ((1 << 256) - 1), target[1] & ((1 << 256) - 1)) if preserve_zero else (0, 0)
            assert (observed[0] & mask, observed[1] & mask) == expected
            digest = sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()
            record['seed_sha256'] = digest
            if digest not in seen:
                seen.add(digest)
                changed = dict(base, values=values)
                finite = [[bitplane_rank(rows(values, cut, point)) for cut in (0, 1)] for point in (2, 3)]
                record['finite_ranks'] = finite
                oracle = AsymmetricPressure(changed, 1)
                scores = []
                for multiplier in (1, .5, .2):
                    try:
                        value = oracle.evaluate(np.r_[natural * multiplier / math.sqrt(2), math.pi / 4])[0] / math.sqrt(2)
                        scores.append(value)
                    except (ArithmeticError, ValueError):
                        report['unresolved_pressure_calls'] += 1
                record['natural_rho'] = min(scores) if scores else None
                retained = dict(changed, seed_sha256=digest, target_index=index, finite_ranks=finite,
                                natural_rho=record['natural_rho'], coefficient_bitplanes=record['coefficient_bitplanes'])
                report['retained_seeds'].append(retained)
                # Keep both favorable finite ranks and favorable natural cost.
                all_retained = report['retained_seeds']
                by_finite = sorted(all_retained, key=lambda r: (sum(map(sum, r['finite_ranks'])), r['natural_rho'] or 1e6))[:32]
                by_natural = sorted(all_retained, key=lambda r: r['natural_rho'] or 1e6)[:32]
                report['retained_seeds'] = list({r['seed_sha256']: r for r in by_finite + by_natural}.values())
        report['cases'].append(record)
        if index % 20 == 0 or not any(remainder):
            report.update(completed=False, elapsed_seconds=time.monotonic() - began,
                          completed_targets=len(report['cases']), unique_seeds=len(seen))
            save(args.output, report)
            print(json.dumps({k: report[k] for k in ('completed_targets', 'feasible_targets', 'unique_seeds', 'elapsed_seconds')}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic() - began, completed_targets=len(report['cases']),
                  unique_seeds=len(seen), target_queue_exhausted=len(report['cases']) == len(tasks))
    save(args.output, report)


if __name__ == '__main__':
    main()
