#!/usr/bin/env python3
"""Close a polynomial delta span under independent constant paired frames.

For each elementary E_ij, the two binary operators L=E_left+E^T_right
and Q=E_left E^T_right span the differences from the F4 transvections.
The physical pairs are x/v and y/u.  All new directions have zero
crossing endpoint.  Exact degree constraints are solved on the full
resulting tensors.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

from dense_coupled_reflection_frames import digest, save
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import degree, valuation
from solve_companion_affine_degree import (
    Span, scaled, coefficient, full_vector, high_vector, unpack, witness_entries,
)


def shift(value, offset):
    return value << offset if offset >= 0 else value >> -offset


def operators():
    mask = (1 << 25) - 1
    masks = {(leg, letter): sum(mask << (25 * index) for index in range(256)
                                if index // 4 ** leg % 4 == letter)
             for leg in range(4) for letter in range(4)}
    result = []
    for pair, (left, right) in enumerate(((0, 3), (1, 2))):
        for target in range(4):
            for source in range(4):
                if target == source:
                    continue
                first = (masks[left, source], 25 * (target - source) * 4 ** left)
                second = (masks[right, target], 25 * (source - target) * 4 ** right)
                for kind in ('linear', 'quadratic'):
                    result.append({'pair': pair, 'target': target, 'source': source,
                                   'kind': kind, 'first': first, 'second': second})
    return result


def image(vector, operator):
    m1, d1 = operator['first']
    m2, d2 = operator['second']
    if operator['kind'] == 'linear':
        return tuple(shift(value & m1, d1) ^ shift(value & m2, d2) for value in vector)
    return tuple(shift(shift(value & m2, d2) & m1, d1) for value in vector)


def insert(rows, vector):
    a, b = vector
    while a or b:
        pivot = (a | b).bit_length() - 1
        value = coefficient((a, b), pivot)
        if pivot not in rows:
            row = scaled((a, b), (0, 1, 3, 2)[value])
            rows[pivot] = row
            return pivot, row
        ca, cb = scaled(rows[pivot], value)
        a, b = a ^ ca, b ^ cb
    return None


def constraint_vector(vector, low_through=-1):
    values = unpack(vector)
    high = high_vector(values)
    if low_through < 0:
        return high
    width = low_through + 1
    mask = (1 << width) - 1
    low = tuple(sum((value[plane] & mask) << (width * index)
                    for index, value in enumerate(values)) for plane in range(2))
    return tuple((a << (256 * width)) | b for a, b in zip(high, low))


def reconstruct(base, vectors, witness):
    a, b = base
    word = witness_entries(witness)
    for index, scalar in word:
        x, y = scaled(vectors[index], scalar)
        a, b = a ^ x, b ^ y
    values = unpack((a, b))
    for index, (x, y) in enumerate(values):
        i, j, k, l = (index // 4 ** leg % 4 for leg in range(4))
        assert (x.bit_count() % 2, y.bit_count() % 2) == (int(i == l and j == k), 0)
    return {'values': values, 'seed_sha256': digest(values),
            'actual_degree': max(degree(value) for value in values),
            'monomial_content': min(valuation(value) for value in values if value != (0, 0)),
            'basis_witness': word, 'witness_terms': len(word)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('basis', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--basis-output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--stop-at-degree16', action='store_true')
    args = parser.parse_args()
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    raw = args.basis.read_bytes()
    original = json.loads(gzip.decompress(raw))
    for name, expected in original['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    rows = {record['pivot']: tuple(int(value, 16) for value in record['bitplanes'])
            for record in original['basis']}
    assert len(rows) == original['dimension']
    for pivot, vector in rows.items():
        assert (vector[0] | vector[1]).bit_length() - 1 == pivot
        assert coefficient(vector, pivot) == 1
    vectors = [vector for _, vector in sorted(rows.items(), reverse=True)]
    base = tensor(4)
    packed_base = full_vector(base['values'])
    target = high_vector(base['values'])
    high = Span(True)
    for index, vector in enumerate(vectors):
        high.insert(high_vector(unpack(vector)), index)
    queue = [packed_base] + list(vectors)
    ops = operators()
    report = {
        'scope': 'exact independent paired-frame span closure and degree witnesses; no pressure certificate',
        'started_utc': started_utc,
        'input_basis_sha256': sha256(raw).hexdigest(),
        'initial_dimension': len(vectors), 'current_dimension': len(vectors),
        'operator_count': len(ops), 'queue_processed': 0, 'images_processed': 0,
        'requested_seconds': args.seconds, 'closure_complete': False,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, 'solve_companion_affine_degree.py',
                                       'reflection_dual_rail_gaussian.py',
                                       'dense_coupled_reflection_frames.py')},
    }
    save(args.output, report)
    cursor = 0
    stop = False
    while cursor < len(queue) and time.monotonic() - started < args.seconds:
        for operator in ops:
            changed = image(queue[cursor], operator)
            report['images_processed'] += 1
            added = insert(rows, changed)
            if added is not None:
                _, vector = added
                index = len(vectors)
                vectors.append(vector)
                queue.append(vector)
                high.insert(high_vector(unpack(vector)), index)
        cursor += 1
        if cursor % 10 == 0 or cursor == 1:
            possible, witness, _ = high.solve(target)
            report.update(queue_processed=cursor, current_dimension=len(vectors),
                          degree16_feasible=possible, elapsed_seconds=time.monotonic() - started)
            if possible and 'degree16_witness' not in report:
                report['degree16_witness'] = reconstruct(packed_base, vectors, witness)
                assert report['degree16_witness']['actual_degree'] == 16
            save(args.output, report)
            print(json.dumps({key: report[key] for key in
                              ('queue_processed', 'current_dimension', 'degree16_feasible',
                               'elapsed_seconds')}), flush=True)
            if possible and args.stop_at_degree16:
                stop = True
                break
    degree_results = []
    for bound in range(16, 24):
        possible, witness, residual = high.solve(target, 256 * (bound - 16))
        result = {'degree_bound': bound, 'feasible': possible,
                  'remainder': [hex(value) for value in residual]}
        if possible:
            result.update(reconstruct(packed_base, vectors, witness))
            assert 16 <= result['actual_degree'] <= bound
        degree_results.append(result)
    protected_results = []
    for low in (0, 1, 2, 3, 4, 7):
        constraints = Span(True)
        for index, vector in enumerate(vectors):
            constraints.insert(constraint_vector(vector, low), index)
        width = 256 * (low + 1)
        protected_target = target[0] << width, target[1] << width
        possible, witness, residual = constraints.solve(protected_target)
        record = {'degree_bound': 16, 'unchanged_low_coefficients_through': low,
                  'feasible': possible, 'remainder': [hex(value) for value in residual]}
        if possible:
            record.update(reconstruct(packed_base, vectors, witness))
            assert record['actual_degree'] == 16
            mask = (1 << (low + 1)) - 1
            assert all((a & mask, b & mask) == (c & mask, d & mask)
                       for (a, b), (c, d) in zip(record['values'], base['values']))
        protected_results.append(record)
    basis = {'scope': 'exact polynomial delta basis after the recorded paired-frame closure prefix',
             'input_basis_sha256': report['input_basis_sha256'],
             'source_sha256': report['source_sha256'],
             'dimension': len(vectors), 'ambient_field_coordinates': 6400,
             'basis': [{'pivot': (a | b).bit_length() - 1, 'bitplanes': [hex(a), hex(b)]}
                       for a, b in vectors],
             'closure_complete': cursor == len(queue)}
    encoded = json.dumps(basis, separators=(',', ':')).encode()
    compressed = gzip.compress(encoded, mtime=0)
    args.basis_output.parent.mkdir(parents=True, exist_ok=True)
    args.basis_output.write_bytes(compressed)
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  closure_complete=cursor == len(queue), queue_processed=cursor,
                  current_dimension=len(vectors), elapsed_seconds=time.monotonic() - started,
                  stopped_at_degree16=stop, degree_results=degree_results,
                  protected_low_jet_results=protected_results,
                  output_basis_sha256=sha256(compressed).hexdigest(),
                  output_basis_bytes=len(compressed))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in
                      ('completed', 'closure_complete', 'current_dimension',
                       'elapsed_seconds', 'output_basis_bytes')}))


if __name__ == '__main__':
    main()
