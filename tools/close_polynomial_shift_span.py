#!/usr/bin/env python3
"""Exact degree-bounded polynomial-shift closure of endpoint-zero deltas.

The working echelon order is degree first.  Its rows of degree at most d
span the entire degree-at-most-d part, so legal multiplication by w can
be tested on those rows alone.  Output uses the established entry-first
packing for the affine Hasse-target solver.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

from reflection_laurent_physical_frames import add
from search_fused_mixed_weighted import atomic_write
from solve_companion_affine_degree import Span, coefficient, scaled, unpack


ENTRIES = 256
DEGREE_BOUND = 24
WIDTH = DEGREE_BOUND + 1
MASK = (1 << WIDTH) - 1
SPREAD = [sum(1 << (ENTRIES * k) for k in range(8) if value >> k & 1)
          for value in range(256)]


def degree_packing(vector):
    output = [0, 0]
    for entry, pair in enumerate(unpack(vector)):
        for plane, value in enumerate(pair):
            expanded = 0
            for start in (0, 8, 16, 24):
                expanded |= SPREAD[(value >> start) & 255] << (ENTRIES * start)
            output[plane] |= expanded << entry
    return tuple(output)


def entry_packing(vector):
    output = [0, 0]
    for plane, value in enumerate(vector):
        remaining = value
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            position = bit.bit_length() - 1
            degree, entry = divmod(position, ENTRIES)
            output[plane] |= 1 << (WIDTH * entry + degree)
    return tuple(output)


def insert(rows, vector):
    a, b = vector
    while a or b:
        pivot = (a | b).bit_length() - 1
        value = coefficient((a, b), pivot)
        if pivot not in rows:
            row = scaled((a, b), (0, 1, 3, 2)[value])
            rows[pivot] = row
            return pivot, row
        first, second = scaled(rows[pivot], value)
        a, b = a ^ first, b ^ second
    return None


def remainder(rows, vector):
    a, b = vector
    while a or b:
        pivot = (a | b).bit_length() - 1
        if pivot not in rows:
            return a, b
        first, second = scaled(rows[pivot], coefficient((a, b), pivot))
        a, b = a ^ first, b ^ second
    return 0, 0


def filtration(rows):
    degrees = [pivot // ENTRIES for pivot in rows]
    return [sum(value <= degree for value in degrees) for degree in range(WIDTH)]


def zero_invariants(vector):
    """Check every endpoint and both paired polynomial contractions."""
    values = unpack(vector)
    assert all(a.bit_count() % 2 == b.bit_count() % 2 == 0 for a, b in values)
    for first, second in ((0, 3), (1, 2)):
        other = [leg for leg in range(4) if leg not in (first, second)]
        for i in range(4):
            for j in range(4):
                letters = [0] * 4
                letters[other[0]], letters[other[1]] = i, j
                total = (0, 0)
                for letter in range(4):
                    letters[first] = letters[second] = letter
                    index = sum(value * 4 ** leg for leg, value in enumerate(letters))
                    total = add(total, values[index])
                assert total == (0, 0)
    return 256, 32


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('basis', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--basis-output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=600)
    args = parser.parse_args()
    raw = args.basis.read_bytes()
    original = json.loads(gzip.decompress(raw))
    assert original['ambient_field_coordinates'] == ENTRIES * WIDTH
    for name, expected in original['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    input_vectors = [tuple(int(value, 16) for value in row['bitplanes']) for row in original['basis']]
    assert len(input_vectors) == original['dimension']
    started = time.monotonic()
    report = {'scope': 'Exact smallest bounded-degree shift-closed vector space reached from the given endpoint-zero delta span; no pressure certificate.',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'input_basis_sha256': sha256(raw).hexdigest(), 'initial_dimension': len(input_vectors),
              'degree_bound': DEGREE_BOUND, 'requested_seconds': args.seconds,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'solve_companion_affine_degree.py',
                                             'reflection_laurent_physical_frames.py')},
              'completed': False, 'closure_complete': False, 'queue_processed': 0}
    rows = {}
    endpoint_checks = trace_checks = packing_checks = 0
    for vector in input_vectors:
        converted = degree_packing(vector)
        assert entry_packing(converted) == vector
        packing_checks += 1
        assert insert(rows, converted) is not None
        e, t = zero_invariants(vector)
        endpoint_checks += e
        trace_checks += t
    report['initial_degree_filtration'] = filtration(rows)
    queue = [vector for pivot, vector in sorted(rows.items()) if pivot // ENTRIES < DEGREE_BOUND]
    cursor = 0
    atomic_write(args.output, report)
    while cursor < len(queue) and time.monotonic() - started < args.seconds:
        vector = queue[cursor]
        changed = tuple(value << ENTRIES for value in vector)
        assert (changed[0] | changed[1]).bit_length() <= ENTRIES * WIDTH
        added = insert(rows, changed)
        if added is not None:
            pivot, reduced = added
            if pivot // ENTRIES < DEGREE_BOUND:
                queue.append(reduced)
        cursor += 1
        if cursor % 64 == 0 or cursor == len(queue):
            report.update(current_dimension=len(rows), queue_processed=cursor, queue_length=len(queue),
                          elapsed_seconds=time.monotonic() - started)
            atomic_write(args.output, report)
            print(json.dumps({key: report[key] for key in
                              ('current_dimension', 'queue_processed', 'queue_length', 'elapsed_seconds')}), flush=True)
    complete = cursor == len(queue)
    closure_controls = 0
    if complete:
        for pivot, vector in rows.items():
            if pivot // ENTRIES < DEGREE_BOUND:
                assert remainder(rows, tuple(value << ENTRIES for value in vector)) == (0, 0)
                closure_controls += 1
    output = Span()
    for vector in rows.values():
        packed = entry_packing(vector)
        e, t = zero_invariants(packed)
        endpoint_checks += e
        trace_checks += t
        assert degree_packing(packed) == vector
        packing_checks += 1
        assert output.insert(packed, 0)
    basis = {'scope': 'Exact degree-bounded polynomial-shift closure prefix of the source delta span.',
             'input_basis_sha256': report['input_basis_sha256'],
             'source_sha256': report['source_sha256'],
             'packing': 'field coefficient at entry*25+degree in two F4 bitplanes',
             'ambient_field_coordinates': ENTRIES * WIDTH, 'dimension': len(output.rows),
             'degree_bound': DEGREE_BOUND, 'closure_complete': complete,
             'basis': [{'pivot': pivot, 'bitplanes': [hex(row[0]), hex(row[1])]}
                       for pivot, row in sorted(output.rows.items(), reverse=True)]}
    encoded = json.dumps(basis, separators=(',', ':')).encode()
    compressed = gzip.compress(encoded, mtime=0)
    args.basis_output.parent.mkdir(parents=True, exist_ok=True)
    args.basis_output.write_bytes(compressed)
    report.update(completed=True, closure_complete=complete, current_dimension=len(rows),
                  queue_processed=cursor, queue_length=len(queue),
                  final_degree_filtration=filtration(rows),
                  endpoint_coefficient_checks=endpoint_checks, polynomial_trace_checks=trace_checks,
                  packing_round_trips=packing_checks, independent_closure_controls=closure_controls,
                  elapsed_seconds=time.monotonic() - started,
                  output_basis_sha256=sha256(compressed).hexdigest(),
                  output_basis_uncompressed_sha256=sha256(encoded).hexdigest(),
                  output_basis_bytes=len(compressed))
    atomic_write(args.output, report)
    print(json.dumps({key: report[key] for key in
                      ('completed', 'closure_complete', 'current_dimension', 'elapsed_seconds', 'output_basis_bytes')}))


if __name__ == '__main__':
    main()
