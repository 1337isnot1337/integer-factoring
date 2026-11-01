#!/usr/bin/env python3
"""Exact degree reduction and full polynomial span of q4 constant companions.

All variables are affine deltas from the original strict crossing.  The
high-coefficient span answers degree16,...,23 simultaneously.  A separate
full polynomial span is retained for further exact constraints.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms, digest, save
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, degree, valuation
from search_dense_coupled_companions import translate_terms


def scaled(pair, coefficient):
    a, b = pair
    if coefficient == 1:
        return a, b
    if coefficient == 2:
        return b, a ^ b
    if coefficient == 3:
        return a ^ b, a
    return 0, 0


def coefficient(pair, position):
    return ((pair[0] >> position) & 1) | (((pair[1] >> position) & 1) << 1)


class Span:
    def __init__(self, witnesses=False):
        self.rows = {}
        self.witnesses = witnesses

    def insert(self, vector, index):
        a, b = vector
        wa, wb = (1 << index, 0) if self.witnesses else (0, 0)
        while a or b:
            pivot = (a | b).bit_length() - 1
            value = coefficient((a, b), pivot)
            if pivot not in self.rows:
                inverse = (0, 1, 3, 2)[value]
                a, b = scaled((a, b), inverse)
                wa, wb = scaled((wa, wb), inverse)
                self.rows[pivot] = (a, b, wa, wb)
                return True
            ra, rb, rwa, rwb = self.rows[pivot]
            ca, cb = scaled((ra, rb), value)
            cwa, cwb = scaled((rwa, rwb), value)
            a, b, wa, wb = a ^ ca, b ^ cb, wa ^ cwa, wb ^ cwb
        return False

    def solve(self, target, cutoff=0):
        a, b = target
        wa = wb = 0
        while a or b:
            pivot = (a | b).bit_length() - 1
            if pivot < cutoff:
                return True, (wa, wb), (a, b)
            if pivot not in self.rows:
                return False, (wa, wb), (a, b)
            value = coefficient((a, b), pivot)
            ra, rb, rwa, rwb = self.rows[pivot]
            ca, cb = scaled((ra, rb), value)
            cwa, cwb = scaled((rwa, rwb), value)
            a, b, wa, wb = a ^ ca, b ^ cb, wa ^ cwa, wb ^ cwb
        return True, (wa, wb), (0, 0)


SPREAD = [sum(1 << (256 * degree) for degree in range(8) if value >> degree & 1)
          for value in range(256)]


def full_vector(values):
    return tuple(sum(value[plane] << (25 * index) for index, value in enumerate(values))
                 for plane in range(2))


def high_vector(values):
    a = b = 0
    for index, (low, high) in enumerate(values):
        a ^= SPREAD[(low >> 17) & 255] << index
        b ^= SPREAD[(high >> 17) & 255] << index
    return a, b


def unpack(pair):
    mask = (1 << 25) - 1
    return [((pair[0] >> (25 * index)) & mask,
             (pair[1] >> (25 * index)) & mask) for index in range(256)]


def witness_entries(pair):
    output = []
    remaining = pair[0] | pair[1]
    while remaining:
        bit = remaining & -remaining
        index = bit.bit_length() - 1
        output.append([index, coefficient(pair, index)])
        remaining ^= bit
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('evidence/dense-coupled-companion-pressure.json'))
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--basis-output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    raw_archive, raw_catalogue = args.archive.read_bytes(), args.catalogue.read_bytes()
    archive, catalogue = json.loads(raw_archive), json.loads(raw_catalogue)
    assert archive['completed'] and archive['complete_companion_screen']
    assert archive['catalogue_sha256'] == sha256(raw_catalogue).hexdigest()
    for name, expected in archive['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    base = tensor(4)
    _, factors = bell_normalize(base)
    base_full, target = full_vector(base['values']), high_vector(base['values'])
    high_span, full_span = Span(True), Span(False)
    full_deltas = []
    terms_cache = {}
    report = {
        'scope': 'exact affine degree feasibility and full polynomial delta span; no pressure certificate',
        'started_utc': started_utc,
        'archive_sha256': sha256(raw_archive).hexdigest(),
        'catalogue_sha256': sha256(raw_catalogue).hexdigest(),
        'input_seed_sha256': digest(base['values']),
        'variables_are_screen_indices': True,
        'tested_degree_bounds': list(range(16, 24)),
        'coordinate_packing': 'full: entry*25+degree; high: (degree-17)*256+entry',
        'records_reconstructed': 0,
        'full_polynomial_span_dimension': 0,
        'high_coefficient_span_dimension': 0,
        'complete_catalogue': False,
        'first_feasible_prefixes': {},
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, 'search_dense_coupled_companions.py',
                                       'dense_coupled_reflection_frames.py',
                                       'reflection_dual_rail_gaussian.py',
                                       'reflection_laurent_physical_frames.py')},
    }
    save(args.output, report)
    for index, row in enumerate(archive['screen']):
        generator_index = row['generator_index']
        if generator_index not in terms_cache:
            # The archive is grouped by generator; discard stale large terms.
            generator = catalogue['records'][generator_index]
            terms_cache = {generator_index: rank_one_terms(base, generator['u'], generator['v'], factors)}
        terms = terms_cache[generator_index]
        changed = transformed_from_terms(base, translate_terms(terms, row['constant']), row['pole'])
        assert digest(changed['values']) == row['seed_sha256']
        delta = [(a ^ c, b ^ d) for (a, b), (c, d) in zip(base['values'], changed['values'])]
        packed = full_vector(delta)
        full_deltas.append(packed)
        high_span.insert(high_vector(delta), index)
        full_span.insert(packed, index)
        if index % 128 == 127 or index + 1 == len(archive['screen']):
            feasible = []
            for bound in report['tested_degree_bounds']:
                possible, _, _ = high_span.solve(target, cutoff=256 * (bound - 16))
                if possible:
                    feasible.append(bound)
                    report['first_feasible_prefixes'].setdefault(str(bound), index + 1)
            report.update(records_reconstructed=index + 1,
                          full_polynomial_span_dimension=len(full_span.rows),
                          high_coefficient_span_dimension=len(high_span.rows),
                          feasible_degree_bounds_so_far=feasible,
                          elapsed_seconds=time.monotonic() - started)
            save(args.output, report)
            print(json.dumps({key: report[key] for key in
                              ('records_reconstructed', 'full_polynomial_span_dimension',
                               'high_coefficient_span_dimension', 'feasible_degree_bounds_so_far',
                               'elapsed_seconds')}), flush=True)
    results = []
    for bound in report['tested_degree_bounds']:
        possible, witness, residual = high_span.solve(target, cutoff=256 * (bound - 16))
        result = {'degree_bound': bound, 'feasible': possible,
                  'target_remainder': [hex(value) for value in residual]}
        if possible:
            word = witness_entries(witness)
            a, b = base_full
            for index, scalar in word:
                da, db = scaled(full_deltas[index], scalar)
                a, b = a ^ da, b ^ db
            values = unpack((a, b))
            actual = max(degree(value) for value in values)
            assert 16 <= actual <= bound
            for index, (low, high) in enumerate(values):
                x, y, u, v = (index // 4 ** leg % 4 for leg in range(4))
                assert (low.bit_count() % 2, high.bit_count() % 2) == (int(x == v and y == u), 0)
            result.update(witness=word, witness_terms=len(word), actual_degree=actual,
                          seed_sha256=digest(values), values=values,
                          monomial_content=min(valuation(value) for value in values if value != (0, 0)))
        results.append(result)
    basis = {
        'scope': 'complete normalized echelon basis of the polynomial affine delta span',
        'archive_sha256': report['archive_sha256'],
        'source_sha256': report['source_sha256'],
        'packing': 'field coefficient at entry*25+degree in two F4 bitplanes',
        'ambient_field_coordinates': 6400,
        'dimension': len(full_span.rows),
        'basis': [{'pivot': pivot, 'bitplanes': [hex(row[0]), hex(row[1])]}
                  for pivot, row in sorted(full_span.rows.items(), reverse=True)],
    }
    encoded = json.dumps(basis, separators=(',', ':')).encode()
    compressed = gzip.compress(encoded, mtime=0)
    args.basis_output.parent.mkdir(parents=True, exist_ok=True)
    args.basis_output.write_bytes(compressed)
    report.update(completed=True, complete_catalogue=True,
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started,
                  degree_results=results,
                  full_span_basis_sha256=sha256(compressed).hexdigest(),
                  full_span_basis_uncompressed_sha256=sha256(encoded).hexdigest(),
                  full_span_basis_bytes=len(compressed),
                  minimum_feasible_degree=min((row['actual_degree'] for row in results if row['feasible']),
                                              default=None))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in
                      ('completed', 'elapsed_seconds', 'full_polynomial_span_dimension',
                       'minimum_feasible_degree', 'full_span_basis_bytes')}))


if __name__ == '__main__':
    main()
