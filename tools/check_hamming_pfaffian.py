#!/usr/bin/env python3
"""Finite profiles and exact small full-span optima for bounded bit weights.

--full repeats the 2^21 native search and raw profiles through width five.
The default replays its stored positive witness and checks the native oracle.
No asymptotic width improvement is inferred from these finite data.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
import tempfile

from check_modular_recombination import (
    candidates, crossing, decode, encode, insert, ranks, replay, select,
)
from check_partition_channels import rank_mod
from hamming_pfaffian import matrices, full_span_candidates

ROOT = Path(__file__).resolve().parents[1]


def source_basis(matrix, dimension):
    original = list(candidates(matrix, dimension, False).items())
    basis = {}
    for i, (vector, terms) in enumerate(original):
        insert(basis, vector, 1 << i)
    return original, list(basis.values())


def write_input(path, width, basis):
    dimension = 1 << width
    with path.open('w') as out:
        out.write(f'{dimension} {len(basis)}\n')
        for vector, mask in basis:
            out.write(' '.join(f'{row:x}' for row in decode(vector, dimension))+'\n')


def at_index(index, basis):
    gray = index ^ (index >> 1)
    vector, mask = 0, 0
    for i, (other, description) in enumerate(basis):
        if (gray >> i) & 1:
            vector ^= other
            mask ^= description
    return vector, mask


def independent_ranks(vector, dimension):
    rows = decode(vector, dimension)
    dense = [[(row >> j) & 1 for j in range(dimension**2)] for row in rows]
    other = [[dense[x*dimension+y][u*dimension+v]
              for y in range(dimension) for v in range(dimension)]
             for x in range(dimension) for u in range(dimension)]
    return rank_mod(dense, 2), rank_mod(other, 2)


def oracle_checks(binary, path, basis, dimension, threshold):
    rng = random.Random(20261030)
    indices = [0, 1, (1 << len(basis))-1]
    indices += [rng.randrange(1 << len(basis)) for _ in range(128)]
    result = subprocess.run([str(binary), str(path), str(threshold), '--classify'],
                            input=''.join(f'{i}\n' for i in indices),
                            capture_output=True, text=True, check=True)
    actual = [int(line) for line in result.stdout.splitlines()]
    expected = []
    for j, index in enumerate(indices):
        vector, mask = at_index(index, basis)
        row, column = ranks(decode(vector, dimension), dimension)
        if j < 12:
            assert (row, column) == independent_ranks(vector, dimension)
        cost = row*column
        expected.append(cost if cost <= threshold else 0)
    assert actual == expected
    return len(indices)


def search(binary, path, width, original, basis, threshold):
    result = subprocess.run([str(binary), str(path), str(threshold)],
                            capture_output=True, text=True, check=True)
    lines = result.stdout.splitlines()
    found = next(line for line in lines if line.startswith('FOUND ')).split()
    cutoff = int(found[2])
    dimension = 1 << width
    pieces = []
    total = 0
    for line in lines:
        if not line.startswith('piece '):
            continue
        vector, mask = at_index(int(line.split()[1]), basis)
        row, column = independent_ranks(vector, dimension)
        assert 0 < row*column <= cutoff
        terms = [term for i, (old, description) in enumerate(original)
                 if (mask >> i) & 1 for term in description]
        pieces.append({'row_rank': row, 'column_rank': column,
                       'terms': [{'shift': degree, 'permutation': list(p)}
                                 for degree, p in terms]})
        total ^= vector
    assert total == encode(crossing(dimension), dimension)
    assert max(p['row_rank']*p['column_rank'] for p in pieces) == cutoff
    return {'width': width, 'span_dimension': len(basis), 'minimum_cutoff': cutoff,
            'native_nonzero_tensors_enumerated': (1 << len(basis))-1,
            'native_threshold': threshold, 'input_sha256': sha256(path.read_bytes()).hexdigest(),
            'search_stdout': result.stdout, 'witness': pieces,
            'witness_piece_count': len(pieces)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = ROOT/'tools/search_hamming_span.c'
    evidence = ROOT/'evidence/hamming-pfaffian.json'
    saved = json.loads(evidence.read_text()) if evidence.exists() else None
    stored_cases = {case['width']: case for case in saved['full_span_cases']} if saved else {}
    raw_cases, full_span = [], []
    with tempfile.TemporaryDirectory(prefix='factoring-hamming-') as name:
        temporary = Path(name)
        binary = temporary/'search'
        subprocess.run(['cc', '-O3', '-std=c11', '-Wall', '-Wextra', '-Werror',
                        '-o', str(binary), str(source)], check=True)
        for width in range(1, 6 if args.full else 4):
            matrix = matrices(width)
            dimension = 1 << width
            profile = [(degree, *ranks(rows, dimension))
                       for degree, rows in sorted(matrix.items())]
            raw_cases.append({'width': width, 'coefficient_pieces': len(profile),
                              'max_product': max(r*c for degree, r, c in profile),
                              'max_width_sum': max((r-1).bit_length()+(c-1).bit_length()
                                                   for degree, r, c in profile),
                              'ranks': profile, 'crossing_entries': dimension**4})
            if width > 3:
                continue
            original, basis = source_basis(matrix, dimension)
            assert len(basis) <= 8*(4*(width+3)+1)
            path = temporary/f'input-{width}.txt'
            write_input(path, width, basis)
            threshold = {1: 4, 2: 50, 3: 360}[width]
            tested = oracle_checks(binary, path, basis, dimension, threshold)
            if width < 3 or args.full:
                case = search(binary, path, width, original, basis, threshold)
            else:
                case = dict(stored_cases[width])
                assert case['input_sha256'] == sha256(path.read_bytes()).hexdigest()
            replay(case, matrix, dimension)
            if width == 2:
                independent = select(dict(full_span_candidates(width)), dimension)
                assert independent['minimum_cutoff_in_this_list'] == case['minimum_cutoff']
            case['native_oracle_comparisons'] = tested
            full_span.append(case)
            print(json.dumps({'width': width, 'span_dimension': len(basis),
                              'full_span_cutoff': case['minimum_cutoff'],
                              'witness_piece_count': case['witness_piece_count']}), flush=True)
        # The bounded-weight extension has the same direct crossing proof.
        weighted_entries = 0
        for width, weights in [(2, (1, 1, 2, 1)), (3, (1, 2, 1, 2, 1))]:
            weighted = matrices(width, weights)
            weighted_entries += (1 << width)**4
            assert len(weighted) <= 4*(2+sum(weights[1:]))+1
    paths = [Path(__file__), source, ROOT/'tools/hamming_pfaffian.py']
    result = {'status': 'effective family and finite ranks; no improved asymptotic exponent',
              'full': args.full, 'raw_profiles': raw_cases, 'full_span_cases': full_span,
              'bounded_weight_crossing_entries': weighted_entries,
              'sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
