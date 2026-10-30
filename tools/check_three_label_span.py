#!/usr/bin/env python3
"""Prepare and independently check the finite native three-label span screen.

The default verifies the native rank classifier and Gray-code interval logic
against independent Python elimination. --verify-log audits a completed full
run; it does not independently repeat all 2^34 classifications. Every cutoff from one through nine is tracked.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
import tempfile

from check_modular_crossing import regular_partitions
from check_modular_recombination import (
    candidates, crossing, decode, encode, insert, profile, ranks, represent,
)
from check_partition_channels import rank_mod

ROOT = Path(__file__).resolve().parents[1]


def full_basis():
    basis = {}
    for vector in candidates(profile(regular_partitions(2, 12)), 3, False):
        insert(basis, vector, 0)
    assert len(basis) == 34
    return [vector for vector, mask in basis.values()]


def write_basis(path, generators):
    path.write_text(''.join(f'{v & ((1 << 64)-1):x} {v >> 64:x}\n' for v in generators))


def independent_cost(vector):
    rows = decode(vector, 3)
    dense = [[(row >> j) & 1 for j in range(9)] for row in rows]
    other = [[dense[x*3+y][u*3+v] for y in range(3) for v in range(3)]
             for x in range(3) for u in range(3)]
    return rank_mod(dense, 2) * rank_mod(other, 2)


def run_native(binary, path, start, end):
    text = subprocess.run([str(binary), str(path), str(start), str(end)],
                          check=True, capture_output=True, text=True).stdout
    return [json.loads(line) for line in text.splitlines()]


def verify_interval(binary, path, generators, start, end):
    records = run_native(binary, path, start, end)
    assert records[0] == {'dimension': len(generators), 'start': start, 'end': end}
    summary = records[-1]
    assert summary['complete'] and summary['checked'] == end-start
    counts = [0]*9
    bases = [{} for _ in range(9)]
    for index in range(start, end):
        gray = index ^ (index >> 1)
        vector = 0
        for i, generator in enumerate(generators):
            if (gray >> i) & 1:
                vector ^= generator
        cost = independent_cost(vector)
        if vector and cost <= 9:
            for cutoff in range(cost, 10):
                counts[cutoff-1] += 1
                insert(bases[cutoff-1], vector, 0)
    assert summary['counts'] == counts
    assert summary['span_ranks'] == [len(basis) for basis in bases]
    for record in records:
        if 'new_basis' in record:
            vector = int(record['new_basis'], 16)
            cutoff = record['cutoff']
            assert vector and independent_cost(vector) <= cutoff
            assert represent(bases[cutoff-1], vector) is not None
    return end-start


def quick(binary, temporary):
    rng = random.Random(20261030)
    path = temporary/'basis.txt'
    cases = []
    # Arbitrary matrices plus deliberately small row ranks exercise both
    # acceptance and rejection near product nine, not just dense negatives.
    for rank in (1, 2, 3, 4, 5, 9):
        for _ in range(40):
            factors = [rng.randrange(1 << 9) for _ in range(rank)]
            rows = []
            for i in range(9):
                row = 0
                for factor in factors:
                    if rng.randrange(2):
                        row ^= factor
                rows.append(row)
            cases.append(encode(rows, 3))
    cases += [0, encode(crossing(3), 3), (1 << 81)-1]
    accepts = 0
    for vector in cases:
        write_basis(path, [vector])
        records = run_native(binary, path, 0, 2)
        cost = independent_cost(vector)
        expected = [int(vector != 0 and cost <= cutoff) for cutoff in range(1, 10)]
        assert records[-1]['counts'] == expected
        accepts += expected[-1]
    assert 0 < accepts < len(cases)
    generators = full_basis()[:12]
    write_basis(path, generators)
    checked = verify_interval(binary, path, generators, 0, 4096)
    checked += verify_interval(binary, path, generators, 341, 2987)
    return {'classifier_cases': len(cases), 'accepted_classifier_cases': accepts,
            'independent_interval_vectors': checked}


def verify_log(path):
    generators = full_basis()
    source_basis = {}
    for generator in generators:
        insert(source_basis, generator, 0)
    records = [json.loads(line) for line in path.read_text().splitlines()]
    assert records[0] == {'dimension': 34, 'start': 0, 'end': 1 << 34}
    summary = records[-1]
    assert summary['complete'] and summary['checked'] == 1 << 34
    low_bases = [{} for _ in range(9)]
    for record in records:
        if 'new_basis' not in record:
            continue
        vector = int(record['new_basis'], 16)
        cutoff = record['cutoff']
        assert represent(source_basis, vector) is not None
        assert vector and independent_cost(vector) <= cutoff
        gray = record['index'] ^ (record['index'] >> 1)
        expected = 0
        for i, generator in enumerate(generators):
            if (gray >> i) & 1:
                expected ^= generator
        assert vector == expected
        assert insert(low_bases[cutoff-1], vector, 0)
    assert [len(basis) for basis in low_bases] == summary['span_ranks']
    summary = dict(summary)
    summary['crossing_in_cutoff_span'] = [represent(
        basis, encode(crossing(3), 3)) is not None for basis in low_bases]
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', type=Path)
    parser.add_argument('--verify-log', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = ROOT/'tools/search_three_label_span.c'
    if args.prepare:
        write_basis(args.prepare, full_basis())
        print(json.dumps({'basis_dimension': 34, 'basis_sha256': sha256(
            args.prepare.read_bytes()).hexdigest()}))
        return
    if args.verify_log:
        result = verify_log(args.verify_log)
    else:
        with tempfile.TemporaryDirectory(prefix='factoring-span-') as name:
            temporary = Path(name)
            binary = temporary/'screen'
            subprocess.run(['cc', '-O3', '-std=c11', '-Wall', '-Wextra', '-Werror',
                            '-o', str(binary), str(source)], check=True)
            result = quick(binary, temporary)
    result['source_sha256'] = sha256(source.read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
