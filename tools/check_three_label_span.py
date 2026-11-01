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
    target = encode(crossing(3), 3)
    summary['crossing_in_cutoff_span'] = [represent(
        basis, target) is not None for basis in low_bases]
    minimum = summary['crossing_in_cutoff_span'].index(True)+1
    summary['minimum_cutoff'] = minimum
    if minimum > 1:
        basis = low_bases[minimum-2]
        def reduced(vector):
            for pivot in sorted(basis, reverse=True):
                if (vector >> pivot) & 1:
                    vector ^= basis[pivot][0]
            return vector
        residual = reduced(target)
        assert residual
        free = residual.bit_length()-1
        functional = sum(((reduced(1 << j) >> free) & 1) << j for j in range(81))
        assert (functional & target).bit_count() % 2 == 1
        assert all((functional & vector).bit_count() % 2 == 0
                   for vector, mask in basis.values())
        summary['separating_functional_below_minimum'] = hex(functional)
    summary['log_sha256'] = sha256(path.read_bytes()).hexdigest()
    summary['audit_checker_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return summary


def verify_small_products():
    """Independent complete classification below four via rank-one realignment."""
    source = {}
    for generator in full_basis():
        insert(source, generator, 0)
    small_ranks = {vector: rank_mod(
        [[(vector >> (3*i+j)) & 1 for j in range(3)] for i in range(3)], 2)
        for vector in range(1, 512)}
    accepted = set()
    tested = 0
    from check_modular_recombination import permute
    for left, left_rank in small_ranks.items():
        for right, right_rank in small_ranks.items():
            if left_rank*right_rank > 3:
                continue
            tested += 1
            vector = sum(right << (9*i) for i in range(9) if (left >> i) & 1)
            if represent(source, vector) is None:
                continue
            accepted.add(vector)
            exchanged = encode(permute(decode(vector, 3), 3, (0, 2, 1, 3)), 3)
            assert represent(source, exchanged) is not None
            accepted.add(exchanged)
    basis = {}
    functional = 0x120000000000000000
    for vector in accepted:
        assert independent_cost(vector) <= 3
        assert (functional & vector).bit_count() % 2 == 0
        insert(basis, vector, 0)
    target = encode(crossing(3), 3)
    assert (functional & target).bit_count() % 2 == 1
    assert represent(basis, target) is None
    assert tested == 47677 and len(accepted) == 489 and len(basis) == 27
    return {'rank_one_outer_products_tested': tested,
            'all_product_below_four_tensors': len(accepted),
            'low_product_span_rank': len(basis),
            'separating_functional': hex(functional)}


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
    result['independent_low_product_classification'] = verify_small_products()
    result['source_sha256'] = sha256(source.read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
