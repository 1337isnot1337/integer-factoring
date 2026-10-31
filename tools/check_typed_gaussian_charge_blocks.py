#!/usr/bin/env python3
"""Replay exact charge-orbit compression for typed Gaussian coefficients.

Builds the native probe, independently compares complete coefficient-span
ranks for three small types, and records larger exact rank tables.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import heapq
import json
from math import factorial, isqrt
from pathlib import Path
import subprocess
import tempfile

from check_schur_gaussian_reduction import rank_binary, seed_data
from check_typed_gaussian_generic_blocks import coefficient_matrices


ROOT = Path(__file__).resolve().parent.parent
NATIVE = ROOT / 'tools/probe_typed_gaussian_charge_blocks.cpp'


def run_probe(binary, seed, counts, all_codes=True):
    command = [str(binary), str(seed), ','.join(map(str, counts)), '20000']
    if all_codes:
        command.append('all')
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    records = [json.loads(line) for line in result.stdout.splitlines()]
    for record in records:
        record.pop('seconds', None)
    summary = records[-1]
    assert summary['done']
    expected = factorial(sum(counts))
    for count in counts:
        expected //= factorial(count)
    assert summary['D'] == expected
    blocks = [record for record in records if 'block' in record]
    for cut in (0, 1):
        assert sum(b['size'] * b['multiplicity'] for b in blocks if b['cut'] == cut) == expected ** 2
    coefficients = [r for r in records if 'degree' in r]
    ranks = [(r['R'], r['C']) for r in records if 'code' in r]
    if all_codes:
        assert ranks[0] == (0, 0)
        assert ranks[-1] == (expected ** 2, expected ** 2)
        assert all(a > 0 and b > 0 for a, b in ranks[1:]), 'coefficient list is dependent'
    return {'counts': counts, 'dimension': expected, 'blocks': blocks,
            'coefficient_ranks': coefficients, 'all_code_rank_pairs': ranks,
            'block_count': len(blocks), 'largest_block': summary['largest_block']}


def independent_full_comparison(record):
    dimension, row, column = coefficient_matrices(tuple(record['counts']))
    assert dimension == record['dimension']
    assert all(not any(row[k]) and not any(column[k]) for k in range(1, len(row), 2))
    degree = len(row) - 1
    codes = 1 << (degree // 2 + 1)
    assert codes == len(record['all_code_rank_pairs'])
    current_r = [0] * dimension ** 2
    current_c = [0] * dimension ** 2
    previous = 0
    for step in range(codes):
        code = step ^ (step >> 1)
        if step:
            k = 2 * ((code ^ previous).bit_length() - 1)
            current_r = [a ^ b for a, b in zip(current_r, row[k])]
            current_c = [a ^ b for a, b in zip(current_c, column[k])]
        assert (rank_binary(current_r), rank_binary(current_c)) == record['all_code_rank_pairs'][code]
        previous = code
    return codes


def shortest_path_certificate(record):
    pairs = record['all_code_rank_pairs']
    weights = [isqrt(a * b) for a, b in pairs]
    count = len(weights)
    distance = [None] * count
    predecessor = [None] * count
    distance[0] = 0
    queue = [(0, 0)]
    while queue:
        total, state = heapq.heappop(queue)
        if distance[state] != total:
            continue
        for code in range(1, count):
            target, cost = state ^ code, total + weights[code]
            if distance[target] is None or cost < distance[target]:
                distance[target] = cost
                predecessor[target] = state, code
                heapq.heappush(queue, (cost, target))
    target = count - 1
    witness = []
    while target:
        target, code = predecessor[target]
        witness.append(code)
    # An exact rational upper bound for the same path, without floating roots.
    scale = 10 ** 6
    upper = 0
    for code in witness:
        a, b = pairs[code]
        square = a * b * scale * scale
        root = isqrt(square)
        upper += root + (root * root < square)
    lower = distance[-1]
    assert lower * scale <= upper
    dimension = record['dimension']
    assert lower ** 12 > dimension ** 17
    return {'field': 'F_2', 'integer_floor_shortest_path_lower': lower,
            'witness_codes': witness, 'rational_upper_numerator': upper,
            'rational_upper_denominator': scale,
            'exact_optimum': lower * scale == upper,
            'above_D_power_17_over_12': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cxx', default='g++')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='typed-charge-') as directory:
        directory = Path(directory)
        binary, seed = directory / 'probe', directory / 'seed.txt'
        subprocess.run([args.cxx, '-O3', '-std=c++17', str(NATIVE), '-o', str(binary)], check=True)
        local = seed_data()['local']
        seed.write_text(''.join(str(entry['matrix'][0][0]) + '\n' for entry in local))
        records = []
        independent_codes = 0
        for counts in ((1, 1, 1, 0), (2, 1, 1, 0), (1, 1, 1, 1),
                       (2, 1, 1, 1), (6, 2, 2, 0)):
            record = run_probe(binary, seed, counts)
            if len(records) < 3:
                checked = independent_full_comparison(record)
                independent_codes += checked
                record['independent_full_matrix_code_checks'] = checked
            record['recombination_certificate'] = shortest_path_certificate(record)
            records.append(record)
        dense = run_probe(binary, seed, (4, 3, 3, 0), all_codes=False)
        # The probe must reject a seed outside its advertised structural scope.
        bad = [entry['matrix'][0][0] for entry in local]
        bad[1] ^= 3  # Preserve the endpoint but violate integer charge.
        invalid = directory / 'invalid-seed.txt'
        invalid.write_text(''.join(str(value) + '\n' for value in bad))
        rejected = subprocess.run([str(binary), str(invalid), '1,1,1,0'],
                                  text=True, capture_output=True)
        assert rejected.returncode != 0
        assert 'charge conservation' in rejected.stderr
    result = {
        'status': 'exact rank-compression controls; no new factoring exponent',
        'independent_full_matrix_code_checks': independent_codes,
        'exhaustive_span_cases': records,
        'dense_coefficient_case': dense,
        'invalid_seed_rejected': True,
        'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ('tools/check_typed_gaussian_charge_blocks.py',
                                       'tools/probe_typed_gaussian_charge_blocks.cpp',
                                       'tools/check_typed_gaussian_generic_blocks.py',
                                       'tools/check_schur_gaussian_reduction.py')},
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
