#!/usr/bin/env python3
"""Construct and screen actual homogeneous codes for overlapping Ward forms."""

import argparse
import hashlib
from itertools import combinations, product
import json
from math import comb, isqrt
from pathlib import Path
import subprocess

from analyze_pfaffian_binary_screen import analyze, unpack_binary
from check_double_endpoint_retracts import bit_indices, native_text, optimize
from check_modular_crossing import rank_two
from explore_complement_edge_pencils import all_length_labels
from explore_symplectic_boundary_codes import coefficient_basis, nullspace, reduced_basis


def subsets(width, degree):
    return [sum(1 << i for i in indices) for indices in combinations(range(width), degree)] if 0 <= degree <= width else []


def overlapping_code(columns, length):
    pairs, number = len(columns), max(columns).bit_length()
    assert rank_two(columns) == number
    labels = all_length_labels(pairs, length)
    physical = [sum(mask << (2 * i) for i, mask in enumerate(word)) for word in labels]
    lookup = {mask: i for i, mask in enumerate(physical)}
    inputs, outputs, sigmas = [], [], []
    total = fixed = 0
    syndrome_records = []
    for sigma in range(1 << pairs):
        singles = bit_indices(sigma)
        even = [i for i in range(pairs) if not sigma >> i & 1]
        if (length - len(singles)) % 2:
            continue
        grade = (length - len(singles)) // 2
        t = len(even)
        if not 0 <= grade - number <= t - 2 * number:
            continue
        weights = [sum(((columns[i] >> a) & 1) << j for j, i in enumerate(even)) for a in range(number)]
        if rank_two(weights) != number or any((a & b).bit_count() % 2 for a, b in product(weights, repeat=2)):
            continue
        states = subsets(t, grade)
        constraints = []
        for weight in weights:
            down = {mask: 0 for mask in subsets(t, grade - 1)}
            up = {mask: 0 for mask in subsets(t, grade + 1)}
            for j, source in enumerate(states):
                for bit in bit_indices(weight):
                    if source >> bit & 1:
                        down[source ^ (1 << bit)] ^= 1 << j
                    else:
                        up[source | (1 << bit)] ^= 1 << j
            constraints += list(down.values()) + list(up.values())
        vectors = nullspace(constraints, len(states))
        expected = comb(t - 2 * number, grade - number)
        assert len(vectors) == expected
        total += (1 << len(singles)) * expected
        if not singles:
            fixed += expected
            continue
        syndrome_records.append([sigma, expected, 1 << (len(singles) - 1)])
        for orientation in range(1 << (len(singles) - 1)):
            base = 1 << (2 * singles[0])
            for j, site in enumerate(singles[1:]):
                base |= 1 << (2 * site + ((orientation >> j) & 1))
            transformed = []
            for vector in vectors:
                value = 0
                for j in bit_indices(vector):
                    occupied = base
                    for a in bit_indices(states[j]):
                        occupied |= 3 << (2 * even[a])
                    value ^= 1 << lookup[occupied]
                transformed.append(value)
            for pivot, value in reduced_basis(transformed):
                inputs.append(value)
                outputs.append(pivot)
                sigmas.append(sigma)
    d = len(inputs)
    assert d == (total - fixed) // 2
    for j, value in enumerate(inputs):
        for i, second in enumerate(inputs):
            assert not (value & second).bit_count() % 2
            assert (value >> outputs[i]) & 1 == int(i == j)
    return {'pairs': pairs, 'length': length, 'labels': d, 'ambient_occupations': labels,
            'input_vectors': inputs, 'output_indices': outputs, 'syndromes': sigmas,
            'input_Gram_rows': [0] * d, 'induced_mode_complement_rows': [0] * d,
            'radical': True, 'half_orientation_selector': 'least single pair uses low mode',
            'form_column_types': columns, 'form_count': number,
            'full_product_image_homogeneous_dimension': total, 'zero_syndrome_dimension': fixed,
            'syndrome_dimension_checks': syndrome_records,
            'maximum_input_terms': max(map(int.bit_count, inputs), default=0)}


def integer_gate(d):
    target, lo, hi = d ** 99, 0, d * d + 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** 70 >= target:
            hi = mid
        else:
            lo = mid
    return hi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--columns', default='1,1,2,2,3,3')
    parser.add_argument('--length', type=int, default=6)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--export-dir', type=Path, default=Path('/tmp/overlapping-Ward-codes'))
    parser.add_argument('--native', type=Path)
    parser.add_argument('--optimize', action='store_true')
    args = parser.parse_args()
    columns = list(map(int, args.columns.split(',')))
    code = overlapping_code(columns, args.length)
    if code['labels'] < 2:
        raise ValueError('at least two labels required')
    data, rows, cols = coefficient_basis(code)
    assert data['common_q_plus_r_power'] >= code['form_count']
    assert data['p_endpoint_zero'] and data['q_endpoint_ranks'] == [0, 0]
    name = 'overlap_' + '_'.join(map(str, columns)) + '_L' + str(args.length)
    result = dict(data, name=name, code=code,
                  exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  basis_row_blocks=[[[hex(v) for v in block] for block in matrix] for matrix in rows],
                  basis_column_blocks=[[[hex(v) for v in block] for block in matrix] for matrix in cols])
    args.export_dir.mkdir(parents=True, exist_ok=True)
    basis = args.export_dir / (name + '.json')
    basis.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    native = native_text(data, rows, cols, integer_gate(code['labels']))
    (args.export_dir / (name + '.txt')).write_text(native)
    if args.optimize:
        result['F2_optimization'] = optimize(data, rows, cols)
    if args.native:
        run = subprocess.run([str(args.native), '--field2'], input=native, text=True, stdout=subprocess.PIPE, check=True)
        screen = args.export_dir / (name + '.jsonl')
        screen.write_text(run.stdout)
        result['F2_screen'] = analyze(basis, screen)
        audited = 0
        for event in map(json.loads, run.stdout.splitlines()):
            if event['event'] != 'eligible':
                continue
            mask = unpack_binary(event['code'], len(rows))
            ranks = []
            for cut in (rows, cols):
                work = [[0] * size for size in data['block_sizes']]
                for i in bit_indices(mask):
                    for target, source in zip(work, cut[i]):
                        for j, value in enumerate(source):
                            target[j] ^= value
                ranks.append(sum(rank_two(block) for block in work))
            assert ranks == [event['row_rank'], event['column_rank']]
            audited += 1
        result['independently_audited_eligible_vectors'] = audited
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('basis_row_blocks', 'basis_column_blocks', 'code', 'F2_optimization')}), flush=True)
    if args.optimize:
        print(json.dumps({k: v for k, v in result['F2_optimization'].items() if k != 'combination_ranks'}), flush=True)


if __name__ == '__main__':
    main()
