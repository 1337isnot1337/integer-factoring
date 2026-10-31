#!/usr/bin/env python3
"""An exact 220-dimensional constant-profile family of q16 crossings.

This certifies an explicit affine family at two prescribed gauges. It does
not assert completeness among all tensorized sections or other gauges.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from check_polynomial_crossing_interpolation import (
    EXPECTED, INFINITY_HEIGHTS, LABELS, ZERO_HEIGHTS, column_elimination,
    flatten_position, gauged_matrix, height, polynomial_matrix_product,
    smith_jet_changes,
)
from check_schur_gaussian_reduction import seed_data


ROOT = Path(__file__).resolve().parent.parent
PRECISION = 32
POLYNOMIAL_STRIDE = 32


def xor_values(values):
    result = 0
    for value in values:
        result ^= value
    return result


def row_basis(rows):
    basis = {}
    for row in rows:
        while row:
            lead = row.bit_length() - 1
            if lead in basis:
                row ^= basis[lead]
            else:
                basis[lead] = row
                break
    return list(basis.values())


def binary_matrix_product(left, right):
    output = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            value ^= right[bit.bit_length() - 1]
            row ^= bit
        output.append(value)
    return output


def kronecker_binary(left, right):
    assert len(left) == len(right) == 16
    return [sum(right[j] << (16 * k) for k in range(16) if left[i] >> k & 1)
            for i in range(16) for j in range(16)]


def seed_and_changes():
    base = [0] * 256
    for entry in seed_data()['local']:
        mask = sum(v << (2 * i) for i, v in enumerate(entry['indices']))
        base[mask] = entry['matrix'][0][0]
    for labels, value in zip(LABELS, base):
        x, y, u, v = labels
        assert value.bit_count() % 2 == int(x == v and y == u)
    changes = []
    for heights, cut, infinity in ((ZERO_HEIGHTS, 0, False),
                                    (INFINITY_HEIGHTS, 1, True)):
        orders, left, right = smith_jet_changes(
            gauged_matrix(base, heights, cut, infinity), PRECISION)
        assert orders == EXPECTED
        changes.append((left, right))
    return base, changes


def section_spaces(changes):
    """Solve the complete finite Laurent windows for twists -2,-1,0,1."""
    row_offsets = []
    offset = 0
    for width in EXPECTED:
        row_offsets.append([offset + j * width for j in range(16)])
        offset += 16 * width
    first_width = offset
    cache = {}

    def encode(mask, shift, end):
        key = mask, shift, end
        if key not in cache:
            row, column = flatten_position(LABELS[mask], end)
            left = changes[end][0]
            cache[key] = sum(
                (((left[i][row] << shift) & ((1 << width) - 1))
                 << row_offsets[i][column])
                for i, width in enumerate(EXPECTED))
        return cache[key]

    spaces, reports, raw = {}, [], {}
    for twist, dimension in ((-2, 1), (-1, 6), (0, 21), (1, 56)):
        variables, columns = [], []
        for mask, labels in enumerate(LABELS):
            h0 = height(labels, ZERO_HEIGHTS)
            hi = height(labels, INFINITY_HEIGHTS) + twist
            for exponent in range(-h0, hi + 1):
                variables.append((mask, exponent))
                columns.append(encode(mask, exponent + h0, 0)
                               | (encode(mask, hi - exponent, 1) << first_width))
        basis, kernel = column_elimination(columns)
        assert len(kernel) == dimension
        for witness in kernel:
            assert xor_values(column for i, column in enumerate(columns)
                              if witness >> i & 1) == 0
        spaces[twist] = [
            [(mask, exponent) for i, (mask, exponent) in enumerate(variables)
             if witness >> i & 1] for witness in kernel]
        raw[twist] = variables, kernel
        reports.append({'twist': twist, 'variables': len(variables),
                        'constraint_rank': len(basis), 'dimension': dimension})

    zero_variables, zero_kernel = raw[0]
    positive_variables, positive_kernel = raw[1]
    lookup = {variable: i for i, variable in enumerate(positive_variables)}
    inherited = [
        sum(1 << lookup[mask, exponent + shift]
            for i, (mask, exponent) in enumerate(zero_variables)
            if witness >> i & 1)
        for witness in zero_kernel for shift in (0, 1)]
    quotient_basis, _ = column_elimination(inherited)
    assert len(quotient_basis) == 36
    complements = []
    for witness in positive_kernel:
        remaining = witness
        while remaining:
            lead = remaining.bit_length() - 1
            if lead in quotient_basis:
                remaining ^= quotient_basis[lead][0]
            else:
                quotient_basis[lead] = remaining, 0
                complements.append(
                    [(mask, exponent)
                     for i, (mask, exponent) in enumerate(positive_variables)
                     if witness >> i & 1])
                break
    assert len(complements) == 20 and len(quotient_basis) == 56
    return spaces, complements, reports


def tensor_mask(left, right):
    """Combine each pair of q4 physical labels as a+4b."""
    return sum((((left >> (2 * i)) & 3)
                | (((right >> (2 * i)) & 3) << 2)) << (4 * i)
               for i in range(4))


def tensor_polynomial(left, right):
    """Sparse coefficients of (1+z) left tensor right."""
    terms = set()
    for mask_a, exponent_a in left:
        for mask_b, exponent_b in right:
            mask = tensor_mask(mask_a, mask_b)
            for extra in (0, 1):
                term = mask, exponent_a + exponent_b + extra
                if term in terms:
                    terms.remove(term)
                else:
                    terms.add(term)
    output = {}
    for mask, exponent in terms:
        assert 0 <= exponent < POLYNOMIAL_STRIDE
        output[mask] = output.get(mask, 0) ^ (1 << exponent)
    result = sorted((mask, value) for mask, value in output.items() if value)
    assert all(value.bit_count() % 2 == 0 for _, value in result)
    return result


def polynomial_span(deltas):
    packed = [sum(value << (POLYNOMIAL_STRIDE * mask)
                  for mask, value in delta) for delta in deltas]
    assert len(row_basis(packed[:252])) == 180
    basis, independent = {}, []
    for index, column in enumerate(packed):
        while column:
            lead = column.bit_length() - 1
            if lead in basis:
                column ^= basis[lead]
            else:
                basis[lead] = column
                independent.append(index)
                break
    assert len(independent) == 220
    return independent


def normalized_residue(terms, twist, end, changes):
    """Constant term after the prescribed row factor is divided out."""
    heights = INFINITY_HEIGHTS if end else ZERO_HEIGHTS
    matrix = [[0] * 16 for _ in range(16)]
    for mask, exponent in terms:
        row, column = flatten_position(LABELS[mask], end)
        h = height(LABELS[mask], heights)
        power = h + twist - exponent if end else h + exponent
        assert power >= 0
        matrix[row][column] ^= 1 << power
    left, right = changes[end]
    transformed = polynomial_matrix_product(
        polynomial_matrix_product(left, matrix, PRECISION), right, PRECISION)
    for i, row in enumerate(transformed):
        assert all(value & ((1 << EXPECTED[i]) - 1) == 0 for value in row)
    return [sum(((value >> EXPECTED[i]) & 1) << j
                for j, value in enumerate(row))
            for i, row in enumerate(transformed)]


def residue_certificate(matrices, independent):
    selected = [matrices[i] for i in independent]
    outgoing = [0] * 256
    for matrix in matrices:
        for i, row in enumerate(matrix):
            outgoing[i] |= row
    alive = (1 << 256) - 1
    layers = []
    while alive:
        leaves = [i for i in range(256)
                  if alive >> i & 1 and not (outgoing[i] & alive)]
        assert leaves, 'the common support has a directed cycle'
        layers.append(leaves)
        alive ^= sum(1 << i for i in leaves)
    layer_of = {index: layer for layer, indices in enumerate(layers)
                for index in indices}
    for matrix in matrices:
        for i, row in enumerate(matrix):
            assert all(layer_of[j] < layer_of[i]
                       for j in range(256) if row >> j & 1)
    assert [len(layer) for layer in layers] == [1, 8, 28, 56, 70, 56, 28, 8, 1]

    image = [1 << i for i in range(256)]
    dimensions = [256]
    for _ in range(5):
        image = row_basis(
            row for matrix in selected
            for row in binary_matrix_product(image, matrix))
        dimensions.append(len(image))
    assert dimensions == [256, 240, 176, 80, 16, 0]
    return {'common_image_dimensions': dimensions,
            'common_strict_triangular_layers': layers,
            'layer_sizes': [len(layer) for layer in layers],
            'all_products_of_five_residues_vanish': True}


def check():
    base, changes = seed_and_changes()
    spaces, complements, section_reports = section_spaces(changes)
    explicit_a = [(mask, exponent + 2)
                  for mask in range(192, 256)
                  for exponent in range(base[mask & 63].bit_length())
                  if base[mask & 63] >> exponent & 1]
    assert spaces[-2][0] == explicit_a
    # Key2 denotes the selected quotient representatives, whose twist is1.
    sections = dict(spaces)
    sections[2] = complements
    descriptors = [(-1, i, 0, j) for i in range(6) for j in range(21)]
    descriptors += [(0, j, -1, i) for i in range(6) for j in range(21)]
    descriptors += [(-2, 0, 2, j) for j in range(20)]
    descriptors += [(2, j, -2, 0) for j in range(20)]
    deltas = [tensor_polynomial(sections[t][i], sections[u][j])
              for t, i, u, j in descriptors]
    independent = polynomial_span(deltas)
    endpoint_reports = []
    for end in (0, 1):
        residues = {
            t: [normalized_residue(section, 1 if t == 2 else t, end, changes)
                for section in elements]
            for t, elements in sections.items()}
        a = residues[-2][0]
        assert not any(binary_matrix_product(a, a))
        for b in residues[2]:
            assert not any(binary_matrix_product(binary_matrix_product(a, b), a))
        matrices = [kronecker_binary(residues[t][i], residues[u][j])
                    for t, i, u, j in descriptors]
        endpoint_reports.append({
            'endpoint': 'infinity' if end else 'zero',
            **residue_certificate(matrices, independent)})

    expected = tuple(sorted(a + b for a in EXPECTED for b in EXPECTED))
    data_hash = sha256(json.dumps(
        {'sections': sections, 'descriptors': descriptors, 'deltas': deltas},
        separators=(',', ':'), sort_keys=True).encode()).hexdigest()
    return {
        'status': 'pass',
        'scope': 'The explicit 220-dimensional affine family at the two prescribed tensor-product gauges.',
        'no_improved_factoring_exponent': True,
        'physical_alphabet': 16,
        'field': 'F2',
        'local_section_spaces': section_reports,
        'inherited_product_generators': 252,
        'inherited_product_rank': 180,
        'new_product_generators': 40,
        'total_generators': len(deltas),
        'total_rank': len(independent),
        'distinct_crossings': str(1 << len(independent)),
        'independent_generator_indices': independent,
        'endpoint_perturbations_vanish_at_one': True,
        'maximum_perturbation_degree': max(value.bit_length() - 1
                                          for delta in deltas for _, value in delta),
        'baseline_one_site_profile': list(EXPECTED),
        'unchanged_two_site_profile': dict(Counter(expected)),
        'endpoint_residue_certificates': endpoint_reports,
        'deterministic_construction_sha256': data_hash,
        'source_sha256': {
            name: sha256((ROOT / name).read_bytes()).hexdigest()
            for name in ('tools/check_tensorized_smith_sections.py',
                         'tools/check_polynomial_crossing_interpolation.py',
                         'tools/check_schur_gaussian_reduction.py')},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
