#!/usr/bin/env python3
"""Exact sharpness controls for the filtered Clifford compression theorem."""

from hashlib import sha256
import json
from math import comb
from pathlib import Path

from clifford_compression import (
    add_columns, bit_rank, clifford_generators, compressed_operator,
    jet_order_multiplicities, multiply_columns, row_jet_coefficients,
    validate_retraction,
)


def product(generators, indices, size):
    answer = tuple(1 << j for j in range(size))
    for index in reversed(indices):
        answer = multiply_columns(generators[index], answer)
    return answer


def commutator(left, right):
    return add_columns(multiply_columns(left, right),
                       multiply_columns(right, left))


def row_control(length, degree, compute_jet=False):
    size = 1 << length
    generators = clifford_generators(length)
    operator = product(generators, tuple(2*i for i in range(degree)), size)
    low_mask = (1 << degree)-1
    kept = [j for j in range(size)
            if j & low_mask or not (j >> degree) & 1]
    embedding = tuple(1 << j for j in kept)
    validate_retraction(length, embedding, embedding)
    assert len(kept) == size-(size >> (degree+1))
    tilted = tuple(reversed(operator))
    assert bit_rank(tilted) == size >> degree
    assert compressed_operator(tilted, embedding, embedding) == 0
    record = {
        'sites': length, 'degree': degree, 'ambient_dimension': size,
        'retained_dimension': len(kept),
        'kernel_operator_rank': size >> degree,
        'strict_threshold_equality': True,
        'row_frame_kernel_verified': True,
    }
    if compute_jet:
        coefficients = row_jet_coefficients(length, embedding, embedding, degree+1)
        counts, ranks = jet_order_multiplicities(coefficients)
        assert counts[:degree] == tuple(comb(2*length, j) for j in range(degree))
        assert counts[degree] < comb(2*length, degree)
        record['row_smith_counts_through_degree'] = counts
        record['full_frame_count_at_degree'] = comb(2*length, degree)
        record['toeplitz_ranks'] = ranks
    return record


def cubic_control():
    length, size = 7, 128
    generators = clifford_generators(length)
    terms = [(1, 3, 8), (5, 7, 9), (1, 5, 10),
             (3, 7, 11), (1, 7, 12), (3, 5, 13)]
    operator = (0,)*size
    for indices in terms:
        operator = add_columns(operator, product(generators, indices, size))
    square = multiply_columns(operator, operator)
    expected_square = product(generators, (1, 3, 5, 7), size)
    assert square == expected_square and any(square)
    first = [commutator(operator, generator) for generator in generators]
    for matrix in first:
        assert not any(multiply_columns(matrix, matrix))
    pairs = 0
    for i, left in enumerate(first):
        for right in first[:i]:
            assert not any(commutator(left, right))
            pairs += 1
    derivative = commutator(square, generators[0])
    assert derivative == product(generators, (3, 5, 7), size)
    assert any(derivative)
    ranks = bit_rank(operator), bit_rank(square), bit_rank(derivative)
    assert ranks == (38, 8, 16)
    return {
        'sites': length, 'ambient_dimension': size,
        'zero_based_generator_terms': terms,
        'operator_rank': ranks[0], 'square_rank': ranks[1],
        'noncentral_square_commutator_rank': ranks[2],
        'first_commutator_squares_checked': len(first),
        'commuting_first_pairs_checked': pairs,
        'rank_refinement_disproved': False,
        'central_square_induction_disproved': True,
    }


def main():
    records = [row_control(length, degree, compute_jet=length <= 4)
               for length in range(2, 7) for degree in range(1, length)]
    root = Path(__file__).resolve().parents[1]
    paths = ['tools/clifford_compression.py', 'tools/check_clifford_threshold_sharpness.py']
    report = {
        'status': 'pass',
        'scope': 'Exact finite controls for two proved characteristic-two statements; no factoring exponent gain.',
        'row_threshold_controls': records,
        'cubic_countercontrol': cubic_control(),
        'source_sha256': {path: sha256((root/path).read_bytes()).hexdigest()
                          for path in paths},
    }
    output = root/'evidence/clifford-threshold-sharpness.json'
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'pass', 'row_threshold_controls': len(records),
                      'row_jets_checked': sum('row_smith_counts_through_degree' in r for r in records),
                      'cubic_ranks': [38, 8, 16], 'output': str(output)}))


if __name__ == '__main__':
    main()
