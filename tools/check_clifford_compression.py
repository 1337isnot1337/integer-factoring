#!/usr/bin/env python3
"""Check Clifford compression identities and exact pressure exclusions."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from random import Random

from clifford_compression import (
    add_columns, bit_rank, clifford_generators, clifford_monomials,
    column_jet_coefficients, compressed_clifford_basis,
    jet_order_multiplicities, multiply_columns, row_jet_coefficients, validate_retraction,
)
from check_block_smith_crossings import reverse_polynomial
from check_integer_shared_minor_floor import (
    add, balanced_polynomial, multiply, positive_certificate,
    profile_polynomial, threshold,
)
from linear_physical_restriction import (
    linear_restriction_matrices, linear_restriction_profiles,
)


ROOT = Path(__file__).resolve().parent.parent


def linear_combination(matrices, mask):
    answer = (0,)*len(matrices[0])
    for index, matrix in enumerate(matrices):
        if mask >> index & 1:
            answer = add_columns(answer, matrix)
    return answer


def graph_retraction(length, dimension, rng, *, coordinate=False):
    size = 1 << length
    if coordinate:
        identity = [1 << i for i in range(dimension)]
        return identity, identity
    complement = size-dimension
    lower = [rng.getrandbits(dimension) for _ in range(complement)]
    right = [rng.getrandbits(complement) for _ in range(dimension)]
    embedding = [((1 << j) | sum(((row >> j) & 1) << (dimension+i)
                                 for i, row in enumerate(lower)))
                 for j in range(dimension)]
    retraction = []
    for j, row in enumerate(right):
        left = 1 << j
        for i, value in enumerate(lower):
            if row >> i & 1:
                left ^= value
        retraction.append(left | (row << dimension))
    permutation = list(range(size))
    rng.shuffle(permutation)

    def change(vector):
        return sum(((vector >> i) & 1) << permutation[i] for i in range(size))

    embedding, retraction = list(map(change, embedding)), list(map(change, retraction))
    validate_retraction(length, embedding, retraction)
    return embedding, retraction


def operator_controls(rng):
    first_layer = commutator_bounds = quadratic_dichotomies = 0
    for length in range(1, 5):
        size = 1 << length
        identity = tuple(1 << i for i in range(size))
        zero = (0,)*size
        generators = clifford_generators(length)
        matrices = (identity,)+generators
        for mask in range(1, 1 << len(matrices)):
            matrix = linear_combination(matrices, mask)
            square = multiply_columns(matrix, matrix)
            assert square in (zero, identity)
            if square == zero:
                assert bit_rank(matrix) == size//2
                index = ((mask >> 1) & -(mask >> 1)).bit_length()-1
                partner = generators[index ^ 1]
                assert add_columns(multiply_columns(matrix, partner),
                                   multiply_columns(partner, matrix)) == identity
            else:
                assert bit_rank(matrix) == size
            first_layer += 1
    for length, degree, samples in ((1, 2, None), (2, 2, None), (3, 2, 768),
                                    (4, 2, 384), (4, 3, 384)):
        all_monomials = clifford_monomials(length)
        masks = [mask for mask in range(len(all_monomials)) if mask.bit_count() <= degree]
        matrices = tuple(all_monomials[mask] for mask in masks)
        choices = range(1, 1 << len(masks)) if samples is None else (
            rng.randrange(1, 1 << len(masks)) for _ in range(samples))
        identity = all_monomials[0]
        generators = clifford_generators(length)
        for choice in choices:
            matrix = linear_combination(matrices, choice)
            selected = max((mask for i, mask in enumerate(masks) if choice >> i & 1),
                           key=int.bit_count)
            derivative = matrix
            for index in range(2*length):
                if selected >> index & 1:
                    partner = generators[index ^ 1]
                    derivative = add_columns(multiply_columns(derivative, partner),
                                             multiply_columns(partner, derivative))
            assert derivative == identity
            assert bit_rank(matrix)*(1 << selected.bit_count()) >= 1 << length
            commutator_bounds += 1
            if degree == 2:
                zero = (0,)*(1 << length)
                square = multiply_columns(matrix, matrix)
                assert square == zero or bit_rank(matrix) >= 1 << (length-1)
                first_squares, brackets = [], []
                for generator in generators:
                    bracket = add_columns(multiply_columns(matrix, generator),
                                          multiply_columns(generator, matrix))
                    bracket_square = multiply_columns(bracket, bracket)
                    assert bracket_square in (zero, identity)
                    first_squares.append(bracket_square)
                    brackets.append(bracket)
                all_singular = all(value == zero for value in first_squares)
                if all_singular:
                    all_singular = all(
                        add_columns(multiply_columns(a, b), multiply_columns(b, a)) == zero
                        for i, a in enumerate(brackets) for b in brackets[i+1:])
                if all_singular:
                    assert square in (zero, identity)
                quadratic_dichotomies += 1
    return {'nonzero_first_layer_operators': first_layer,
            'iterated_commutator_rank_controls': commutator_bounds,
            'quadratic_rank_dichotomies_checked': quadratic_dichotomies}


def restriction_controls(rng):
    records = []
    coefficient_comparisons = profile_comparisons = 0
    parity_controls = []
    for length in range(1, 5):
        size = 1 << length
        for dimension in range(size//2+1, size+1):
            for instance in range(3):
                embedding, retraction = graph_retraction(length, dimension, rng,
                                                        coordinate=instance == 0)
                compressed = compressed_clifford_basis(length, embedding, retraction)
                first = [compressed[0]]+[compressed[1 << i] for i in range(2*length)]
                assert bit_rank(first) == 2*length+1
                jet = column_jet_coefficients(length, embedding, retraction, 2)
                counts, ranks = jet_order_multiplicities(jet)
                assert counts == (1, 2*length)
                records.append({'sites': length, 'physical_dimension': dimension,
                                'instance': instance, 'column_initial_counts': counts,
                                'toeplitz_ranks': ranks})
                if instance == 0 or (length == 4 and dimension in (9, 12, 15) and instance == 1):
                    full = linear_restriction_matrices(length, embedding, retraction)
                    exact_jet = column_jet_coefficients(length, embedding, retraction, length+1)
                    for degree, rows in enumerate(exact_jet):
                        expected = [sum(((value >> degree) & 1) << j for j, value in enumerate(row))
                                    for row in full['column']]
                        assert rows == expected
                        coefficient_comparisons += 1
                    reversed_row = [[reverse_polynomial(value, length) for value in line]
                                    for line in full['row']]
                    for degree, rows in enumerate(row_jet_coefficients(
                            length, embedding, retraction, length+2)):
                        expected = [sum(((value >> degree) & 1) << j for j, value in enumerate(row))
                                    for row in reversed_row]
                        assert rows == expected
                        coefficient_comparisons += 1
                    assert all(full['column'][i][j] == full['column'][j][i]
                               for i in range(dimension**2) for j in range(dimension**2))
                    assert all(full['column'][dimension*x+u][dimension*x+u] == int(x == u)
                               for x in range(dimension) for u in range(dimension))
                    profiles = linear_restriction_profiles(length, embedding, retraction)
                    column_counts = Counter(profiles['column_zero_orders'])
                    assert column_counts[0] == 1 and column_counts[1] == 2*length
                    assert all(count % 2 == 0 for order, count in column_counts.items() if order % 2)
                    assert sum(profiles['row_infinity_orders']) % 2 == (length*dimension**2) % 2
                    profile_comparisons += 1
                    parity_controls.append({'sites': length, 'physical_dimension': dimension,
                                            'instance': instance,
                                            'column_counts': sorted(column_counts.items()),
                                            'row_reversal_order_sum': sum(profiles['row_infinity_orders'])})
    higher = []
    for length, dimension, degree in ((4, 13, 2), (4, 14, 2), (4, 15, 2), (4, 16, 3),
                                       (5, 25, 2), (5, 28, 2), (5, 31, 3), (5, 32, 3)):
        size = 1 << length
        assert ((degree == 2 and 4*dimension > 3*size)
                or 2*(size-dimension)*(1 << degree) < size)
        for instance in range(2):
            embedding, retraction = graph_retraction(length, dimension, rng,
                                                    coordinate=instance == 0)
            compressed = compressed_clifford_basis(length, embedding, retraction)
            low = [value for mask, value in enumerate(compressed) if mask.bit_count() <= degree]
            assert bit_rank(low) == sum(comb(2*length, j) for j in range(degree+1))
            jet = column_jet_coefficients(length, embedding, retraction, degree+1)
            counts, ranks = jet_order_multiplicities(jet)
            assert counts == tuple(comb(2*length, j) for j in range(degree+1))
            row_degree = max(j for j in range(degree+1)
                             if 2*(size-dimension)*(1 << j) < size)
            row_counts, row_ranks = jet_order_multiplicities(
                row_jet_coefficients(length, embedding, retraction, row_degree+1))
            assert row_counts == tuple(comb(2*length, j) for j in range(row_degree+1))
            higher.append({'sites': length, 'physical_dimension': dimension,
                           'instance': instance, 'last_resolved_order': degree,
                           'column_initial_counts': counts, 'toeplitz_ranks': ranks,
                           'last_resolved_row_order': row_degree,
                           'row_initial_counts': row_counts, 'row_toeplitz_ranks': row_ranks})
    half = []
    for length in range(2, 5):
        dimension = 1 << (length-1)
        embedding = retraction = [1 << j for j in range(dimension)]
        counts, ranks = jet_order_multiplicities(
            column_jet_coefficients(length, embedding, retraction, 2))
        assert counts == (1, 2*(length-1))
        half.append({'sites': length, 'physical_dimension': dimension,
                     'column_initial_counts': counts, 'toeplitz_ranks': ranks})
    three_quarters = []
    for length in range(2, 5):
        embedding = retraction = [1 << j for j in range(1 << length) if j & 3 != 3]
        compressed = compressed_clifford_basis(length, embedding, retraction)
        assert compressed[5] == 0  # The operator e_1 e_2 vanishes on its kernel.
        counts, ranks = jet_order_multiplicities(
            column_jet_coefficients(length, embedding, retraction, 3))
        assert counts == (1, 2*length, comb(2*length, 2)-2)
        three_quarters.append({'sites': length, 'physical_dimension': len(embedding),
                               'column_initial_counts': counts, 'toeplitz_ranks': ranks,
                               'vanishing_nonzero_operator_mask': 5})
    return {'first_layer_retractions': len(records), 'first_layer_records': records,
            'full_tensor_coefficient_matrices_compared': coefficient_comparisons,
            'complete_smith_profiles_compared': profile_comparisons,
            'complete_profile_parity_controls': parity_controls,
            'higher_layer_records': higher, 'half_dimension_countercontrols': half,
            'three_quarter_dimension_countercontrols': three_quarters}


def pressure_case(length, dimension, last_order, row_last_order=None):
    counts = [comb(2*length, j) for j in range(last_order+1)]
    fixed, spent = sum(counts), sum(j*c for j, c in enumerate(counts))
    size = dimension**2
    assert fixed <= size
    remaining = dimension
    prefix = 0
    for order, count in enumerate(counts):
        take = min(remaining, count)
        prefix += order*take
        remaining -= take
    prefix += (last_order+1)*remaining
    budget = dimension*length-prefix
    assert prefix >= length and 0 <= budget <= (dimension-1)*length
    # Verify the exponents controlling every permitted prefix-budget step.
    for previous in range(budget):
        assert (size*length-previous-1)//(size-dimension) >= previous//(dimension-1)
    shared_row = profile_polynomial(dimension, length, budget)
    row_counts = None
    if row_last_order is None:
        row = shared_row
    else:
        row_counts = [comb(2*length, j) for j in range(row_last_order+1)]
        row_fixed = sum(row_counts)
        row_spent = sum(j*c for j, c in enumerate(row_counts))
        if row_fixed == size:
            assert row_spent == size*length
            row = row_counts
        else:
            row = add(row_counts, balanced_polynomial(size-row_fixed, size*length-row_spent))
    if fixed == size:
        assert spent == size*length
        column = counts
    else:
        column = add(counts, balanced_polynomial(size-fixed, size*length-spent))
    bound = threshold(dimension)
    polynomial = [bound.denominator*value for value in multiply(row, column)]
    polynomial[length] -= bound.numerator
    leaves = positive_certificate(polynomial)
    return {'sites': length, 'physical_dimension': dimension,
            'last_proved_column_order': last_order, 'column_initial_counts': counts,
            'last_proved_row_order': row_last_order, 'row_initial_counts': row_counts,
            'column_prefix_lower_bound': prefix, 'row_prefix_budget': budget,
            'row_bound_used': 'shared minor' if row_last_order is None else 'forced row multiplicities',
            'row_lower_polynomial': row, 'column_lower_polynomial': column,
            'pressure_threshold': str(bound), 'positive_polynomial': polynomial,
            'bernstein_leaves': leaves}


def run():
    rng = Random(20261101)
    operators = operator_controls(rng)
    restrictions = restriction_controls(rng)
    ranges = {1: range(2, 3), 2: range(3, 5), 3: range(5, 9), 4: range(9, 17),
              5: range(25, 33), 6: range(49, 65), 7: (97, 98)+tuple(range(113, 129))}
    cases = []
    for length, dimensions in ranges.items():
        size = 1 << length
        for dimension in dimensions:
            row_degree = max(j for j in range(2*length+1)
                             if (1 << (j+1))*(size-dimension) < size)
            column_degree = max(row_degree, 1, 2 if 4*dimension > 3*size else 0)
            cases.append((length, dimension, column_degree, row_degree if row_degree else None))
    certificates = [pressure_case(*case) for case in cases]
    assert 2*577**2 > 816**2
    sources = ('tools/clifford_compression.py', 'tools/check_clifford_compression.py',
               'tools/linear_physical_restriction.py', 'tools/check_integer_shared_minor_floor.py')
    return {'status': 'pass',
            'scope': 'Constant elementary-tensor physical retractions in characteristic two, without additional gauges; finite controls over F2.',
            'first_layer_theorem': 'q>2^(t-1) implies column multiplicities m0=1,m1=2t',
            'second_layer_theorem': 'q>3*2^(t-2) implies m0=1,m1=2t,m2=binom(2t,2)',
            'higher_layer_theorem': '2(2^t-q)<2^(t-d) implies mj=binom(2t,j) for j<=d',
            'row_layer_theorem': 'The same generic threshold forces both row-reversal and column orders through d.',
            'near_full_asymptotic_floor': 'codimension n^(beta+o(1)) gives liminf rho>=h2(min((1-beta)/2,1/4))',
            'certified_rho_floor': 'rho>577/816>1/sqrt(2)',
            'pressure_cases': len(certificates), 'operator_controls': operators,
            'restriction_controls': restrictions, 'pressure_certificates': certificates,
            'floating_arithmetic_used': False, 'factoring_exponent_improved': False,
            'source_sha256': {path: sha256((ROOT/path).read_bytes()).hexdigest() for path in sources}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'evidence/clifford-compression.json')
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('status', 'first_layer_theorem', 'second_layer_theorem',
                                                  'higher_layer_theorem', 'pressure_cases',
                                                  'operator_controls', 'certified_rho_floor')}, indent=2))
