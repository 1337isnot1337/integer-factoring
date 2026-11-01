#!/usr/bin/env python3
"""Exact controls for the all-type elementary pressure floor over F_2."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb, factorial
from pathlib import Path

import check_block_smith_crossings as smith
from check_cancellation_module_automata import pmul


def one_plus_powers(bound):
    values = [1]
    for _ in range(bound):
        values.append(pmul(values[-1], 3))
    return values


def word_masks(length, weight):
    return tuple(sum(1 << i for i in selected)
                 for selected in combinations(range(length), weight))


def multiplicities(length, weight):
    numerator = factorial(length)
    row = [numerator // (factorial(weight-j)*factorial(length-weight-j)*factorial(2*j))
           for j in range(weight+1)]
    column = [numerator // (factorial(j)**2*factorial(length-2*j))
              for j in range(weight+1)]
    return row, column


def typed_matrices(length, weight):
    words = word_masks(length, weight)
    q = len(words)
    size = q*q
    row = [[0]*size for _ in range(size)]
    column = [[0]*size for _ in range(size)]
    powers = one_plus_powers(length)
    row_groups = Counter((x & y, x ^ y) for x, y in product(words, repeat=2))
    column_groups = Counter((x ^ u, x & (x ^ u)) for x, u in product(words, repeat=2))
    row_counts = Counter((difference.bit_count()//2) for _, difference in row_groups)
    column_counts = Counter((difference.bit_count()//2) for difference, _ in column_groups)
    expected_row, expected_column = multiplicities(length, weight)
    assert [row_counts[j] for j in range(weight+1)] == expected_row
    assert [column_counts[j] for j in range(weight+1)] == expected_column
    for (_, difference), block_size in row_groups.items():
        j = difference.bit_count()//2
        assert block_size == comb(2*j, j)
    for (difference, _), block_size in column_groups.items():
        j = difference.bit_count()//2
        assert block_size == comb(length-2*j, weight-j)
    entries = 0
    for xi, yi, ui, vi in product(range(q), repeat=4):
        x, y, u, v = words[xi], words[yi], words[ui], words[vi]
        value = 1
        for bit in range(length):
            xb, yb, ub, vb = ((word >> bit) & 1 for word in (x, y, u, v))
            straight = int(xb == ub and yb == vb)
            crossing = int(xb == vb and yb == ub)
            value = pmul(value, straight | ((straight ^ crossing) << 1))
        assert value.bit_count() % 2 == int(x == v and y == u)
        difference = x ^ y
        if (x & y, difference) == (u & v, u ^ v):
            j = difference.bit_count()//2
            distance = (x & difference & ~u).bit_count()
            predicted_row = powers[2*j-2*distance] << (2*distance)
        else:
            predicted_row = 0
        assert predicted_row == value
        difference = x ^ u
        if ((y ^ v) == difference and (v & difference) == (x & difference)):
            j = difference.bit_count()//2
            distance = (x & ~difference & ~y).bit_count()
            predicted_column = powers[2*distance] << (2*j)
        else:
            predicted_column = 0
        assert predicted_column == value
        row[q*xi+yi][q*ui+vi] = value
        column[q*xi+ui][q*yi+vi] = value
        entries += 1
    assert max(value.bit_length()-1 for line in row for value in line) == 2*weight
    return row, column, entries


def weighted_sum(orders, scale, degree=0):
    return sum((scale**(degree-order) for order in orders), Fraction(0))


def check_small_type(length, weight):
    row, column, entries = typed_matrices(length, weight)
    degree = 2*weight
    q = comb(length, weight)
    profiles = []
    for matrix in (row, column):
        zero, _ = smith.smith_profile(matrix, degree)
        reversed_matrix = [[smith.reverse_polynomial(value, degree) for value in line]
                           for line in matrix]
        infinity, _ = smith.smith_profile(reversed_matrix, degree)
        profiles.append((zero, infinity))
    assert profiles[0][0] == (0,)*(q*q)
    if length == 2*weight:
        assert profiles[0][1] == profiles[1][0]
        assert profiles[1][1] == (0,)*(q*q)
        assert all(left+right == 4*weight
                   for left, right in zip(profiles[0][1], reversed(profiles[0][1])))
    row_counts, column_counts = multiplicities(length, weight)
    pressure_checks = []
    for scale in (Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3)):
        actual = weighted_sum(profiles[0][1], scale, degree)*weighted_sum(profiles[1][0], scale)
        block_floor = (sum(count*scale**(2*j) for j, count in enumerate(row_counts))
                       *sum(count*scale**(-2*j) for j, count in enumerate(column_counts)))
        opposite = weighted_sum(profiles[1][1], scale, degree)*weighted_sum(profiles[0][0], scale)
        assert actual >= block_floor >= q**3
        assert opposite >= q**4
        pressure_checks.append({
            'scale': str(scale), 'actual_product': str(actual),
            'block_floor': str(block_floor), 'uniform_floor': q**3,
        })
    return {
        'length': length, 'weight': weight, 'alphabet': q,
        'degree': degree, 'tensor_entries_checked': entries,
        'row_infinity_counts': sorted(Counter(profiles[0][1]).items()),
        'column_zero_counts': sorted(Counter(profiles[1][0]).items()),
        'pressure_checks': pressure_checks,
    }


def middle_matrix(rank):
    words = word_masks(2*rank, rank)
    powers = one_plus_powers(rank)
    return [[powers[(left & right).bit_count()] for right in words] for left in words]


def partial_smith_orders(matrix, precision):
    """Find every order below precision; a zero remainder stays unclassified.

    This is exact elimination over F_2[t]/(t^precision), with a small
    multiplication table. The caller uses the proved middle-slice duality
    to complete the profile, rather than treating a zero jet as infinity.
    """
    modulus = 1 << precision
    mask = modulus-1
    table = [[pmul(left, right) & mask for right in range(modulus)]
             for left in range(modulus)]
    inverse = {left: next(right for right in range(1, modulus, 2)
                          if table[left][right] == 1)
               for left in range(1, modulus, 2)}
    valuation = [precision] + [(value & -value).bit_length()-1
                               for value in range(1, modulus)]
    work = [[value & mask for value in row] for row in matrix]
    dimension = len(work)
    orders = []
    for pivot in range(dimension):
        order = precision
        selected_row = selected_column = None
        for i in range(pivot, dimension):
            for j in range(pivot, dimension):
                candidate = valuation[work[i][j]]
                if candidate < order:
                    order, selected_row, selected_column = candidate, i, j
                if order == 0:
                    break
            if order == 0:
                break
        if selected_row is None:
            break
        orders.append(order)
        work[pivot], work[selected_row] = work[selected_row], work[pivot]
        if selected_column != pivot:
            for row in work:
                row[pivot], row[selected_column] = row[selected_column], row[pivot]
        pivot_row = work[pivot]
        times = table[inverse[pivot_row[pivot] >> order]]
        for j in range(pivot, dimension):
            pivot_row[j] = times[pivot_row[j]]
        assert pivot_row[pivot] == 1 << order
        support = [j for j in range(pivot+1, dimension) if pivot_row[j]]
        for i in range(pivot+1, dimension):
            row = work[i]
            factor = row[pivot] >> order
            if factor:
                times = table[factor]
                for j in support:
                    row[j] ^= times[pivot_row[j]]
            row[pivot] = 0
    assert orders == sorted(orders)
    return orders


def middle_profile(rank):
    matrix = middle_matrix(rank)
    if rank == 0:
        return {'rank': 0, 'dimension': 1, 'profile_counts': [[0, 1]]}
    low = partial_smith_orders(matrix, rank)
    dimension = len(matrix)
    assert 2*len(low) <= dimension
    profile = low + [rank]*(dimension-2*len(low)) + [2*rank-order for order in reversed(low)]
    assert sum(profile) == rank*dimension
    if rank <= 4:
        independent, _ = smith.smith_profile(matrix, rank)
        assert tuple(profile) == independent
    expected_low = Counter({i: comb(2*rank-1, i)-(comb(2*rank-1, i-1) if i else 0)
                            for i in range(rank)})
    assert Counter(low) == expected_low
    return {
        'rank': rank, 'dimension': dimension,
        'visible_orders_below_rank': len(low),
        'profile_counts': sorted(Counter(profile).items()),
        'finite_multiplicity_pattern_matches': True,
        'completion_uses_proved_square_identity': True,
    }


def coefficient_product(factors, selected_degree):
    coefficients = [1] + [0]*selected_degree
    for factor in factors:
        updated = [0]*(selected_degree+1)
        for i, value in enumerate(coefficients):
            if value:
                for j, coefficient in enumerate(factor):
                    if i+j <= selected_degree and coefficient:
                        updated[i+j] ^= pmul(value, coefficient)
        coefficients = updated
    return coefficients[selected_degree]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-johnson', type=int, default=6)
    parser.add_argument('--output', default='evidence/direct-type-elementary-floor.json')
    args = parser.parse_args()
    if not 1 <= args.max_johnson <= 6:
        raise ValueError('The finite replay supports middle ranks one through six.')
    records = [check_small_type(length, weight)
               for length in range(2, 6) for weight in range(1, length//2+1)]
    integer_cases = 0
    for length in range(2, 129):
        for weight in range(1, length//2+1):
            q = comb(length, weight)
            row, column = multiplicities(length, weight)
            p = [comb(weight, j)*comb(length-weight, j) for j in range(weight+1)]
            assert sum(p) == q
            for j, value in enumerate(p):
                assert row[j]*column[j] == q*value*comb(length, 2*j)
                assert comb(length, 2*j) >= value
                assert row[j]*column[j] >= q*value*value
            assert sum(row[j]*comb(2*j, j) for j in range(weight+1)) == q*q
            assert sum(column[j]*comb(length-2*j, weight-j)
                       for j in range(weight+1)) == q*q
            integer_cases += 1

    coefficient_cases = 0
    for rank in range(33):
        for intersection in range(rank+1):
            original = ([[1, 5]]*intersection
                        +[[1, 3]]*(2*(rank-intersection))
                        +[[1, 1]]*intersection)
            simplified = ([[1, 4, 5]]*intersection
                          +[[1, 0, 5]]*(rank-intersection))
            actual = coefficient_product(original, rank)
            assert actual == coefficient_product(simplified, rank)
            assert actual == ((1 << (2*rank)) if intersection == rank else 0)
            coefficient_cases += 1
    matrix_square_entries = 0
    for rank in range(4):
        matrix = middle_matrix(rank)
        dimension = len(matrix)
        for i, j in product(range(dimension), repeat=2):
            value = 0
            for h in range(dimension):
                value ^= pmul(matrix[i][h], matrix[h][j])
            assert value == ((1 << (2*rank)) if i == j else 0)
            matrix_square_entries += 1
    profiles = [middle_profile(rank) for rank in range(args.max_johnson+1)]
    root = Path(__file__).resolve().parents[1]
    report = {
        'status': 'pass',
        'scope': 'Uniform floor for full common binary types of the elementary seed; no factoring lower bound.',
        'proved_pressure_floor': '3/4',
        'half_threshold_separation': '9/16 > 1/2',
        'typed_tensor_cases': records,
        'tensor_entries_checked': sum(record['tensor_entries_checked'] for record in records),
        'integer_binomial_cases': integer_cases,
        'square_identity_coefficient_cases': coefficient_cases,
        'square_identity_matrix_entries': matrix_square_entries,
        'middle_slice_profiles': profiles,
        'general_middle_multiplicity_formula_proved': False,
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in (
        'status', 'proved_pressure_floor', 'tensor_entries_checked',
        'integer_binomial_cases', 'square_identity_coefficient_cases',
        'square_identity_matrix_entries',
    )}, indent=2))


if __name__ == '__main__':
    main()
