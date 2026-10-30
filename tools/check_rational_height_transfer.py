#!/usr/bin/env python3
"""Exact finite checks for the rational degree/height transfer.

The asymptotic argument is in research/rational-height-transfer.md. These
checks do not construct a faster crossing family or run a fast factorer.
"""

from fractions import Fraction
import json
from math import prod
from random import Random


def polynomial_add(left, right):
    result = left.copy()
    for monomial, value in right.items():
        result[monomial] = result.get(monomial, 0) + value
        if not result[monomial]:
            del result[monomial]
    return result


def polynomial_mul(left, right):
    result = {}
    for a, av in left.items():
        for b, bv in right.items():
            monomial = tuple(x+y for x, y in zip(a, b))
            result[monomial] = result.get(monomial, 0) + av*bv
    return {m: v for m, v in result.items() if v}


def check_circuit(scalars, operations, expand=True):
    """Every scalar input, including zero and one, has syntactic degree one."""
    scalars = [Fraction(value) for value in scalars]
    count = len(scalars)
    values = scalars.copy()
    degrees = [1]*count
    polynomials = []
    for i in range(count):
        monomial = tuple(int(i == j) for j in range(count))
        polynomials.append({monomial: 1})
    for operation, a, b in operations:
        assert 0 <= a < len(values) and 0 <= b < len(values)
        if operation == 'add':
            values.append(values[a]+values[b])
            degrees.append(max(degrees[a], degrees[b]))
            if expand:
                polynomials.append(polynomial_add(polynomials[a], polynomials[b]))
        elif operation == 'mul':
            values.append(values[a]*values[b])
            degrees.append(degrees[a]+degrees[b])
            if expand:
                polynomials.append(polynomial_mul(polynomials[a], polynomials[b]))
        else:
            raise ValueError(operation)
    size = len(values)
    degree = max(degrees)
    input_bits = max(max(abs(value.numerator).bit_length(), value.denominator.bit_length(), 1)
                     for value in scalars)
    height_bound = 1+degree*((size+1)*input_bits+size)
    common = prod(value.denominator for value in scalars)
    for index, (value, local_degree) in enumerate(zip(values, degrees), start=1):
        assert max(abs(value.numerator).bit_length(), value.denominator.bit_length()) <= height_bound
        assert pow(common, local_degree, value.denominator) == 0
        if expand:
            polynomial = polynomials[index-1]
            assert all(sum(monomial) <= local_degree for monomial in polynomial)
            norm = sum(abs(coefficient) for coefficient in polynomial.values())
            assert norm <= 2**(local_degree*index)
            reconstructed = sum(Fraction(coefficient)*prod(x**e for x, e in zip(scalars, monomial))
                                for monomial, coefficient in polynomial.items())
            assert reconstructed == value
    return {
        'nodes': size,
        'syntactic_degree': degree,
        'input_bits': input_bits,
        'largest_output_bits': max(max(abs(x.numerator).bit_length(), x.denominator.bit_length())
                                   for x in values),
        'proved_height_bound': height_bound,
    }


def schedule_checks():
    stopping_cases = local_cases = label_cases = 0
    for x0 in (18, 19, 32):
        logs = [x0]
        for j in range(12):
            logs.append(2*logs[-1]+6*j+2*x0)
        for ell in range(13):
            assert 2*sum(logs[1:ell+1]) <= 4*logs[ell]
            assert 2*ell <= logs[ell]
            for stop in range(ell+1):
                degree_log = logs[stop]+sum(2+2*x for x in logs[stop+1:ell+1])
                assert degree_log <= 6*logs[ell]
                stopping_cases += 1
            # rho >= 1/2 bounds rho^(-ell) by 2^ell. This integer
            # comparison verifies the displayed label-width envelope.
            assert ell-2*logs[ell] <= -2*x0
            label_cases += 1
    for g in range(2, 65):
        # A scalar-one product accumulator is charged as a degree-one leaf.
        raw_degree = g*g*(g+g*g+1)+1
        assert raw_degree <= 4*g**4
        local_cases += 1
    return stopping_cases, local_cases, label_cases


def matmul(left, right):
    if not left:
        return []
    columns = list(zip(*right))
    return [[sum(a*b for a, b in zip(row, column)) for column in columns] for row in left]


def rank_factor(matrix):
    """Exact pivot columns with coefficient rows found by rational elimination."""
    rows = len(matrix)
    columns = len(matrix[0]) if rows else 0
    reduced = [[Fraction(value) for value in row] for row in matrix]
    pivot_columns = []
    rank = 0
    for column in range(columns):
        pivot = next((r for r in range(rank, rows) if reduced[r][column]), None)
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        lead = reduced[rank][column]
        reduced[rank] = [value/lead for value in reduced[rank]]
        for r in range(rows):
            if r == rank or not reduced[r][column]:
                continue
            coefficient = reduced[r][column]
            reduced[r] = [a-coefficient*b for a, b in zip(reduced[r], reduced[rank])]
        pivot_columns.append(column)
        rank += 1
        if rank == rows:
            break
    left = [[Fraction(row[j]) for j in pivot_columns] for row in matrix]
    right = reduced[:rank]
    if rank:
        assert matmul(left, right) == matrix
    else:
        assert all(not value for row in matrix for value in row)
    return rank, left, right


def factorization_checks():
    rng = Random(8022026)
    cases = entries = 0
    for dimension in range(1, 8):
        for rank in range(dimension+1):
            if rank:
                left = [[Fraction(int(i == j)) if i < rank else
                         Fraction(rng.randrange(-20, 21), rng.randrange(1, 17))
                         for j in range(rank)] for i in range(dimension)]
                right = [[Fraction(int(i == j)) if j < rank else
                          Fraction(rng.randrange(-20, 21), rng.randrange(1, 17))
                          for j in range(dimension)] for i in range(rank)]
                matrix = matmul(left, right)
            else:
                matrix = [[Fraction(0)]*dimension for _ in range(dimension)]
            actual, factor_left, factor_right = rank_factor(matrix)
            assert actual == rank
            cases += 1
            entries += dimension**2
            if actual:
                corrupted = [row.copy() for row in factor_right]
                corrupted[0][0] += 1
                assert matmul(factor_left, corrupted) != matrix
    return cases, entries


def verify():
    rng = Random(7302026)
    records = []
    for case in range(18):
        scalars = [1, 0, Fraction(-3-case, 5+2*case), Fraction(7+case, 11+case)]
        degrees = [1]*len(scalars)
        operations = []
        for _ in range(24):
            a, b = rng.randrange(len(degrees)), rng.randrange(len(degrees))
            operation = 'mul' if rng.randrange(2) and degrees[a]+degrees[b] <= 12 else 'add'
            operations.append((operation, a, b))
            degrees.append(degrees[a]+degrees[b] if operation == 'mul' else max(degrees[a], degrees[b]))
        records.append(check_circuit(scalars, operations))

    # Shared gates repeatedly double the coefficient norm without increasing
    # degree. This checks the S term in the height bound.
    doubling = [('add', j, j) for j in range(32)]
    records.append(check_circuit([Fraction(1, 17)], doubling))

    # Repeated squaring is the size-only and degree-zero-constant control.
    squaring = [('add', 0, 0)]
    for _ in range(12):
        current = len(squaring)
        squaring.append(('mul', current, current))
    square_record = check_circuit([1], squaring, expand=False)
    assert square_record['syntactic_degree'] == 2**12
    assert square_record['largest_output_bits'] == 2**12+1
    assert square_record['largest_output_bits'] > square_record['nodes']**2
    records.append(square_record)

    stopping, local, labels = schedule_checks()
    factor_cases, factor_entries = factorization_checks()
    return {
        'scope': 'finite exact degree, rational-height, and rank-factor checks; no improved crossing profile',
        'rational_circuits': records,
        'schedule_stopping_cases': stopping,
        'grouped_query_degree_cases': local,
        'label_envelope_cases': labels,
        'exact_rank_factorizations': factor_cases,
        'reconstructed_matrix_entries': factor_entries,
        'negative_controls': [
            'repeated squaring defeats a size-only or degree-zero-constant height estimate',
            'perturbing a nonzero factor changes the reconstructed matrix',
        ],
        'new_factoring_exponent_claimed': False,
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
