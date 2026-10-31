#!/usr/bin/env python3
"""Complete closed-point Smith profiles of the finite fused crossings.

Exact F16 Berlekamp factorization identifies every singular place. Minor
valuations give the profiles; multiplication modulo f^L independently
checks nonlinear places without choosing roots in a splitting extension.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import lcm
from pathlib import Path

from fused_six_vertex_crossings import construct, power, pdiv, gcd, pscale
from f16_identity_operator_profiles import MUL, INV, add, multiply
from check_fused_multipoint_profiles import determinant, cut_blocks, field_rank


def polynomial_modulo(value, modulus):
    return pdiv(value, modulus)[1]


def multiply_modulo(first, second, modulus):
    return polynomial_modulo(multiply(first, second), modulus)


def power_modulo(value, exponent, modulus):
    result = {0: 1}
    while exponent:
        if exponent & 1:
            result = multiply_modulo(result, value, modulus)
        value = multiply_modulo(value, value, modulus)
        exponent //= 2
    return result


def canonical(polynomial):
    return tuple(polynomial.get(k, 0) for k in range(max(polynomial) + 1))


def from_canonical(coefficients):
    return {k: value for k, value in enumerate(coefficients) if value}


def field_kernel(matrix):
    work = [row[:] for row in matrix]
    columns = len(work[0])
    pivots = []
    r = 0
    for column in range(columns):
        pivot = next((i for i in range(r, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        scalar = INV[work[r][column]]
        work[r] = [MUL[scalar][value] for value in work[r]]
        for i in range(len(work)):
            if i == r or not work[i][column]:
                continue
            scalar = work[i][column]
            work[i] = [a ^ MUL[scalar][b] for a, b in zip(work[i], work[r])]
        pivots.append(column)
        r += 1
        if r == len(work):
            break
    result = []
    for column in range(columns):
        if column in pivots:
            continue
        vector = [0] * columns
        vector[column] = 1
        for i, pivot in enumerate(pivots):
            vector[pivot] = work[i][column]
        for row in matrix:
            total = 0
            for a, b in zip(row, vector):
                total ^= MUL[a][b]
            assert total == 0
        result.append({i: value for i, value in enumerate(vector) if value})
    return result


def berlekamp_basis(polynomial):
    degree = max(polynomial)
    x16 = power_modulo({1: 1}, 16, polynomial)
    powers = [{0: 1}]
    for _ in range(1, degree):
        powers.append(multiply_modulo(powers[-1], x16, polynomial))
    matrix = [[powers[j].get(i, 0) ^ int(i == j) for j in range(degree)]
              for i in range(degree)]
    return field_kernel(matrix)


def irreducible_factors(polynomial):
    if max(polynomial) == 0:
        return []
    polynomial = pscale(polynomial, INV[polynomial[max(polynomial)]])
    derivative = {k - 1: value for k, value in polynomial.items() if k & 1}
    if not derivative:
        root = {k // 2: power(value, 8) for k, value in polynomial.items()}
        return irreducible_factors(root) * 2
    common = gcd(polynomial, derivative)
    if max(common) > 0:
        quotient, remainder = pdiv(polynomial, common)
        assert not remainder
        return irreducible_factors(common) + irreducible_factors(quotient)
    basis = berlekamp_basis(polynomial)
    if len(basis) == 1:
        return [polynomial]
    parts = [polynomial]
    for basis_element in basis:
        for constant in range(16):
            splitter = add(basis_element, {0: constant}) if constant else basis_element
            updated = []
            for part in parts:
                common = gcd(part, splitter)
                if common and 0 < max(common) < max(part):
                    quotient, remainder = pdiv(part, common)
                    assert not remainder
                    updated.extend((common, quotient))
                else:
                    updated.append(part)
            parts = updated
        if len(parts) == len(basis):
            break
    assert len(parts) == len(basis)
    return parts


@lru_cache(maxsize=None)
def verify_irreducible(coefficients):
    polynomial = from_canonical(coefficients)
    degree = max(polynomial)
    x = polynomial_modulo({1: 1}, polynomial)
    current = x
    for k in range(1, degree + 1):
        current = power_modulo(current, 16, polynomial)
        if k <= degree // 2:
            assert max(gcd(polynomial, add(current, x))) == 0
    assert current == x
    return True


def factor_valuation(polynomial, irreducible):
    order = 0
    while polynomial:
        quotient, remainder = pdiv(polynomial, irreducible)
        if remainder:
            return order
        polynomial = quotient
        order += 1
    raise ValueError('Zero polynomials have no finite valuation.')


def quotient_rank(matrix, irreducible, length):
    modulus = {0: 1}
    for _ in range(length):
        modulus = multiply(modulus, irreducible)
    dimension = max(modulus)
    size = len(matrix)
    result = [[0] * (size * dimension) for _ in range(size * dimension)]
    for a in range(size):
        for b in range(size):
            value = matrix[a][b]
            for j in range(dimension):
                remainder = polynomial_modulo({k + j: scalar for k, scalar in value.items()}, modulus)
                for i, scalar in remainder.items():
                    result[a * dimension + i][b * dimension + j] = scalar
    return field_rank(result)


def check_instance(m, theta):
    seed = construct(m, theta)
    assert seed['invariant'] and seed['endpoint_identity']
    alphabet, degree = m + 1, seed['degree']
    blocks = [cut_blocks(seed['matrix'], alphabet, cut) for cut in (0, 1)]
    factors = {}
    all_minors = []
    determinant_products = nonzero_minors = 0
    for cut in blocks:
        cut_minors = []
        for block in cut:
            size = len(block)
            family = []
            for k in range(1, size + 1):
                values = []
                for rows in combinations(range(size), k):
                    for columns in combinations(range(size), k):
                        value = determinant([[block[i][j] for j in columns] for i in rows])
                        if value:
                            values.append(value)
                assert values
                family.append(values)
                nonzero_minors += len(values)
            determinant_value = family[-1][0]
            product_value = {0: 1}
            for factor in irreducible_factors(determinant_value):
                key = canonical(factor)
                assert verify_irreducible(key)
                factors[key] = factor
                product_value = multiply(product_value, factor)
            assert pscale(product_value, determinant_value[max(determinant_value)]) == determinant_value
            determinant_products += 1
            cut_minors.append(family)
        all_minors.append(cut_minors)
    closed_points = []
    quotient_checks = 0
    for key, factor in sorted(factors.items(), key=lambda item: (len(item[0]), item[0])):
        profiles = []
        residue_degree = max(factor)
        for cut_index, cut in enumerate(all_minors):
            orders = []
            for block_index, family in enumerate(cut):
                previous = 0
                block_orders = []
                for values in family:
                    current = min(factor_valuation(value, factor) for value in values)
                    block_orders.append(current - previous)
                    previous = current
                assert block_orders == sorted(block_orders)
                if residue_degree > 1:
                    for length in range(1, max(block_orders) + 2):
                        actual = quotient_rank(blocks[cut_index][block_index], factor, length)
                        expected = residue_degree * sum(max(0, length - order) for order in block_orders)
                        assert actual == expected
                        quotient_checks += 1
                orders.extend(block_orders)
            profiles.append(sorted(orders))
        closed_points.append({'factor': list(key), 'residue_degree': residue_degree,
                              'cuts': profiles})
    infinity_profiles = []
    for cut in all_minors:
        orders = []
        for family in cut:
            previous = 0
            for k, values in enumerate(family, 1):
                current = min(degree * k - max(value) for value in values)
                orders.append(current - previous)
                previous = current
        infinity_profiles.append(sorted(orders))
    expanded_profiles = []
    for row in closed_points:
        for conjugate in range(row['residue_degree']):
            expanded_profiles.append({'point': f"root:{row['factor']}:{conjugate}",
                                      'cuts': row['cuts']})
    expanded_profiles.append({'point': None, 'cuts': infinity_profiles})
    totals = [sum(sum(row['cuts'][cut]) for row in expanded_profiles) for cut in (0, 1)]
    assert totals == [alphabet * alphabet * degree] * 2
    return {'m': m, 'theta': theta, 'alphabet': alphabet, 'degree': degree,
            'splitting_extension_degree_over_f16': lcm(*(row['residue_degree'] for row in closed_points)),
            'closed_points': closed_points, 'profiles': expanded_profiles,
            'controls': {'exact_determinant_products': determinant_products,
                         'nonzero_minors': nonzero_minors,
                         'nonlinear_quotient_ranks': quotient_checks,
                         'complete_determinant_degree_ledger': totals}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [check_instance(m, theta) for m, theta in ((2, 2), (3, 4), (4, 4), (5, 4))]
    directory = Path(__file__).parent
    sources = ['check_fused_closed_profiles.py', 'check_fused_multipoint_profiles.py',
               'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py',
               'f4_lattice_linear.py']
    totals = {'instances': len(cases),
              'closed_points': sum(len(case['closed_points']) for case in cases),
              'split_places_including_infinity': sum(len(case['profiles']) for case in cases),
              'exact_determinant_products': sum(case['controls']['exact_determinant_products'] for case in cases),
              'nonzero_minors': sum(case['controls']['nonzero_minors'] for case in cases),
              'nonlinear_quotient_ranks': sum(case['controls']['nonlinear_quotient_ranks'] for case in cases)}
    result = {'schema_version': 1,
              'claim': 'Complete finite closed-point profiles; no improved factoring exponent.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest()
                                for name in sources}, 'cases': cases, 'totals': totals}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(totals, sort_keys=True))


if __name__ == '__main__':
    main()
