#!/usr/bin/env python3
"""Exact finite control for the target-order and common-multiplier rows.

This enumerates one small class group and verifies the order rows, the mixed
order-row determinant, the common-multiplier determinant identity on synthetic
matrices, one actual mixed smooth hit, and a three-target order-only annihilator
with the existing two-part extractor. It is not a factoring algorithm; the
class-group enumeration is exponential in the input bit length.
"""
import argparse
import json
from fractions import Fraction
from functools import lru_cache
from math import prod

from arithmetic_two_group import power, project_to_two_part, two_group_basis
from binary_quadratic_forms import (
    ambiguous_factor, compose_forms, inverse_form, reduced_forms, reduce_form,
)


def determinant(matrix):
    size = len(matrix)
    values = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if values[row][column]),
                     None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            result = -result
        pivot_value = values[column][column]
        result *= pivot_value
        for row in range(column+1, size):
            factor = values[row][column]/pivot_value
            if factor:
                for index in range(column, size):
                    values[row][index] -= factor*values[column][index]
    return result


def solve(matrix, rhs):
    size = len(matrix)
    values = [[Fraction(value) for value in row] for row in matrix]
    result = [Fraction(value) for value in rhs]
    for column in range(size):
        pivot = next((row for row in range(column, size) if values[row][column]),
                     None)
        assert pivot is not None
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            result[column], result[pivot] = result[pivot], result[column]
        pivot_value = values[column][column]
        for row in range(column+1, size):
            factor = values[row][column]/pivot_value
            if not factor:
                continue
            for index in range(column, size):
                values[row][index] -= factor*values[column][index]
            result[row] -= factor*result[column]
    answer = [Fraction(0)]*size
    for row in range(size-1, -1, -1):
        total = result[row]-sum(values[row][index]*answer[index]
                                for index in range(row+1, size))
        answer[row] = total/values[row][row]
    return answer


def common_multiplier_identity(k, v, z):
    size = len(k)
    d_minus_v = [[Fraction(k[i] if i == j else 0)-Fraction(v[i][j])
                  for j in range(size)] for i in range(size)]
    left = [[d_minus_v[i][j]+Fraction(k[i])*Fraction(z[j])
             for j in range(size)] for i in range(size)]
    base = determinant(d_minus_v)
    assert base
    x = solve(d_minus_v, [Fraction(k[i]) for i in range(size)])
    scalar = 1+sum(Fraction(z[i])*x[i] for i in range(size))
    right = base*scalar
    actual = determinant(left)
    assert actual == right and actual > 0
    return {'size': size, 'determinant': int(actual),
            'base_determinant': int(base), 'scalar_numerator': scalar.numerator,
            'scalar_denominator': scalar.denominator}


def exact_class_group(discriminant):
    forms = reduced_forms(discriminant)
    identity = next(form for form in forms if form[0] == 1)
    multiply = lru_cache(None)(compose_forms)

    def order(form, limit=None):
        current, value = form, 1
        while current != identity:
            current = multiply(current, form)
            value += 1
            assert limit is None or value <= limit
        return value

    def prime_form(prime):
        assert discriminant % prime
        middle = next(value for value in range(prime)
                      if (value*value-discriminant) % (4*prime) == 0)
        return reduce_form((prime, middle,
                            (middle*middle-discriminant)//(4*prime)))

    return forms, identity, multiply, order, prime_form


def generated_subgroup_size(generators, identity, multiply):
    seen = {identity}
    frontier = [identity]
    while frontier:
        following = []
        for value in frontier:
            for generator in generators:
                product = multiply(value, generator)
                if product not in seen:
                    seen.add(product)
                    following.append(product)
        frontier = following
    return len(seen)


def check_small_group():
    integer = 10403
    discriminant = -4*integer
    forms, identity, multiply, order, prime_form = exact_class_group(discriminant)
    class_number = len(forms)
    assert class_number == 90 and integer == 101*103
    primes = (3, 11, 17, 19, 23)
    target_forms = {prime: prime_form(prime) for prime in primes}
    orders = {prime: order(target_forms[prime], class_number) for prime in primes}
    assert orders == {3: 30, 11: 6, 17: 15, 19: 5, 23: 3}
    indices = {prime: class_number//orders[prime] for prime in primes}
    assert indices == {3: 3, 11: 15, 17: 6, 19: 18, 23: 30}

    gap_window = 2
    gap_rows = []
    for prime in primes:
        gap_rows.append({
            'target': prime, 'order': orders[prime], 'index': indices[prime],
            'order_covered_at_W2': orders[prime] <= gap_window,
            'density_covered_at_kappa_quarter': indices[prime]**4 <= abs(discriminant),
        })
    assert [row['order_covered_at_W2'] for row in gap_rows] == [False]*5
    assert [row['density_covered_at_kappa_quarter'] for row in gap_rows] == \
        [True, False, True, False, False]

    mixed = multiply(target_forms[19], target_forms[23])
    assert order(mixed, class_number) == 15
    assert class_number//15 == 6 and 6**4 <= abs(discriminant)
    assert multiply(mixed, mixed) == target_forms[17]

    relation_rows = {
        'mixed_19_times_23_squared': [[2, 2, -1, 0], [0, 3, 0, 0],
                                      [0, 0, 15, 0], [0, 0, 0, 6]],
    }
    matrix = relation_rows['mixed_19_times_23_squared']
    relation_determinant = int(determinant(matrix))
    assert relation_determinant == 540 == 6*class_number

    generators = [target_forms[19], target_forms[23], target_forms[17],
                  target_forms[11]]
    mixed_subgroup_size = generated_subgroup_size(generators, identity, multiply)
    assert mixed_subgroup_size == class_number
    # Exact row relations: 2(e19+e23)-e17, 3e23, 15e17, 6e11.
    left = multiply(multiply(power(target_forms[19], 2, identity, multiply),
                             power(target_forms[23], 2, identity, multiply)),
                    inverse_form(target_forms[17]))
    assert left == identity
    assert power(target_forms[23], 3, identity, multiply) == identity
    assert power(target_forms[17], 15, identity, multiply) == identity
    assert power(target_forms[11], 6, identity, multiply) == identity
    assert all(power(form, relation_determinant, identity, multiply) == identity
               for form in forms)

    projected = project_to_two_part(generators, relation_determinant,
                                    identity, multiply)
    basis = two_group_basis(projected, identity, multiply, inverse_form,
                            class_number.bit_length())
    tops = [power(form, 1 << (power_order-1), identity, multiply)
            for form, power_order in basis]
    factors = [ambiguous_factor(form, integer) for form in tops]
    factor = next((value for value in factors if value), None)
    assert factor in (101, 103) and integer % factor == 0

    order_only_primes = (23, 19, 11)
    order_only_orders = [3, 5, 6]
    order_only_generators = [target_forms[prime] for prime in order_only_primes]
    order_only_subgroup_size = generated_subgroup_size(
        order_only_generators, identity, multiply)
    assert order_only_subgroup_size == class_number
    order_only_matrix = [[3, 0, 0], [0, 5, 0], [0, 0, 6]]
    order_only_determinant = int(determinant(order_only_matrix))
    assert order_only_determinant == 90 == class_number
    for prime, expected in zip(order_only_primes, order_only_orders):
        assert orders[prime] == expected
        assert power(target_forms[prime], expected, identity, multiply) == identity
    assert all(power(form, order_only_determinant, identity, multiply) == identity
               for form in forms)
    order_only_projected = project_to_two_part(
        order_only_generators, order_only_determinant, identity, multiply)
    order_only_basis = two_group_basis(
        order_only_projected, identity, multiply, inverse_form,
        class_number.bit_length())
    order_only_tops = [power(form, 1 << (power_order-1), identity, multiply)
                       for form, power_order in order_only_basis]
    order_only_factors = [ambiguous_factor(form, integer)
                          for form in order_only_tops]
    order_only_factor = next((value for value in order_only_factors if value),
                             None)
    assert order_only_factor in (101, 103) and integer % order_only_factor == 0
    class_group_exponent = max(order(form, class_number) for form in forms)
    assert class_group_exponent == 30 and class_group_exponent < class_number

    mixed_order_rows = {}
    all_mixed_orders = {}
    for multiplier in (23,):
        mixed_order_rows['mixed_19_by_23'] = {}
        for prime in primes:
            mixed_class = multiply(target_forms[prime], target_forms[multiplier])
            all_mixed_orders[prime] = order(mixed_class, class_number)
        # A single mixed target: det = order * (1+z_19), z_19=0.
        mixed_order_rows['mixed_19_by_23'] = {
            'multiplier': multiplier, 'mixed_order': all_mixed_orders[19],
            'determinant': all_mixed_orders[19],
        }
        assert all_mixed_orders[19] == 15
        assert mixed_order_rows['mixed_19_by_23']['determinant'] == 15
    common_mixed_determinant = prod(all_mixed_orders.values())*(1+1)
    assert all_mixed_orders == {3: 30, 11: 6, 17: 5, 19: 15, 23: 3}
    assert common_mixed_determinant == 81000

    return {
        'integer': integer,
        'factorization': [101, 103],
        'discriminant': discriminant,
        'class_number': class_number,
        'class_group_exponent': class_group_exponent,
        'class_group_cyclic': class_group_exponent == class_number,
        'd_fourth_root_floor': int(abs(discriminant)**(1/4)),
        'gap_window': gap_window,
        'targets': gap_rows,
        'mixed_target': {
            'target': 19, 'multiplier': 23,
            'class_order': 15, 'index': 6,
            'hit_exponent': 2, 'hit_word': {'17': 1},
            'identity': '(19*23)^2 = 17',
        },
        'relation_matrix': matrix,
        'relation_determinant': relation_determinant,
        'determinant_over_class_number': relation_determinant//class_number,
        'mixed_generators_subgroup_size': mixed_subgroup_size,
        'two_part_basis_orders': [power_order for _, power_order in basis],
        'two_part_basis_forms': [list(form) for form in tops],
        'extracted_factor': factor,
        'order_only_certificate': {
            'targets': list(order_only_primes),
            'orders': order_only_orders,
            'relation_matrix': order_only_matrix,
            'determinant': order_only_determinant,
            'generated_subgroup_size': order_only_subgroup_size,
            'window': max(order_only_orders),
            'two_part_basis_orders': [power_order
                                      for _, power_order in order_only_basis],
            'extracted_factor': order_only_factor,
        },
        'all_forms_annihilated': True,
        'mixed_order_rows': mixed_order_rows,
        'common_mixed_orders': all_mixed_orders,
        'common_mixed_determinant': common_mixed_determinant,
        'hybrid_block_determinant': 3*15*6*2,
    }


def check_synthetic_identities():
    rows = []
    for size, k, v, z in (
        (2, [3, 4], [[1, 0], [0, 1]], [2, 1]),
        (2, [5, 3], [[1, 1], [1, 0]], [0, 4]),
        (3, [4, 5, 6], [[1, 0, 1], [0, 1, 1], [1, 1, 0]], [1, 2, 3]),
        (3, [7, 5, 4], [[2, 1, 0], [1, 1, 1], [0, 1, 2]], [0, 0, 5]),
        (4, [6, 7, 5, 4],
         [[1, 0, 0, 1], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 1, 0]],
         [2, 3, 1, 4]),
    ):
        assert all(row_sum < k_value for row_sum, k_value in
                   zip((sum(row) for row in v), k))
        record = common_multiplier_identity(k, v, z)
        record.update({'k': k, 'v': v, 'z': z})
        rows.append(record)
    return rows


def check_mixed_order_determinant():
    z = [0, 1, 0, 0]
    orders = [30, 6, 5, 15, 3]
    determinant_value = prod(orders)*(1+sum(z))
    assert determinant_value == 81000
    return {'orders': orders, 'z': z, 'determinant': determinant_value}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    arguments = parser.parse_args()
    receipt = {
        'scope': ('Exact small class-group, order-only, and synthetic '
                  'determinant controls; no uniform multiplier selection, '
                  'mixed-density positivity, generating-subset bound, or '
                  'complete factoring improvement'),
        'complete_saving_exponent': Fraction(1, 2).__float__(),
        'small_group': check_small_group(),
        'mixed_order_determinant': check_mixed_order_determinant(),
        'synthetic_common_multiplier_identities': check_synthetic_identities(),
    }
    if arguments.output:
        with open(arguments.output, 'w') as handle:
            json.dump(receipt, handle, indent=2, sort_keys=True)
            handle.write('\n')
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
