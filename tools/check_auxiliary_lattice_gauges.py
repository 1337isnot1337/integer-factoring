#!/usr/bin/env python3
"""Exact canonical-lattice and local-minimum controls for auxiliary gauges."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random

import check_smith_profile_transfer as basic


def canonical_forms(depth):
    result = []
    for h in range(depth+1):
        for a in range(h+1):
            b = h-a
            for p in range(1 << b):
                if a == 0 or b == 0 or p & 1:
                    result.append((a, b, p))
    return result


def primitive(a, b, p):
    order = min(a, b, basic.valuation(p) if p else a+b+1)
    return a-order, b-order, p >> order


def neighbors(form):
    a, b, p = form
    return [primitive(a, b+1, p),
            primitive(a, b+1, p ^ (1 << b)),
            primitive(a+1, b, (p << 1) & ((1 << b)-1))]


def gauge_matrix(form):
    a, b, p = form
    return [[1 << a, p], [0, 1 << b]]


def product(left, right):
    return [[basic.multiply(left[a][0], right[0][b]) ^
             basic.multiply(left[a][1], right[1][b])
             for b in range(2)] for a in range(2)]


def determinant(matrix):
    return basic.multiply(matrix[0][0], matrix[1][1]) ^ basic.multiply(matrix[0][1], matrix[1][0])


def adjugate(matrix):
    return [[matrix[1][1], matrix[0][1]], [matrix[1][0], matrix[0][0]]]


def reverse(polynomial, degree):
    return sum(1 << (degree-k) for k in range(degree+1) if polynomial >> k & 1)


def flatten(local, cut, side):
    dimension = len(next(iter(local.values())))
    matrix = [[0]*4 for _ in range(4*dimension*dimension)]
    for (x, y, u, v), entries in local.items():
        left, right = (2*x+y, 2*u+v) if cut == 0 else (2*x+u, 2*y+v)
        for a, b in itertools.product(range(dimension), repeat=2):
            row, column = (left, right) if side == 0 else (right, left)
            matrix[dimension*dimension*row+dimension*a+b][column] = entries[a][b]
    return matrix


def adaptive_orders(matrix, degree):
    maximum = len(matrix[0])*degree+1
    precision = min(maximum, degree+1)
    while True:
        try:
            return basic.smith_orders(matrix, degree, precision=precision)
        except AssertionError:
            if precision >= maximum:
                raise
            precision = min(maximum, 2*precision)


def profiles(local, degree, gauge, audit=False):
    det = determinant(gauge)
    assert det
    shift = basic.valuation(det)
    inverse_numerator = adjugate(gauge)
    out = [None]*8
    audits = 0
    for infinity in (False, True):
        changed = {}
        for labels, entries in local.items():
            if infinity:
                entries = [[reverse(value, degree) for value in row] for row in entries]
            changed[labels] = product(product(gauge, entries), inverse_numerator)
        bound = max(value.bit_length()-1 for entries in changed.values()
                    for row in entries for value in row)
        for cut, side in itertools.product(range(2), repeat=2):
            matrix = flatten(changed, cut, side)
            orders = adaptive_orders(matrix, bound)
            if audit:
                assert orders == basic.smith_orders(matrix, bound)
                precision = max(orders)+1
                assert basic.toeplitz_rank(matrix, precision) == sum(precision-order for order in orders)
                audits += 2
            index = 4*cut+2*side+int(infinity)
            out[index] = [degree+shift-order for order in orders] if infinity else [order-shift for order in orders]
    return out, audits


def check_enumeration(depth=6):
    forms = canonical_forms(depth)
    counts = [sum(a+b == h for a, b, p in forms) for h in range(depth+1)]
    assert counts == [1]+[3*2**(h-1) for h in range(1, depth+1)]
    assert len(forms) == 1+3*(2**depth-1)
    signatures = set()
    membership_checks = neighbor_checks = 0
    for form in forms:
        a, b, p = form
        matrix = gauge_matrix(form)
        assert determinant(matrix) == 1 << (a+b)
        assert product(matrix, adjugate(matrix)) == [[1 << (a+b), 0], [0, 1 << (a+b)]]
        signature = 0
        for f, g in itertools.product(range(1 << depth), repeat=2):
            first = not (f & ((1 << a)-1))
            second = not ((g ^ basic.multiply(f >> a, p)) & ((1 << b)-1))
            if first and second:
                signature |= 1 << ((f << depth) | g)
            membership_checks += 1
        assert signature.bit_count() == 1 << (2*depth-a-b)
        assert signature not in signatures
        signatures.add(signature)
        adjacent = neighbors(form)
        assert len(set(adjacent)) == 3
        levels = [x+y for x, y, p in adjacent]
        if a+b:
            assert sorted(levels) == [a+b-1, a+b+1, a+b+1]
        else:
            assert levels == [1, 1, 1]
        for other in adjacent:
            assert form in neighbors(other)
            neighbor_checks += 1
    return {'maximum_depth': depth, 'classes_by_depth': counts,
            'total_classes': len(forms), 'finite_lattice_membership_checks': membership_checks,
            'bidirectional_neighbor_checks': neighbor_checks}


def check_unit_invariance():
    generator = random.Random(51719)
    forms = canonical_forms(4)
    families = [('elementary', *basic.elementary(2)),
                ('ordered pairing', *basic.pairing((1, 6, 9), 0)),
                ('nilpotent loop', *basic.loop(13, 18))]
    comparisons = audits = 0
    units = [[[1, 2], [0, 1]], [[1, 0], [4, 1]], [[3, 0], [0, 1]]]
    for name, local, degree in families:
        for form in generator.sample(forms, 6):
            gauge = gauge_matrix(form)
            expected, checked = profiles(local, degree, gauge, audit=True)
            audits += checked
            for unit in units:
                assert basic.valuation(determinant(unit)) == 0
                changed, _ = profiles(local, degree, product(unit, gauge))
                assert changed == expected
                comparisons += 8
            scaled = [[entry << 2 for entry in row] for row in gauge]
            changed, _ = profiles(local, degree, scaled)
            assert changed == expected
            comparisons += 8
    return {'spectral_invariance_comparisons': comparisons,
            'full_precision_and_toeplitz_checks': audits}


def check_convexity_and_certificates():
    local, degree = basic.loop(13, 18)
    cache = {}
    def get(form):
        if form not in cache:
            cache[form] = profiles(local, degree, gauge_matrix(form))[0]
        return cache[form]
    def exponents(form, index):
        values = get(form)[index]
        return values if index % 2 else [-value for value in values]
    def cost(form, index, s):
        return sum(s**value for value in exponents(form, index))
    partial_sum_checks = rational_convexity_checks = 0
    for form in canonical_forms(3):
        adjacent = neighbors(form)
        for first, second in itertools.combinations(adjacent, 2):
            for index in range(8):
                middle = sorted(exponents(form, index), reverse=True)
                left = sorted(exponents(first, index), reverse=True)
                right = sorted(exponents(second, index), reverse=True)
                for k in range(1, 5):
                    assert 2*sum(middle[:k]) <= sum(left[:k])+sum(right[:k])
                    partial_sum_checks += 1
                for s in (Fraction(3, 2), Fraction(2), Fraction(3)):
                    assert 2*cost(form, index, s) <= cost(first, index, s)+cost(second, index, s)
                    rational_convexity_checks += 1
    certificates = []
    exhaustive_checks = 0
    for index in range(8):
        s = Fraction(3)
        current = (0, 0, 0)
        path = [current]
        for _ in range(64):
            candidates = neighbors(current)
            following = min(candidates, key=lambda form: cost(form, index, s))
            if cost(following, index, s) >= cost(current, index, s):
                break
            current = following
            path.append(current)
        else:
            raise AssertionError('No finite neighbor certificate found in this control case')
        optimum = cost(current, index, s)
        assert all(cost(other, index, s) >= optimum for other in neighbors(current))
        for other in canonical_forms(6):
            assert cost(other, index, s) >= optimum
            exhaustive_checks += 1
        certificates.append({'profile_index': index, 's': str(s),
                             'gauge': current, 'value': str(optimum),
                             'neighbor_values': [str(cost(other, index, s)) for other in neighbors(current)],
                             'strict_descent_steps': len(path)-1})
    return {'partial_sum_convexity_checks': partial_sum_checks,
            'exact_rational_path_convexity_checks': rational_convexity_checks,
            'neighbor_certificates': certificates,
            'larger_depth_comparison_checks': exhaustive_checks,
            'distinct_gauges_evaluated': len(cache)}


def scalar_profiles(local, degree):
    result = [None]*8
    for infinity in (False, True):
        changed = {labels: [[reverse(value[0][0], degree) if infinity else value[0][0]]]
                   for labels, value in local.items()}
        for cut, side in itertools.product(range(2), repeat=2):
            orders = basic.smith_orders(flatten(changed, cut, side), degree)
            result[4*cut+2*side+int(infinity)] = [degree-value for value in orders] if infinity else orders
    return result


def check_common_projections():
    component_checks = rational_checks = trace_checks = 0
    def compare(full, scalar):
        nonlocal component_checks, rational_checks
        for index in range(8):
            for first, second in zip(full[index], scalar[index]):
                assert first >= second if index % 2 else first <= second
                component_checks += 1
            sign = 1 if index % 2 else -1
            assert sum(Fraction(3)**(sign*x) for x in full[index]) >= sum(
                Fraction(3)**(sign*x) for x in scalar[index])
            rational_checks += 1
    local, degree = basic.loop(13, 18)
    shadow = {labels: [[entries[0][0]]] for labels, entries in local.items()}
    scalar = scalar_profiles(shadow, degree)
    for form in ((0, 0, 0), (0, 3, 0), (2, 0, 0)):
        full, _ = profiles(local, degree, gauge_matrix(form))
        compare(full, scalar)
    # Dimension three is invertible in F_2, so normalized trace is trace.
    # The nonconstant unit conjugation changes diagonal coordinates while
    # preserving this common scalar projection.
    local, _ = basic.elementary(3)
    for (x, y, u, v), entries in local.items():
        if x != u:
            entries[1][0] ^= 6
        if y != v:
            entries[2][1] ^= 6
        if x == v and y != u:
            entries[0][2] ^= 6
    degree = 2
    shadow = {labels: [[entries[0][0] ^ entries[1][1] ^ entries[2][2]]]
              for labels, entries in local.items()}
    scalar = scalar_profiles(shadow, degree)
    gauge = [[1, 2, 0], [0, 1, 2], [0, 0, 1]]
    inverse = [[1, 2, 4], [0, 1, 2], [0, 0, 1]]
    assert basic.matrix_product(gauge, inverse) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    full = [None]*8
    for infinity in (False, True):
        changed = {}
        for labels, entries in local.items():
            entries = [[reverse(value, degree) if infinity else value for value in row]
                       for row in entries]
            changed[labels] = basic.matrix_product(basic.matrix_product(gauge, entries), inverse)
            assert changed[labels][0][0] ^ changed[labels][1][1] ^ changed[labels][2][2] == (
                entries[0][0] ^ entries[1][1] ^ entries[2][2])
            trace_checks += 1
        bound = max(value.bit_length()-1 for entries in changed.values() for row in entries for value in row)
        for cut, side in itertools.product(range(2), repeat=2):
            orders = basic.smith_orders(flatten(changed, cut, side), bound)
            full[4*cut+2*side+int(infinity)] = [degree-value for value in orders] if infinity else orders
    compare(full, scalar)
    assert (1 ^ 1) == 0  # Trace of I_2 cannot be normalized to one in F_2.
    return {'componentwise_smith_projection_checks': component_checks,
            'rational_laurent_sum_comparisons': rational_checks,
            'nonconstant_gauge_trace_invariance_checks': trace_checks,
            'characteristic_dividing_dimension_negative_control': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'status': 'passed', 'field_cardinality': 2,
              'enumeration': check_enumeration(),
              'unit_and_scalar_invariance': check_unit_invariance(),
              'convexity_and_certificates': check_convexity_and_certificates(),
              'common_scalar_projections': check_common_projections(),
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'base_smith_checker_sha256': hashlib.sha256(Path(basic.__file__).read_bytes()).hexdigest(),
              'scope': 'Exact finite controls for lattice-gauge enumeration and fixed-s neighbor certificates; no termination, ramified-gauge exclusion, or improved factoring bound.'}
    print(json.dumps(report, indent=2))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
