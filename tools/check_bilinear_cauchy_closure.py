#!/usr/bin/env python3
"""Exact two-stage Cauchy contraction over a subfield and its Fermat torus.

The input polynomial is a Frobenius-binomial product of formal degree
strictly below the subfield size. The algorithm does not enumerate either
summation domain; finite enumerations below are independent test oracles.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random

from check_frobenius_rational_trace import field_of_degree
from fermat_torus_binomial_contraction import torus_binomial_sum


COUNTS = Counter()


def inverse_zero(field, value):
    return field.inverse(value) if value else 0


def factor_value(field, factor, value):
    a, b, exponent = factor
    return a ^ field.mul(b, field.power(value, 1 << (exponent % field.degree)))


def polynomial_value_derivative(field, factors, value):
    values = [factor_value(field, factor, value) for factor in factors]
    prefix = [1]
    for term in values:
        prefix.append(field.mul(prefix[-1], term))
    suffix, derivative = 1, 0
    for index in reversed(range(len(factors))):
        _, slope, exponent = factors[index]
        if exponent == 0:
            derivative ^= field.mul(slope, field.mul(prefix[index], suffix))
        suffix = field.mul(suffix, values[index])
    return prefix[-1], derivative


def inverse_factors(field, denominators):
    """Positive-power inverse-or-zero: bits 0 and 2,...,d."""
    result = []
    for a, b, exponent in denominators:
        first, second = a, b
        for power in range(field.degree + 1):
            if power != 1:
                result.append((first, second, exponent + power))
            first, second = field.mul(first, first), field.mul(second, second)
    assert len(result) == field.degree * len(denominators)
    return result


def outer_value(field, w, numerators, denominators, value, monomial=0, scalar=1):
    modulus = (1 << w) + 1
    result = field.mul(scalar, field.power(value, monomial % modulus))
    for factor in numerators:
        result = field.mul(result, factor_value(field, factor, value))
    for factor in denominators:
        result = field.mul(result, inverse_zero(field, factor_value(field, factor, value)))
    return result


def cauchy_torus_sum(field, w, polynomial, gamma, numerators, denominators,
                     monomial=0, scalar=1):
    """Sum W(v)*sum P(u)*iota(1+gamma*u*v) in polynomially many operations."""
    if w < 1 or field.degree % (2*w):
        raise ValueError('The explicit field must contain F_(2^(2w)).')
    if any(not isinstance(e, int) or e < 0 for _, _, e in polynomial):
        raise ValueError('Input polynomial exponents must be nonnegative integers.')
    size = 1 << w
    formal_degree = sum(1 << e for _, _, e in polynomial)
    if formal_degree >= size:
        raise ValueError('The formal input degree must be strictly below 2^w.')
    outer_factors = list(numerators) + inverse_factors(field, denominators)
    if gamma == 0:
        top = 1 if formal_degree == size-1 else 0
        for _, slope, _ in polynomial:
            top = field.mul(top, slope)
        result = torus_binomial_sum(w, outer_factors, multiply=field.mul,
                                    monomial=monomial, scalar=field.mul(scalar, top))
        return result, {'gamma_zero': True, 'exceptional_point_in_torus': False,
                        'expanded_factors': len(outer_factors)}
    inverse_gamma = field.inverse(gamma)
    gamma_b = field.power(gamma, size)
    substituted = [(a, field.mul(b, field.power(inverse_gamma, 1 << e)), e+w)
                   for a, b, e in polynomial]
    factors = outer_factors + substituted + inverse_factors(field, [(gamma_b, gamma, 1)])
    main = torus_binomial_sum(w, factors, multiply=field.mul, monomial=monomial,
                             scalar=field.mul(scalar, gamma_b))
    # The potential pole is unique because squaring is a field automorphism.
    square = field.mul(gamma_b, inverse_gamma)
    exceptional = field.power(square, 1 << (field.degree-1))
    in_torus = field.power(exceptional, size+1) == 1
    correction = 0
    if in_torus:
        point = field.inverse(field.mul(gamma, exceptional))
        assert field.power(point, size) == point and point
        p_value, derivative = polynomial_value_derivative(field, polynomial, point)
        if size == 2:
            derivative ^= p_value
        inner = field.mul(point, derivative)
        correction = field.mul(inner, outer_value(field, w, numerators, denominators,
                                                  exceptional, monomial, scalar))
    return main ^ correction, {'gamma_zero': False, 'exceptional_point_in_torus': in_torus,
        'exceptional_point': exceptional if in_torus else None,
        'correction': correction, 'expanded_factors': len(factors)}


def cauchy_torus_circuit(field, w, polynomial, gamma, numerators, denominators,
                         monomial=0, scalar=1):
    """A division-free, coefficient-independent version of the same formula."""
    if w < 1 or field.degree % (2*w):
        raise ValueError('The explicit field must contain F_(2^(2w)).')
    size = 1 << w
    formal_degree = sum(1 << e for _, _, e in polynomial)
    if formal_degree >= size:
        raise ValueError('The formal input degree must be strictly below 2^w.')
    iota = lambda x: field.power(x, (1 << (field.degree+1))-3)
    indicator_zero = lambda x: 1 ^ field.power(x, (1 << field.degree)-1)
    inverse_gamma, gamma_b = iota(gamma), field.power(gamma, size)
    outer_factors = list(numerators) + inverse_factors(field, denominators)
    substituted = [(a, field.mul(b, field.power(inverse_gamma, 1 << e)), e+w)
                   for a, b, e in polynomial]
    factors = outer_factors + substituted + inverse_factors(field, [(gamma_b, gamma, 1)])
    main = torus_binomial_sum(w, factors, multiply=field.mul, monomial=monomial,
                             scalar=field.mul(scalar, gamma_b))
    exceptional = field.power(field.mul(gamma_b, inverse_gamma), 1 << (field.degree-1))
    membership = indicator_zero(1 ^ field.power(exceptional, size+1))
    point = iota(field.mul(gamma, exceptional))
    p_value, derivative = polynomial_value_derivative(field, polynomial, point)
    if size == 2:  # A dimension-dependent choice, not a coefficient branch.
        derivative ^= p_value
    outer = field.mul(scalar, field.power(exceptional, monomial % (size+1)))
    for factor in numerators:
        outer = field.mul(outer, factor_value(field, factor, exceptional))
    for factor in denominators:
        outer = field.mul(outer, iota(factor_value(field, factor, exceptional)))
    correction = field.mul(membership, field.mul(outer, field.mul(point, derivative)))
    top = 1 if formal_degree == size-1 else 0
    for _, slope, _ in polynomial:
        top = field.mul(top, slope)
    uncoupled = torus_binomial_sum(w, outer_factors, multiply=field.mul, monomial=monomial,
                                   scalar=field.mul(scalar, top))
    uncoupled = field.mul(indicator_zero(gamma), uncoupled)
    return main ^ correction ^ uncoupled


def complete_reference(field, w, polynomial, gamma, numerators, denominators,
                       monomial=0, scalar=1):
    size = 1 << w
    subfield = [u for u in range(field.size) if field.power(u, size) == u]
    torus = [v for v in range(1, field.size) if field.power(v, size+1) == 1]
    assert len(subfield) == size and len(torus) == size+1
    result = 0
    for v in torus:
        inner = 0
        for u in subfield:
            p_value, _ = polynomial_value_derivative(field, polynomial, u)
            denominator = 1 ^ field.mul(gamma, field.mul(u, v))
            if denominator == 0:
                COUNTS['zero_cauchy_denominator_points'] += 1
            inner ^= field.mul(p_value, inverse_zero(field, denominator))
        result ^= field.mul(inner, outer_value(field, w, numerators, denominators,
                                               v, monomial, scalar))
        COUNTS['outer_reference_points'] += 1
        COUNTS['joint_reference_points'] += len(subfield)
    return result


def multiply_linear(field, polynomial, constant, slope):
    result = [0] * (len(polynomial)+1)
    for index, coefficient in enumerate(polynomial):
        result[index] ^= field.mul(coefficient, constant)
        result[index+1] ^= field.mul(coefficient, slope)
    return result


def formal_pole_controls():
    for w in range(1, 5):
        field = field_of_degree(2*w)
        size = 1 << w
        nonzero = [u for u in range(1, field.size) if field.power(u, size) == u]
        denominator = [1]
        for u in nonzero:
            denominator = multiply_linear(field, denominator, 1, u)
        assert denominator == [1]+[0]*(size-2)+[1]
        numerator = list(denominator)
        for omitted in nonzero:
            quotient = [1]
            for u in nonzero:
                if u != omitted:
                    quotient = multiply_linear(field, quotient, 1, u)
            for index, value in enumerate(quotient):
                numerator[index] ^= value
        assert numerator == [0]*(size-1)+[1]
        COUNTS['exact_formal_rational_identities'] += 1
        COUNTS['distinct_simple_poles'] += len(nonzero)
    # A formal degree b polynomial is excluded even though its field values
    # agree with a lower-degree polynomial. Ignoring this changes the trace.
    field, w = field_of_degree(4), 2
    outside = next(v for v in range(1, field.size) if field.power(v, 4) != v)
    subfield = [u for u in range(field.size) if field.power(u, 4) == u]
    actual = 0
    for u in subfield:
        actual ^= field.mul(field.power(u, 4), inverse_zero(field, 1 ^ field.mul(u, outside)))
    z = field.inverse(outside)
    naive = field.mul(field.power(outside, 4), field.mul(field.power(z, 4),
                      inverse_zero(field, field.power(outside, 4) ^ outside)))
    assert actual != naive
    try:
        cauchy_torus_sum(field, w, [(0, 1, w)], 1, [], [])
    except ValueError:
        COUNTS['degree_bound_rejections'] += 1
    else:
        raise AssertionError('Invalid input degree was accepted.')
    COUNTS['unrestricted_degree_formula_counterexamples'] += 1
    # The smallest subfield needs the extra finite-part term.
    field = field_of_degree(2)
    actual = inverse_zero(field, 1) ^ inverse_zero(field, 0)
    invalid = field.mul(1, inverse_zero(field, 0))
    assert actual == 1 and invalid == 0
    COUNTS['missing_b2_correction_counterexamples'] += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(0xCA0C4E2026)
    models = []
    for w, degree in [(1,2), (1,4), (2,4), (2,8), (3,6), (4,8)]:
        field = field_of_degree(degree)
        models.append({'subfield_bits': w, 'coefficient_field_bits': degree,
                       'irreducible_polynomial': field.polynomial})
        for trial in range(24):
            # Distinct binary places guarantee formal degree <=b-1.
            exponents = [e for e in range(w) if rng.randrange(2)]
            if trial % 6 == 0:
                exponents = list(range(w))
            polynomial = [(rng.randrange(field.size), rng.randrange(field.size), e) for e in exponents]
            if trial % 7 == 0:
                polynomial = [(0, 1, 0)]
            gamma = 0 if trial % 8 == 0 else rng.randrange(1, field.size)
            numerators = [(rng.randrange(field.size), rng.randrange(field.size), rng.randrange(6*w))
                          for _ in range(trial % 3)]
            denominators = [(rng.randrange(field.size), rng.randrange(field.size), rng.randrange(6*w))
                            for _ in range(trial % 2)]
            if trial % 10 == 1:
                denominators = [(1, 1, 0)]  # Actual outer pole at v=1.
            if trial % 10 == 3:
                denominators = [(0, 0, 0)]  # Entirely zero denominator.
            monomial = rng.randrange(-4*((1<<w)+1), 4*((1<<w)+1)+1)
            scalar = rng.randrange(field.size)
            observed, details = cauchy_torus_sum(field, w, polynomial, gamma, numerators,
                                                 denominators, monomial, scalar)
            circuit = cauchy_torus_circuit(field, w, polynomial, gamma, numerators,
                                           denominators, monomial, scalar)
            expected = complete_reference(field, w, polynomial, gamma, numerators,
                                           denominators, monomial, scalar)
            assert observed == circuit == expected
            COUNTS['complete_two_stage_sums'] += 1
            COUNTS['gamma_zero_cases'] += int(details['gamma_zero'])
            COUNTS['exceptional_points_in_torus'] += int(details['exceptional_point_in_torus'])
            COUNTS['nonzero_point_corrections'] += int(details.get('correction', 0) != 0)
            COUNTS['absent_exceptional_points'] += int(not details['gamma_zero'] and not details['exceptional_point_in_torus'])
            COUNTS['b2_complete_sums'] += int(w == 1)
    formal_pole_controls()
    class CountedField:
        def __init__(self, original):
            self.original, self.degree, self.size = original, original.degree, original.size
            self.gates = 0
        def mul(self, a, b):
            self.gates += 1
            return self.original.mul(a, b)
        def power(self, a, exponent):
            result = 1
            while exponent:
                if exponent & 1:
                    result = self.mul(result, a)
                a = self.mul(a, a)
                exponent >>= 1
            return result
    schedules = []
    for fill in (0, 1, 2, 7):
        field = CountedField(field_of_degree(4))
        cauchy_torus_circuit(field, 2, [(fill, fill, 0), (fill, fill, 1)], fill,
                             [(fill, fill, 0), (fill, fill, 3)], [(fill, fill, 2)], 1, fill)
        schedules.append(field.gates)
    assert len(set(schedules)) == 1
    COUNTS['coefficient_independent_gate_schedules'] = len(schedules)
    COUNTS['multiplication_gates_in_each_schedule'] = schedules[0]
    sources = ['check_bilinear_cauchy_closure.py', 'fermat_torus_binomial_contraction.py',
               'check_frobenius_rational_trace.py', 'check_carry_fourier_elimination.py']
    report = {'schema_version':1, 'scope':'Exact restricted two-variable closure and explicit excluded-case controls; no integer multiplication phase or factoring exponent claim.',
        'counts':dict(sorted(COUNTS.items())), 'models':models,
        'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources}}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['counts'],sort_keys=True))


if __name__ == '__main__':
    main()
