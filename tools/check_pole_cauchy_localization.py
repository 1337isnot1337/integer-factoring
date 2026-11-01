#!/usr/bin/env python3
"""Exact controls for the pole moment quotient's localization obstruction."""

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb

from interval_unit_jets import QuadraticJets
from pole_cauchy_circuit_moments import compressed_inverse_context


def pair_add(x, y, modulus):
    return ((x[0] + y[0]) % modulus, (x[1] + y[1]) % modulus)


def evaluate(coefficients, value, ring):
    result = (0, 0)
    for coefficient in reversed(coefficients):
        result = pair_add(ring.pair_product(result, value), coefficient, ring.modulus)
    return result


def valuation(value, cap):
    return cap if value == 0 else (abs(value) & -abs(value)).bit_length() - 1


def pair_valuation(value, cap):
    return min(valuation(value[0], cap), valuation(value[1], cap))


def literal_packet(corrections, degree, ring):
    """Independent root-by-root multiplication in the (Y-1) basis."""
    elementary = [(0, 0)] * (degree + 1)
    elementary[0] = (1, 0)
    for index, correction in enumerate(corrections):
        for j in range(min(index + 1, degree), 0, -1):
            term = ring.pair_product(elementary[j - 1], correction)
            elementary[j] = pair_add(elementary[j], term, ring.modulus)
    polynomial = [(0, 0)] * (degree + 1)
    for j, value in enumerate(elementary):
        for r in range(degree - j + 1):
            scale = (-1) ** (degree - r) * comb(degree - j, r)
            term = (scale * value[0], scale * value[1])
            polynomial[r] = pair_add(polynomial[r], term, ring.modulus)
    return tuple(polynomial), elementary


def physical_trace(number, first_a, first_b, count, displacement, precision):
    """Literal pole sum; no matrix inverse or packet identities are used."""
    ring = QuadraticJets(precision, 0)
    result = (0, 0)
    factors = []
    for i in range(count):
        a = first_a + 2 * i
        for j in range(count):
            b = first_b + 2 * j
            difference = number - a * b
            if difference == 0:
                factors.append([a, b])
            order = min(displacement, valuation(difference, displacement))
            remaining = displacement - order
            if remaining >= precision:
                continue
            denominator = ring.constant(difference // (1 << order), 1 << remaining)
            term = ring.multiply(ring.constant(0, 1 << remaining), ring.inverse(denominator))[0]
            result = pair_add(result, term, ring.modulus)
    z = ring.constant(number, 1 << displacement)
    correction = ring.multiply(
        ring.constant(0, comb(count, 2) * (1 << displacement)), ring.inverse(z))[0]
    result = pair_add(result, (-correction[0], -correction[1]), ring.modulus)
    return result, factors


def verify_finite_packets():
    controls = []
    coefficient_comparisons = 0
    newton_bounds = 0
    for s, precision, first_a, first_b, number in (
        (2, 2, 1, 1, 9), (3, 3, 3, 5, 45), (3, 4, 7, 1, 105),
        (4, 5, 1, 9, 225), (4, 8, 5, 3, 399),
        (5, 12, 1, 1, 105), (5, 16, 11, 7, 729),
        (6, 2, 1, 1, 10201), (6, 3, 3, 131, 13837),
        (6, 16, 1, 1, 10201), (6, 17, 3, 131, 13837), (6, 30, 5, 9, 6553),
        (7, 30, 1, 1, 15), (7, 50, 13, 3, 10201),
    ):
        count = 1 << s
        maximum_a = first_a + 2 * (count - 1)
        maximum_b = first_b + 2 * (count - 1)
        horizon = max(abs(number - first_a * first_b),
                      abs(number - maximum_a * maximum_b)).bit_length()
        displacement = max(s + 1, horizon, 1)
        context = compressed_inverse_context(
            number, [(first_a, count)], [(first_b, count)], displacement, precision)
        degree = context['quotient_degree']
        ring = QuadraticJets(precision, 0)
        z = (number % ring.modulus, (1 << displacement) % ring.modulus)
        left_corrections = []
        right_corrections = []
        for i in range(count):
            a = first_a + 2 * i
            reciprocal = ring.inverse(ring.constant(a))[0]
            root = ring.pair_product(z, reciprocal)
            left_corrections.append(((root[0] - 1) % ring.modulus, root[1]))
            right_corrections.append(((first_b + 2 * i - 1) % ring.modulus, 0))
        for key, corrections in (
            ('column_packet', right_corrections),
            ('reciprocal_row_packet', left_corrections),
        ):
            literal, elementary = literal_packet(corrections, degree, ring)
            assert literal == context[key]
            coefficient_comparisons += len(literal)
            for j in range(1, degree + 1):
                bound = min(precision, j + s - 1 - valuation(j, precision))
                assert pair_valuation(elementary[j], precision) >= bound
                newton_bounds += 1
            assert evaluate(context[key], (1, 0), ring) == (0, 0)
        assert evaluate(context['row_quotient'], z, ring) == (0, 0)
        lower_multiplicity = min(degree, max(0, s - 1 - (precision - 1).bit_length() + 1))
        assert 0 <= lower_multiplicity <= degree
        for corrections in (left_corrections, right_corrections):
            _, elementary = literal_packet(corrections, degree, ring)
            assert all(elementary[j] == (0, 0)
                       for j in range(degree - lower_multiplicity + 1, degree + 1))
        control = dict(s=s, matrix_dimension=count, precision=precision,
                       displacement=displacement, quotient_degree=degree,
                       first_row=first_a, first_column=first_b, number=number,
                       common_packet_power_lower_bound=lower_multiplicity,
                       evaluation_certificate=[0, 0])
        if count == 64 and number in (10201, 13837) and precision > displacement:
            residue, factors = physical_trace(number, first_a, first_b, count,
                                              displacement, precision - displacement)
            assert factors == [[101, 101 if number == 10201 else 137]]
            assert residue == (1, 0)
            assert 2 * displacement >= precision
            control.update(physical_trace_absolute_precision=precision - displacement,
                           physical_trace_residue=list(residue), factor_pairs=factors,
                           zero_polynomial_passes_squared_residual=True)
        controls.append(control)
    return controls, coefficient_comparisons, newton_bounds


def verify_cancellation():
    ring = QuadraticJets(12, 0)
    z = (7, 8)
    inverse_five = ring.inverse(ring.constant(5))[0]
    x = ring.pair_product(z, inverse_five)
    displacement = pair_add(z, tuple(-v for v in ring.pair_product(x, (3, 0))), ring.modulus)
    scaled_inverse = ring.pair_product((5, 0), ring.inverse(ring.constant(*z))[0])
    assert ring.pair_product(displacement, scaled_inverse) == (2, 0)
    return dict(precision=12, z=list(z), reduced_pole=list(displacement),
                scaled_inverse=list(scaled_inverse), product=[2, 0],
                scope='Exact algebraic cancellation control, not an interval instance')


def poly_trim(polynomial):
    while polynomial and not polynomial[-1]:
        polynomial.pop()
    return polynomial


def poly_add(left, right):
    result = [Fraction(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i] += value
    return poly_trim(result)


def poly_scale(polynomial, value):
    return poly_trim([coefficient * value for coefficient in polynomial])


def poly_multiply(left, right):
    if not left or not right:
        return []
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return poly_trim(result)


def poly_divide(dividend, divisor):
    remainder = dividend[:]
    quotient = [Fraction(0)] * max(0, len(dividend) - len(divisor) + 1)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        factor = remainder[-1] / divisor[-1]
        quotient[shift] = factor
        remainder = poly_add(remainder, [Fraction(0)] * shift + poly_scale(divisor, -factor))
    return poly_trim(quotient), remainder


def poly_extended_gcd(left, right):
    a, b = left[:], right[:]
    s0, s1 = [Fraction(1)], []
    t0, t1 = [], [Fraction(1)]
    while b:
        quotient, remainder = poly_divide(a, b)
        a, b = b, remainder
        s0, s1 = s1, poly_add(s0, poly_scale(poly_multiply(quotient, s1), -1))
        t0, t1 = t1, poly_add(t0, poly_scale(poly_multiply(quotient, t1), -1))
    scale = Fraction(1) / a[-1]
    return poly_scale(a, scale), poly_scale(s0, scale), poly_scale(t0, scale)


def divide_at_one(polynomial, modulus):
    result = [(0, 0)] * (len(polynomial) - 1)
    total = (0, 0)
    for i in range(len(polynomial) - 1, 0, -1):
        total = pair_add(total, polynomial[i], modulus)
        result[i - 1] = total
    assert pair_add(polynomial[0], total, modulus) == (0, 0)
    return result


def verify_canonical_lift_failure():
    number, count, displacement, precision = 13837, 64, 15, 17
    context = compressed_inverse_context(number, [(3, count)], [(131, count)],
                                         displacement, precision)
    modulus = 1 << precision
    ring = QuadraticJets(precision, 0)
    f = list(context['column_packet'])
    g = list(context['reciprocal_row_packet'])
    cancelled = 0
    while len(f) > 1 and evaluate(f, (1, 0), ring) == evaluate(g, (1, 0), ring) == (0, 0):
        f = divide_at_one(f, modulus)
        g = divide_at_one(g, modulus)
        cancelled += 1
    assert cancelled == 4
    assert all(b == 0 for a, b in f + g)
    f = [Fraction(a if a < modulus // 2 else a - modulus) for a, b in f]
    g = [Fraction(a if a < modulus // 2 else a - modulus) for a, b in g]
    assert list(map(int, f)) == [-4863, 21684, -22462, 44004, -22545, -56728,
                                -5220, -34712, 16111, 27876, -49086, -12364, 1]
    assert list(map(int, g)) == [17409, -57932, -55998, -284, 56815, -21400,
                                -19044, 52840, 25071, -5660, 18754, 22196, 1]
    common = [Fraction((-1) ** (cancelled - j) * comb(cancelled, j))
              for j in range(cancelled + 1)]
    for lifted, key in ((f, 'column_packet'), (g, 'reciprocal_row_packet')):
        reconstructed = poly_multiply(common, lifted)
        assert tuple((int(value) % modulus, 0) for value in reconstructed) == context[key]
    divisor, inverse, other = poly_extended_gcd(f, g)
    assert divisor == [Fraction(1)]
    assert poly_add(poly_multiply(f, inverse), poly_multiply(g, other)) == [Fraction(1)]
    for polynomial in (f, g):
        derivative = [i * polynomial[i] for i in range(1, len(polynomial))]
        assert poly_extended_gcd(polynomial, derivative)[0] == [Fraction(1)]
    derivative = [i * f[i] for i in range(1, len(f))]
    _, h = poly_divide(poly_multiply([Fraction(0)] + derivative, inverse), g)
    degree = len(g) - 1
    powers = [Fraction(degree)]
    for j in range(1, degree):
        powers.append(-sum(g[degree - i] * powers[j - i] for i in range(1, j))
                      - j * g[degree - j])
    trace_from_newton = sum(h[i] * powers[i] for i in range(len(h)))
    trace_from_columns = Fraction(0)
    for j in range(degree):
        _, column = poly_divide([Fraction(0)] * j + h, g)
        trace_from_columns += column[j] if j < len(column) else 0
    assert trace_from_newton == trace_from_columns
    coefficient = trace_from_columns - comb(degree, 2)
    assert coefficient == Fraction(2274178272512455570345056001700682549607049910138,
                                   210931526610747294728550895089720540989500093909)
    order = valuation(coefficient.numerator, precision) - valuation(coefficient.denominator, precision)
    assert order == 1
    actual, factors = physical_trace(number, 3, 131, count, displacement, 2)
    assert actual == (1, 0) and factors == [[101, 137]]
    return dict(number=number, matrix_dimension=count, displacement=displacement,
                precision=precision, cancelled_common_power=cancelled, virtual_dimension=degree,
                balanced_F_coefficients=list(map(int, f)), balanced_G_coefficients=list(map(int, g)),
                exact_lifts_coprime=True, exact_lifts_squarefree=True,
                polynomial_moment_agreement_precision=precision,
                virtual_trace_coefficient_numerator=coefficient.numerator,
                virtual_trace_coefficient_denominator=coefficient.denominator,
                virtual_trace_valuation=displacement + order,
                actual_trace_residue_mod_4=list(actual), virtual_trace_residue_mod_4=[0, 0],
                normalized_actual_trace_residue_mod_2power_P=[1 << displacement, 0],
                normalized_virtual_trace_residue_mod_2power_P=[0, 0])


def verify_invisible_functionals():
    controls = []
    for s in (2, 3, 4, 5):
        count = 1 << s
        order = count - 1
        precision = 2 * count - 2 - s
        modulus = 1 << precision
        weights = [(-1) ** i * comb(order, i) for i in range(count)]
        for degree in range(2 * count + 1):
            moment = sum(weights[i] * pow(1 + 2 * i, degree, modulus) for i in range(count))
            assert moment % modulus == 0
        for location in (0, count // 2, count - 1):
            assert weights[location] % 2 == 1
        controls.append(dict(matrix_dimension=count, invisible_precision=precision,
                             literal_monomial_checks=2 * count + 1,
                             pole_selector_parities=[1, 1, 1]))
    return controls


def build_receipt():
    controls, coefficient_comparisons, newton_bounds = verify_finite_packets()
    s = 40
    count = 1 << s
    number = (count - 3) ** 2
    displacement = 2 * s + 2
    precision = displacement + 2
    context = compressed_inverse_context(number, [(1, count)], [(1, count)],
                                         displacement, precision)
    ring = QuadraticJets(precision, 0)
    z = (number % ring.modulus, (1 << displacement) % ring.modulus)
    assert evaluate(context['column_packet'], (1, 0), ring) == (0, 0)
    assert evaluate(context['reciprocal_row_packet'], (1, 0), ring) == (0, 0)
    assert evaluate(context['row_quotient'], z, ring) == (0, 0)
    large = dict(s=s, matrix_dimension=count, number=number, displacement=displacement,
                 precision=precision, quotient_degree=context['quotient_degree'],
                 evaluation_certificate=[0, 0], enumerated_nodes=0,
                 scope='Compressed execution; the complete-residue theorem supplies its certificate')
    receipt = dict(schema='pole-cauchy-localization-v1',
                   claim='Direct scaled inversion is obstructed in the stated moment quotient',
                   factoring_exponent_improvement=False, finite_controls=controls,
                   literal_packet_coefficient_comparisons=coefficient_comparisons,
                   literal_newton_valuation_checks=newton_bounds,
                   large_compressed_control=large, cancellation_control=verify_cancellation(),
                   canonical_lift_failure=verify_canonical_lift_failure(),
                   invisible_functional_controls=verify_invisible_functionals(),
                   scalar_only_saturation_status='Not ruled out and not constructed')
    encoded = json.dumps(receipt, sort_keys=True, separators=(',', ':')).encode()
    receipt['controls_sha256'] = hashlib.sha256(encoded).hexdigest()
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    receipt = build_receipt()
    with open(args.output, 'w', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(dict(status='passed',
                          finite_controls=len(receipt['finite_controls']),
                          coefficient_comparisons=receipt['literal_packet_coefficient_comparisons'],
                          controls_sha256=receipt['controls_sha256']), sort_keys=True))


if __name__ == '__main__':
    main()
