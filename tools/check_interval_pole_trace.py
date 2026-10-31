#!/usr/bin/env python3
"""Independent exact pole-filter and global-recursion controls."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import random

from interval_orbit_gcd import interval_factor_gcd
from interval_pole_trace import BilinearRectangle, pole_count, pole_filter_bits
from interval_unit_jets import QuadraticJets


def reduce_fraction(value, modulus):
    assert value.denominator % 2
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def literal_filter(delta, scale_bits, modulus):
    # Exact inversion in Q(w), independently of valuation/parity recursion.
    scale = 1 << scale_bits
    norm = delta * delta - delta * scale + scale * scale
    return (reduce_fraction(Fraction(scale * scale, norm), modulus),
            reduce_fraction(Fraction(scale * delta, norm), modulus))


def polynomial_value(ring, coefficients):
    result = ring.constant()
    for coefficient in reversed(coefficients):
        result = ring.add(ring.multiply(result, ring.argument()), ring.constant(coefficient))
    return result


def literal_trace(a, b, c, d, left, right):
    return sum((Fraction(1, -1 - (a * u * v + b * u + c * v + d))
                for u in range(*left) for v in range(*right)), Fraction())


def check():
    rng = random.Random(20261039)
    filter_controls = weighted_controls = tree_controls = child_controls = boundary_controls = 0
    positive_rectangles = positive_weighted_controls = 0
    max_nodes = max_bits = 0
    for _ in range(384):
        number = 2 * rng.randrange(1, 1024) + 1
        left = tuple(sorted(rng.choices(range(1, 33), k=2)))
        right = tuple(sorted(rng.choices(range(1, 33), k=2)))
        precision, degree = rng.randrange(1, 14), rng.randrange(4)
        ring = QuadraticJets(precision, degree)
        scale_bits = pole_filter_bits(number, left, right, precision)
        trace, weighted = ring.constant(), ring.constant()
        expected = sum(a * b == number for a in range(*left) for b in range(*right))
        for a in range(*left):
            weight = ring.inverse(ring.add(ring.argument(), ring.constant(-a)))
            for b in range(*right):
                value = literal_filter(number - a * b, scale_bits, ring.modulus)
                assert value == (int(a * b == number), 0)
                filter_controls += 1
                term = ring.constant(*value)
                trace = ring.add(trace, term)
                weighted = ring.add(weighted, ring.multiply(term, weight))
        assert trace == ring.constant(expected)
        g = interval_factor_gcd(number, left, right)
        derivative = tuple(i * coefficient for i, coefficient in enumerate(g) if i)
        assert weighted == ring.multiply(polynomial_value(ring, derivative),
                                         ring.inverse(polynomial_value(ring, g)))
        weighted_controls += 1
        count, stats = pole_count(number, left, right)
        assert count == expected == len(g) - 1
        positive_rectangles += bool(count)
        odd_length = min(sum(a % 2 for a in range(*left)),
                         sum(b % 2 for b in range(*right)))
        assert stats['nodes'] <= max(1, 4 * odd_length - 1)
        assert stats['max_depth'] <= max(0, odd_length - 1).bit_length()
        tree_controls += 1
        max_nodes = max(max_nodes, stats['nodes'])
        max_bits = max(max_bits, stats['max_coefficient_bits'])

    # Compare every child's zero set with the parent, not just its coefficients.
    for _ in range(128):
        state = BilinearRectangle(2 * rng.randrange(1, 9), 2 * rng.randrange(8) + 1,
                                  2 * rng.randrange(8) + 1, rng.randrange(-200, 50),
                                  (0, 12), (0, 12))
        parent_zeros = {(u, v) for u in range(12) for v in range(12)
                        if state.evaluate(u, v) == 0}
        actual = set()
        for e, child in enumerate(state.children()):
            f = (state.b * e + state.d) % 2
            for u in range(*child.left):
                for v in range(*child.right):
                    assert state.evaluate(e + 2 * u, f + 2 * v) == 2 * child.evaluate(u, v)
                    if child.evaluate(u, v) == 0:
                        actual.add((e + 2 * u, f + 2 * v))
                    child_controls += 1
        assert actual == parent_zeros

        # Close the two unit interval shifts by changing only bilinear parameters.
        a, b, c, d = (rng.randrange(1, 9), rng.randrange(1, 9),
                      rng.randrange(1, 9), rng.randrange(20))
        left, right = (2, 8), (1, 7)
        base = literal_trace(a, b, c, d, left, right)
        assert literal_trace(a, b, c + a, d + b, left, right) - base == (
            literal_trace(a, b, c, d, (left[1], left[1] + 1), right)
            - literal_trace(a, b, c, d, (left[0], left[0] + 1), right))
        assert literal_trace(a, b + a, c, d + c, left, right) - base == (
            literal_trace(a, b, c, d, left, (right[1], right[1] + 1))
            - literal_trace(a, b, c, d, left, (right[0], right[0] + 1)))
        boundary_controls += 2

    # Integer-semiprime and square controls, including narrow offset rectangles.
    for a in range(3, 50, 2):
        for b in range(a, 50, 2):
            number = a * b
            left, right = (max(1, a - 5), a + 7), (max(1, b - 4), b + 8)
            actual, stats = pole_count(number, left, right)
            assert actual == sum(number % u == 0 and right[0] <= number // u < right[1]
                                 for u in range(*left))
            tree_controls += 1
            if a in (3, 5, 7):
                ring = QuadraticJets(11, 3)
                displacement = pole_filter_bits(number, left, right, ring.precision)
                weighted = ring.constant()
                for u in range(*left):
                    weight = ring.inverse(ring.add(ring.argument(), ring.constant(-u)))
                    for v in range(*right):
                        value = literal_filter(number - u * v, displacement, ring.modulus)
                        weighted = ring.add(weighted, ring.multiply(ring.constant(*value), weight))
                g = interval_factor_gcd(number, left, right)
                derivative = tuple(i * coefficient for i, coefficient in enumerate(g) if i)
                assert len(g) > 1
                assert weighted == ring.multiply(polynomial_value(ring, derivative),
                                                 ring.inverse(polynomial_value(ring, g)))
                positive_weighted_controls += 1

    # The old modular-product family is globally coprime and prunes at the root.
    huge = []
    for bits in (12, 80, 200):
        a, b, precision = (1 << bits) + 1, (1 << (bits + 1)) + 1, 32
        number, left, right = a * b, (1, a), (1, b)
        count, stats = pole_count(number, left, right)
        assert count == 0 and stats['nodes'] == stats['enclosure_prunes'] == 1
        fake_order = min((a - 1) // 2, (b - 1) // 2) - precision + 1
        assert fake_order > 0
        true_count, true_stats = pole_count(number, (1, a + 1), (1, b + 1))
        assert true_count == true_stats['nodes'] == true_stats['corner_solves'] == 1
        displacement = pole_filter_bits(number, (1, a + 1), (1, b + 1), precision)
        for u, v in ((1, 1), (a - 1, b - 1), (a, b)):
            assert literal_filter(number - u * v, displacement, 1 << precision) == (int(u * v == number), 0)
        huge.append({'factor_bit_scale': bits, 'number_bits': number.bit_length(),
                     'raw_modular_fake_factor_order_at_least': fake_order,
                     'true_pair_count': count,
                     'displacement_bits': pole_filter_bits(number, left, right, precision),
                     'tree_nodes': stats['nodes'], 'included_corner_pair_count': true_count,
                     'included_corner_tree_nodes': true_stats['nodes']})
    return {'scope': 'Exact global pole filter and enumerative bilinear closure; no improved factoring bound.',
            'seed': 20261039, 'literal_quadratic_filter_controls': filter_controls,
            'weighted_logarithmic_derivative_controls': weighted_controls,
            'positive_random_rectangles': positive_rectangles,
            'additional_positive_weighted_controls': positive_weighted_controls,
            'global_count_controls': tree_controls, 'literal_child_pair_controls': child_controls,
            'boundary_translation_controls': boundary_controls,
            'random_max_tree_nodes': max_nodes, 'random_max_coefficient_bits': max_bits,
            'huge_artificial_factor_controls': huge,
            'unresolved': 'Evaluate or merge the normalized global resolvent states without the enumerative tree cost.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
