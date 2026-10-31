#!/usr/bin/env python3
"""Exact compact certificates for the complete homogeneous-filter obstruction.

The accompanying ordinary proof covers every power-of-two interval length,
every charge, and every positive scale.  This checker proves its remaining
compact inequalities with outward integer arithmetic, audits each complete
partition, and preserves a counterexample to an unnecessary stronger claim.
"""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from charged_gaussian_floor import Dyadic, split_rectangle


F = Fraction


def interval(arithmetic, left, right):
    return (arithmetic.rational(left.numerator, left.denominator)[0],
            arithmetic.rational(right.numerator, right.denominator)[1])


def absolute_interval(value):
    if value[0] <= 0 <= value[1]:
        return 0, max(-value[0], value[1])
    return (-value[1], -value[0]) if value[1] < 0 else value


def ground(n, charge):
    """Exact coefficient of t after reciprocal cancellation, for J>=0."""
    return (n * n + sum(max(F(0), charge - (2 * i + 1)) for i in range(n))
            - sum(abs(4 * i - 2 * n + 2 + charge) for i in range(n)))


def ground_interval(arithmetic, n, left, right):
    # All breakpoints are integers.  A continuous piecewise affine function
    # takes its extrema at an endpoint or a breakpoint of the interval.
    points = [left, right]
    points.extend(F(value) for value in range(int(left) + 1, int(right) + 1))
    values = [ground(n, point) for point in points]
    return interval(arithmetic, min(values), max(values))


def cancelled_pressure(arithmetic, n, rectangle):
    lo_t, hi_t, lo_j, hi_j = rectangle
    scale = interval(arithmetic, lo_t, hi_t)
    charge = interval(arithmetic, lo_j, hi_j)
    result = arithmetic.multiply(ground_interval(arithmetic, n, lo_j, hi_j), scale)
    for i in range(n):
        odd = arithmetic.rational(2 * i + 1)
        energies = (absolute_interval(arithmetic.subtract(odd, charge)),
                    arithmetic.add(odd, charge))
        for energy in energies:
            result = arithmetic.add(result, arithmetic.softplus(
                arithmetic.multiply(energy, scale)))
        result = arithmetic.add(result, arithmetic.multiply(
            arithmetic.rational(2), arithmetic.softplus(arithmetic.multiply(odd, scale))))
    physical = arithmetic.zero
    physical_scale = arithmetic.divide(scale, arithmetic.root_two)
    for i in range(n):
        energy = absolute_interval(arithmetic.add(
            arithmetic.rational(4 * i - 2 * n + 2), charge))
        physical = arithmetic.add(physical, arithmetic.softplus(
            arithmetic.multiply(energy, physical_scale)))
    return arithmetic.subtract(result, arithmetic.multiply(
        arithmetic.multiply(arithmetic.rational(2), arithmetic.root_two), physical))


def support_floor_interval(arithmetic, n, rectangle):
    scale = interval(arithmetic, rectangle[0], rectangle[1])
    floor = arithmetic.multiply(arithmetic.rational(n), arithmetic.softplus(
        arithmetic.multiply(arithmetic.rational(2 * n), scale)))
    return arithmetic.subtract(cancelled_pressure(arithmetic, n, rectangle), floor)


def compact_gap_interval(arithmetic, domain):
    """2f(u)+f(2u)-2sqrt(2)f(sqrt(2)u)-1/10, on [0,2]."""
    argument = interval(arithmetic, *domain)
    positive = arithmetic.add(
        arithmetic.multiply(arithmetic.rational(2), arithmetic.softplus(argument)),
        arithmetic.softplus(arithmetic.multiply(arithmetic.rational(2), argument)))
    negative = arithmetic.multiply(arithmetic.multiply(arithmetic.rational(2),
                                                       arithmetic.root_two),
        arithmetic.softplus(arithmetic.multiply(arithmetic.root_two, argument)))
    return arithmetic.subtract(arithmetic.subtract(positive, negative),
                               arithmetic.rational(1, 10))


def split_interval(domain):
    left, right = domain
    middle = (left + right) / 2
    return ('0', (left, middle)), ('1', (middle, right))


def construct_cover(function, root, split):
    pending = [('', root)]
    leaves = []
    evaluated = 0
    while pending:
        path, domain = pending.pop()
        lower = function(domain)[0]
        evaluated += 1
        if lower > 0:
            leaves.append([path, str(lower)])
        else:
            assert len(path) < 36, (path, domain, lower)
            pending.extend((path + symbol, child) for symbol, child in split(domain))
        assert evaluated < 200000
    return {
        'domain': [str(value) for value in root],
        'evaluated_nodes': evaluated,
        'positive_leaves': len(leaves),
        'minimum_lower_numerator': str(min(int(value) for _, value in leaves)),
        'cover': leaves,
    }


def audit_cover(function, record, split):
    """Replay the partition, with no reliance on its generation traversal."""
    root = tuple(F(value) for value in record['domain'])
    entries = dict(record['cover'])
    assert len(entries) == len(record['cover']), 'Duplicate leaf path'
    proper_prefixes = {path[:length] for path in entries for length in range(len(path))}
    visited = set()
    nodes = 0

    def visit(path, domain):
        nonlocal nodes
        nodes += 1
        if path in entries:
            assert path not in proper_prefixes, 'Overlapping leaf paths'
            assert int(entries[path]) == function(domain)[0] > 0, 'Incorrect lower bound'
            visited.add(path)
            return
        assert path in proper_prefixes, 'Missing portion of the claimed cover'
        for symbol, child in split(domain):
            visit(path + symbol, child)

    visit('', root)
    assert visited == set(entries)
    assert nodes == record['evaluated_nodes']
    assert len(visited) == record['positive_leaves']
    return nodes


def exact_profile_controls(arithmetic):
    profiles = overlaps = 0
    for n in (1, 2, 3, 4, 8, 16, 32):
        for numerator in range(8 * n + 1):
            charge = F(numerator, 4)
            zero = [F(2 * i + 1) + sign * charge
                    for i in range(n) for sign in (-1, 1)]
            infinity = [F(2 * i + 1 - 2 * n) for i in range(2 * n)]
            physical = [F(4 * i - 2 * n + 2) + sign * charge
                        for i in range(n) for sign in (-1, 1)]
            expected_ground = (-sum(min(x, 0) for x in zero + infinity)
                               + sum(min(x, 0) for x in physical))
            assert expected_ground == ground(n, charge)
            assert Counter(map(abs, infinity)) == Counter(
                F(2 * i + 1) for i in range(n) for _ in range(2))
            assert Counter(map(abs, physical)) == Counter(
                abs(F(4 * i - 2 * n + 2) + charge)
                for i in range(n) for _ in range(2))
            assert sum(abs(4 * i - 2 * n + 2) for i in range(n)) <= n * n
            profiles += 1
        for charge in (F(0), F(1, 3), F(2 * n - 1), F(2 * n)):
            for scale in (F(1, 128 * n), F(1, 2 * n), F(3, 2 * n)):
                original = arithmetic.finite_pressure(
                    n, interval(arithmetic, scale, scale), interval(arithmetic, charge, charge))
                cancelled = cancelled_pressure(arithmetic, n, (scale, scale, charge, charge))
                assert max(original[0], cancelled[0]) <= min(original[1], cancelled[1])
                overlaps += 1
    return {'exact_signed_profile_controls': profiles,
            'independent_complete_expression_overlap_controls': overlaps}


def analytic_controls():
    controls = []

    def require(name, condition):
        assert condition, name
        controls.append(name)

    require('sqrt2 below 10/7', F(10, 7) ** 2 > 2)
    require('small-scale pressure floor', F(2, 21) - F(8, 128) > 0)
    require('edge kernel above 1/16', F(8, 21) - F(5, 16) > F(1, 16))
    require('e^3 above 20', sum(F(3) ** k / factorial(k) for k in range(9)) > 20)
    require('e^4 above 50', sum(F(4) ** k / factorial(k) for k in range(10)) > 50)
    require('root scale below 3/(2n) for n<=16', F(20 * 3, 2 * 16 + 3) > 1)
    require('root scale below 3/n for n<=256', F(400 * 3, 256 + 3) > 1)
    require('root scale below 3/n at n512', F(400 * 3, 512 + 3) > 1)
    require('bounded edge scale at n512', F(512, 30) - F(1, 20) > 15)
    require('bounded edge scale at n1024', F(1024, 40) - F(1, 20) > 20)
    require('large edge scalar exponential', 2500 > 128 * 4)
    require('large edge length threshold', F(3 * 1024, 256 * 11) > F(3, 4))
    require('power-of-two threshold monotonicity', all(2 * (k + 1) > k + 2
                                                      for k in range(10, 100)))
    require('tail charge window contained in five units', 2 + 2 * F(10, 7) < 5)
    require('tail alpha above three halves', 2 - F(5, 11) > F(3, 2))
    return controls


def stronger_floor_negative_control(arithmetic):
    """The sum floor is false in the interior: n4, J0, t3/5."""
    scale = arithmetic.rational(3, 5)
    difference = arithmetic.finite_pressure(4, scale, arithmetic.zero)
    for i in range(4):
        difference = arithmetic.subtract(difference, arithmetic.softplus(
            arithmetic.multiply(arithmetic.rational(4 * i + 2), scale)))
    assert difference[1] < 0
    return {
        'rejected_claim': 'G_n(J,t) > sum_{i<n} f((4*i+2)*t) for every charge',
        'n': 4, 'charge': '0', 'scale': '3/5',
        'difference_interval_numerators': [str(value) for value in difference],
        'upper_endpoint_is_strictly_negative': True,
    }


def main():
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=directory.parent / 'evidence' /
                        'global-homogeneous-boolean-filter-floor.json')
    args = parser.parse_args()
    arithmetic = Dyadic(80)
    scalar_function = lambda domain: compact_gap_interval(arithmetic, domain)
    scalar = construct_cover(scalar_function, (F(0), F(2)), split_interval)
    audited = audit_cover(scalar_function, scalar, split_interval)
    bases = []
    for n in (1, 2, 4, 8, 16, 32, 64, 128, 256):
        rectangle = (F(1, 128 * n), F(3, 2 * n) if n <= 16 else F(3, n),
                     F(0 if n <= 16 else 2 * n - 5), F(2 * n))
        function = lambda domain, n=n: support_floor_interval(arithmetic, n, domain)
        record = construct_cover(function, rectangle, split_rectangle)
        audited += audit_cover(function, record, split_rectangle)
        record['n'] = n
        bases.append(record)

    for mutated in (dict(scalar, cover=scalar['cover'][1:]),
                    dict(scalar, cover=[[path, str(int(value) + (index == 0))]
                                       for index, (path, value) in enumerate(scalar['cover'])])):
        try:
            audit_cover(scalar_function, mutated, split_interval)
        except AssertionError:
            pass
        else:
            raise AssertionError('Malformed certificate was accepted')

    names = ('charged_gaussian_floor.py', Path(__file__).name)
    evidence = {
        'status': 'Every compact support-floor and scalar certificate passed exact outward arithmetic.',
        'scope': 'All charges and scales for homogeneous binary one-level Boolean occupation filters, at every admissible power-of-two length; no factoring exponent improvement.',
        'source_baseline': 'c2e951fa1f391edd97ba587d15ec13d3d6c8892d',
        'dyadic_precision': 80, 'common_denominator': str(arithmetic.scale),
        'positive_cover_leaves': scalar['positive_leaves'] + sum(row['positive_leaves'] for row in bases),
        'independently_audited_tree_nodes': audited,
        'rejected_certificate_negative_controls': 2,
        'analytic_rational_controls': analytic_controls(),
        **exact_profile_controls(arithmetic),
        'stronger_floor_negative_control': stronger_floor_negative_control(arithmetic),
        'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
        'scalar_gap': scalar, 'finite_bases': bases,
        'proof_dependencies': [
            'homogeneous-boolean-filter-bounds.md: exact ratio and large-scale threshold',
            'homogeneous-boolean-charge-bands.md: support loss, exterior floor, and interior band',
            'global-charged-gaussian-floor.md: all-charge Gaussian floor and scalar edge kernel',
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, separators=(',', ':')) + '\n')
    summary = {key: value for key, value in evidence.items() if key not in ('scalar_gap', 'finite_bases')}
    summary['scalar_gap'] = {key: value for key, value in scalar.items() if key != 'cover'}
    summary['finite_bases'] = [{key: value for key, value in row.items() if key != 'cover'} for row in bases]
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
