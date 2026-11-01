#!/usr/bin/env python3
"""Exact pressure certificate for one pure correction at every tensor length."""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

from check_geometric_low_flag_obstruction import verify as verify_geometric
from check_pure_two_copy_flag_floor import profile_and_flags


def ceiling_root(value, degree):
    """Integer Newton iteration; both final inequalities are checked exactly."""
    assert value > 0 and degree > 0
    current = 1 << ((value.bit_length()+degree-1)//degree)
    while True:
        following = ((degree-1)*current+value//pow(current, degree-1))//degree
        if following >= current:
            break
        current = following
    while pow(current, degree) > value:
        current -= 1
    while pow(current+1, degree) <= value:
        current += 1
    answer = current+int(pow(current, degree) < value)
    assert pow(answer-1, degree) < value <= pow(answer, degree)
    return answer


def pressure_certificate(coefficients, m=82, n=113):
    """Prove A(s)^n >= Q(s^(n/(2m)))^(2m) on the whole unit interval."""
    a = list(coefficients); a[10] -= 1
    assert all(c >= 0 for c in a) and a[0] == a[1] == 1
    assert len(a) == 21 and n > m > 0
    # On s<=1/D, use u=s^(n/m)<=s, Q(s^(n/(2m)))<=1+4u,
    # and the exact derivative bound (1+4u)^(2m)<=1+16m*u.
    D = 1024
    assert (D+4)**(2*m-1) < 2*D**(2*m-1)
    assert n**m * D**(n-m) > (16*m)**m
    root_depth = 48

    def a_numerator(left, depth):
        return sum(c*pow(left, i)*(1 << (depth*(20-i))) for i, c in enumerate(a) if c)

    def q_numerator(upper):
        return sum(pow(upper, i)*(1 << (root_depth*(10-i))) for i in (0, 2, 6, 8, 10))

    pending = [(1, depth) for depth in range(1, 11)]
    leaves = []; trials = 0; started = time.monotonic()
    while pending:
        left, depth = pending.pop()
        shift = 2*m*root_depth-n*depth
        assert shift >= 0 and depth <= 40
        upper = ceiling_root(pow(left+1, n) << shift, 2*m)
        # upper/2^root_depth is an upper bound for the rational power of
        # the interval's right endpoint. This comparison is the proof gate.
        assert (pow(upper, 2*m) << (n*depth)) >= (pow(left+1, n) << (2*m*root_depth))
        lhs, rhs = pow(a_numerator(left, depth), n), pow(q_numerator(upper), 2*m)
        denominator_difference = 20*m*root_depth-20*depth*n
        if denominator_difference >= 0:
            lhs <<= denominator_difference
        else:
            rhs <<= -denominator_difference
        trials += 1
        if lhs > rhs:
            leaves.append([left, depth, upper])
        else:
            assert depth < 35, (left, depth, 'No exact interval certificate at this rate')
            pending.extend(((2*left, depth+1), (2*left+1, depth+1)))
        if trials % 500 == 0:
            print(json.dumps({'interval_trials': trials, 'accepted_intervals': len(leaves),
                              'elapsed_seconds': time.monotonic()-started}), file=sys.stderr, flush=True)
    leaves.sort(key=lambda row: Fraction(row[0], 1 << row[1]))
    endpoint = Fraction(1, D)
    for left, depth, _ in leaves:
        assert Fraction(left, 1 << depth) == endpoint
        endpoint = Fraction(left+1, 1 << depth)
    assert endpoint == 1
    assert 2*m*m > n*n
    # The lower partition does not support an unjustified 3/4 floor.
    t = Fraction(9, 10)
    lower = sum(c*t**(6*i) for i, c in enumerate(a))
    physical = sum(t**(4*i) for i in (0, 2, 6, 8, 10))
    assert lower**4 < physical**6
    return {'rate_floor': [m, n], 'partition_lower_bound_coefficients': a,
            'small_interval_denominator': D,
            'small_interval_derivative_bound': 16*m,
            'small_interval_exact_inequalities_verified': True,
            'dyadic_root_upper_bound_depth': root_depth, 'dyadic_intervals': leaves,
            'interval_count': len(leaves), 'maximum_interval_depth': max(row[1] for row in leaves),
            'strictly_above_one_over_sqrt_two': True, 'overclaimed_three_quarters_rejected': True}


def verify():
    geometry = verify_geometric()
    data = profile_and_flags()
    pressure = pressure_certificate(data['relative_profile_polynomial'])
    # Exact finite controls for the factorization used in the all-k proof.
    controls = 0
    for denominator in (2, 3, 5, 8, 13):
        for numerator in range(denominator+1):
            s = Fraction(numerator, denominator)
            p = sum(c*s**i for i, c in enumerate(data['relative_profile_polynomial']))
            delta = s**10; a = p-delta
            assert a >= 1
            for copies in (1, 2, 3, 8, 32):
                difference = delta*sum(p**(copies-1-j)*a**j for j in range(copies))
                assert difference == p**copies-a**copies and difference >= delta
                controls += 1
    root = Path(__file__).resolve().parent
    names = (Path(__file__).name, 'check_geometric_low_flag_obstruction.py',
             'check_pure_two_copy_flag_floor.py')
    return {'status': 'Exact all-tensor-length one-pure-term pressure floor certified',
            'scope': 'Every positive tensor length, fixed natural physical heights, any characteristic-two coefficient field, arbitrary rational physically pure correction with a strict polynomial endpoint.',
            'reference_digest': data['reference_digest'], 'geometric_rank_one_checks': geometry,
            'relative_profile_polynomial': data['relative_profile_polynomial'],
            'pressure_certificate': pressure, 'finite_all_copy_factorization_controls': controls,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names},
            'geometric_certificate_sha256': sha256((root.parent/'evidence/geometric-low-flag-obstruction.json').read_bytes()).hexdigest(),
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
