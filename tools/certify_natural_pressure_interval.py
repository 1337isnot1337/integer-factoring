#!/usr/bin/env python3
"""Exact outward pressure arithmetic on a frozen complete natural profile.

This tool does not reconstruct the input profile's polynomial matrix.
Every certificate decision uses integers or fractions, never floating logs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import time


def _ceil(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def power_interval(base: Fraction, exponent: int, bits: int) -> tuple[int, int]:
    """Bounds for base**exponent, in units of 2**(-bits)."""
    if not 0 < base <= 1 or type(exponent) is not int or exponent < 0:
        raise ValueError('Use a base in (0,1] and a nonnegative integer exponent.')
    if type(bits) is not int or bits < 2:
        raise ValueError('Use at least two precision bits.')
    scale = 1 << bits
    lower = base.numerator * scale // base.denominator
    upper = _ceil(base * scale)
    result_lower = result_upper = scale
    while exponent:
        if exponent & 1:
            result_lower = result_lower * lower // scale
            result_upper = (result_upper * upper + scale - 1) // scale
        exponent >>= 1
        if exponent:
            lower = lower * lower // scale
            upper = (upper * upper + scale - 1) // scale
    return result_lower, result_upper


def polynomial_interval(exponents, base: Fraction, bits: int) -> tuple[int, int]:
    counts = Counter(exponents)
    if not counts:
        raise ValueError('The positive polynomial must have at least one term.')
    lower = upper = 0
    for exponent, multiplicity in counts.items():
        lo, hi = power_interval(base, exponent, bits)
        lower += multiplicity * lo
        upper += multiplicity * hi
    return lower, upper


def _small_log_interval(t: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    """Bounds for log((1+t)/(1-t)), for 0 <= t <= 1/3."""
    if not 0 <= t <= Fraction(1, 3) or type(terms) is not int or terms < 1:
        raise ValueError('Use t in [0,1/3] and a positive term count.')
    square = t * t
    power = t
    total = Fraction(0)
    for j in range(terms):
        total += power / (2 * j + 1)
        power *= square
    lower = 2 * total
    tail = 2 * power / ((2 * terms + 1) * (1 - square))
    return lower, lower + tail


@lru_cache(None)
def _log_two(terms: int) -> tuple[Fraction, Fraction]:
    return _small_log_interval(Fraction(1, 3), terms)


def log_interval(value: Fraction, bits: int, terms: int) -> tuple[int, int]:
    """Outward dyadic bounds for the natural logarithm of a positive rational."""
    if value <= 0 or type(bits) is not int or bits < 2:
        raise ValueError('Use a positive rational and at least two precision bits.')
    n, d = value.numerator, value.denominator
    k = n.bit_length() - d.bit_length()
    below = n < (d << k) if k >= 0 else (n << (-k)) < d
    if below:
        k -= 1
    reduced = value / (1 << k) if k >= 0 else value * (1 << (-k))
    assert 1 <= reduced < 2
    t = (reduced - 1) / (reduced + 1)
    lo, hi = _small_log_interval(t, terms)
    log2_lo, log2_hi = _log_two(terms)
    if k >= 0:
        lo, hi = lo + k * log2_lo, hi + k * log2_hi
    else:
        lo, hi = lo + k * log2_hi, hi + k * log2_lo
    scale = 1 << bits
    scaled_lower = lo * scale
    return scaled_lower.numerator // scaled_lower.denominator, _ceil(hi * scale)


def certify(orders, energies, rho: Fraction, base: Fraction, *, bits=192, terms=96):
    """Certify B log(P) < A log(Q), A/B=2*rho, P=sum z^(A*sigma)."""
    orders, energies = tuple(orders), tuple(energies)
    if rho <= 0 or not 0 < base < 1:
        raise ValueError('Use a positive bound and a rational base strictly between zero and one.')
    if (not orders or not energies or 0 not in orders or 0 not in energies
            or any(type(x) is not int or x < 0 for x in orders + energies)):
        raise ValueError('Use nonnegative integer orders and energies, each containing zero.')
    support = 2 * rho
    a, b = support.numerator, support.denominator
    p = polynomial_interval([a * x for x in orders], base, bits)
    q = polynomial_interval([b * x for x in energies], base, bits)
    scale = 1 << bits
    assert p[0] >= scale and q[0] >= scale
    log_p = (log_interval(Fraction(p[0], scale), bits, terms)[0],
             log_interval(Fraction(p[1], scale), bits, terms)[1])
    log_q = (log_interval(Fraction(q[0], scale), bits, terms)[0],
             log_interval(Fraction(q[1], scale), bits, terms)[1])
    gap = (a * log_q[0] - b * log_p[1], a * log_q[1] - b * log_p[0])
    return {
        'rho_strict_bound': [rho.numerator, rho.denominator],
        'support_parameter': [a, b], 'rational_base': [base.numerator, base.denominator],
        'precision_bits': bits, 'logarithm_terms': terms,
        'scale_denominator_power_of_two': bits,
        'polynomial_p_dyadic_interval': list(p), 'polynomial_q_dyadic_interval': list(q),
        'log_p_dyadic_interval': list(log_p), 'log_q_dyadic_interval': list(log_q),
        'a_log_q_minus_b_log_p_dyadic_interval': list(gap),
        'profile_channel_count': len(orders), 'physical_alphabet_size': len(energies),
        'passed': gap[0] > 0,
    }


def load_profile(path: Path):
    data = json.loads(path.read_bytes())
    profile = data.get('profile', data)
    from hot_code_charge_profiles import hot_energy, validate_code
    x, y = profile['x_letters'], profile['y_letters']
    validate_code(x, profile['modes'])
    validate_code(y, profile['modes'])
    if x != y or profile['energies_x'] != profile['energies_y']:
        raise ValueError('This interface requires identical paired coordinate codes.')
    if profile['energies_x'] != list(map(hot_energy, x)):
        raise ValueError('The saved energies do not match the code.')
    orders = profile['orders_w']
    rebuilt = sorted(2 * value for block in profile['blocks'] for value in block['orders_t'])
    if (len(orders) != len(x)**2 or orders != rebuilt
            or any(len(block['orders_t']) != block['size'] for block in profile['blocks'])):
        raise ValueError('The saved blocks do not contain the complete profile.')
    return profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--rho', type=int, nargs=2, required=True)
    parser.add_argument('--base', type=int, nargs=2, required=True)
    parser.add_argument('--bits', type=int, default=192)
    parser.add_argument('--terms', type=int, default=96)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    began = time.monotonic()
    profile = load_profile(args.input)
    result = certify(profile['orders_w'], profile['energies_x'], Fraction(*args.rho),
                     Fraction(*args.base), bits=args.bits, terms=args.terms)
    result.update({
        'scope': 'Exact arithmetic on a frozen complete natural profile. '
                 'The profile construction is separately certified; no floating decision is used.',
        'input_path': str(args.input), 'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'profile_sha256': sha256(json.dumps(profile, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        'elapsed_seconds': time.monotonic() - began,
    })
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': result['passed'], 'rho': result['rho_strict_bound'],
                      'positive_gap_numerator': result['a_log_q_minus_b_log_p_dyadic_interval'][0],
                      'elapsed_seconds': result['elapsed_seconds']}))
    if not result['passed']:
        raise SystemExit('The requested strict pressure certificate did not pass.')


if __name__ == '__main__':
    main()
