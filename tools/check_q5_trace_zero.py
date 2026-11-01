#!/usr/bin/env python3
"""Exact level-five specialization; a hard branch explicitly queries tau mod 4."""

from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from lambert_theta_sections import AdditivePhaseAlgebra
from tate_coordinate_compiler import (
    compile_coordinate, delta_composition, embed_coefficient, hensel_coefficients,
)

ROOT = Path(__file__).resolve().parents[1]


def strip_fives(integer):
    while integer % 5 == 0:
        integer //= 5
    return integer


def q5_value(integer, divisor_count_mod_four):
    """Native A_5 value at an odd original index, with at most one charged query.

    Basis: c^2+c+1=0, c maps to Z^2+Z^3. The callback must return tau(M)
    modulo four. No fast unrestricted callback is implemented here.
    """
    if not isinstance(integer, int) or integer < 1 or integer % 2 == 0:
        raise ValueError("an odd positive original index is required")
    reduced = strip_fives(integer)
    root = isqrt(reduced)
    if root*root == reduced:
        return 2 ^ int(root % 5 in (1, 4))
    if reduced % 5 in (1, 4):
        return 0
    residue = divisor_count_mod_four(reduced)
    if residue not in (0, 2):
        raise ValueError("a nonsquare has even divisor count modulo four")
    return residue//2


def square_bit(integer):
    return int(integer > 0 and isqrt(integer)**2 == integer)


def unary_square_value(integer):
    return square_bit(integer) ^ (square_bit(integer//5) if integer % 5 == 0 else 0)


def quadratic_parity(integer, a, b):
    """Parity of positive solutions a*x^2+b*y^2=integer, by charged enumeration."""
    parity = steps = 0
    for y in range(1, isqrt(max(0, (integer-a)//b))+1):
        steps += 1
        remainder = integer-b*y*y
        if remainder % a == 0:
            parity ^= square_bit(remainder//a)
    return parity, steps


def q5_quadratic_tate_coefficient(integer):
    """Exact original-nome value, using four positive binary quadratic forms.

    For odd input this is the desired signed divisor coefficient. Even
    input gives the Tate coefficient, namely the signed odd-divisor data.
    """
    if not isinstance(integer, int) or integer < 1:
        raise ValueError("a positive original coefficient index is required")
    first = unary_square_value(integer)
    second = unary_square_value(integer//2) if integer % 2 == 0 else 0
    cubic = steps = 0
    for a, b in ((1, 2), (1, 10), (5, 2), (5, 10)):
        bit, count = quadratic_parity(integer, a, b)
        cubic ^= bit
        steps += count
    return (first ^ second ^ cubic) | ((first ^ second) << 1), steps


def literal_divisors(integer):
    values = set()
    for d in range(1, isqrt(integer)+1):
        if integer % d == 0:
            values.update((d, integer//d))
    return values


def literal_signed_value(integer, plan):
    residues = [0]*5
    for d in literal_divisors(integer):
        residues[d % 5] ^= 1
        residues[-d % 5] ^= 1
    algebra = AdditivePhaseAlgebra(5)
    for value in range(4):
        if algebra.signed_residues(embed_coefficient(value, plan)) == tuple(residues):
            return value
    raise AssertionError("literal signed vector left the two-dimensional image")


def is_prime_small(integer):
    return integer >= 2 and all(integer % d for d in range(2, isqrt(integer)+1))


def main():
    plan = compile_coordinate(5, intrinsic=True)
    assert plan['nodal_modulus'] == 0b111
    assert embed_coefficient(2, plan) == (1 << 2) ^ (1 << 3)
    queries = squares = cheap_zero = squarefree_controls = quadratic_steps = 0
    for integer in range(1, 16384, 2):
        called = []

        def provider(reduced):
            called.append(reduced)
            assert reduced % 5 and isqrt(reduced)**2 != reduced
            return len(literal_divisors(reduced)) % 4

        value = q5_value(integer, provider)
        assert value == literal_signed_value(integer, plan)
        quadratic, steps = q5_quadratic_tate_coefficient(integer)
        assert quadratic == value
        quadratic_steps += steps
        reduced = strip_fives(integer)
        square = isqrt(reduced)**2 == reduced
        squares += square
        cheap_zero += not square and reduced % 5 in (1, 4)
        queries += len(called)
        assert len(called) <= 1
        if all(reduced % (d*d) for d in range(2, isqrt(reduced)+1)):
            promised = 3 if reduced == 1 else int(is_prime_small(reduced) and reduced % 5 in (2, 3))
            assert promised == value
            squarefree_controls += 1
    coefficients = hensel_coefficients(plan, 64)
    composed = delta_composition(coefficients, 64)
    for integer in range(1, 65):
        quadratic, _ = q5_quadratic_tate_coefficient(integer)
        assert composed[integer] == quadratic
        if integer % 2:
            assert composed[integer] == q5_value(integer, lambda n: len(literal_divisors(n)) % 4)

    primes = [p for p in range(3, 200, 2) if p != 5 and is_prime_small(p)]
    semiprimes = 0
    for i, p in enumerate(primes):
        for q in primes[i+1:]:
            assert q5_value(p*q, lambda n: 0) == 0
            assert q5_quadratic_tate_coefficient(p*q)[0] == 0
            semiprimes += 1

    large_cheap = []
    for bits in (63, 64, 255, 256, 511, 512):
        base = (1 << bits)+1

        def forbid_query(_):
            raise AssertionError("the advertised cheap case queried divisor counts")

        integer = 5**17 * base*base
        reduced_root = isqrt(strip_fives(integer))
        expected = 2 ^ int(reduced_root % 5 in (1, 4))
        assert q5_value(integer, forbid_query) == expected
        large_cheap.append({"bits": integer.bit_length(), "case": "square after removing fives", "value": expected})
        integer *= 21
        assert q5_value(integer, forbid_query) == 0
        large_cheap.append({"bits": integer.bit_length(), "case": "nonsquare quadratic residue", "value": 0})

    supplied_hard = []
    for p in (3, 7):
        for exponent in (1, 3, 5, 7):
            integer = p**exponent * 11**200 * 19**100
            tau = (exponent+1)*201*101
            expected = int(exponent % 4 == 1)
            called = []

            def supplied_provider(reduced):
                called.append(reduced)
                assert reduced == integer
                return tau % 4

            assert q5_value(integer, supplied_provider) == expected
            assert called == [integer]
            supplied_hard.append({"bits": integer.bit_length(), "odd_prime": p,
                                  "odd_prime_exponent": exponent,
                                  "supplied_tau_mod_four": tau % 4, "value": expected})
    contrast = []
    for integer in (33, 363):
        value = q5_value(integer, lambda n: len(literal_divisors(n)) % 4)
        contrast.append({"integer": integer, "residue_mod_5": integer % 5,
                         "is_prime": is_prime_small(integer),
                         "is_square": isqrt(integer)**2 == integer,
                         "tau_mod_four": len(literal_divisors(integer)) % 4,
                         "value": value})
    assert [x['value'] for x in contrast] == [0, 1]
    result = {
        "scope": "Exact level-five remaining coordinate; unrestricted tau-mod-four engine is not implemented",
        "literal_odd_indices": 8192,
        "square_cases_without_queries": squares,
        "nonsquare_residue_cases_without_queries": cheap_zero,
        "supplied_small_tau_queries": queries,
        "quadratic_form_coefficients": 8192,
        "charged_quadratic_enumeration_steps": quadratic_steps,
        "squarefree_promise_controls": squarefree_controls,
        "tate_delta_prefix_controls": 64,
        "distinct_squarefree_semiprime_controls": semiprimes,
        "large_cheap_controls": large_cheap,
        "large_supplied_factor_controls": supplied_hard,
        "same_small_metadata_contrast": contrast,
        "source_sha256": {name: sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/check_q5_trace_zero.py', 'tools/tate_coordinate_compiler.py')},
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
