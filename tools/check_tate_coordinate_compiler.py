#!/usr/bin/env python3
"""Tate/Delta composition, division identities, and internal coefficient controls."""

from hashlib import sha256
import json
from math import comb, isqrt
from pathlib import Path
from random import Random

from cyclotomic_character_decoder import modular_inverse, polynomial_product
from lambert_theta_sections import AdditivePhaseAlgebra
from tate_coordinate_compiler import (
    algebraic_coefficient, compile_coordinate, delta_composition,
    division_polynomial, embed_coefficient, evaluate_psi, hensel_coefficients,
    nodal_trace, plan_algebra, q3_algebraic_coefficient,
    trace_one_element, trace_tate_coefficient,
)

ROOT = Path(__file__).resolve().parents[1]


def pascal_mask(exponent):
    """(1+z)^exponent modulo two, independently of the division recurrence."""
    result, shift = 1, 1
    while exponent:
        if exponent & 1:
            result ^= result << shift
        exponent >>= 1
        shift <<= 1
    return result


def check_nodal_identity():
    cases = []
    for index in range(1, 66):
        polynomial = division_polynomial(index).get(0, 0)
        degree = polynomial.bit_length()-1
        numerator = 0
        for i in range(degree+1):
            if polynomial >> i & 1:
                numerator ^= pascal_mask(2*(degree-i)) << i
        left = polynomial_product(numerator, pascal_mask(index*index))
        right = polynomial_product((1 << index) ^ 1, pascal_mask(2*degree))
        right <<= index*(index-1)//2
        assert left == right
        cases.append({"index": index, "nodal_x_degree": degree})
    return cases


def literal_diagonal(index, plan):
    """Dense bivariate formal division, independent of binary diagonal steps."""
    algebra = plan_algebra(plan)
    inverse = modular_inverse(plan['denominator'][0, 0], algebra.modulus)
    values = {}
    for i in range(index+1):
        for j in range(index+1):
            value = plan['numerator'].get((i, j), 0)
            for (k, ell), b in plan['denominator'].items():
                if k == ell == 0 or k > i or ell > j:
                    continue
                value ^= algebra.product(b, values[i-k, j-ell])
            values[i, j] = algebra.product(value, inverse)
    return values[index, index]


def literal_tate_coefficients(prime, limit, constant):
    """Direct Lambert coefficients with only odd divisors, no Hensel or Delta."""
    algebra = AdditivePhaseAlgebra(prime)
    result = [constant]+[0]*limit
    for divisor in range(1, limit+1, 2):
        value = algebra.phase(divisor) ^ algebra.phase(-divisor)
        for integer in range(divisor, limit+1, divisor):
            result[integer] ^= value
    return result


def check_plans():
    records = []
    for prime in (3, 5, 7, 11, 13, 17, 19, 23, 31):
        plan = compile_coordinate(prime)
        intrinsic = compile_coordinate(prime, intrinsic=True)
        algebra = AdditivePhaseAlgebra(prime)
        assert not plan['equation'].get((0, 0), 0)
        assert algebra.product(plan['derivative'], plan['derivative_inverse']) == 1
        assert plan['equation'] == {position: embed_coefficient(value, intrinsic)
                                    for position, value in intrinsic['equation'].items()}
        eta = trace_one_element(intrinsic)
        assert nodal_trace(eta, intrinsic) == 1
        assert nodal_trace(intrinsic['constant'], intrinsic) == ((prime*prime-1)//8) % 2
        records.append({"prime": prime, "constant": plan['constant'],
                        "intrinsic_coefficient_degree": intrinsic['coefficient_modulus'].bit_length()-1,
                        "nodal_modulus": hex(intrinsic['nodal_modulus']),
                        "trace_one_element": eta,
                        "division_terms": sum(x.bit_count() for x in plan['psi'].values()),
                        "equation_terms": len(plan['equation']),
                        "numerator_terms": len(plan['numerator']),
                        "denominator_terms": len(plan['denominator']),
                        "degree_v": plan['degree_v'], "degree_w": plan['degree_w'],
                        "common_linear_space_bound_in_v": plan['linear_space_bound']})
    return records


def check_prefixes():
    records = []
    for prime, limit in ((3, 64), (5, 64), (7, 48), (11, 24), (13, 16)):
        plan = compile_coordinate(prime, intrinsic=True)
        full_plan = compile_coordinate(prime)
        coefficients = hensel_coefficients(plan, limit)
        algebra = plan_algebra(plan)
        assert not any(evaluate_psi(plan['psi'], coefficients, limit, algebra))
        products = 0
        for index in range(limit+1):
            result = algebraic_coefficient(index, plan)
            assert result['value'] == coefficients[index]
            assert all(step['child_terms'] <= plan['linear_space_bound'] for step in result['steps'])
            products += result['ring_products']
            if index <= 5:
                expected = plan['constant'] if index == 0 else literal_diagonal(index, plan)
                assert result['value'] == expected
            if index <= 8:
                assert algebraic_coefficient(index, full_plan)['value'] == embed_coefficient(result['value'], plan)
        composed = delta_composition(coefficients, limit)
        assert [embed_coefficient(value, plan) for value in composed] == literal_tate_coefficients(
            prime, limit, embed_coefficient(plan['constant'], plan))
        trace = [nodal_trace(value, plan) for value in composed]
        for n in range(limit+1):
            artin_schreier = trace[n] ^ (trace[n//2] if n % 2 == 0 else 0)
            delta = int(n > 0 and n % 2 == 1 and isqrt(n)**2 == n)
            delta_q = int(n > 0 and n % prime == 0 and (n//prime) % 2 == 1
                          and isqrt(n//prime)**2 == n//prime)
            assert artin_schreier == delta ^ delta_q
            if n:
                assert trace[n] == trace_tate_coefficient(n, prime)
                eta = trace_one_element(plan)
                remainder = composed[n] ^ plan_algebra(plan).product(eta, trace[n])
                assert nodal_trace(remainder, plan) == 0
        records.append({"prime": prime, "prefix_limit_in_v": limit,
                        "coefficients_compared": limit+1,
                        "dense_diagonal_limit": 5,
                        "composed_coefficients_compared": limit+1,
                        "odd_original_indices": (limit+1)//2,
                        "ring_products_in_compiler_queries": products})
    return records


def check_large_internal_indices():
    plan = compile_coordinate(3, intrinsic=True)
    for n in range(1, 129):
        assert q3_algebraic_coefficient(n) == (comb(4*n-2, n-1)//n) % 2
    rng = Random(20261101)
    indices = [0, 1, 2, 3, 4, 5]
    for bits in (32, 64, 128, 256, 512, 1024):
        indices += [1 << (bits-1), (1 << (bits-1))+1, (1 << bits)-1,
                    (1 << (bits-1)) | rng.getrandbits(bits-1)]
    records = []
    for index in indices:
        result = algebraic_coefficient(index, plan)
        expected = q3_algebraic_coefficient(index)
        assert result['value'] == expected
        assert len(result['steps']) == index.bit_length()
        records.append({"internal_index": str(index), "internal_index_bits": index.bit_length(),
                        "value": expected, "steps": len(result['steps']),
                        "ring_products": result['ring_products']})
    return {"small_exact_binomial_checks": 128, "large_internal_controls": records}


def check_large_original_trace_indices():
    records = []
    for prime in (3, 5, 7, 11, 13, 17, 19, 23, 31):
        for bits in (64, 256, 512):
            square = ((1 << bits)+1)**2
            extra_prime = 11 if prime == 3 else 3
            for extra in (1, extra_prime):
                integer = 2**17 * prime**23 * square * extra
                value = int(extra == 1)
                assert trace_tate_coefficient(integer, prime) == value
                records.append({"prime": prime, "original_index_bits": integer.bit_length(),
                                "extra_prime": extra, "trace_value": value})
    return records


def main():
    result = {
        "scope": "Algebraic coefficients in v and finite Tate/Delta compositions; no fast original-index engine",
        "nodal_rational_identities": check_nodal_identity(),
        "plans": check_plans(),
        "prefixes": check_prefixes(),
        "internal_index_controls": check_large_internal_indices(),
        "large_original_trace_controls": check_large_original_trace_indices(),
        "source_sha256": {name: sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/tate_coordinate_compiler.py',
                                       'tools/check_tate_coordinate_compiler.py',
                                       'tools/lambert_theta_sections.py',
                                       'tools/cyclotomic_character_decoder.py')},
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
