#!/usr/bin/env python3
"""Exact phase-seven parameter, odd-index Sturm controls and charged evaluator.

Finite prefix providers below are explicitly charged. They are not fast
original-index coefficient oracles. The polynomial identity is checked
without truncation; the theta exclusion uses a finite separating covector.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import isqrt
from pathlib import Path

from cyclotomic_character_decoder import modular_product, polynomial_product
from lambert_theta_sections import AdditivePhaseAlgebra
from tate_coordinate_compiler import (
    binary_add, binary_product, compile_coordinate, division_polynomial,
    embed_coefficient,
)

ROOT = Path(__file__).resolve().parents[1]


def binary_power(polynomial, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, polynomial)
        exponent //= 2
        if exponent:
            polynomial = polynomial_product(polynomial, polynomial)
    return result


def divide_monic(polynomial, divisor):
    """Polynomials in Y, with binary polynomial coefficients in u."""
    remainder = dict(polynomial)
    quotient = {}
    degree = max(divisor)
    assert divisor[degree] == 1
    while remainder and max(remainder) >= degree:
        shift = max(remainder)-degree
        factor = remainder[max(remainder)]
        quotient[shift] = quotient.get(shift, 0) ^ factor
        remainder = binary_add(remainder, {
            j+shift: polynomial_product(factor, value)
            for j, value in divisor.items()})
    return quotient, remainder


def parameter_identity():
    # D=1+u+u^2, K=Y^3+u*Y^2+(1+u)*Y+1.
    psi = division_polynomial(7)
    terms = [(i, j) for j, value in psi.items()
             for i in range(value.bit_length()) if value >> i & 1]
    denominator_power = max(i+4*j for i, j in terms)
    cleared = {}
    for i, j in terms:
        term = binary_power(7, denominator_power-i-4*j) << j
        cleared = binary_add(cleared, {i: term})
    kernel = {3: 1, 2: 2, 1: 3, 0: 1}
    quotient, remainder = divide_monic(cleared, kernel)
    assert not remainder
    assert binary_product(quotient, kernel) == cleared
    # N(Y)=Y^3+Y+1 and E(Y)=Y^2+Y are coprime: N(0)=N(1)=1.
    assert (0b1011 & 1) == 1 and ((0b1011).bit_count() % 2) == 1
    return {
        "cleared_D_power": denominator_power,
        "cleared_terms": sum(v.bit_count() for v in cleared.values()),
        "kernel_Y_degree": 3,
        "quotient_Y_degree": max(quotient),
        "quotient_terms": sum(v.bit_count() for v in quotient.values()),
        "quotient_u_masks_by_Y_degree": {str(j): hex(v) for j, v in sorted(quotient.items())},
        "remainder": 0,
    }


def kummer_identity():
    # Work over F4=F2[w]/(w^2+w+1), with coefficients in u.
    # (u+w^2)(Y+w)^3+(u+w)(Y+w^2)^3 = K(Y,u).
    def product(left, right):
        result = [0]*(len(left)+len(right)-1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                result[i+j] ^= modular_product(a, b, 7)
        return result

    def add(left, right):
        return [(left[i] if i < len(left) else 0) ^
                (right[i] if i < len(right) else 0)
                for i in range(max(len(left), len(right)))]

    first = [[1], [3], [2], [1]]  # coefficients of (Y+w)^3
    second = [[1], [2], [3], [1]]
    result = [add(product(a, [3, 1]), product(b, [2, 1]))
              for a, b in zip(first, second)]
    assert result == [[1, 0], [1, 1], [0, 1], [1, 0]]
    return {"coefficient_field_modulus": "w^2+w+1", "exact_cleared_identity": True}


def product_truncated(left, right, cut):
    if left.bit_count() > right.bit_count():
        left, right = right, left
    return polynomial_product(left, right) & cut


def frobenius_series(series, power, cut):
    return sum(1 << (power*i)
               for i in range(min(series.bit_length(), (cut.bit_length()-1)//power+1))
               if series >> i & 1)


def literal_coordinate_planes(prime, limit):
    """Finite Lambert sieve: enumerate odd divisors up to limit."""
    plan = compile_coordinate(prime, intrinsic=True)
    phase = AdditivePhaseAlgebra(prime)
    width = (prime-1)//2
    inverse = {embed_coefficient(v, plan): v for v in range(1 << width)}
    weights = {r: inverse[phase.phase(r) ^ phase.phase(-r)] for r in range(prime)}
    values = [plan['constant']]+[0]*limit
    for divisor in range(1, limit+1, 2):
        weight = weights[divisor % prime]
        for integer in range(divisor, limit+1, divisor):
            values[integer] ^= weight
    return [sum(((v >> i) & 1) << n for n, v in enumerate(values))
            for i in range(width)]


def field_series_product(left, right, cut):
    """Three binary planes, in F8=F2[c]/(c^3+c+1)."""
    result = [0, 0, 0]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            coefficient = modular_product(1 << i, 1 << j, 11)
            value = product_truncated(a, b, cut)
            for k in range(3):
                if coefficient >> k & 1:
                    result[k] ^= value
    return result


def partition_parity_prefix(limit):
    """b_7(n): O(limit^2) scalar XOR operations, O(limit) bits stored."""
    values = bytearray(limit+1)
    values[0] = 1
    updates = 0
    for part in range(1, limit+1):
        if part % 7:
            for integer in range(part, limit+1):
                values[integer] ^= values[integer-part]
                updates += 1
    return values, updates


def check_original_parameter(limit):
    cut = (1 << (limit+1))-1
    partitions, updates = partition_parity_prefix((limit-1)//4)
    u = sum(int(value) << (4*n+1) for n, value in enumerate(partitions))
    u2 = frobenius_series(u, 2, cut)
    D = 1 ^ u ^ u2
    delta = sum(1 << (n*n) for n in range(1, isqrt(limit)+1, 2))
    D4 = 1 ^ frobenius_series(u, 4, cut) ^ frobenius_series(u, 8, cut)
    assert u == product_truncated(delta, D4, cut)
    X = literal_coordinate_planes(7, limit)
    Y = [product_truncated(a, D, cut) for a in X]
    Y2 = field_series_product(Y, Y, cut)
    Y3 = field_series_product(Y, Y2, cut)
    residual = [Y3[i] ^ product_truncated(u, Y2[i], cut) ^
                product_truncated(1 ^ u, Y[i], cut) ^ int(i == 0)
                for i in range(3)]
    assert residual == [0, 0, 0]
    # Tr(1)=1 and Tr(c)=Tr(c^2)=0 in this basis.
    assert Y[0] == u
    unary = sum(1 << (n*n) for n in range(1, isqrt(limit)+1))
    unary ^= sum(1 << (7*n*n) for n in range(1, isqrt(limit//7)+1))
    assert X[0] == unary ^ frobenius_series(unary, 2, cut)
    return {"original_coefficients_including_constant": limit+1,
            "partition_coefficients": len(partitions),
            "partition_scalar_XOR_updates": updates,
            "cubic_and_parameter_residuals": [0, 0, 0],
            "trace_identity": True}


def square_residue_series(prime, limit):
    variables = []
    names = []
    for scale in (1, prime):
        for residue in range((prime+1)//2):
            value = sum(1 << (scale*n*n)
                        for n in range(1, isqrt(limit//scale)+1)
                        if min(n % prime, (-n) % prime) == residue)
            names.append({"scale": scale, "unsigned_residue": residue})
            variables.append(value)
    return names, variables


def theta_span_certificate(limit=4096, degree=6):
    cut = (1 << (limit+1))-1
    names, variables = square_residue_series(7, limit)
    polynomials = {(): 1}
    columns = []
    basis = {}
    for total_degree in range(degree+1):
        for monomial in combinations_with_replacement(range(len(variables)), total_degree):
            value = (1 if not monomial else product_truncated(
                polynomials[monomial[:-1]], variables[monomial[-1]], cut))
            polynomials[monomial] = value
            columns.append(value)
            while value:
                pivot = (value & -value).bit_length()-1
                if pivot in basis:
                    value ^= basis[pivot]
                else:
                    basis[pivot] = value
                    break
    targets = literal_coordinate_planes(7, limit)
    residual = targets[1]
    for pivot in sorted(basis):
        if residual >> pivot & 1:
            residual ^= basis[pivot]
    assert residual
    first_obstruction = (residual & -residual).bit_length()-1
    covector = 1 << first_obstruction
    for pivot in sorted(basis, reverse=True):
        if (covector & basis[pivot]).bit_count() % 2:
            covector ^= 1 << pivot
    assert all((covector & column).bit_count() % 2 == 0 for column in columns)
    target_value = sum(((covector & t).bit_count() % 2) << i
                       for i, t in enumerate(targets))
    assert target_value == 6
    return {"scope": "Only polynomials of total degree at most six in these eight series",
            "prefix_limit": limit, "total_degree": degree, "variables": names,
            "columns": len(columns), "binary_rank": len(basis),
            "covector_indices": [n for n in range(limit+1) if covector >> n & 1],
            "largest_covector_index": first_obstruction,
            "monomial_pairings_nonzero": 0,
            "target_pairing_native_F8": target_value}


# Variables here are a=Theta_1, b=Theta_2, d=Theta_3, all at scale one.
# The ordinary modular/Sturm proof is in research/q7-cubic-parameter.md.
Q7_ODD_POLYNOMIALS = (
    ((0,), (0, 1), (0, 2), (0, 0, 0), (0, 0, 2), (0, 1, 1),
     (0, 1, 2), (1, 1, 1), (0, 0, 0, 1), (0, 0, 0, 2),
     (0, 1, 1, 2), (0, 1, 2, 2), (1, 2, 2, 2)),
    ((2,), (0, 2), (1, 2), (0, 0, 0), (0, 0, 2), (0, 1, 2),
     (1, 2, 2), (2, 2, 2), (0, 0, 1, 2), (0, 1, 1, 1),
     (0, 1, 1, 2), (0, 2, 2, 2), (1, 2, 2, 2)),
)


def check_odd_identity(limit):
    cut = (1 << (limit+1))-1
    odd_mask = sum(1 << n for n in range(1, limit+1, 2))
    _, variables = square_residue_series(7, limit)
    theta = variables[1:4]
    V = 0
    for value in variables:
        V ^= value
    targets = literal_coordinate_planes(7, limit)
    assert not (targets[0] ^ V) & odd_mask
    records = []
    for coordinate, monomials in enumerate(Q7_ODD_POLYNOMIALS, 1):
        value = V
        for monomial in monomials:
            term = 1
            for variable in monomial:
                term = product_truncated(term, theta[variable], cut)
            value ^= term
        assert not (targets[coordinate] ^ value) & odd_mask
        records.append({"native_coordinate": coordinate,
                        "polynomial_monomials": [list(m) for m in monomials]})
    return {"status": "Exact finite congruence covering the ordinary modular/Sturm proof bound",
            "largest_index": limit, "odd_indices": (limit+1)//2,
            "training_largest_index": 4096, "polynomial_total_degree": 4,
            "trace_checked": True, "candidates": records}



def sturm_lift_data():
    # Native coefficient duals are 1,c^2,c under the F8/F2 trace.
    def trace(value):
        squared = modular_product(value, value, 11)
        return value ^ squared ^ modular_product(squared, squared, 11)

    for coordinate, dual in enumerate((1, 4, 2)):
        assert [trace(modular_product(dual, 1 << i, 11))
                for i in range(3)] == [int(i == coordinate) for i in range(3)]
    level = 4*(4*7*7)
    index = level*level*3//4*48//49
    bound = 2*index//12
    assert (level, index, bound) == (784, 451584, 75264)
    assert 131072 >= bound
    # An independent integer Eisenstein-prefix control, before reducing mod 2.
    limit = 4096
    sigma = [0]*(limit+1)
    for divisor in range(1, limit+1):
        for n in range(divisor, limit+1, divisor):
            sigma[n] += divisor
    for n in range(1, limit+1):
        value = sigma[n] - (3*sigma[n//2] if n % 2 == 0 else 0)
        value += 2*sigma[n//4] if n % 4 == 0 else 0
        assert value % 2 == int(n % 2 == 1 and isqrt(n)**2 == n)
    return {"proof": "research/q7-cubic-parameter.md, section 5; ordinary modularity and Sturm hypotheses",
            "weight": 2, "theta_base_level": 196,
            "subgroup_used_for_bound": "Gamma_1(784)",
            "SL2_index": index, "Sturm_bound_inclusive": bound,
            "required_odd_coefficients": bound//2,
            "checked_odd_coefficients": 65536,
            "native_trace_duals": [1, 4, 2],
            "integer_Eisenstein_control_coefficients": limit}


def normalized_factors(monomial):
    # Frobenius: Theta_a^m = product over set bits j of Theta_a(t^(2^j)).
    return tuple((variable+1, 1 << bit)
                 for variable, multiplicity in sorted(Counter(monomial).items())
                 for bit in range(multiplicity.bit_length())
                 if multiplicity >> bit & 1)


def constrained_quadratic_parity(integer, factors):
    """Charged positive representations with unsigned residues modulo seven."""
    residue, scale = factors[0]
    if len(factors) == 1:
        if integer <= 0 or integer % scale:
            return 0, 1
        root = isqrt(integer//scale)
        return int(root*root*scale == integer and
                   min(root % 7, (-root) % 7) == residue), 1
    minimum_rest = sum(r*r*w for r, w in factors[1:])
    bound = isqrt(max(0, (integer-minimum_rest)//scale))
    parity = steps = 0
    for root in range(1, bound+1):
        steps += 1
        if min(root % 7, (-root) % 7) == residue:
            child, charged = constrained_quadratic_parity(
                integer-scale*root*root, factors[1:])
            parity ^= child
            steps += charged
    return parity, steps


def q7_odd_original_coefficient(integer):
    """All-index odd-input formula; O(N poly(log N)) explicit enumerator."""
    if not isinstance(integer, int) or integer < 1 or integer % 2 == 0:
        raise ValueError("a positive odd original index is required")
    unary = int(isqrt(integer)**2 == integer)
    if integer % 7 == 0:
        unary ^= int(isqrt(integer//7)**2 == integer//7)
    value = 7*unary  # 1+c+c^2 in the native basis
    cache = {}
    steps = 0
    for coordinate, monomials in enumerate(Q7_ODD_POLYNOMIALS, 1):
        bit_value = 0
        for monomial in monomials:
            factors = normalized_factors(monomial)
            if factors not in cache:
                cache[factors], charged = constrained_quadratic_parity(integer, factors)
                steps += charged
            bit_value ^= cache[factors]
        value ^= bit_value << coordinate
    return {"value": value, "charged_steps": steps,
            "distinct_quadratic_predicates": len(cache)}


def check_original_evaluator():
    controls = sorted(set(range(1, 128, 2)) | set(range(129, 1024, 46)))
    target = literal_coordinate_planes(7, max(controls))
    charged = 0
    for integer in controls:
        result = q7_odd_original_coefficient(integer)
        expected = sum(((plane >> integer) & 1) << i for i, plane in enumerate(target))
        assert result['value'] == expected
        charged += result['charged_steps']
    factors = {normalized_factors(m) for P in Q7_ODD_POLYNOMIALS for m in P}
    assert max(map(len, factors)) == 3
    return {"literal_odd_controls": len(controls), "largest_control_index": max(controls),
            "charged_enumeration_steps": charged,
            "distinct_quadratic_predicates": len(factors),
            "maximum_positive_variables_per_predicate": 3,
            "time_bound": "O(N poly(log(N+1)))", "workspace_bound": "O(log(N+1)) bits"}


def q5_calibration(limit):
    cut = (1 << (limit+1))-1
    V = sum(1 << (n*n) for n in range(1, isqrt(limit)+1))
    V ^= sum(1 << (5*n*n) for n in range(1, isqrt(limit//5)+1))
    V2 = frobenius_series(V, 2, cut)
    V3 = product_truncated(V, V2, cut)
    assert literal_coordinate_planes(5, limit) == [V ^ V2 ^ V3, 1 ^ V ^ V2]
    return {"phase": 5, "degree": 3, "checked_coefficients_including_constant": limit+1}


def main():
    result = {
        "scope": "Exact cubic parameter and odd-index theta engine with an ordinary Sturm proof; exponential evaluator; no faster factoring bound",
        "polynomial_identity": parameter_identity(),
        "kummer_identity": kummer_identity(),
        "original_parameter_controls": check_original_parameter(8193),
        "theta_span_exclusion": theta_span_certificate(),
        "positive_calibration": q5_calibration(16384),
        "odd_index_identity_controls": check_odd_identity(131072),
        "modular_lift_data": sturm_lift_data(),
        "charged_original_evaluator": check_original_evaluator(),
        "source_sha256": {name: sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/check_q7_cubic_parameter.py',
                                       'tools/tate_coordinate_compiler.py',
                                       'tools/lambert_theta_sections.py',
                                       'tools/cyclotomic_character_decoder.py',
                                       'research/q7-cubic-parameter.md')},
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
