"""Construct and certify bounded-degree theta programs for odd original indices.

This is a partial compiler: it can certify a program, exclude its specified
finite module, or report failed validation. It does not prove that every
phase admits a program. Query evaluation enumerates quadratic representations.
"""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, isqrt

from cyclotomic_character_decoder import polynomial_product
from lambert_theta_sections import AdditivePhaseAlgebra
from tate_coordinate_compiler import compile_coordinate, embed_coefficient


def require_odd_prime(prime):
    if not isinstance(prime, int) or prime < 3 or prime % 2 == 0 or any(
            prime % d == 0 for d in range(3, isqrt(prime)+1, 2)):
        raise ValueError("an odd prime phase is required")


def sturm_bound(prime):
    require_odd_prime(prime)
    # H=Gamma_0(16*q^2) intersect Gamma_1(q), SL2 index=24*q*(q^2-1).
    return 4*prime*(prime*prime-1)


def product_truncated(left, right, cut):
    if left.bit_count() > right.bit_count():
        left, right = right, left
    return polynomial_product(left, right) & cut


def native_phase_weights(prime):
    """Polynomial-size inverse of the intrinsic embedding, not 2^h enumeration."""
    require_odd_prime(prime)
    plan = compile_coordinate(prime, intrinsic=True)
    basis = {}
    for i in range((prime-1)//2):
        value, native = embed_coefficient(1 << i, plan), 1 << i
        while value:
            pivot = value.bit_length()-1
            if pivot in basis:
                column, decoder = basis[pivot]
                value ^= column
                native ^= decoder
            else:
                basis[pivot] = value, native
                break
        assert value, "intrinsic embedding must be injective"
    phase = AdditivePhaseAlgebra(prime)
    weights = []
    for residue in range(prime):
        value = phase.phase(residue) ^ phase.phase(-residue)
        native = 0
        for pivot in sorted(basis, reverse=True):
            if value >> pivot & 1:
                column, decoder = basis[pivot]
                value ^= column
                native ^= decoder
        assert not value
        weights.append(native)
    return plan['constant'], weights


def literal_planes(prime, limit):
    """Finite odd-divisor sieve, with no exponential embedding lookup."""
    constant, weights = native_phase_weights(prime)
    values = [constant]+[0]*limit
    for divisor in range(1, limit+1, 2):
        weight = weights[divisor % prime]
        for index in range(divisor, limit+1, divisor):
            values[index] ^= weight
    return [sum(((v >> i) & 1) << n for n, v in enumerate(values))
            for i in range((prime-1)//2)]


def theta_data(prime, limit):
    # Each variable is the integral full coset sum n=+a modulo q,
    # equivalently positive n in the unsigned pair +a,-a.
    variables = [sum(1 << (n*n) for n in range(1, isqrt(limit)+1)
                     if min(n % prime, (-n) % prime) == a)
                 for a in range(1, (prime+1)//2)]
    trace = sum(1 << (n*n) for n in range(1, isqrt(limit)+1))
    trace ^= sum(1 << (prime*n*n) for n in range(1, isqrt(limit//prime)+1))
    return variables, trace


def fit_program(prime, limit):
    """Exact F2 fit for degree <=4 plus the scalar odd-square trace column."""
    require_odd_prime(prime)
    width = (prime-1)//2
    variables, trace = theta_data(prime, limit)
    cut = (1 << (limit+1))-1
    odd = sum(1 << n for n in range(1, limit+1, 2))
    columns = []
    monomials = []
    polynomials = {(): 1}
    basis = {}

    def insert(monomial, polynomial):
        relation = 1 << len(columns)
        value = polynomial & odd
        columns.append(value)
        monomials.append(monomial)
        while value:
            pivot = (value & -value).bit_length()-1
            if pivot in basis:
                column, decoder = basis[pivot]
                value ^= column
                relation ^= decoder
            else:
                basis[pivot] = value, relation
                break

    insert(None, trace)
    for degree in range(5):
        for monomial in combinations_with_replacement(range(width), degree):
            value = (1 if not monomial else product_truncated(
                polynomials[monomial[:-1]], variables[monomial[-1]], cut))
            polynomials[monomial] = value
            insert(monomial, value)
    assert len(columns) == 1+comb(width+4, 4)
    targets = literal_planes(prime, limit)
    coefficients = [0]*len(columns)
    missing = []
    witnesses = []
    for coordinate, target in enumerate(targets):
        residual = target & odd
        solution = 0
        for pivot in sorted(basis):
            if residual >> pivot & 1:
                column, decoder = basis[pivot]
                residual ^= column
                solution ^= decoder
        if residual:
            first = (residual & -residual).bit_length()-1
            covector = 1 << first
            for pivot in sorted(basis, reverse=True):
                if (covector & basis[pivot][0]).bit_count() % 2:
                    covector ^= 1 << pivot
            assert all((covector & column).bit_count() % 2 == 0 for column in columns)
            target_pairing = sum(((covector & t).bit_count() % 2) << i
                                 for i, t in enumerate(targets))
            assert target_pairing >> coordinate & 1
            missing.append(coordinate)
            witnesses.append({"coordinate": coordinate, "target_pairing": target_pairing,
                              "indices": [n for n in range(limit+1) if covector >> n & 1]})
        else:
            for j in range(len(columns)):
                if solution >> j & 1:
                    coefficients[j] ^= 1 << coordinate
    program = [{"monomial": None if m is None else list(m), "native_coefficient": a}
               for m, a in zip(monomials, coefficients) if a]
    return {"phase": prime, "fit_limit": limit, "columns": len(columns),
            "binary_rank": len(basis), "missing_coordinates": missing,
            "separating_covectors": witnesses, "program": program}


def validate_program(result, limit):
    prime = result['phase']
    variables, trace = theta_data(prime, limit)
    cut = (1 << (limit+1))-1
    values = [0]*((prime-1)//2)
    for instruction in result['program']:
        monomial = instruction['monomial']
        if monomial is None:
            value = trace
        else:
            value = 1
            for variable in monomial:
                value = product_truncated(value, variables[variable], cut)
        coefficient = instruction['native_coefficient']
        for i in range(len(values)):
            if coefficient >> i & 1:
                values[i] ^= value
    targets = literal_planes(prime, limit)
    odd = sum(1 << n for n in range(1, limit+1, 2))
    failures = []
    for coordinate, (target, value) in enumerate(zip(targets, values)):
        residual = (target ^ value) & odd
        if residual:
            first = (residual & -residual).bit_length()-1
            failures.append({"coordinate": coordinate, "first_failed_index": first})
    return {"checked_limit": limit, "odd_coefficients_per_coordinate": (limit+1)//2,
            "failures": failures}


def compile_program(prime):
    bound = sturm_bound(prime)
    columns = 1+comb((prime-1)//2+4, 4)
    training = min(bound, max(4096, 2*columns))
    result = fit_program(prime, training)
    result['Sturm_bound_inclusive'] = bound
    if result['missing_coordinates']:
        result['status'] = 'bounded_module_excluded'
        return result
    result['validation'] = validate_program(result, bound+128)
    result['status'] = ('validation_failed' if result['validation']['failures']
                        else 'certified_by_ordinary_modular_Sturm_proof')
    return result


def normalized_factors(monomial):
    return tuple((variable+1, 1 << bit)
                 for variable, multiplicity in sorted(Counter(monomial).items())
                 for bit in range(multiplicity.bit_length())
                 if multiplicity >> bit & 1)


def constrained_parity(integer, factors, prime):
    residue, scale = factors[0]
    if len(factors) == 1:
        if integer <= 0 or integer % scale:
            return 0, 1
        root = isqrt(integer//scale)
        return int(root*root*scale == integer and
                   min(root % prime, (-root) % prime) == residue), 1
    minimum_rest = sum(r*r*w for r, w in factors[1:])
    bound = isqrt(max(0, (integer-minimum_rest)//scale))
    parity = steps = 0
    for root in range(1, bound+1):
        steps += 1
        if min(root % prime, (-root) % prime) == residue:
            bit, charged = constrained_parity(integer-scale*root*root, factors[1:], prime)
            parity ^= bit
            steps += charged
    return parity, steps


def original_coefficient(integer, result):
    """Only accepts certified programs; enumerative cost is explicitly returned."""
    if result.get('status') != 'certified_by_ordinary_modular_Sturm_proof':
        raise ValueError("a certified full-coordinate program is required")
    if not isinstance(integer, int) or integer < 1 or integer % 2 == 0:
        raise ValueError("a positive odd original index is required")
    prime = result['phase']
    value = steps = 0
    for instruction in result['program']:
        monomial = instruction['monomial']
        if monomial is None:
            bit = int(isqrt(integer)**2 == integer)
            if integer % prime == 0:
                bit ^= int(isqrt(integer//prime)**2 == integer//prime)
            steps += 2
        elif not monomial:
            bit = 0
        else:
            bit, charged = constrained_parity(integer, normalized_factors(monomial), prime)
            steps += charged
        if bit:
            value ^= instruction['native_coefficient']
    return {"value": value, "charged_steps": steps}
