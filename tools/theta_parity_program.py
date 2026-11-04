"""Certify parity-refined theta programs for odd original indices.

The original theta compiler splits the nonzero residues modulo an odd prime q
into unsigned pairs +a,-a.  At phases q=1 mod 4 its specified degree-four
module misses one native coordinate.  This variant keeps the two parity
classes modulo 2q separately.  Each variable is a positive coset theta sum
n = +-a or n = +-(q-a) modulo 2q, so the parity information is retained.

A successful program is certified by a conservative ordinary modular/Sturm
argument: the homogenized weight-two forms lie on
Gamma_0(64 q^2) intersect Gamma_1(8 q), whose SL2 index is
384 q (q^2-1); the inclusive weight-two Sturm bound is one sixth of that,
64 q (q^2-1).  Evaluation still enumerates quadratic representations and is
exponential in the input bit length.  No factoring exponent is improved.
"""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, isqrt

from theta_coordinate_program import (
    literal_planes, product_truncated, require_odd_prime,
)


def sturm_bound(prime):
    require_odd_prime(prime)
    # H=Gamma_0(64*q^2) intersect Gamma_1(8*q), SL2 index=384*q*(q^2-1).
    return 64*prime*(prime*prime-1)


def variable_residues(prime, variable):
    """The two residues modulo 2q represented by one parity variable."""
    require_odd_prime(prime)
    half = (prime-1)//2
    if not 0 <= variable < 2*half:
        raise ValueError("parity variable index out of range")
    a = variable//2+1
    if variable % 2 == 0:
        return frozenset((a % (2*prime), (-a) % (2*prime)))
    return frozenset(((prime-a) % (2*prime), (prime+a) % (2*prime)))


def theta_data(prime, limit):
    require_odd_prime(prime)
    width = prime-1
    variables = []
    for variable in range(width):
        residues = variable_residues(prime, variable)
        variables.append(sum(1 << (n*n) for n in range(1, isqrt(limit)+1)
                             if n % (2*prime) in residues))
    trace = sum(1 << (n*n) for n in range(1, isqrt(limit)+1))
    trace ^= sum(1 << (prime*n*n) for n in range(1, isqrt(limit//prime)+1))
    return variables, trace


def fit_program(prime, limit):
    """Exact F2 fit in the parity-refined degree <=4 module."""
    require_odd_prime(prime)
    variables, trace = theta_data(prime, limit)
    width = len(variables)
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
    return {"phase": prime, "fit_limit": limit, "variables": width,
            "columns": len(columns), "binary_rank": len(basis),
            "missing_coordinates": missing, "separating_covectors": witnesses,
            "program": program}


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
    width = prime-1
    columns = 1+comb(width+4, 4)
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


def normalized_factors(prime, monomial):
    return tuple((variable_residues(prime, variable), 1 << bit)
                 for variable, multiplicity in sorted(Counter(monomial).items())
                 for bit in range(multiplicity.bit_length())
                 if multiplicity >> bit & 1)


def constrained_parity(integer, factors, modulus):
    residues, scale = factors[0]
    if len(factors) == 1:
        if integer <= 0 or integer % scale:
            return 0, 1
        root = isqrt(integer//scale)
        return int(root*root*scale == integer and root % modulus in residues), 1
    minimum_rest = sum(min(residues)**2*scale for residues, scale in factors[1:])
    bound = isqrt(max(0, (integer-minimum_rest)//scale))
    parity = steps = 0
    for root in range(1, bound+1):
        steps += 1
        if root % modulus in residues:
            bit, charged = constrained_parity(integer-scale*root*root,
                                               factors[1:], modulus)
            parity ^= bit
            steps += charged
    return parity, steps


def original_coefficient(integer, result):
    """Evaluate a certified parity program; the enumerative cost is returned."""
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
            bit, charged = constrained_parity(
                integer, normalized_factors(prime, monomial), 2*prime)
            steps += charged
        if bit:
            value ^= instruction['native_coefficient']
    return {"value": value, "charged_steps": steps}
