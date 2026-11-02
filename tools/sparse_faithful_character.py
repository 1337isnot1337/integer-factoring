"""Sparse decoding from one supplied faithful character; no coefficient engine.

The full nonreduced power-basis observation is retained. Berlekamp splitting
acquires a small auxiliary field; consecutive moments and a Vandermonde
recurrence decode each interleaved sparse binary vector.
"""
from functools import cache
from math import isqrt, prod

from cyclotomic_character_decoder import (
    modular_compose, modular_inverse, modular_power, modular_product,
    parameters, polynomial_remainder,
)


def least_prime_factor(value):
    divisor = 3
    while divisor*divisor <= value:
        if value % divisor == 0:
            return divisor
        divisor += 2
    return value


def polynomial_gcd(left, right):
    while right:
        left, right = right, polynomial_remainder(left, right)
    return left


def binary_kernel(rows, columns):
    rows, pivots, rank = list(rows), [], 0
    for column in range(columns):
        pivot = next((j for j in range(rank, len(rows)) if rows[j] >> column & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for j in range(len(rows)):
            if j != rank and rows[j] >> column & 1:
                rows[j] ^= rows[rank]
        pivots.append(column)
        rank += 1
    return [((1 << free) | sum(((row >> free) & 1) << pivot
                              for row, pivot in zip(rows, pivots)))
            for free in range(columns) if free not in pivots]


@cache
def primitive_component(phi):
    """A deterministic irreducible factor via the classical Berlekamp kernel.

    phi is squarefree here. Kernel elements take values 0/1 on its field
    factors, and a kernel basis separates any two distinct factors.
    """
    dimension = phi.bit_length()-1
    columns = [modular_power(2, 2*j, phi) ^ (1 << j) for j in range(dimension)]
    rows = [sum(((column >> i) & 1) << j for j, column in enumerate(columns))
            for i in range(dimension)]
    factor = phi
    for value in binary_kernel(rows, dimension):
        part = polynomial_gcd(value, factor)
        if 1 < part < factor:
            factor = part
    return factor


def field_solve(matrix, vector, modulus):
    """Solve a nonsingular square system, or return None on singular input."""
    size = len(vector)
    rows = [list(row)+[value] for row, value in zip(matrix, vector)]
    for column in range(size):
        pivot = next((j for j in range(column, size) if rows[j][column]), None)
        if pivot is None:
            return None
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = modular_inverse(rows[column][column], modulus)
        rows[column] = [modular_product(x, inverse, modulus) for x in rows[column]]
        for j in range(size):
            if j != column and rows[j][column]:
                multiple = rows[j][column]
                rows[j] = [x ^ modular_product(multiple, y, modulus)
                           for x, y in zip(rows[j], rows[column])]
    return [row[-1] for row in rows]


def decode_odd_sparse(observed, odd_order, phi, sparsity):
    if not observed:
        return 0
    field = primitive_component(phi)
    moments = [0]+[modular_compose(observed, modular_power(2, j, field), field)
                   for j in range(1, 2*sparsity+1)]
    recurrence = None
    for size in range(sparsity, 0, -1):
        matrix = [[moments[i+j+1] for j in range(size)] for i in range(size)]
        coefficients = field_solve(matrix, moments[size+1:2*size+1], field)
        if coefficients is not None:
            recurrence = coefficients+[1]
            break
    if recurrence is None:
        raise ValueError('Observation has no certified sparse preimage')
    result, point = 0, 1
    for exponent in range(odd_order):
        value = 0
        for coefficient in reversed(recurrence):
            value = modular_product(value, point, field) ^ coefficient
        if value == 0:
            result |= 1 << exponent
        point = modular_product(point, 2, field)
    if result.bit_count() > sparsity or polynomial_remainder(result, phi) != observed:
        raise ValueError('Observation has no certified sparse binary preimage')
    return result


def decode_sparse_faithful(order, observed, sparsity):
    """Decode P of degree<order/2, weight<=s from P(T) mod Phi_order.

    For odd part u>1, require 2s<least_prime_factor(u). For u=1 the
    full observation already has the complete target dimension.
    """
    e, u, phi, modulus = parameters(order)
    if (not isinstance(sparsity, int) or sparsity < 0
            or not isinstance(observed, int) or observed < 0
            or observed.bit_length() >= modulus.bit_length()):
        raise ValueError('Canonical observation and nonnegative sparsity required')
    if u > 1 and 2*sparsity >= least_prime_factor(u):
        raise ValueError('The faithful sparse-distance certificate is insufficient')
    if u == 1:
        result = observed
    else:
        result = 0
        for residue in range(e):
            block = sum(((observed >> (residue+e*j)) & 1) << j
                        for j in range(phi.bit_length()-1))
            recovered = decode_odd_sparse(block, u, phi, sparsity)
            result |= sum(((recovered >> j) & 1) << (residue+e*j) for j in range(u))
    if result.bit_count() > sparsity or polynomial_remainder(result, modulus) != observed:
        raise ValueError('Observation has no certified sparse preimage')
    return result


def is_prime(value):
    return value >= 2 and all(value % d for d in range(2, isqrt(value)+1))


def acquire_rough_prime_levels(sparsity, count):
    """Deterministic CRT/least-prime acquisition, with exp(O(s))*poly(count) cost.

    Its uniform bound uses classical Linnik; no numerical Linnik constant is
    assumed. This intentionally simple reference tests candidates by division.
    """
    if sparsity < 1 or count < 1:
        raise ValueError('Positive sparsity and prime count required')
    forbidden = [p for p in range(3, 2*sparsity+1, 2) if is_prime(p)]
    odd_product = prod(forbidden)
    base_modulus = 4*odd_product
    base_residue = 2+odd_product*pow(odd_product, -1, 4)
    auxiliary = max(2*sparsity+1, count+1)
    while not is_prime(auxiliary):
        auxiliary += 1
    modulus = base_modulus*auxiliary
    records = []
    for residue in range(1, count+1):
        start = base_residue+base_modulus*((residue-base_residue)*pow(base_modulus, -1, auxiliary) % auxiliary)
        candidate, tested = start, 1
        while not is_prime(candidate):
            candidate += modulus
            tested += 1
        records.append(dict(auxiliary_residue=residue, progression_residue=start,
                            prime=candidate, tested_candidates=tested))
    return dict(sparsity=sparsity, count=count, forbidden_odd_primes=forbidden,
                base_modulus=base_modulus, auxiliary_prime=auxiliary,
                progression_modulus=modulus, levels=records)
