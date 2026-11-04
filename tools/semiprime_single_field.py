"""Decode supplied field coefficients under the distinct-semiprime promise.

This module does not compute divisor coefficients from the composite index.
"""
from functools import cache
from math import gcd, isqrt

from cyclotomic_character_decoder import modular_power, modular_product, parameters
from odd_character_residues import primitive_generator
from sparse_faithful_character import is_prime, primitive_component


@cache
def component_data(prime):
    """Acquire a primitive odd-order field and its signed-residue labels."""
    if prime < 7 or prime % 4 != 3 or not is_prime(prime):
        raise ValueError('A prime level q>=7, q=3 modulo 4, is required')
    _, order, phi, _ = parameters(prime-1)
    field = primitive_component(phi)
    generator = primitive_generator(prime)
    logs = {pow(generator, j, prime): j % order for j in range(prime-1)}
    powers = tuple(modular_power(2, j, field) for j in range(order))
    return order, field, generator, logs, powers


def squared_factor_residues(integer, prime, observed):
    """Return the unordered squared-factor residues from ONE field element.

    Correctness promises integer=p*r for distinct primes, prime not dividing
    integer, and observed=sum_(d|integer) omega**log_g(d) in the acquired field.
    Equal character roots are retained with multiplicity two.
    """
    _, field, generator, logs, powers = component_data(prime)
    if integer < 2 or gcd(integer, prime) != 1:
        raise ValueError('A positive coprime composite index is required')
    if not isinstance(observed, int) or observed < 0 or observed.bit_length() >= field.bit_length():
        raise ValueError('A canonical field coefficient is required')
    product = powers[logs[integer % prime]]
    total = observed ^ 1 ^ product
    roots = [j for j, value in enumerate(powers)
             if modular_product(value, value, field)
             ^ modular_product(total, value, field) ^ product == 0]
    if len(roots) not in (1, 2):
        raise ValueError('Coefficient has no semiprime preimage in the root subgroup')
    if len(roots) == 1:
        roots *= 2
    return tuple(sorted(pow(generator, 2*j, prime) for j in roots))


def recover_semiprime(integer, supplied):
    """Recover a promised distinct semiprime from supplied (prime, coefficient).

    Prime levels must be distinct and coprime to the input, and their product
    must exceed integer**2+1. No divisor or original-index oracle is supplied
    by this routine.
    """
    value, modulus, seen = 0, 1, set()
    for prime, observed in supplied:
        if prime in seen:
            raise ValueError('Distinct prime levels are required')
        seen.add(prime)
        residues = squared_factor_residues(integer, prime, observed)
        target = sum(residues) % prime
        value += modulus * ((target-value) * pow(modulus, -1, prime) % prime)
        modulus *= prime
    if modulus <= integer**2+1:
        raise ValueError('The supplied CRT product does not certify the integer lift')
    total = isqrt(value+2*integer)
    if total*total != value+2*integer:
        raise ValueError('The lifted factor-sum square is inconsistent')
    discriminant = total*total-4*integer
    if discriminant < 0:
        raise ValueError('The lifted quadratic has no integer factors')
    root = isqrt(discriminant)
    if root*root != discriminant or (total-root) % 2:
        raise ValueError('The lifted quadratic has no integer factors')
    left, right = (total-root)//2, (total+root)//2
    if left < 2 or left*right != integer:
        raise ValueError('Recovered factors do not certify the input')
    return left, right
