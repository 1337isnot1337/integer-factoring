"""Recover divisors from fixed-characteristic global residue polynomials.

This specializes Garg--Schost sparse interpolation.  The fast residue
evaluator is an explicit missing input, not implemented by this module.
"""

from check_fixed_characteristic_transfer import first_primes
from divisor_moment_recovery import smallest_integer_root


def support_polynomial(support, modulus):
    polynomial = [1]
    for residue in sorted(support):
        result = [0]*(len(polynomial)+1)
        for index, coefficient in enumerate(polynomial):
            result[index] = (result[index]-residue*coefficient) % modulus
            result[index+1] = (result[index+1]+coefficient) % modulus
        polynomial = result
    return polynomial


def recover_from_global_residues(integer, characteristic, divisor_bound, oracle):
    """Use A_q=sum_(d|N) Z^(d mod q) over F_p for small prime q.

    The callback returns the dense q-coefficient representative.  Its exact
    evaluation cost must be supplied separately.  D>=tau(N) is a promise.
    """
    if integer < 2 or characteristic < 2 or divisor_bound < 2:
        raise ValueError("N>=2, a fixed prime characteristic and D>=tau(N) are required")
    bits, bound = integer.bit_length(), divisor_bound
    bad_bound = bits*bound*(bound-1)//2
    query_count = bad_bound+bound*(bits+1)+3
    primes = first_primes(query_count)
    reconstruction_threshold = 1 << (bound*(bits+1)+1)
    rank, modulus, coefficients, retained = 0, 1, [], 0
    support_histogram = {}
    for prime in primes:
        values = tuple(oracle(prime))
        if len(values) != prime or any(not isinstance(value, int) or not 0 <= value < characteristic
                                      for value in values):
            raise ValueError("the oracle must return the canonical dense residue polynomial")
        support = [index for index, value in enumerate(values) if value]
        size = len(support)
        support_histogram[size] = support_histogram.get(size, 0)+1
        if size > bound:
            raise ValueError("the residue support exceeds the supplied divisor bound")
        if size < rank or size == 0:
            continue
        if size > rank:
            rank, modulus, coefficients, retained = size, 1, [0]*size, 0
        if modulus > reconstruction_threshold:
            continue
        polynomial = support_polynomial(support, prime)
        inverse = pow(modulus, -1, prime)
        coefficients = [left+modulus*((right-left)*inverse % prime)
                        for left, right in zip(coefficients, polynomial[:-1])]
        modulus *= prime
        retained += 1
    if rank < 2 or modulus <= reconstruction_threshold:
        raise ValueError("insufficient collision-free residue information")
    polynomial = [value if 2*value <= modulus else value-modulus for value in coefficients]+[1]
    quotient = [0]*rank
    quotient[-1] = 1
    for index in range(rank-2, -1, -1):
        quotient[index] = polynomial[index+1]+quotient[index+1]
    if polynomial[0]+quotient[0]:
        raise ValueError("the recovered divisor polynomial omits the root one")
    factor = smallest_integer_root(quotient, integer)
    if integer % factor:
        raise ValueError("the recovered least node is not a divisor")
    return {"smallest_prime_factor": factor, "divisor_count": rank,
            "divisor_polynomial": polynomial, "query_count": query_count,
            "largest_cyclic_degree": primes[-1], "retained_prime_count": retained,
            "crt_modulus_bits": modulus.bit_length(),
            "residue_support_histogram": dict(sorted(support_histogram.items()))}
