"""Exact finite checks for the fixed-characteristic factoring transfer.

The finite prefix tests use a trial-division divisor oracle. They verify the
transfer driver, not the hypothetical fast weighted evaluator.
"""

from functools import lru_cache
from itertools import combinations, product
from math import comb, isqrt
import hashlib
import json
from pathlib import Path


@lru_cache(None)
def first_primes(count):
    """Sieve with doubling; the proof gives the ceiling 16*(count+1)^2."""
    assert count >= 1
    limit = max(16, 2 * count)
    ceiling = 16 * (count + 1) ** 2
    while True:
        assert limit <= ceiling
        sieve = bytearray(b"\x01") * (limit + 1)
        sieve[0:2] = b"\x00\x00"
        for prime in range(2, isqrt(limit) + 1):
            if sieve[prime]:
                sieve[prime * prime:limit + 1:prime] = (
                    b"\x00" * ((limit - prime * prime) // prime + 1))
        found = tuple(i for i in range(2, limit + 1) if sieve[i])
        if len(found) >= count:
            return found[:count]
        limit = min(2 * limit, ceiling)


def quotient(poly, characteristic, degree):
    """Sparse coordinates of the canonical residue modulo Z^degree-1."""
    assert characteristic >= 2 and degree >= 1
    result = {}
    for exponent, coefficient in poly.items():
        assert exponent >= 0
        index = exponent % degree
        result[index] = (result.get(index, 0) + coefficient) % characteristic
    return {e: c for e, c in result.items() if c}


def sparse_nonzero(poly, characteristic, n, support_bound):
    assert n >= 1 and support_bound >= 1
    reduced = {e: c % characteristic for e, c in poly.items()
               if c % characteristic}
    assert len(reduced) <= support_bound
    assert all(0 <= e < 2 ** n for e in reduced)
    primes = first_primes(n * (support_bound - 1) + 1)
    return any(quotient(reduced, characteristic, q) for q in primes)


def divisor_sparsity_exponent(n):
    """The computable cap 2^u, without constructing that large integer."""
    assert n >= 1
    k = n.bit_length() - 1
    return n if k < 18 else (3 * n + k - 1) // k


def recover_min_factor(number, characteristic, support_bound, oracle):
    """Recover from an exact weighted quotient oracle; every call is charged."""
    assert number >= 2 and support_bound >= 1
    n = number.bit_length()
    primes = first_primes(n * (support_bound - 1) + 1)
    prefix = 0
    evaluations = 0
    for k in range(n - 1, -1, -1):
        child = 2 * prefix
        found = False
        for q in primes:
            residue = dict(oracle(k, child, q))
            evaluations += 1
            if (1 >> k) == child:
                index = 1 % q
                residue[index] = (residue.get(index, 0) - 1) % characteristic
            if any(value % characteristic for value in residue.values()):
                found = True
                break
        prefix = child if found else child + 1
    return prefix, evaluations


def finite_divisors(number):
    """Slow exact oracle support used only in finite verification."""
    return tuple(d for d in range(1, number + 1) if number % d == 0)


def main():
    stats = {}
    polynomial_cases = 0
    isolation_checks = 0
    for characteristic, exponents, max_terms in ((2, 10, 4), (3, 8, 3), (5, 6, 3)):
        n = (exponents - 1).bit_length()
        assert not sparse_nonzero({}, characteristic, n, 1)
        for terms in range(1, max_terms + 1):
            primes = first_primes(n * (terms - 1) + 1)
            for support in combinations(range(exponents), terms):
                e0 = support[0]
                isolated = [q for q in primes
                            if all((e - e0) % q for e in support[1:])]
                assert isolated
                isolation_checks += 1
                for coefficients in product(range(1, characteristic), repeat=terms):
                    poly = dict(zip(support, coefficients))
                    assert sparse_nonzero(poly, characteristic, n, terms)
                    polynomial_cases += 1
    stats["nonzero_sparse_polynomials"] = polynomial_cases
    stats["support_isolation_checks"] = isolation_checks

    # Monic quotient arithmetic does not require a reduced ring or q != p.
    quotient_checks = 0
    for characteristic in (2, 3, 5):
        for q in (2, 3, 5, 7):
            for a in range(17):
                for b in range(17):
                    assert quotient({a + b: 1}, characteristic, q) == {
                        ((a % q) + (b % q)) % q: 1}
                    quotient_checks += 1
    stats["monomial_product_reductions"] = quotient_checks

    # Negative controls: parity, one quotient, or base-field values lose support.
    assert (1 + 1) % 2 == 0
    assert quotient({1: 1, 3: 1}, 2, 2) == {}
    assert quotient({1: 1, 3: 1}, 2, 3)
    assert all((x * x + x) % 2 == 0 for x in range(2))
    assert quotient({1: 1, 2: 1}, 2, 3)
    stats["negative_controls"] = [
        "unweighted parity of two terms is zero",
        "one bad exponent modulus misses two terms",
        "all base-field values miss a nonzero polynomial",
    ]

    for count in range(1, 257):
        primes = first_primes(count)
        assert len(primes) == count
        assert primes[-1] <= 16 * (count + 1) ** 2
        assert 8 * (count + 1) ** 2 >= count * (4 + 2 * count)
    for m in range(1, 129):
        central = comb(2 * m, m)
        assert central >= 2 ** m
        for prime in first_primes(2 * m):
            if prime > 2 * m:
                break
            value, valuation = central, 0
            while value % prime == 0:
                value //= prime
                valuation += 1
            assert prime ** valuation <= 2 * m
    stats["prime_bound_cases"] = 256
    stats["central_binomial_cases"] = 128

    for k in range(18, 1025):
        assert k * (k + 2) <= 2 ** (k // 2)
        for n in (2 ** k, 2 ** k + 1, 2 ** (k + 1) - 1):
            assert 2 ** ((k + 1) // 2) * (k + 2) * k <= n
            u = divisor_sparsity_exponent(n)
            assert u * k >= 3 * n
            assert (u - 1) * k < 3 * n
    stats["divisor_bound_integer_cases"] = 3 * (1025 - 18)

    prefix_inputs = prefix_evaluations = 0
    for characteristic in (2, 3):
        for number in range(2, 513):
            divisors = finite_divisors(number)
            # A separate elementary cap shortens these small tests:
            # pairing d with N/d gives tau(N) <= 2*floor(sqrt(N)).
            support_bound = 2 * isqrt(number)
            assert len(divisors) <= support_bound

            def oracle(k, prefix, q):
                return quotient({d: 1 for d in divisors if d >> k == prefix},
                                characteristic, q)

            found, calls = recover_min_factor(
                number, characteristic, support_bound, oracle)
            assert found == divisors[1]
            prefix_inputs += 1
            prefix_evaluations += calls
    stats["prefix_inputs"] = prefix_inputs
    stats["weighted_quotient_calls"] = prefix_evaluations
    stats["status"] = "all finite exact checks passed; no faster evaluator tested"
    stats["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
