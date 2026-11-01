"""Prime-form roots with deterministic Tonelli--Shanks construction.

The helpers require a known prime. The factor-base constructor supplies
primes by a sieve. Polynomial worst-case nonresidue search uses the
pinned analytic input, as detailed in the companion research note.
"""
from dataclasses import dataclass
from math import isqrt

from binary_quadratic_forms import inverse_form, reduce_form


def primes_through(bound):
    if not isinstance(bound, int) or bound < 0:
        raise ValueError("The prime cutoff must be a nonnegative integer.")
    if bound < 2:
        return []
    sieve = bytearray(b'\x01') * (bound + 1)
    sieve[0:2] = b'\x00\x00'
    for divisor in range(2, isqrt(bound) + 1):
        if sieve[divisor]:
            start = divisor * divisor
            count = (bound - start) // divisor + 1
            sieve[start:bound + 1:divisor] = b'\x00' * count
    return [value for value in range(2, bound + 1) if sieve[value]]


def _bump(stats, field, count=1):
    if stats is not None:
        stats[field] = stats.get(field, 0) + count


def nonresidue_scan_bound(prime):
    """The explicit scan bound; primality is an input promise."""
    if not isinstance(prime, int) or prime < 3 or prime % 2 == 0:
        raise ValueError("The input must be a known odd prime.")
    return min(prime - 1, (10000 * (prime.bit_length() + 1)) ** 32)


def least_quadratic_nonresidue(prime, stats=None):
    """Return the least integer nonresidue, for a known odd prime."""
    bound = nonresidue_scan_bound(prime)
    for candidate in range(2, bound + 1):
        _bump(stats, 'nonresidue_tests')
        if pow(candidate, (prime - 1) // 2, prime) == prime - 1:
            if stats is not None:
                stats['largest_nonresidue'] = max(
                    stats.get('largest_nonresidue', 0), candidate)
            return candidate
    raise ValueError("The known-prime/analytic scan-bound contract failed.")


def modular_square_root(value, prime, stats=None):
    """Return a root in [0,p), or None for a nonsquare, for known prime p."""
    if not isinstance(prime, int) or prime < 2:
        raise ValueError("The modulus must be a known prime.")
    value %= prime
    if prime == 2 or value == 0:
        return value
    if prime % 2 == 0:
        raise ValueError("An even modulus larger than two is not prime.")
    _bump(stats, 'legendre_tests')
    if pow(value, (prime - 1) // 2, prime) != 1:
        return None
    if prime % 4 == 3:
        _bump(stats, 'direct_exponent_roots')
        return pow(value, (prime + 1) // 4, prime)

    odd = prime - 1
    valuation = 0
    while odd % 2 == 0:
        odd //= 2
        valuation += 1
    nonresidue = least_quadratic_nonresidue(prime, stats)
    correction = pow(nonresidue, odd, prime)
    root = pow(value, (odd + 1) // 2, prime)
    remainder = pow(value, odd, prime)
    order_bits = valuation
    while remainder != 1:
        index = 0
        power = remainder
        while power != 1 and index < order_bits:
            power = power * power % prime
            index += 1
            _bump(stats, 'order_test_squarings')
        if index >= order_bits:
            raise ValueError("The Tonelli--Shanks known-prime invariant failed.")
        adjustment = pow(correction, 1 << (order_bits - index - 1), prime)
        adjustment_squared = adjustment * adjustment % prime
        root = root * adjustment % prime
        remainder = remainder * adjustment_squared % prime
        correction = adjustment_squared
        order_bits = index
        _bump(stats, 'tonelli_rounds')
    assert root * root % prime == value
    return root


@dataclass(frozen=True)
class PrimeFormRecord:
    prime: int
    root: int
    positive: tuple
    negative: tuple


def build_prime_form_factor_base(discriminant, cutoff):
    """Construct all unramified split prime-form records up to cutoff."""
    if (not isinstance(discriminant, int) or discriminant >= 0
            or discriminant % 4 not in (0, 1)):
        raise ValueError("A negative quadratic-order discriminant is required.")
    records = []
    excluded = []
    stats = {}
    primes = primes_through(cutoff)
    for prime in primes:
        if discriminant % prime == 0:
            excluded.append(prime)
            continue
        if prime == 2:
            if discriminant % 8 != 1:
                continue
            root = 1
        else:
            root = modular_square_root(discriminant, prime, stats)
            if root is None:
                continue
            if root % 2 != discriminant % 2:
                root += prime
        assert 0 <= root < 2 * prime
        assert (root * root - discriminant) % (4 * prime) == 0
        form = reduce_form((prime, root, (root * root - discriminant) // (4 * prime)))
        records.append(PrimeFormRecord(prime, root, form, inverse_form(form)))
    stats.update({'primes_tested': len(primes), 'split_prime_records': len(records),
                  'excluded_discriminant_primes': len(excluded)})
    return {'records': records, 'excluded_primes': excluded, 'stats': stats}
