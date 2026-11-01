"""Exact prime-form sample constructor for the fixed-strip research note.

The finite constructor accepts its prime cutoff as input.  It does not encode
an uncomputed numerical constant from the analytic theorem, certify arbitrary
cutoffs, or enumerate the class group.  The companion checker performs a
separate all-class enumeration only for small controls.
"""
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt

from binary_quadratic_forms import reduce_form


@lru_cache(None)
def primes_below(cutoff):
    """Deliberately elementary deterministic prime enumeration."""
    return tuple(n for n in range(3, cutoff, 2)
                 if all(n % d for d in range(3, isqrt(n) + 1, 2)))


def _log_series_bounds(ratio, terms):
    """Bounds for log((1+ratio)/(1-ratio)), with 0 <= ratio <= 1/3."""
    assert 0 <= ratio <= Fraction(1, 3) and terms >= 1
    square = ratio * ratio
    power = ratio
    lower = Fraction(0)
    for j in range(terms):
        lower += 2 * power / (2 * j + 1)
        power *= square
    # Every later denominator is at least the first omitted denominator.
    tail = 2 * power / ((2 * terms + 1) * (1 - square))
    return lower, lower + tail


@lru_cache(None)
def log_rational_bounds(integer, precision_bits):
    """Exact rational bounds with width at most 2**(-precision_bits).

    No floating logarithms or assumptions about a real-number floor occur.
    """
    assert integer >= 1 and precision_bits >= 1
    exponent = integer.bit_length() - 1
    terms = 1
    # A conservative common geometric-tail bound for both series.
    while (exponent + 1) * 9 * (1 << precision_bits) > 4 * 3 ** (2 * terms + 1):
        terms += 1
    log_two_low, log_two_high = _log_series_bounds(Fraction(1, 3), terms)
    scale = 1 << exponent
    ratio = Fraction(integer - scale, integer + scale)
    unit_low, unit_high = _log_series_bounds(ratio, terms)
    low = exponent * log_two_low + unit_low
    high = exponent * log_two_high + unit_high
    assert high - low <= Fraction(1, 1 << precision_bits)
    return low, high


@lru_cache(None)
def log_lower_dyadic(integer, precision_bits):
    """ell satisfies 0 <= log(integer)-ell/2**k <= 2**(1-k)."""
    low, high = log_rational_bounds(integer, precision_bits)
    scale = 1 << precision_bits
    scaled = low * scale
    numerator = scaled.numerator // scaled.denominator
    rounded = Fraction(numerator, scale)
    assert rounded <= low <= high
    assert high - rounded <= Fraction(2, scale)
    return numerator


def prime_form_sample(discriminant, cutoff, precision_bits=12):
    """Return exact integer masses and source records for both orientations.

    ``cutoff`` is x in the note; primes equal to x have zero weight and are
    omitted.  Primes dividing 2*discriminant are excluded without factoring
    the discriminant.  The dictionary may combine several source primes.
    """
    assert discriminant < 0 and discriminant % 4 in (0, 1)
    assert cutoff >= 4 and precision_bits >= 1
    masses = {}
    records = []
    for prime in primes_below(cutoff):
        if gcd(prime, discriminant) != 1:
            continue
        roots = [middle for middle in range(discriminant % 2, 2 * prime, 2)
                 if (middle * middle - discriminant) % (4 * prime) == 0]
        assert len(roots) in (0, 2)
        if not roots:
            continue
        log_numerator = log_lower_dyadic(prime, precision_bits)
        mass = log_numerator * prime * (cutoff - prime)
        assert mass > 0
        for middle in roots:
            raw = (prime, middle, (middle * middle - discriminant) // (4 * prime))
            assert gcd(gcd(*raw[:2]), raw[2]) == 1
            form = reduce_form(raw)
            masses[form] = masses.get(form, 0) + mass
            records.append((prime, middle, form, mass))
    return masses, records


def uniform_distance_certificate(forms, masses):
    """Exact l1 distance from uniform; it bounds every nontrivial bias.

    This function receives the full class list.  Obtaining that list is not
    part of the polynomial-cost sample construction or its analytic proof.
    """
    assert forms and set(masses) <= set(forms)
    total = sum(masses.values())
    assert total > 0
    order = len(forms)
    return Fraction(sum(abs(order * masses.get(form, 0) - total) for form in forms),
                    order * total)
