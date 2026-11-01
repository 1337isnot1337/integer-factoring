"""Explicit coarse integer cutoffs under the pinned Hecke zero-free input."""
from fractions import Fraction


def ceil_log_two(value):
    """ceil(log2(value)) for positive rational value, with integer arithmetic."""
    value = Fraction(value)
    assert value > 0
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if exponent >= 0:
        if value.numerator > value.denominator << exponent:
            exponent += 1
    elif value.numerator << (-exponent) > value.denominator:
        exponent += 1
    return exponent


def generator_cutoff(discriminant):
    assert discriminant < 0 and discriminant % 4 in (0, 1)
    bits = (-discriminant).bit_length()
    return (10000 * (bits + 1)) ** 32


def sample_cutoff_bits(discriminant_bits, epsilon):
    """Compute cutoff and exact error envelope, never enumerate the primes."""
    epsilon = Fraction(epsilon)
    assert discriminant_bits >= 2 and 0 < epsilon <= Fraction(1, 2)
    density_bits = ceil_log_two(1 / epsilon)
    parameter = discriminant_bits + density_bits + 2
    ratio = 10000 * parameter / epsilon
    base = (ratio.numerator + ratio.denominator - 1) // ratio.denominator
    cutoff = base ** 32
    error = (Fraction(100 * (discriminant_bits + 1), base)
             + Fraction(16, base ** 15)
             + Fraction(16 * (discriminant_bits + 1), base ** 31))
    assert error <= epsilon / 32
    unrounded_bias = error / (Fraction(1, 6) - error)
    weight_bits = ceil_log_two(48 / epsilon)
    rounding_distance = Fraction(24, 1 << weight_bits)
    assert unrounded_bias < epsilon / 4
    assert rounding_distance <= epsilon / 2
    assert unrounded_bias + rounding_distance < epsilon
    return {
        'base': base, 'cutoff': cutoff, 'weight_bits': weight_bits,
        'normalized_error': error, 'unrounded_bias': unrounded_bias,
        'rounding_distance': rounding_distance,
    }
