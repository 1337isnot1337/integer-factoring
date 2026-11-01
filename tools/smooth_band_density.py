"""Integer parameters for the multiband smooth-form existence theorem.

This module only computes bounds. It neither enumerates the prime base nor
finds a represented norm. Its analytic dependency is stated in the note.
"""
from math import factorial, isqrt


def largest_band_coefficient(bands):
    """Return an index, maximum coefficient, and full coefficient list.

    The polynomial is (1+z+...+z^(bands-1))^bands. Sliding-window
    convolution uses O(bands^3) integer additions and O(bands^2) entries.
    """
    if not isinstance(bands, int) or bands < 1:
        raise ValueError("The number of bands must be a positive integer.")
    coefficients = [1]
    for _ in range(bands):
        new = []
        window = 0
        for index in range(len(coefficients) + bands - 1):
            if index < len(coefficients):
                window += coefficients[index]
            if 0 <= index - bands < len(coefficients):
                window -= coefficients[index - bands]
            new.append(window)
        coefficients = new
    maximum = max(coefficients)
    return coefficients.index(maximum), maximum, coefficients


def balanced_band_parameters(discriminant_bits):
    """Compute explicit balanced parameters, or return failed eligibility.

    Large integers are returned only for eligible inputs. They have
    O(sqrt(b log b)) bits asymptotically. In particular this function does
    not construct the radius 2^r, whose bit length is linear in b.
    """
    b = discriminant_bits
    if not isinstance(b, int) or b < 2:
        raise ValueError("The discriminant bit length must be at least two.")
    log_bits = (b - 1).bit_length()
    radicand = 8 * b * log_bits
    smooth_bits = isqrt(radicand)
    if smooth_bits * smooth_bits < radicand:
        smooth_bits += 1
    radius_bits = (b + 1) // 2 - 3
    result = {
        'discriminant_bits': b,
        'ceil_log_two_bits': log_bits,
        'smoothness_exponent': smooth_bits,
        'radius_exponent': radius_bits,
        'eligible': False,
        'conditions': {'u_at_least_four': radius_bits >= 4 * smooth_bits},
    }
    if not result['conditions']['u_at_least_four']:
        return result

    q = 0
    while radius_bits > smooth_bits << q:
        q += 1
    denominator = q * smooth_bits
    s = ((q + 1) * radius_bits + denominator - 1) // denominator
    # A bit-length comparison is slightly stronger than comparison with
    # the power of two. It avoids constructing that power at this stage.
    result['conditions'].update({
        'wide_prime_bands': smooth_bits >= 4 * s * q,
        'analytic_shell_floor': (
            (1000000 * (b + 1)) ** 128).bit_length() <= smooth_bits,
        'distinct_prime_budget': (
            (512 * s * s * q) ** 4).bit_length() <= smooth_bits,
        'final_prime_budget': (
            (64 * s * smooth_bits) ** 8).bit_length() <= 5 * smooth_bits,
    })
    result.update({'eta_denominator': q, 'band_count': s})
    if not all(result['conditions'].values()):
        return result

    density_denominator = (256 * (b + 2) * (s + 1) * s * s * smooth_bits
                           * (256 * q) ** s * factorial(s))
    epsilon_bits = (4 * density_denominator - 1).bit_length()
    sample_base = 40000 * (b + epsilon_bits + 2) * density_denominator
    sample_cutoff = sample_base ** 32
    smooth_cutoff = 1 << smooth_bits
    largest_cutoff = max(sample_cutoff, smooth_cutoff)
    multiplicity = 2 * s + 3
    result.update({
        'eligible': True,
        'density_denominator': density_denominator,
        'bias_reciprocal_bits': epsilon_bits,
        'sample_base': sample_base,
        'sample_cutoff': sample_cutoff,
        'smooth_cutoff': smooth_cutoff,
        'largest_prime_cutoff': largest_cutoff,
        'norm_multiplicity': multiplicity,
        'relation_diagonal': multiplicity + 1,
    })
    return result


def parameter_summary(parameters):
    """Return a JSON-safe summary without enormous integers or float logs."""
    result = {key: value for key, value in parameters.items() if key in {
        'discriminant_bits', 'ceil_log_two_bits', 'smoothness_exponent',
        'radius_exponent', 'eligible', 'conditions', 'eta_denominator',
        'band_count', 'norm_multiplicity', 'relation_diagonal',
    }}
    if not parameters['eligible']:
        return result
    result.update({
        'density_denominator_bits': parameters['density_denominator'].bit_length(),
        'sample_cutoff_bits': parameters['sample_cutoff'].bit_length(),
        'largest_prime_cutoff_bits': parameters['largest_prime_cutoff'].bit_length(),
        'guaranteed_norm_bits': (parameters['sample_cutoff'].bit_length()
                                 + 2 * parameters['radius_exponent']),
        'accepted_norm_bits_upper_bound': (
            parameters['norm_multiplicity']
            * parameters['largest_prime_cutoff'].bit_length()),
    })
    return result
