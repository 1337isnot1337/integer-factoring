"""Exact existence parameters for powered prime samples and short norms."""
from math import isqrt


def ceil_log_two_integer(value):
    if not isinstance(value, int) or value <= 0:
        raise ValueError("A positive integer is required.")
    return (value - 1).bit_length()


def powered_prime_parameters(discriminant_bits, width=None):
    """Return integer parameters, with explicit eligibility flags.

    The default width balances factor-base setup and one norm decoder.
    A supplied width must be a positive multiple of 32. No prime sample,
    radius, target norm, or class-group enumeration is materialized.
    """
    b = discriminant_bits
    if not isinstance(b, int) or b < 2:
        raise ValueError("The discriminant bit length must be at least two.")
    if width is None:
        w = isqrt(b // 2048)
        if 2048 * w * w < b:
            w += 1
        width = 32 * w
    elif not isinstance(width, int) or width <= 0 or width % 32:
        raise ValueError("The width must be a positive integer multiple of 32.")
    else:
        w = width // 32
    radius = (b + 1) // 2 - 3
    factors = radius // width
    padding_bits = ceil_log_two_integer(10000 * (b + w + 2))
    bias_bits = w - padding_bits
    result = {
        'discriminant_bits': b, 'smoothness_exponent': width,
        'root_base_exponent': w, 'radius_exponent': radius,
        'smooth_form_multiplicity': factors, 'sample_bias_bits': bias_bits,
        'eligible': False,
        'conditions': {
            'positive_smooth_multiplicity': factors >= 1,
            'positive_sample_bias_bits': bias_bits >= 1,
            'analytic_shell_floor': (
                (1000000 * (b + 1)) ** 32).bit_length() <= width,
            'enough_distinct_primes': (64 * width * max(factors, 1)).bit_length() <= width,
        },
    }
    if not all(result['conditions'].values()):
        return result

    radius_remainder = (b + 1) // 2 - factors * width
    density_denominator = ((b + 2) * (32 * width * factors) ** factors
                           << radius_remainder)
    needed_bits = ceil_log_two_integer(4 * density_denominator)
    walk_steps = (needed_bits + bias_bits - 1) // bias_bits
    multiplicity = 2 * factors + walk_steps
    sample_base = 10000 * (b + bias_bits + 2) * (1 << bias_bits)
    assert sample_base <= 1 << w
    result.update({
        'eligible': True, 'density_denominator': density_denominator,
        'density_bias_bits': needed_bits, 'walk_steps': walk_steps,
        'sample_base': sample_base, 'sample_cutoff': sample_base ** 32,
        'norm_multiplicity': multiplicity, 'relation_diagonal': multiplicity + 1,
        'guaranteed_norm_log_two_bound': 2 * radius + walk_steps * width,
        'accepted_norm_log_two_bound': multiplicity * width,
        'orientation_table_exponent': (multiplicity + 1) // 2,
        'batch_orientation_exponent': width + (multiplicity + 1) // 2,
    })
    return result


def powered_parameter_summary(parameters):
    excluded = {'density_denominator', 'sample_base', 'sample_cutoff'}
    result = {key: value for key, value in parameters.items() if key not in excluded}
    if parameters['eligible']:
        result['density_denominator_bits'] = parameters['density_denominator'].bit_length()
        result['sample_cutoff_bits'] = parameters['sample_cutoff'].bit_length()
    return result
