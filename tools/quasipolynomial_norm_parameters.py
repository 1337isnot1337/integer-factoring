"""Exact conservative powered-prime parameters for a quasipolynomial base."""

from powered_prime_norm_parameters import ceil_log_two_integer


def divisor_table_exponent(height_exponent):
    """An integer E with sqrt((B+1)*tau(n)) <= 2^E for n <= 2^B.

The elementary split at U=floor(B/ceil(log2(B+1))^3) is used when
U>=2; a direct tau(n)<=2^B bound handles the finite small branch.
This bound does not include polynomial sorting and arithmetic factors.
"""
    b = height_exponent
    if not isinstance(b, int) or b < 1:
        raise ValueError('The height exponent must be a positive integer.')
    q = ceil_log_two_integer(b+1)
    cutoff = b//q**3
    if cutoff >= 2:
        lower_log = cutoff.bit_length()-1
        divisor_exponent = cutoff*q+(b+lower_log-1)//lower_log
    else:
        divisor_exponent = b
    return {
        'small_prime_split': cutoff,
        'divisor_count_exponent_bound': divisor_exponent,
        'table_size_exponent_bound': (divisor_exponent+q+1)//2,
    }


def quasipolynomial_norm_parameters(discriminant_bits):
    b = discriminant_bits
    if not isinstance(b, int) or b < 2:
        raise ValueError('The discriminant bit length must be at least two.')
    q = ceil_log_two_integer(b)
    w, v = q*q, 32*q*q
    r = (b+1)//2-3
    k = r//v
    a = ceil_log_two_integer(10000*(b+w+2))
    h = w-a
    conditions = {
        'positive_smooth_multiplicity': k >= 1,
        'positive_sample_bias_bits': h >= 1,
        'analytic_shell_floor': ceil_log_two_integer(1000000*(b+1)) <= w,
        'enough_distinct_primes': ceil_log_two_integer(64*v*max(k,1)) <= v,
    }
    result = {
        'discriminant_bits': b, 'logarithm_ceiling': q,
        'root_base_exponent': w, 'smoothness_exponent': v,
        'radius_exponent': r, 'smooth_form_multiplicity': k,
        'sample_bias_bits': h, 'conditions': conditions,
        'eligible': all(conditions.values()),
    }
    if not result['eligible']:
        return result
    residual = (b+1)//2-k*v
    density_upper_bits = (2+ceil_log_two_integer(b+2)+residual
                          +k*ceil_log_two_integer(32*v*k))
    t = (density_upper_bits+h-1)//h
    multiplicity = 2*k+t
    height = v*multiplicity
    result.update({
        'density_bias_bits_upper_bound': density_upper_bits,
        'walk_steps': t, 'norm_multiplicity': multiplicity,
        'accepted_norm_log_two_bound': height,
        'ordered_walk_enumeration_exponent_bound': (v+1)*t,
        **divisor_table_exponent(height),
    })
    return result
