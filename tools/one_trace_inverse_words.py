"""One-input, half-width native factoring predicates; no domain contractor.

Reuses the Witt trace and the explicitly charged exact-halving gates. The
cofactor is a division-free polynomial inverse, followed by an interval test.
"""
from math import prod

from witt_trace_word_circuits import extract_bits, word_query_parameters


def inverse_query_parameters(target, degree=None):
    if target < 3 or target % 2 == 0:
        raise ValueError('The inverse query requires an odd target at least three')
    bits = target.bit_length()
    candidate_bits = (bits+1)//2
    degree = 4*candidate_bits+1 if degree is None else degree
    certificate = word_query_parameters(candidate_bits, degree)
    count_precision = degree+1
    inverse_bits = bits+candidate_bits
    starting_precision = count_precision+candidate_bits+inverse_bits
    inverse_precision = starting_precision-candidate_bits
    return dict(bits=bits, candidate_bits=candidate_bits, degree=degree,
                genus=certificate['genus'], hasse_gap=certificate['hasse_gap'],
                count_precision=count_precision, inverse_bits=inverse_bits,
                starting_precision=starting_precision,
                inverse_precision=inverse_precision,
                inverse_terms=(inverse_precision+2)//3,
                exact_halvings=candidate_bits+inverse_bits)


def odd_inverse_polynomial(value, precision):
    """x sum_(j=0)^(ceil(h/3)-1) (1-x*x)^j, also defined on even x.

    For odd x, x*x=1 mod 8, so this is its inverse modulo 2^h.
    The caller masks even inputs; this is not a total unit-inverse oracle.
    """
    if precision < 1:
        raise ValueError('Positive inverse precision required')
    modulus = 1 << precision
    ratio = (1-value*value) % modulus
    answer = 0
    for _ in range((precision+2)//3):
        answer = (answer*ratio+1) % modulus
    return value*answer % modulus


def inverse_factor_prefix_value(target, prefix, scalar_trace, degree=None):
    """Actual Horner/known-division pointwise predicate, not a fast sum.

    Prefixes describe the candidate's ceil(bits(N)/2)-bit word, MSB first.
    The scalar_trace must be supplied at the plan's starting precision.
    """
    plan = inverse_query_parameters(target, degree)
    a, b = plan['candidate_bits'], plan['bits']
    if len(prefix) > a or any(x not in '01' for x in prefix):
        raise ValueError('Bad candidate prefix')
    digits, remaining = extract_bits(scalar_trace, plan['starting_precision'], a)
    x = sum(bit << j for j, bit in enumerate(digits))
    assert remaining == plan['inverse_precision']
    inverse = odd_inverse_polynomial(x, remaining)
    cofactor, final_precision = extract_bits(target*inverse, remaining, plan['inverse_bits'])
    assert final_precision == plan['count_precision']
    odd_proper_candidate = digits[0]*(1-prod(1-bit for bit in digits[1:]))
    cofactor_in_range = prod(1-bit for bit in cofactor[b:])
    cofactor_nontrivial = 1-prod(1-bit for bit in cofactor[1:b])
    prefix_match = prod(digits[a-1-i] if bit == '1' else 1-digits[a-1-i]
                        for i, bit in enumerate(prefix))
    return odd_proper_candidate*cofactor_in_range*cofactor_nontrivial*prefix_match


def smallest_factor_from_inverse_counts(target, count_query):
    """Conditional consumer; count_query(N,prefix) must supply the exact count.

    This module supplies no implementation of that missing count primitive.
    """
    if target < 2:
        raise ValueError('Factoring target must be at least two')
    if target % 2 == 0:
        return 2, 0
    queries = 1
    if not count_query(target, ''):
        return target, queries
    prefix = ''
    for _ in range((target.bit_length()+1)//2):
        queries += 1
        prefix += '0' if count_query(target, prefix+'0') else '1'
    return int(prefix, 2), queries
