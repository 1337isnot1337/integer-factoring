"""Pointwise native-word circuits with explicitly charged exact halving.

This is NOT a contraction algorithm. The original polynomial carry API
does not implement these precision-losing division gates. Ring/Frobenius
arithmetic and H_k retain the precision backend's preceding credit.
"""
from math import prod
from functools import cache

from teichmuller_cyclic_carries import parity_lift_coefficients


@cache
def parity_coefficients(precision):
    return tuple(parity_lift_coefficients(precision))


def word_query_parameters(bits, degree):
    if bits < 1 or degree < 1 or degree % 2 != 1:
        raise ValueError('Positive word size and odd coefficient degree required')
    genus = ((1 << (2*bits))-3*(1 << bits)+2)//6
    if (1 << degree) <= 4*genus*genus:
        raise ValueError('The stated trace-surjectivity certificate is insufficient')
    count_precision = 2*degree+1
    zero_bits = 2*bits+1
    return dict(bits=bits, degree=degree, genus=genus,
                hasse_gap=(1 << degree)-4*genus*genus,
                count_precision=count_precision, zero_bits=zero_bits,
                starting_precision=count_precision+bits+zero_bits)


def trace_value(ring, point):
    """Actual lifted trace; callers supply a Teichmuller domain point."""
    answer, conjugate = 0, point
    for _ in range(ring.degree):
        answer = ring.add(answer, conjugate)
        conjugate = ring.frobenius(conjugate)
    assert not any(ring.coefficients(answer)[1:])
    return answer


def parity_bit(value, precision):
    """Evaluate the acquired integer H_k, rather than reading a bit."""
    if precision < 1:
        raise ValueError('Positive scalar precision required')
    answer, modulus = 0, 1 << precision
    for coefficient in reversed(parity_coefficients(precision)):
        answer = (answer*value+coefficient) % modulus
    assert answer in (0, 1)
    return answer


def extract_bits(value, precision, number):
    """Horner parity gates and known exact divisions, losing one bit each."""
    if not 0 <= number < precision:
        raise ValueError('Insufficient precision for the stated halving depth')
    value %= 1 << precision
    bits = []
    for _ in range(number):
        bit = parity_bit(value, precision)
        numerator = value-bit
        assert numerator % 2 == 0
        value = numerator//2
        precision -= 1
        bits.append(bit)
    return tuple(bits), precision


def factor_prefix_value(target, prefix, left_trace, right_trace, degree):
    """Pointwise factor/prefix predicate, 0 or 1 at the charged count precision.

    A prefix is written most-significant bit first. This procedure never
    finds or sums domain points and makes no fast-evaluation claim for sums.
    """
    if target < 2:
        raise ValueError('Factoring target must be at least two')
    b = target.bit_length()
    if len(prefix) > b or any(x not in '01' for x in prefix):
        raise ValueError('Bad binary prefix')
    plan = word_query_parameters(b, degree)
    precision = plan['starting_precision']
    left, remaining = extract_bits(left_trace, precision, b)
    right, right_remaining = extract_bits(right_trace, precision, b)
    assert remaining == right_remaining == plan['count_precision']+plan['zero_bits']
    x = sum(bit << j for j, bit in enumerate(left))
    y = sum(bit << j for j, bit in enumerate(right))
    low_bits, final_precision = extract_bits(x*y-target, remaining, plan['zero_bits'])
    assert final_precision == plan['count_precision']
    product_zero = prod(1-bit for bit in low_bits)
    left_nontrivial = 1-prod(1-bit for bit in left[1:])
    right_nontrivial = 1-prod(1-bit for bit in right[1:])
    prefix_match = prod(left[b-1-i] if digit == '1' else 1-left[b-1-i]
                        for i, digit in enumerate(prefix))
    return product_zero*left_nontrivial*right_nontrivial*prefix_match
