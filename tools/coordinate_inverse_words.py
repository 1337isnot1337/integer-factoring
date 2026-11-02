"""Bijective odd-candidate words from trace-dual field coordinates.

The ring model is supplied with a basic irreducible polynomial. These are
pointwise circuits and a conditional count consumer, not a fast domain sum.
The parity primitives and predecessor's no-wrap proof retain source credit.
"""
from math import prod

from teichmuller_cyclic_carries import UnramifiedRing
from witt_trace_word_circuits import parity_bit, trace_value


def coordinate_query_parameters(target, trial_bits=1):
    """Larger trial_bits requires completed trial division through 2^t.

    This function does not certify roughness. The consumer below performs
    the pretest before requesting a count at the smaller output precision.
    """
    if target < 5 or target % 2 == 0 or trial_bits < 1:
        raise ValueError('Odd target >=5 and positive trial width required')
    b = target.bit_length()
    a = (b+1)//2
    degree = a-1
    count_precision = min(degree, (b+trial_bits-1)//trial_bits)
    inverse_bits = a+b
    return dict(bits=b, candidate_bits=a, degree=degree,
                trial_bits=trial_bits, count_precision=count_precision,
                inverse_bits=inverse_bits,
                starting_precision=count_precision, exact_halvings=0,
                subtract_steps=inverse_bits, coordinate_traces=degree)


def trace_dual_basis(polynomial):
    """Acquire the trace-dual of 1,z,...,z^(d-1), in polynomial time.

    Returns binary-packed residue-field elements. Irreducibility is the
    same supplied-model assumption as UnramifiedRing, not searched here.
    """
    field = UnramifiedRing(1, polynomial)
    d = field.degree
    rows = [sum(trace_value(field, field.mul(1 << i, 1 << j)) << j
                for j in range(d)) | (1 << (d+i)) for i in range(d)]
    for j in range(d):
        pivot = next(i for i in range(j, d) if rows[i] >> j & 1)
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(d):
            if i != j and rows[i] >> j & 1:
                rows[i] ^= rows[j]
    inverse = [row >> d for row in rows]
    return tuple(sum(((inverse[i] >> j) & 1) << i for i in range(d))
                 for j in range(d))


def coordinate_bits(ring, dual_lifts, point):
    """Actual lifted traces and acquired integer parity polynomials.

    point is a Teichmuller input; dual_lifts are its model's lifted dual
    basis. No residue-word bit access is used to produce the output bits.
    """
    if len(dual_lifts) != ring.degree:
        raise ValueError('One lifted dual element per field coordinate required')
    bits = []
    for c in dual_lifts:
        conjugate, value = ring.mul(c,point), 0
        for _ in range(ring.degree):
            value = ring.add(value,conjugate)
            conjugate = ring.mul(conjugate,conjugate)
        assert not any(ring.coefficients(value)[1:])
        bits.append(parity_bit(value,ring.precision))
    return tuple(bits)


def cofactor_bit_circuit(target, digits, width, precision):
    """N/x modulo 2^width from odd x's supplied low-to-high bits.

    Every register stays a Boolean scalar via integer polynomial gates.
    Only the public target's bits are read; no variable bit extraction,
    division, comparison or residue-field enumeration occurs here.
    """
    modulus = 1 << precision
    residual = [(target >> j) & 1 for j in range(width)]
    answer = []
    def xor(a,b):
        return (a+b-2*a*b) % modulus
    for i in range(width):
        bit = residual[i]
        answer.append(bit)
        borrow = 0
        for j in range(i,width):
            left = residual[j]
            right = bit*digits[j-i] if j-i < len(digits) else 0
            residual[j] = xor(xor(left,right),borrow)
            borrow = ((1-left)*(right+borrow-right*borrow)
                      +left*right*borrow) % modulus
        assert residual[i] == 0
    return tuple(answer)


def coordinate_factor_prefix_value(target, prefix, ring, dual_lifts, point,
                                   trial_bits=1):
    plan = coordinate_query_parameters(target, trial_bits)
    d, b = plan['degree'], plan['bits']
    if ring.degree != d or ring.precision != plan['starting_precision']:
        raise ValueError('Ring model must match the complete query plan')
    if len(prefix) > d or any(c not in '01' for c in prefix):
        raise ValueError('Bad odd-candidate prefix')
    digits = coordinate_bits(ring, dual_lifts, point)
    cofactor = cofactor_bit_circuit(target,(1,)+digits,plan['inverse_bits'],
                                    ring.precision)
    candidate_nontrivial = 1-prod(1-bit for bit in digits)
    cofactor_in_range = prod(1-bit for bit in cofactor[b:])
    cofactor_nontrivial = 1-prod(1-bit for bit in cofactor[1:b])
    prefix_match = prod(digits[d-1-i] if bit == '1' else 1-digits[d-1-i]
                        for i, bit in enumerate(prefix))
    return candidate_nontrivial*cofactor_in_range*cofactor_nontrivial*prefix_match


def smallest_factor_from_coordinate_counts(target, count_query, trial_bits=1):
    """Charged pretest plus conditional exact-count consumer.

    count_query(N,prefix,t) must return the unweighted count, after this
    pretest. It is not implemented by this module. Prefixes describe the
    d-bit high part w of the odd candidate x=1+2w, MSB first.
    """
    if target < 2 or trial_bits < 1:
        raise ValueError('Target >=2 and positive trial width required')
    if target % 2 == 0:
        return 2, 0
    if target == 3:
        return 3, 0
    for divisor in range(3, min(1 << trial_bits, target)+1, 2):
        if target % divisor == 0:
            return divisor, 0
    plan = coordinate_query_parameters(target, trial_bits)
    if not count_query(target, '', trial_bits):
        return target, 1
    prefix = ''
    for _ in range(plan['degree']):
        prefix += '0' if count_query(target, prefix+'0', trial_bits) else '1'
    return 1+2*int(prefix, 2), plan['degree']+1
