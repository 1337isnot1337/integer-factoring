"""One-bit native factoring queries with the existing polynomial hash family.

No parity contractor is implemented. The control field/hash acquisition
retains check_bit_tagged_witness_recovery.py's implementation and credit.
"""
from math import prod

from check_bit_tagged_witness_recovery import field_for_bound
from coordinate_inverse_words import coordinate_bits, cofactor_bit_circuit


def native_parity_state(target, ring, dual_lifts, point):
    """Actual pure-field coordinate bits and unprefixed factor flag."""
    b = target.bit_length()
    a = (b+1)//2
    if target < 5 or target % 2 == 0 or ring.precision != 1 or ring.degree != a-1:
        raise ValueError('Odd target >=5 and matching precision-one field required')
    digits = coordinate_bits(ring,dual_lifts,point)
    cofactor = cofactor_bit_circuit(target,(1,)+digits,a+b,1)
    flag = ((1-prod(1-bit for bit in digits))
            *prod(1-bit for bit in cofactor[b:])
            *(1-prod(1-bit for bit in cofactor[1:b])))
    return digits,flag


def masked_parity_value(digits, flag, basis, control_degree, bucket, tag=None):
    """Direct polynomial Boolean mask; control-field values are constants.

    No embedding of the control field into the coefficient field is used.
    A tag multiplies by one coordinate bit. Sum this scalar over the full
    native field to obtain the required parity-query answer.
    """
    if len(basis) != len(digits) or not 0 <= bucket < (1 << control_degree):
        raise ValueError('Hash dimensions or bucket invalid')
    if tag is not None and not 0 <= tag < len(digits):
        raise ValueError('Invalid coordinate tag')
    value = flag
    for j in range(control_degree):
        bit = sum(digit for digit,c in zip(digits,basis) if c >> j & 1) % 2
        value = value*(1+bit+((bucket >> j) & 1)) % 2
    if tag is not None:
        value *= digits[tag]
    return value


def smallest_factor_from_native_parities(target, parity_query, trial_bits=1):
    """Complete conditional recovery with a charged pretest.

    parity_query(N,basis,h,bucket,tag) supplies one scalar parity of the
    preceding pointwise composition. The function is an assumed provider,
    not implemented by this module. Prefix counts are not requested.
    """
    if target < 2 or trial_bits < 1:
        raise ValueError('Target >=2 and positive trial width required')
    if target % 2 == 0:
        return 2,0
    if target == 3:
        return 3,0
    for divisor in range(3,min(1 << trial_bits,target)+1,2):
        if target % divisor == 0:
            return divisor,0
    b = target.bit_length()
    d = (b+1)//2-1
    bound = 1 << min(d,(b+trial_bits-1)//trial_bits)
    control = field_for_bound(d,bound)
    best,queries = target,0
    for point in range(control.size):
        basis = control.hash_basis(point,d)
        for bucket in range(control.size):
            scalar = parity_query(target,basis,control.degree,bucket,None)
            queries += 1
            if scalar != 1:
                continue
            word = 0
            for i in range(d):
                bit = parity_query(target,basis,control.degree,bucket,i)
                queries += 1
                if bit not in (0,1):
                    raise ValueError('Parity provider must return a bit')
                word |= bit << i
            candidate = 1+2*word
            if 1 < candidate < best and target % candidate == 0:
                best = candidate
    return best,queries
