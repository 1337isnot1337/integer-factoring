"""Recover the smallest factor from supplied odd divisor power sums.

The exact coefficient/moment provider is an unresolved algorithmic input.
This module implements deterministic recovery after those moments arrive.
"""

from fractions import Fraction
from math import gcd, isqrt, lcm


def positive_hankel_rank(moments, bound):
    matrix = [[Fraction(moments[i+j]) for j in range(bound)] for i in range(bound)]
    for column in range(bound):
        pivot = matrix[column][column]
        if pivot == 0:
            if any(matrix[i][j] for i in range(column, bound) for j in range(column, bound)):
                raise ValueError("the supplied moments do not give a positive Hankel matrix")
            return column
        if pivot < 0:
            raise ValueError("the supplied moments do not give a positive Hankel matrix")
        for row in range(column+1, bound):
            factor = matrix[row][column]/pivot
            for j in range(column+1, bound):
                matrix[row][j] -= factor*matrix[column][j]
            matrix[row][column] = 0
    return bound


def solve_exact(matrix, right):
    size = len(right)
    rows = [[Fraction(value) for value in row]+[Fraction(value)]
            for row, value in zip(matrix, right)]
    for column in range(size):
        pivot = rows[column][column]
        if not pivot:
            raise ValueError("the moment system is singular")
        for row in range(column+1, size):
            factor = rows[row][column]/pivot
            for j in range(column, size+1):
                rows[row][j] -= factor*rows[column][j]
    solution = [Fraction(0)]*size
    for row in range(size-1, -1, -1):
        solution[row] = (rows[row][-1]-sum(rows[row][j]*solution[j]
                                         for j in range(row+1, size)))/rows[row][row]
    return solution


def primitive_integer_polynomial(coefficients):
    denominator = lcm(*(value.denominator for value in coefficients))
    integers = [int(value*denominator) for value in coefficients]
    content = gcd(*integers)
    return [value//content for value in integers] if content else []


def polynomial_remainder(left, right):
    output = list(map(Fraction, left))
    while output and len(output) >= len(right):
        factor = output[-1]/right[-1]
        offset = len(output)-len(right)
        for i, coefficient in enumerate(right):
            output[offset+i] -= factor*coefficient
        while output and not output[-1]:
            output.pop()
    return output


def sturm_chain(polynomial):
    derivative = [i*value for i, value in enumerate(polynomial)][1:]
    chain = [list(polynomial), derivative]
    while chain[-1]:
        remainder = polynomial_remainder(chain[-2], chain[-1])
        if not remainder:
            break
        chain.append(primitive_integer_polynomial([-value for value in remainder]))
    return chain


def sign_at_half_integer(polynomial, odd_numerator):
    value, denominator_power = polynomial[-1], 2
    for coefficient in reversed(polynomial[:-1]):
        value = value*odd_numerator+coefficient*denominator_power
        denominator_power *= 2
    return (value > 0)-(value < 0)


def sturm_variations(chain, odd_numerator):
    signs = [sign_at_half_integer(polynomial, odd_numerator) for polynomial in chain]
    signs = [sign for sign in signs if sign]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def smallest_integer_root(polynomial, upper_bound):
    """Isolate the smallest root in [2,upper_bound] using half-integer cuts."""
    chain = sturm_chain(polynomial)
    initial = sturm_variations(chain, 3)
    if initial == sturm_variations(chain, 2*upper_bound+1):
        raise ValueError("no root in the supplied interval")
    lower, upper = 2, upper_bound
    while lower < upper:
        middle = (lower+upper)//2
        if initial > sturm_variations(chain, 2*middle+1):
            upper = middle
        else:
            lower = middle+1
    value = 0
    for coefficient in reversed(polynomial):
        value = value*lower+coefficient
    if value:
        raise ValueError("the isolated root is not an integer")
    return lower


def factor_from_odd_moments(integer, divisor_count_bound, moment):
    """Use a_j=sum_(d|N) d^(2j+3), 0<=j<2D, where D>=tau(N).

    A callback must supply these exact values; finite checkers compute them
    from an explicitly charged divisor list. No fast provider is supplied.
    """
    if integer < 2 or divisor_count_bound < 2:
        raise ValueError("N>=2 and a divisor-count bound at least two are required")
    bound = divisor_count_bound
    moments = [moment(j) for j in range(2*bound)]
    if any(not isinstance(value, int) or value <= 0 for value in moments):
        raise ValueError("positive integer moments are required")
    rank = positive_hankel_rank(moments, bound)
    if rank < 2:
        raise ValueError("the divisor measure of N>=2 has at least two nodes")
    matrix = [[moments[i+j] for j in range(rank)] for i in range(rank)]
    coefficients = solve_exact(matrix, [-moments[i+rank] for i in range(rank)])
    if any(value.denominator != 1 for value in coefficients):
        raise ValueError("the divisor-node recurrence must be integral")
    polynomial = [int(value) for value in coefficients]+[1]
    for shift in range(2*bound-rank):
        if sum(coefficient*moments[shift+i] for i, coefficient in enumerate(polynomial)):
            raise ValueError("the supplied moment sequence fails its recurrence")
    # Remove the known divisor node 1; the next node is p^2 for the least p|N.
    quotient = [0]*rank
    quotient[-1] = polynomial[-1]
    for index in range(rank-2, -1, -1):
        quotient[index] = polynomial[index+1]+quotient[index+1]
    if polynomial[0]+quotient[0]:
        raise ValueError("the divisor measure must include the node one")
    node = smallest_integer_root(quotient, integer*integer)
    factor = isqrt(node)
    if factor*factor != node or integer % factor:
        raise ValueError("the least nontrivial node is not a squared divisor of N")
    return {"smallest_prime_factor": factor, "divisor_count": rank,
            "moment_queries": 2*bound, "largest_eisenstein_weight": 4*bound+2,
            "node_polynomial": polynomial,
            "maximum_moment_bits": max(value.bit_length() for value in moments),
            "maximum_polynomial_coefficient_bits": max(abs(value).bit_length() for value in polynomial)}


def factor_with_odd_moments(integer, moment):
    """Trial divide to 2^ceil(sqrt(bits(N))), then use a uniform rank bound."""
    if integer < 2:
        raise ValueError("N>=2 is required")
    bits = integer.bit_length()
    threshold_bits = isqrt(bits-1)+1
    trial_bound = 1 << threshold_bits
    for candidate in range(2, min(trial_bound, integer)+1):
        if integer % candidate == 0:
            return {"smallest_prime_factor": candidate, "moment_queries": 0,
                    "trial_bound": trial_bound, "stage": "trial_division"}
    divisor_bound = 1 << ((bits+threshold_bits-1)//threshold_bits)
    result = factor_from_odd_moments(integer, divisor_bound, moment)
    result.update({"trial_bound": trial_bound, "stage": "supplied_moments"})
    return result
