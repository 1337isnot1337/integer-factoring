"""Collect affine trace indicators in a square-zero coefficient algebra."""
from joint_frobenius_carries import joint_trace
from trace_constraint_contraction import independent_constraints, trace_packets


class SquareZeroAlgebra:
    """E[z_0,...,z_(r-1)]/(z_i^2), packed in the squarefree monomial basis.

    This is NOT a field. Only addition (XOR), multiplication and nonnegative
    powering are used. Inverses of its arbitrary elements are not provided.
    """
    def __init__(self, field, rank):
        self.field = field
        self.dimension = 1 << rank
        self.bits = field.degree
        self.mask = field.size - 1

    def coefficients(self, value):
        return [(i, (value >> (i * self.bits)) & self.mask)
                for i in range(self.dimension) if (value >> (i * self.bits)) & self.mask]

    def mul(self, left, right):
        result = 0
        right_terms = self.coefficients(right)
        for i, a in self.coefficients(left):
            for j, b in right_terms:
                if not i & j:
                    result ^= self.field.mul(a, b) << ((i | j) * self.bits)
        return result

    def power(self, value, exponent):
        if exponent < 0:
            raise ValueError('Only nonnegative powers are defined')
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            value = self.mul(value, value)
            exponent >>= 1
        return result


def square_zero_constrained_trace(field, width, t, variables, numerators,
                                  denominators, constraints, nonzero=False):
    if width < 1 or t < 1 or field.degree % width or width % t or variables < 1:
        raise ValueError('Incompatible dimensions')
    independent = independent_constraints(field, width, variables, constraints)
    if independent is None:
        return 0, {'inconsistent': True, 'rank': None, 'algebra_dimension': 0}
    rank = len(independent)
    algebra = SquareZeroAlgebra(field, rank)
    factors = list(numerators)
    zero = (0,) * (variables * (width // t))
    for k, (coefficients, bit) in enumerate(independent):
        shift = field.degree * (1 << k)
        for row, terms in trace_packets(field, width, t, coefficients, bit):
            factors.append((row, [(1, zero)] + [(c << shift, h) for c, h in terms]))
    # Only the ORIGINAL E-valued denominator factors are expanded by the
    # inherited positive-power inverse identity. The new factors are numerators.
    packed = joint_trace(algebra, field.degree, width, t, variables,
                         factors, denominators, nonzero)
    answer = (packed >> (field.degree * ((1 << rank) - 1))) & (field.size - 1)
    return answer, {'inconsistent': False, 'rank': rank,
                    'algebra_dimension': 1 << rank, 'indicator_factors': len(factors) - len(numerators)}
