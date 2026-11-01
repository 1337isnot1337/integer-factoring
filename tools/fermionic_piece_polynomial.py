#!/usr/bin/env python3
"""Complete characteristic-two crossing polynomials by a Wick Pfaffian.

This module computes finite tensors. It does not certify an asymptotic rank
bound. See research/fermionic-piece-polynomial.md for the ordinary proof.
The underlying normally ordered operators are from anthropics/formal-math,
revision 16a5ce41304fcdaea4e0be7c132e1c361dc26da5 (Appendix B, Apache-2.0).
"""

from functools import lru_cache


def regular_two(parts: tuple[int, ...]) -> bool:
    return len(set(parts)) == len(parts) and all(part > 0 and part % 2 for part in parts)


@lru_cache(None)
def binomial_polynomial(exponent: int) -> int:
    """Bit j is [z^j](1+z)^exponent over F_2."""
    if exponent < 0:
        raise ValueError("negative polynomial exponent")
    if not exponent:
        return 1
    bit = exponent & -exponent
    previous = binomial_polynomial(exponent - bit)
    return previous ^ (previous << bit)


@lru_cache(maxsize=100000)
def polynomial_product(left: int, right: int) -> int:
    """Carry-free multiplication of two F_2[z] bitsets."""
    if left.bit_count() > right.bit_count():
        left, right = right, left
    result = 0
    while left:
        bit = left & -left
        result ^= right << (bit.bit_length() - 1)
        left ^= bit
    return result


def creation_coefficient(label: int, mode: int) -> int:
    return int((label - 1) & ~(mode - 1) == 0)


def annihilation_coefficient(label: int, mode: int) -> int:
    return int((label - 1) & mode == 0)


@lru_cache(maxsize=100000)
def monomial_piece_polynomial(
    operator_a: tuple[int, ...], operator_b: tuple[int, ...],
    output: tuple[int, ...], input_: tuple[int, ...],
) -> tuple[int, int]:
    """Return (bitset, offset) for sum_s [output]Q_s(a,b)input * z^s.

    Every argument must be a distinct-odd monomial. The last three have
    a common weight h. The first may have any weight A; the generalized
    locality parameter is A+h. Bit j represents exponent offset+j.
    """
    if not all(regular_two(parts) for parts in
               (operator_a, operator_b, output, input_)):
        raise ValueError("use explicit distinct-odd monomial arguments")
    if not sum(operator_b) == sum(output) == sum(input_):
        raise ValueError("operator_b, input, and output must have one weight")
    locality = sum(operator_a) + sum(operator_b)
    vertices = ([('U', k) for k in output] + [('B', k) for k in operator_b]
                + [('A', k) for k in operator_a] + [('Y', k) for k in input_])
    count = len(vertices)
    if count % 2:
        return 0, -locality
    edges = [[0] * count for _ in range(count)]
    for i, (first_group, first_mode) in enumerate(vertices):
        for j in range(i + 1, count):
            second_group, second_mode = vertices[j]
            groups = first_group, second_group
            value = 0
            if groups == ('U', 'B'):
                if creation_coefficient(second_mode, first_mode):
                    value = binomial_polynomial(second_mode)
            elif groups == ('U', 'A'):
                if creation_coefficient(second_mode, first_mode):
                    value = binomial_polynomial(second_mode) << first_mode
            elif groups == ('U', 'Y'):
                if first_mode == second_mode:
                    value = 1 << second_mode
            elif groups == ('B', 'A'):
                if annihilation_coefficient(first_mode, second_mode):
                    value = 1 << second_mode
            elif groups == ('B', 'Y'):
                if annihilation_coefficient(first_mode, second_mode):
                    value = binomial_polynomial(first_mode) << second_mode
            elif groups == ('A', 'Y'):
                if annihilation_coefficient(first_mode, second_mode):
                    value = binomial_polynomial(first_mode)
            edges[i][j] = value

    @lru_cache(None)
    def pfaffian(remaining: int) -> int:
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        candidates = rest
        value = 0
        while candidates:
            bit = candidates & -candidates
            other = bit.bit_length() - 1
            candidates ^= bit
            if edges[first][other]:
                value ^= polynomial_product(edges[first][other], pfaffian(rest ^ bit))
        return value

    return pfaffian((1 << count) - 1), -locality


def linear_piece_polynomial(operator_a, operator_b, output, input_) -> tuple[int, int]:
    """Extend the first operator argument linearly; coefficients reduce mod2."""
    terms = [monomial_piece_polynomial(a, operator_b, output, input_)
             for a, coefficient in operator_a if coefficient % 2]
    if not terms:
        return 0, 0
    offset = min(shift for _, shift in terms)
    result = 0
    for polynomial, shift in terms:
        result ^= polynomial << (shift - offset)
    return result, offset


def piece_matrices(labels, duals) -> dict[int, list[int]]:
    """Rows (x,y), columns (v,u); each matrix row is a bitset."""
    if len(labels) != len(duals):
        raise ValueError("one dual polynomial is required per label")
    dimension = len(labels)
    matrices = {}
    for x in range(dimension):
        for y_index, y in enumerate(labels):
            row = x * dimension + y_index
            for v_index, v in enumerate(labels):
                for u_index, u in enumerate(labels):
                    polynomial, offset = linear_piece_polynomial(duals[x], v, u, y)
                    while polynomial:
                        bit = polynomial & -polynomial
                        shift = offset + bit.bit_length() - 1
                        polynomial ^= bit
                        if shift not in matrices:
                            matrices[shift] = [0] * (dimension * dimension)
                        matrices[shift][row] ^= 1 << (v_index * dimension + u_index)
    return matrices


def reflected_duals(labels, modes: int):
    """Bit-complement monomial duals inside modes 1,3,...,2*modes-1.

    modes must be a power of two. Labels have a common weight. Their
    particle counts may differ; the reflected weights then differ too.
    """
    if modes < 1 or modes & (modes - 1):
        raise ValueError("modes must be a power of two")
    if not labels or len({sum(label) for label in labels}) != 1:
        raise ValueError("labels must have one weight")
    if not all(regular_two(label) and all(part < 2 * modes for part in label)
               for label in labels):
        raise ValueError("labels must be distinct odd parts in the specified box")
    return tuple(((tuple(sorted(2 * modes - part for part in label)), 1),)
                 for label in labels)
