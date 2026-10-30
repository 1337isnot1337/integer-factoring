"""Exact quadratic-code coefficient queries for the 23-term crossing family.

Research interface only: no all-size rank-profile improvement is asserted.
Quadratic Boolean functions are sets of square-free monomial masks.
"""

from functools import lru_cache
from itertools import combinations

from hamming_tensor_power import local_terms


def quadratic_monomials(variables):
    return (0,) + tuple(1 << j for j in range(variables)) + tuple(
        (1 << i) | (1 << j) for i, j in combinations(range(variables), 2)
    )


def label_polynomial(variables, width, label):
    monomials = quadratic_monomials(variables)
    assert 1 <= width <= len(monomials) and 0 <= label < 1 << width
    return frozenset(m for j, m in enumerate(monomials[:width]) if label >> j & 1)


def value(polynomial, point):
    return sum(point & monomial == monomial for monomial in polynomial) & 1


def truth_word(variables, polynomial):
    return sum(value(polynomial, point) << point for point in range(1 << variables))


def polynomial_product(left, right):
    answer = set()
    for a in left:
        for b in right:
            monomial = a | b
            if monomial in answer:
                answer.remove(monomial)
            else:
                answer.add(monomial)
    return frozenset(answer)


def gauss_sum(variables, polynomial):
    """Sum (-1)^f exactly, eliminating a quadratic variable pair each step."""
    current = set(polynomial)
    assert all(m < 1 << variables and m.bit_count() <= 2 for m in current)
    pairs = 0
    while True:
        edge = next((m for m in current if m.bit_count() == 2), None)
        if edge is None:
            if any(current):
                return 0
            return (-1 if 0 in current else 1) * (1 << (variables - pairs))
        bit_i = edge & -edge
        bit_j = edge ^ bit_i
        left, right, remainder = set(), set(), set()
        for monomial in current:
            if monomial == edge:
                continue
            if monomial & bit_i:
                left.add(monomial ^ bit_i)
            elif monomial & bit_j:
                right.add(monomial ^ bit_j)
            else:
                remainder.add(monomial)
        current = remainder.symmetric_difference(polynomial_product(left, right))
        pairs += 1


def joint_histogram(variables, polynomials):
    """Sixteen counts via Walsh inversion; no truth table is constructed."""
    assert len(polynomials) == 4
    character_sums = []
    for mask in range(16):
        polynomial = frozenset()
        for j, other in enumerate(polynomials):
            if mask >> (3 - j) & 1:
                polynomial = polynomial.symmetric_difference(other)
        character_sums.append(gauss_sum(variables, polynomial))
    counts = []
    for pattern in range(16):
        total = sum((-1 if (mask & pattern).bit_count() & 1 else 1) * scalar
                    for mask, scalar in enumerate(character_sums))
        assert total % 16 == 0
        counts.append(total // 16)
    assert min(counts) >= 0 and sum(counts) == 1 << variables
    return tuple(counts)


def z_oneplus_factors(polynomial):
    """Recognize the exact source local factor z^a(1+z)^b."""
    assert polynomial > 0
    degree_z = (polynomial & -polynomial).bit_length() - 1
    polynomial >>= degree_z
    degree_oneplus = 0
    while polynomial != 1:
        assert polynomial.bit_count() % 2 == 0
        quotient = 0
        while polynomial.bit_length() >= 2:
            shift = polynomial.bit_length() - 2
            quotient ^= 1 << shift
            polynomial ^= 3 << shift
        assert polynomial == 0
        polynomial = quotient
        degree_oneplus += 1
    return degree_z, degree_oneplus


@lru_cache(maxsize=1)
def factored_local_terms():
    answer = []
    for local, scalar in local_terms():
        a, b = z_oneplus_factors(scalar)
        c, d = z_oneplus_factors(local[0])
        answer.append((a + c, b + d, tuple(
            None if not entry else z_oneplus_factors(entry) for entry in local
        )))
    assert len(answer) == 23
    return tuple(answer)


def binomial_profile_from_histogram(histogram):
    """At most 23 pairs (a,b), representing sum z^a(1+z)^b over F_2."""
    assert len(histogram) == 16 and min(histogram) >= 0
    answer = set()
    for low_a, low_b, local in factored_local_terms():
        a, b = low_a, low_b
        for count, factors in zip(histogram, local):
            if count:
                if factors is None:
                    break
                a += count * factors[0]
                b += count * factors[1]
        else:
            if (a, b) in answer:
                answer.remove((a, b))
            else:
                answer.add((a, b))
    return tuple(sorted(answer))


def binomial_profile(variables, polynomials):
    return binomial_profile_from_histogram(joint_histogram(variables, polynomials))


def coefficient(profile, degree):
    assert degree >= 0
    return sum(degree >= a and ((degree - a) & ~b) == 0
               for a, b in profile) & 1


def expanded_profile(profile):
    """Small-instance verifier only; production coefficient queries use Lucas."""
    result = 0
    for a, b in profile:
        current = 1
        bit = 0
        while b:
            if b & 1:
                current ^= current << (1 << bit)
            b >>= 1
            bit += 1
        result ^= current << a
    return result


def alternating_rank(variables, polynomial):
    rows = [0] * variables
    for monomial in polynomial:
        if monomial.bit_count() == 2:
            bit = monomial & -monomial
            other = monomial ^ bit
            i, j = bit.bit_length() - 1, other.bit_length() - 1
            rows[i] ^= other
            rows[j] ^= bit
    return binary_rank(rows)


def binary_rank(vectors):
    basis = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in basis:
                vector ^= basis[pivot]
            else:
                basis[pivot] = vector
                break
    return len(basis)


def quadratic_annihilator_dimension(variables, polynomial):
    monomials = quadratic_monomials(variables)
    word = truth_word(variables, polynomial)
    images = [word & truth_word(variables, frozenset([m])) for m in monomials]
    return len(monomials) - binary_rank(images)
