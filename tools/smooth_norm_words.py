"""Convert supplied smooth represented norms directly to factor-base words.

This module does not search for smooth represented norms. Its optional
word replay uses the repository's small-discriminant composition controls,
which are not the asymptotically efficient form-arithmetic implementation.
"""

from dataclasses import dataclass
from math import gcd

from arithmetic_two_group import power
from binary_quadratic_forms import compose_forms, inverse_form, reduce_form
from bounded_smooth_norm import bezout, smooth_factors


@dataclass(frozen=True)
class PrimeTerm:
    prime: int
    exponent: int
    form: tuple


@dataclass(frozen=True)
class NormWord:
    normalized_vector: tuple
    input_norm: int
    primitive_norm: int
    middle: int
    trailing: int
    terms: tuple

    @property
    def multiplicity(self):
        return sum(abs(term.exponent) for term in self.terms)


def norm_to_word(form, x, y, smooth_bound, multiplicity_budget=None):
    """Return a direct word, or None when a supplied witness is inadmissible.

    Nonzero vectors are divided by their gcd first. Smoothness,
    coprimality and the multiplicity budget refer to the primitive norm.
    The stored prime orientation uses a middle coefficient in [0, p].
    No divisor split, square-root search or class-group search is used.
    """
    a, b, c = form
    discriminant = b * b - 4 * a * c
    if (a <= 0 or discriminant >= 0 or gcd(gcd(a, b), c) != 1
            or smooth_bound < 2):
        raise ValueError("expected a primitive positive form and smooth_bound >= 2")
    if multiplicity_budget is not None and multiplicity_budget < 0:
        raise ValueError("multiplicity budget must be nonnegative")
    divisor = gcd(x, y)
    if divisor == 0:
        return None
    input_norm = a * x * x + b * x * y + c * y * y
    x, y = x // divisor, y // divisor
    norm = a * x * x + b * x * y + c * y * y
    if gcd(norm, discriminant) != 1:
        return None
    factors = smooth_factors(norm, smooth_bound)
    if factors is None:
        return None
    if multiplicity_budget is not None and len(factors) > multiplicity_budget:
        return None

    v, negative_u = bezout(x, y)
    u = -negative_u
    assert x * v - y * u == 1
    middle = 2 * a * x * u + b * (x * v + y * u) + 2 * c * y * v
    trailing = a * u * u + b * u * v + c * v * v
    assert middle * middle - 4 * norm * trailing == discriminant

    terms = []
    for prime in sorted(set(factors)):
        oriented_root = middle % (2 * prime)
        canonical_root = min(oriented_root, 2 * prime - oriented_root)
        sign = 1 if oriented_root == canonical_root else -1
        assert (canonical_root * canonical_root - discriminant) % (4 * prime) == 0
        prime_form = reduce_form((
            prime,
            canonical_root,
            (canonical_root * canonical_root - discriminant) // (4 * prime),
        ))
        terms.append(PrimeTerm(prime, sign * factors.count(prime), prime_form))

    return NormWord((x, y), input_norm, norm, middle, trailing, tuple(terms))


def word_value(word, discriminant):
    """Replay a direct word using exact finite-control group arithmetic."""
    parity = discriminant % 2
    identity = reduce_form((1, parity, (parity - discriminant) // 4))
    result = identity
    for term in word.terms:
        value = term.form if term.exponent >= 0 else inverse_form(term.form)
        result = compose_forms(
            result, power(value, abs(term.exponent), identity, compose_forms)
        )
    return result


def exponent_vector(word, factor_base_primes):
    """Encode a word in a caller-supplied prime-indexed factor base."""
    positions = {prime: index for index, prime in enumerate(factor_base_primes)}
    if len(positions) != len(factor_base_primes):
        raise ValueError("factor-base prime coordinates must be distinct")
    result = [0] * len(factor_base_primes)
    for term in word.terms:
        if term.prime not in positions:
            raise ValueError("word uses a prime outside the factor base")
        result[positions[term.prime]] = term.exponent
    return tuple(result)


def dominant_relation(index, sample_word, target_word, word_budget, sample_budget):
    """Build a certified dominant row from two supplied exponent vectors."""
    if len(sample_word) != len(target_word) or not 0 <= index < len(sample_word):
        raise ValueError("invalid relation coordinates")
    if word_budget < 0 or sample_budget < 0:
        raise ValueError("budgets must be nonnegative")
    if sum(map(abs, sample_word)) > sample_budget:
        raise ValueError("sample word exceeds its budget")
    if sum(map(abs, target_word)) > word_budget:
        raise ValueError("target word exceeds its budget")
    dominant = word_budget + sample_budget + 1
    row = [left - right for left, right in zip(sample_word, target_word)]
    row[index] += dominant
    assert abs(row[index]) > sum(abs(value) for j, value in enumerate(row) if j != index)
    return tuple(row)
