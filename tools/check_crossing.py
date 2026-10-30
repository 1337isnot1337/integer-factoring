"""Exact finite certificates for the crossing-coefficient compression.

This checks algebraic identities and channel encodings. It is not an
implementation or a timing benchmark of a fast factoring algorithm.
Only the Python standard library is required.
"""

from __future__ import annotations

from itertools import product
from math import comb
import random
import unittest


def sign(p: int, q: int) -> int:
    """The sign in upstream Cross.lean, with bit zero first."""
    p1, p2, p3 = ((p >> j) & 1 for j in range(3))
    q1, q2, q3 = ((q >> j) & 1 for j in range(3))
    parity = (
        p1 * q1 + p1 * q2 + p1 * q3 + p2 * q2 + p2 * q3
        + p3 * q3 + q1 * p2 * p3 + p1 * q2 * p3 + p1 * p2 * q3
    ) % 2
    return 1 - 2 * parity


def psi(k: int, p: int, q: int) -> int:
    return sign(p, q) if p ^ q == k else 0


def crossing_r(x: int, y: int, u: int, v: int) -> int:
    return sum(psi(k, x, y) * psi(k, v, u) for k in range(8))


def crossing_c(x: int, y: int, u: int, v: int) -> int:
    return sum(psi(k, x, u) * psi(k, v, y) for k in range(8))


def residual_left(k: int, a: int, x: int, y: int) -> int:
    if x ^ y != k:
        return 0
    return int(x == a) - int(x == 0) * sign(a, a ^ k)


def residual_right(k: int, a: int, u: int, v: int) -> int:
    if v ^ u != k:
        return 0
    return 8 * int(a == v) - sign(a, a ^ k) * sign(v, u)


def residual(x: int, y: int, u: int, v: int) -> int:
    return sum(
        residual_left(k, a, x, y) * residual_right(k, a, u, v)
        for k in range(8) for a in range(1, 8)
    )


def original_coefficients(blocks: tuple[tuple[int, int, int, int], ...]) -> list[int]:
    """Compute prod(R+T C) by independent polynomial convolution."""
    coefficients = [1]
    for block in blocks:
        r, c = crossing_r(*block), crossing_c(*block)
        result = [0] * (len(coefficients) + 1)
        for i, value in enumerate(coefficients):
            result[i] += r * value
            result[i + 1] += c * value
        coefficients = result
    return coefficients


def expanded_coefficient(blocks: tuple[tuple[int, int, int, int], ...], i: int) -> int:
    """The new subset formula, using the 56-channel residual identity."""
    t = len(blocks)
    total = 0
    for mask in range(1 << t):
        j = mask.bit_count()
        if j > i:
            continue
        term = comb(t - j, i - j) * 4 ** (t - i) * (-3) ** (i - j)
        for s, block in enumerate(blocks):
            term *= residual(*block) if mask >> s & 1 else crossing_r(*block)
        total += term
    return total


def encoded_channel_sum(blocks: tuple[tuple[int, int, int, int], ...], i: int) -> int:
    """Actually sum the full 16^t * 7^i padded channel encoding.

    This intentionally slow verifier is restricted to tiny test instances.
    """
    t = len(blocks)
    total = 0
    for mask in range(1 << t):
        j = mask.bit_count()
        if j > i:
            continue
        weight = comb(t - j, i - j) * 4 ** (t - i) * (-3) ** (i - j)
        for classes in product(range(8), repeat=t):
            for labels in product(range(1, 8), repeat=i):
                if any(a != 1 for a in labels[j:]):
                    continue
                left, right, selected = weight, 1, 0
                for s, (x, y, u, v) in enumerate(blocks):
                    if mask >> s & 1:
                        a = labels[selected]
                        selected += 1
                        left *= residual_left(classes[s], a, x, y)
                        right *= residual_right(classes[s], a, u, v)
                    else:
                        left *= psi(classes[s], x, y)
                        right *= psi(classes[s], v, u)
                total += left * right
    return total


def rank_mod(matrix: list[list[int]], modulus: int) -> int:
    a = [[x % modulus for x in row] for row in matrix]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, modulus)
        a[rank] = [x * inv % modulus for x in a[rank]]
        for r in range(rank + 1, len(a)):
            factor = a[r][col]
            if factor:
                a[r] = [(x - factor * y) % modulus for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def exact_rank_formula(t: int, i: int) -> int:
    return sum(comb(t, j) * 8 ** (t - j) * 56 ** j for j in range(i + 1))


def ceil_log2(value: int) -> int:
    if value < 1:
        raise ValueError("value must be positive")
    return (value - 1).bit_length()


class CrossingCompressionTests(unittest.TestCase):
    def test_all_base_entries(self) -> None:
        for x, y, u, v in product(range(8), repeat=4):
            r, c = crossing_r(x, y, u, v), crossing_c(x, y, u, v)
            e = int(x == v and y == u)
            s = residual(x, y, u, v)
            self.assertEqual(r + c, 2 * e)
            self.assertEqual(s, 8 * e - r)
            self.assertEqual(s, 3 * r + 4 * c)

    def test_integer_channel_ranks(self) -> None:
        pairs = list(product(range(8), repeat=2))
        # Columns are ordered (v,u), so E is the identity matrix.
        for function, expected in [(crossing_r, 8), (residual, 56), (crossing_c, 64)]:
            matrix = [[function(x, y, u, v) for v, u in pairs] for x, y in pairs]
            self.assertEqual(rank_mod(matrix, 101), expected)

    def test_scaled_coefficients(self) -> None:
        rng = random.Random(5817)
        for t in range(7):
            for _ in range(32):
                blocks = []
                for _s in range(t):
                    x, y, u = (rng.randrange(8) for _ in range(3))
                    blocks.append((x, y, u, x ^ y ^ u))
                blocks = tuple(blocks)
                expected = [4 ** t * value for value in original_coefficients(blocks)]
                actual = [expanded_coefficient(blocks, i) for i in range(t + 1)]
                self.assertEqual(actual, expected)
                self.assertEqual(sum(actual), 8 ** t * int(all(x == v and y == u for x, y, u, v in blocks)))
                self.assertTrue(all(abs(value) <= 8 ** t for value in actual))

    def test_full_row_and_column_channel_encodings(self) -> None:
        rng = random.Random(7453)
        for t in range(3):
            for _ in range(8):
                blocks = []
                for _s in range(t):
                    x, y, u = (rng.randrange(8) for _ in range(3))
                    blocks.append((x, y, u, x ^ y ^ u))
                blocks = tuple(blocks)
                expected = original_coefficients(blocks)
                transposed = tuple((x, u, y, v) for x, y, u, v in blocks)
                for i in range(t + 1):
                    self.assertEqual(encoded_channel_sum(blocks, i), 4 ** t * expected[i])
                    self.assertEqual(encoded_channel_sum(transposed, t - i), 4 ** t * expected[i])

    def test_rank_width_envelopes(self) -> None:
        self.assertLess(7 ** 26, 2 ** 73)
        for t in range(65):
            for i in range(t + 1):
                r = exact_rank_formula(t, i)
                other = exact_rank_formula(t, t - i)
                self.assertLessEqual(r, 16 ** t * 7 ** i)
                self.assertLessEqual(r * other, 1792 ** t)
                widths = 8 * t + ceil_log2(7 ** i) + ceil_log2(7 ** (t - i))
                self.assertLessEqual(26 * widths, 281 * t + 52)


if __name__ == "__main__":
    unittest.main(verbosity=2)
