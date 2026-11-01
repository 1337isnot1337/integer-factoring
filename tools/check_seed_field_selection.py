#!/usr/bin/env python3
"""Exact finite checks for input-wise field selection and rank specialization.

This is a certificate-mechanism check, not a better crossing seed or a
performance test of the asymptotic factoring algorithm.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path
import random


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def polynomial_remainder(a, b, p):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    while len(a) >= len(b):
        c = a[-1] * pow(b[-1], -1, p) % p
        shift = len(a) - len(b)
        for j, coefficient in enumerate(b):
            a[shift + j] = (a[shift + j] - c * coefficient) % p
        while a and not a[-1]:
            a.pop()
    return a


@dataclass(frozen=True)
class Field:
    p: int
    f: tuple[int, ...]

    def __post_init__(self):
        if not prime(self.p) or len(self.f) < 2 or self.f[-1] != 1:
            raise ValueError("invalid field presentation")
        if any(c < 0 or c >= self.p for c in self.f):
            raise ValueError("noncanonical polynomial")
        for d in range(1, self.e // 2 + 1):
            for coefficients in product(range(self.p), repeat=d):
                if not polynomial_remainder(self.f, coefficients + (1,), self.p):
                    raise ValueError("reducible field polynomial")

    @property
    def e(self):
        return len(self.f) - 1

    @property
    def q(self):
        return self.p**self.e

    def digits(self, a):
        result = []
        for _ in range(self.e):
            a, digit = divmod(a, self.p)
            result.append(digit)
        return result

    def encode(self, coefficients):
        return sum(c * self.p**i for i, c in enumerate(coefficients))

    def add(self, a, b):
        return self.encode([(x + y) % self.p for x, y in zip(self.digits(a), self.digits(b))])

    def neg(self, a):
        return self.encode([(-x) % self.p for x in self.digits(a)])

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        coefficients = [0] * (2 * self.e - 1)
        for i, x in enumerate(self.digits(a)):
            for j, y in enumerate(self.digits(b)):
                coefficients[i + j] = (coefficients[i + j] + x * y) % self.p
        return self.encode(polynomial_remainder(coefficients, self.f, self.p))

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        return next(b for b in range(1, self.q) if self.mul(a, b) == 1)


def rank(matrix, field):
    rows = [row[:] for row in matrix]
    answer = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(answer, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[answer], rows[pivot] = rows[pivot], rows[answer]
        inverse = field.inv(rows[answer][column])
        rows[answer] = [field.mul(x, inverse) for x in rows[answer]]
        for i in range(answer + 1, len(rows)):
            multiplier = rows[i][column]
            if multiplier:
                rows[i] = [field.sub(x, field.mul(multiplier, y)) for x, y in zip(rows[i], rows[answer])]
        answer += 1
    return answer


def rational_rank(matrix):
    rows = [[Fraction(x) for x in row] for row in matrix]
    answer = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(answer, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[answer], rows[pivot] = rows[pivot], rows[answer]
        inverse = 1 / rows[answer][column]
        rows[answer] = [x * inverse for x in rows[answer]]
        for i in range(answer + 1, len(rows)):
            multiplier = rows[i][column]
            rows[i] = [x - multiplier * y for x, y in zip(rows[i], rows[answer])]
        answer += 1
    return answer


def crossing(a):
    d = 1 << a
    return [int(x == v and y == u) for x, y, u, v in product(range(d), repeat=4)]


def cut_matrices(values, a):
    d = 1 << a
    row = [[0] * d**2 for _ in range(d**2)]
    col = [[0] * d**2 for _ in range(d**2)]
    for value, (x, y, u, v) in zip(values, product(range(d), repeat=4)):
        row[x * d + y][u * d + v] = value
        col[x * d + u][y * d + v] = value
    return row, col


def ceil_log2(x):
    return (x - 1).bit_length()


def validate_seed(a, pieces, field):
    if not isinstance(a, int) or a < 1 or not pieces:
        raise ValueError("invalid seed dimensions")
    expected = crossing(a)
    total = [0] * len(expected)
    ranks = []
    for piece in pieces:
        if len(piece) != len(expected) or any(not isinstance(x, int) or x < 0 or x >= field.q for x in piece):
            raise ValueError("malformed tensor")
        total = [field.add(x, y) for x, y in zip(total, piece)]
        if any(piece):
            r, c = (rank(matrix, field) for matrix in cut_matrices(piece, a))
            ranks.append((r, c))
    if total != expected or not ranks:
        raise ValueError("not a scale-one crossing identity")
    w = max(ceil_log2(r) + ceil_log2(c) for r, c in ranks)
    quality = Fraction(w + 2 * ceil_log2(len(ranks)), 4 * a)
    assert quality >= Fraction(1, 2)
    return quality, ranks


def cap(n):
    return max(1, n.bit_length().bit_length() - 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fields = [Field(2, (0, 1)), Field(3, (0, 1)), Field(2, (1, 1, 1)),
              Field(5, (0, 1)), Field(2, (1, 1, 0, 1)), Field(3, (1, 0, 1))]
    axioms = 0
    for f in fields:
        for x in range(f.q):
            assert f.add(x, f.neg(x)) == 0 and f.mul(x, 1) == x
            if x:
                assert f.mul(x, f.inv(x)) == 1
            for y in range(f.q):
                assert f.add(x, y) == f.add(y, x) and f.mul(x, y) == f.mul(y, x)
                for z in range(f.q):
                    assert f.add(f.add(x, y), z) == f.add(x, f.add(y, z))
                    assert f.mul(f.mul(x, y), z) == f.mul(x, f.mul(y, z))
                    assert f.mul(x, f.add(y, z)) == f.add(f.mul(x, y), f.mul(x, z))
                    axioms += 1

    rejected = 0
    for p, polynomial in [(4, (0, 1)), (2, (0, 0, 1)), (3, (2, 0, 1)),
                          (2, (1, 0, 1)), (2, (1, 2, 1)), (2, (1,))]:
        try:
            Field(p, polynomial)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("bad field accepted")

    rng = random.Random(20261030)
    extensions = projections = 0
    for base, extension in [(fields[0], fields[2]), (fields[0], fields[4]), (fields[1], fields[5])]:
        for rows in range(1, 6):
            for cols in range(1, 6):
                for _ in range(4):
                    a = [[rng.randrange(base.q) for _ in range(cols)] for _ in range(rows)]
                    assert rank(a, base) == rank(a, extension)
                    extensions += 1
                    a = [[rng.randrange(extension.q) for _ in range(cols)] for _ in range(rows)]
                    # Constant coefficient is base-linear and fixes 1 even
                    # when the unnormalized field trace annihilates 1.
                    projected = [[x % base.p for x in row] for row in a]
                    assert rank(projected, base) <= extension.e * rank(a, extension)
                    projections += 1
    f4 = fields[2]
    rank_one = [[1, 2], [2, f4.mul(2, 2)]]
    assert rank(rank_one, f4) == 1
    assert rank([[x % 2 for x in row] for row in rank_one], fields[0]) == 2
    projections += 1

    specialization = rank_drops = 0
    prime_fields = [Field(p, (0, 1)) for p in (2, 3, 5, 7, 11)]
    for _ in range(100):
        inner = rng.randrange(1, 5)
        u = [[rng.randrange(-5, 6) for _ in range(inner)] for _ in range(4)]
        v = [[rng.randrange(-5, 6) for _ in range(4)] for _ in range(inner)]
        a = [[sum(u[i][k] * v[k][j] for k in range(inner)) for j in range(4)] for i in range(4)]
        r = rational_rank(a)
        preserved = False
        for f in prime_fields:
            reduced = rank([[x % f.p for x in row] for row in a], f)
            assert reduced <= r
            rank_drops += reduced < r
            preserved |= reduced == r
            specialization += 1
        assert preserved
    assert rank_drops > 0  # Inverse-minor localization has real content.

    events = []
    seed_checks = 0
    for i, field in enumerate(fields):
        a = 1 if i < 4 else 2
        e = crossing(a)
        quality, ranks = validate_seed(a, [e], field)
        assert quality == 1 and ranks == [(1 << (2 * a),) * 2]
        events.append((3 + 17 * i, a, 1, field, quality))
        seed_checks += 1
        random_piece = [rng.randrange(field.q) for _ in e]
        remainder = [field.sub(x, y) for x, y in zip(e, random_piece)]
        validate_seed(a, [random_piece, remainder], field)
        seed_checks += 1
        bad = e[:]
        bad[0] = field.add(bad[0], 1)
        for pieces in ([bad], [e[:-1]], [[field.q] + e[1:]]):
            try:
                validate_seed(a, pieces, field)
            except ValueError:
                rejected += 1
            else:
                raise AssertionError("bad seed accepted")

    selections = field_changes = 0
    last = None
    fallback = (0, 1, 1, fields[0], Fraction(1))
    for n in range(2, 4097):
        b = cap(n)
        eligible = [fallback] + [event for event in events if event[0] <= n and
                    max(event[1], event[2], ceil_log2(event[3].q)) <= b]
        chosen = min(eligible, key=lambda event: (event[4], -event[3].q, -event[1], event[0]))
        assert chosen[4] == min(event[4] for event in eligible)
        identity = (chosen[3].p, chosen[3].f)
        if last is not None:
            field_changes += identity != last
        last = identity
        selections += 1
    assert field_changes >= 3

    caps = 0
    ns = list(range(2, 65537))
    for exponent in range(2, 4097):
        ns.extend(((1 << exponent) - 1, 1 << exponent, (1 << exponent) + 1))
    for n in ns:
        length = n.bit_length()
        b = cap(n)
        assert (1 << b) <= length
        assert b * b * (1 << (4 * b)) <= length**6
        caps += 1
    entropy_envelopes = 0
    for m in range(1, 513):
        n = 1 << (20 * m)
        lam = 1 << (9 * m)
        b = cap(n)
        assert b * lam.bit_length() <= lam
        assert 4 * b + 2 <= lam
        # Exact representative expanded widths n^(9/20).
        assert lam**20 == n**9
        entropy_envelopes += 1

    result = {
        "status": "passed",
        "scope": "Finite field arithmetic, exact seed validation, selection and resource caps; no improved factoring exponent.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "field_orders": [field.q for field in fields],
        "field_axiom_triples": axioms,
        "invalid_fields_and_seeds_rejected": rejected,
        "extension_rank_equalities": extensions,
        "base_field_projection_rank_bounds": projections,
        "integer_matrix_reductions": specialization,
        "bad_prime_rank_drops_observed": rank_drops,
        "valid_crossing_seed_checks": seed_checks,
        "delayed_bounded_selections": selections,
        "selected_field_changes": field_changes,
        "integer_cap_checks": caps,
        "expanded_width_entropy_envelopes": entropy_envelopes,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
