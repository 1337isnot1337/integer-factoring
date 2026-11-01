#!/usr/bin/env python3
"""Exact small replay of the all-word Weyl algebra theorem; no dumps needed."""

import argparse
import hashlib
import json
from pathlib import Path
import random
from time import monotonic


def labels(h, variables):
    if variables == 0:
        return [()] if h == 0 else []
    w = 1 << (variables-1)
    return [low+(a,) for a in range(h//w+1)
            for low in labels(h-a*w, variables-1)]


def multiply(A, B):
    answer = []
    for row in A:
        value = 0
        while row:
            bit = row & -row
            value ^= B[bit.bit_length()-1]
            row ^= bit
        answer.append(value)
    return answer


def add(A, B):
    return [a ^ b for a, b in zip(A, B)]


def rank(A):
    pivots = {}
    for row in A:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def other_cut(A, d):
    answer = [0] * (d*d)
    for a in range(d):
        for y in range(d):
            row = A[a*d+y]
            while row:
                bit = row & -row
                beta, u = divmod(bit.bit_length()-1, d)
                answer[a*d+u] |= 1 << (beta*d+y)
                row ^= bit
    return answer


def normal_right(coefficients, letter):
    p, q, r, s = coefficients
    if letter == "N":
        return r << 1, p ^ (s << 1), 0, r
    return 0, 0, p, q


def algebra_dimension(generators, size):
    pivots, basis = {}, []

    def insert(A):
        value = sum(row << (i*size) for i, row in enumerate(A))
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                mask = (1 << size)-1
                basis.append([(value >> (i*size)) & mask for i in range(size)])
                return
            value ^= pivots[pivot]

    insert([1 << i for i in range(size)])
    position = 0
    while position < len(basis):
        for generator in generators:
            insert(multiply(basis[position], generator))
        position += 1
    return len(basis)


def calculate(h):
    basis = labels(h, h.bit_length())
    ids = {x: i for i, x in enumerate(basis)}
    d, size = len(basis), len(basis)**2
    A = [0] * size
    for a, alpha in enumerate(basis):
        for y, yy in enumerate(basis):
            for beta, bb in enumerate(basis):
                u = tuple(b+z-x for x, z, b in zip(alpha, yy, bb))
                if u in ids and not any(x & z for x, z in zip(alpha, u)):
                    A[a*d+y] |= 1 << (beta*d+ids[u])
    B = other_cut(A, d)
    I = [1 << i for i in range(size)]
    N = add(A, I)
    zero = [0] * size
    assert multiply(A, A) == I
    assert multiply(N, N) == multiply(B, B) == zero

    sums = {}
    for alpha in basis:
        for yy in basis:
            if not any(a & y for a, y in zip(alpha, yy)):
                sigma = tuple(a+y for a, y in zip(alpha, yy))
                if sigma not in sums:
                    sums[sigma] = len(sums)
    k = len(sums)
    U, V = [0] * size, [0] * k
    for a, alpha in enumerate(basis):
        for y, yy in enumerate(basis):
            sigma = tuple(x+z for x, z in zip(alpha, yy))
            if sigma in sums:
                i = sums[sigma]
                V[i] |= 1 << (a*d+y)
                if not any(x & z for x, z in zip(alpha, yy)):
                    U[a*d+y] |= 1 << i
    assert multiply(U, V) == B
    assert not any(multiply(V, U))
    assert rank(U) == rank(V) == rank(B) == k
    NU, VN = multiply(N, U), multiply(V, N)
    assert all(U[a*d+a] == NU[a*d+a] == 0 for a in range(d))
    X = [a | (b << k) for a, b in zip(NU, U)]
    Y = V + VN
    C = add(multiply(N, B), multiply(B, N))
    D = multiply(V, NU)
    assert multiply(X, Y) == C
    assert multiply(Y, X) == D + [row << k for row in D]
    assert multiply(C, N) == multiply(N, C)
    assert multiply(C, B) == multiply(B, C)

    powers = [I]
    for _ in range(5):
        powers.append(multiply(powers[-1], C))
    NB = multiply(N, B)

    def evaluate(coefficients):
        result = zero[:]
        for polynomial, generator in zip(coefficients, (I, N, B, NB)):
            value = zero[:]
            while polynomial:
                bit = polynomial & -polynomial
                value = add(value, powers[bit.bit_length()-1])
                polynomial ^= bit
            result = add(result, multiply(value, generator))
        return result

    words = [("", I, (1, 0, 0, 0))]
    ideal_words = []
    for _ in range(7):
        following = []
        for word, actual, coefficients in words:
            for letter, generator in (("N", N), ("B", B)):
                next_word = word+letter
                actual_next = multiply(actual, generator)
                normal = normal_right(coefficients, letter)
                assert actual_next == evaluate(normal)
                following.append((next_word, actual_next, normal))
                if "B" in next_word:
                    ideal_words.append(actual_next)
        words = following

    top = ids[tuple(int(j == len(basis[0])-1) for j in range(len(basis[0])))]
    rng = random.Random(1337+h)
    affine_checks = 0
    for _ in range(128):
        K = zero[:]
        for _ in range(6):
            K = add(K, rng.choice(ideal_words))
        assert rank(K) <= 2*k
        assert all(K[a*d+a] == 0 for a in range(d))
        for b in (0, 1):
            T = add(add(I, N if b else zero), K)
            row_rank, col_rank = rank(T), rank(other_cut(T, d))
            assert row_rank >= size-2*k
            expected = sum((b ^ int(beta == top)) << (beta*d+beta)
                           for beta in range(d))
            assert T[top*d+top] == expected
            if b:
                assert col_rank >= d-1
            else:
                perturbation_rank = rank(K)
                assert all(perturbation_rank >= g*(d-col_rank*g)
                           for g in range(1, d+1))
            affine_checks += 1

    return {"energy": h, "labels": d, "core_dimension": k,
            "rank_central_commutator": rank(C), "rank_small_core": rank(D),
            "word_algebra_dimension": algebra_dimension((N, B), size),
            "normal_form_word_checks": 254,
            "affine_coset_and_quotient_budget_checks": affine_checks,
            "square_zero_and_involution": True,
            "exact_core_factorizations": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    result = {"source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "checked_shells": [calculate(h) for h in (2, 4, 8)],
              "scope": "Finite exact replay of the written all-energy constant-generator word-algebra theorem. Parameterized coefficient algebras are outside that theorem.",
              "elapsed_seconds": round(monotonic()-started, 3)}
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
