#!/usr/bin/env python3
"""Exact finite checks for the rank-sensitive digit-box Fourier note.

Uses no network access or third-party packages. The determinant under
test is composition/trace based; the reference is Gaussian elimination.
The full large-input factoring evaluator is not implemented here.
"""

from __future__ import annotations

import argparse
from itertools import combinations, product
import json
from pathlib import Path
import random


def det_gauss(matrix, q):
    a = [list(row) for row in matrix]
    n = len(a)
    value = 1
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j][i] % q), None)
        if pivot is None:
            return 0
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            value = -value
        d = a[i][i] % q
        value = value * d % q
        inv = pow(d, -1, q)
        for j in range(i + 1, n):
            multiplier = a[j][i] * inv % q
            for k in range(i, n):
                a[j][k] = (a[j][k] - multiplier * a[i][k]) % q
    return value % q


def power_entry(entry, n, k, i, j, q):
    if k == 0:
        return int(i == j)
    if k == 1:
        return entry(i, j) % q
    left = k // 2
    return sum(
        power_entry(entry, n, left, i, a, q)
        * power_entry(entry, n, k - left, a, j, q)
        for a in range(n)
    ) % q


def det_composition(entry, n, q):
    if q <= n:
        raise ValueError("Trace-composition denominators need q > order")
    if n == 0:
        return 1
    total = 0
    for mask in range(1 << (n - 1)):
        parts = []
        last = 0
        for cut in range(1, n):
            if mask >> (cut - 1) & 1:
                parts.append(cut - last)
                last = cut
        parts.append(n - last)
        term = 1
        for k in parts:
            trace = sum(power_entry(entry, n, k, i, i, q) for i in range(n)) % q
            term = term * trace * pow(k, -1, q) % q
        for j in range(1, len(parts) + 1):
            term = term * pow(j, -1, q) % q
        total = (total + (-1 if len(parts) % 2 else 1) * term) % q
    return (-total if n % 2 else total) % q


def rank_minor(matrix, columns, q):
    rows = len(matrix)
    for r in range(min(rows, columns), 0, -1):
        for ii in combinations(range(rows), r):
            for jj in combinations(range(columns), r):
                d = det_composition(lambda a, b: matrix[ii[a]][jj[b]], r, q)
                if d:
                    return r, ii, jj
    return 0, (), ()


class RootAlgebra:
    def __init__(self, q, modulus):
        self.q = q
        self.modulus = modulus
        self.zero = (0,) * (q - 1)
        self.one = (1,) + (0,) * (q - 2)

    def add(self, a, b):
        return tuple((x + y) % self.modulus for x, y in zip(a, b))

    def scale(self, a, c):
        return tuple(x * c % self.modulus for x in a)

    def phase(self, exponent):
        r = exponent % self.q
        if r == self.q - 1:
            return (self.modulus - 1,) * (self.q - 1)
        return tuple(int(i == r) for i in range(self.q - 1))

    def mul(self, a, b):
        cyc = [0] * self.q
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                cyc[(i + j) % self.q] += x * y
        return tuple((cyc[i] - cyc[-1]) % self.modulus for i in range(self.q - 1))


def weight(x, weights, modulus):
    answer = 1
    for i, v in enumerate(x):
        answer = answer * weights[i][v] % modulus
    return answer


def affine_direct(c, b, domains, weights, alpha, ring):
    result = ring.zero
    for x in product(*domains):
        if any(sum(a * z for a, z in zip(row, x)) % ring.q != rhs for row, rhs in zip(c, b)):
            continue
        result = ring.add(result, ring.scale(ring.phase(sum(a * z for a, z in zip(alpha, x))), weight(x, weights, ring.modulus)))
    return result


def affine_nullity(c, b, domains, weights, alpha, ring):
    s = len(domains)
    r, ii, jj = rank_minor(c, s, ring.q)
    augmented = [list(row) + [rhs] for row, rhs in zip(c, b)]
    if rank_minor(augmented, s + 1, ring.q)[0] != r:
        return ring.zero
    free = [i for i in range(s) if i not in jj]
    pivot = [[c[i][j] for j in jj] for i in ii]
    denominator = det_composition(lambda a, z: pivot[a][z], r, ring.q)
    result = ring.zero
    for values in product(*(domains[i] for i in free)):
        x = [0] * s
        for i, value in zip(free, values):
            x[i] = value
        rhs = [(b[i] - sum(c[i][j] * x[j] for j in free)) % ring.q for i in ii]
        for k, j in enumerate(jj):
            numerator = det_composition(lambda a, z: rhs[a] if z == k else pivot[a][z], r, ring.q)
            x[j] = numerator * pow(denominator, -1, ring.q) % ring.q
        if any(x[i] not in domains[i] for i in range(s)):
            continue
        assert all(sum(a * z for a, z in zip(row, x)) % ring.q == v for row, v in zip(c, b))
        result = ring.add(result, ring.scale(ring.phase(sum(a * z for a, z in zip(alpha, x))), weight(x, weights, ring.modulus)))
    return result


def affine_rank(c, b, domains, weights, alpha, ring):
    s = len(domains)
    r, ii, _ = rank_minor(c, s, ring.q)
    if rank_minor([list(row) + [rhs] for row, rhs in zip(c, b)], s + 1, ring.q)[0] != r:
        return ring.zero
    result = ring.zero
    for v in product(range(ring.q), repeat=r):
        term = ring.phase(-sum(v[j] * b[ii[j]] for j in range(r)))
        for i, domain in enumerate(domains):
            frequency = alpha[i] + sum(c[ii[j]][i] * v[j] for j in range(r))
            one = ring.zero
            for digit in domain:
                one = ring.add(one, ring.scale(ring.phase(digit * frequency), weights[i][digit]))
            term = ring.mul(term, one)
        result = ring.add(result, term)
    return ring.scale(result, pow(pow(ring.q, r, ring.modulus), -1, ring.modulus))


def bilinear_direct(jmatrix, alpha, beta, domains, weights, ring):
    answer = ring.zero
    t = len(beta)
    for x in product(*domains):
        w = weight(x, weights, ring.modulus)
        for y in product(range(ring.q), repeat=t):
            exponent = sum(a * b for a, b in zip(alpha, x)) + sum(a * b for a, b in zip(beta, y))
            exponent += sum(x[i] * jmatrix[i][j] * y[j] for i in range(len(x)) for j in range(t))
            answer = ring.add(answer, ring.scale(ring.phase(exponent), w))
    return answer


def hankel(gamma, m):
    return [[gamma[i + j] for j in range(m)] for i in range(m)]


def run_checks():
    rng = random.Random(20261030)
    determinant_cases = 0
    for q in (3, 5, 7):
        for n in range(min(q, 5)):
            matrices = [()] if n == 0 else ([v[i * n:(i + 1) * n] for i in range(n)] for v in product(range(q), repeat=n * n)) if n <= 2 else ([[rng.randrange(q) for _ in range(n)] for _ in range(n)] for _ in range(32))
            for matrix in matrices:
                actual = det_composition(lambda i, j: matrix[i][j], n, q)
                assert actual == det_gauss(matrix, q), (q, matrix, actual)
                determinant_cases += 1
    try:
        det_composition(lambda i, j: int(i == j), 3, 3)
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid characteristic was silently accepted")

    affine_cases = 0
    bilinear_cases = 0
    prefix_formula_cases = 0
    for q, modulus in ((5, 2), (5, 3), (5, 4), (7, 2)):
        ring = RootAlgebra(q, modulus)
        assert ring.add(ring.one, sum_phase(ring, range(1, q))) == ring.zero
        for s, t in ((1, 1), (2, 1), (2, 2), (3, 2)):
            for case in range(12):
                domains = [tuple(range(rng.choice((1, 2, 4)))) for _ in range(s)]
                weights = [{v: rng.randrange(modulus) for v in domain} for domain in domains]
                jmatrix = [[0 if case == 0 else rng.randrange(q) for _ in range(t)] for _ in range(s)]
                alpha = [rng.randrange(q) for _ in range(s)]
                beta = [0 if case == 0 else rng.randrange(q) for _ in range(t)]
                c = [[jmatrix[i][j] for i in range(s)] for j in range(t)]
                b = [(-v) % q for v in beta]
                expected = affine_direct(c, b, domains, weights, alpha, ring)
                assert affine_nullity(c, b, domains, weights, alpha, ring) == expected
                assert affine_rank(c, b, domains, weights, alpha, ring) == expected
                affine_cases += 1
                if s <= 2:
                    assert bilinear_direct(jmatrix, alpha, beta, domains, weights, ring) == ring.scale(expected, pow(q, t, modulus))
                    bilinear_cases += 1
        for free_bits in range(3):
            for fixed in range(2 ** (2 - free_bits)):
                base = fixed << free_bits
                bits = [rng.randrange(modulus) for _ in range(free_bits)]
                for frequency in range(q):
                    direct = ring.zero
                    for digit in range(base, base + (1 << free_bits)):
                        w = 1
                        for bit in range(free_bits):
                            if digit >> bit & 1:
                                w = w * bits[bit] % modulus
                        direct = ring.add(direct, ring.scale(ring.phase(frequency * digit), w))
                    formula = ring.phase(frequency * base)
                    for bit, w in enumerate(bits):
                        formula = ring.mul(formula, ring.add(ring.one, ring.scale(ring.phase(frequency * (1 << bit)), w)))
                    assert direct == formula
                    prefix_formula_cases += 1

    hankel_cases = []
    completion_cases = 0
    carry_cases = []
    for q, m in ((3, 2), (5, 2), (5, 3), (7, 2)):
        total = q ** (2 * m - 1)
        full = sum(det_gauss(hankel(gamma, m), q) != 0 for gamma in product(range(q), repeat=2 * m - 1))
        assert full >= (q - m) * q ** (2 * m - 2)
        hankel_cases.append({"q": q, "m": m, "full_rank": full, "total": total, "proved_lower_bound": (q - m) * q ** (2 * m - 2)})
        for initial in product(range(q), repeat=m):
            if not any(initial):
                continue
            j = max(i for i, x in enumerate(initial) if x)
            gamma = list(initial) + [0] * (m - 1)
            if j < m - 1:
                gamma[m + j] = 1
            assert det_gauss(hankel(gamma, m), q)
            completion_cases += 1
        ring = RootAlgebra(q, 2 if q != 2 else 3)
        base = 2
        for value in range(1, base ** m):
            legal = [c for c in product(range(q), repeat=m) if sum(c[i] * base ** i for i in range(m)) == value]
            nonzero_initial = []
            for initial in product(range(q), repeat=m):
                fourier = sum_phase(ring, (sum(a * b for a, b in zip(initial, c)) for c in legal))
                if any(initial) and fourier != ring.zero:
                    nonzero_initial.append(initial)
            assert nonzero_initial
            chosen = nonzero_initial[0]
            completions = sum(det_gauss(hankel(chosen + tail, m), q) != 0 for tail in product(range(q), repeat=m - 1))
            assert completions >= (q - m) * q ** (m - 2)
            carry_cases.append({"q": q, "m": m, "N": value, "legal_lower_coefficients": len(legal), "chosen_nonzero_frequency": chosen, "full_rank_extensions": completions, "proved_lower_bound": (q - m) * q ** (m - 2)})

    return {"status": "passed", "scope": "Exact finite transcription checks; no new factoring exponent or Lean proof.", "determinant_cases": determinant_cases, "affine_cases": affine_cases, "bilinear_cases": bilinear_cases, "prefix_formula_cases": prefix_formula_cases, "hankel_cases": hankel_cases, "nonzero_initial_completion_cases": completion_cases, "carry_support_cases": carry_cases}


def sum_phase(ring, exponents):
    value = ring.zero
    for exponent in exponents:
        value = ring.add(value, ring.phase(exponent))
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
