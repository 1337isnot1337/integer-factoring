#!/usr/bin/env python3
"""Exact spherical rank formulas for the full Catalan path crossing.

The Jack calculation uses only the standard library.  NumPy is needed for
the independent F_101 eight-site matrix-rank control.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import permutations, product
import json
from math import comb, factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def partitions(n: int, largest: int | None = None, length: int = 4):
    if n == 0:
        yield ()
        return
    if length == 0:
        return
    if largest is None or largest > n:
        largest = n
    for first in range(largest, 0, -1):
        for rest in partitions(n - first, first, length - 1):
            yield (first,) + rest


def pad(partition: tuple[int, ...]) -> tuple[int, int, int, int]:
    return partition + (0,) * (4 - len(partition))


def dominates(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    a, b = pad(left), pad(right)
    return all(sum(a[:i]) >= sum(b[:i]) for i in range(1, 5))


def jack_eigenvalue(partition: tuple[int, ...]) -> int:
    """Eigenvalue for the alpha=2, four-variable Jack operator."""
    lam = pad(partition)
    return sum(lam[i] * (lam[i] - 1 + 3 - i) for i in range(4))


def jack_operator_column(partition: tuple[int, ...]) -> dict[tuple[int, ...], Fraction]:
    """Expand D_2 m_partition in the monomial-symmetric basis.

    D_2 = sum_i x_i^2 d_i^2 + sum_{i<j}
          (x_i^2 d_i - x_j^2 d_j)/(x_i-x_j).
    Pairing exponent-swapped monomials makes each quotient polynomial.
    """
    lam = pad(partition)
    totals: defaultdict[tuple[int, ...], Fraction] = defaultdict(Fraction)
    diagonal_key = tuple(x for x in sorted(lam, reverse=True) if x)
    totals[diagonal_key] += sum(a * (a - 1) for a in lam)

    for alpha in set(permutations(lam)):
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = alpha[i], alpha[j]
                if a == b:
                    key = tuple(x for x in sorted(alpha, reverse=True) if x)
                    totals[key] += a
                elif a > b:
                    gap = a - b
                    for k in range(gap + 1):
                        coefficient = a if k == 0 or k == gap else gap
                        beta = list(alpha)
                        beta[i], beta[j] = b + gap - k, b + k
                        key = tuple(x for x in sorted(beta, reverse=True) if x)
                        totals[key] += coefficient

    result = {}
    for key, total in totals.items():
        orbit_size = len(set(permutations(pad(key))))
        value = total / orbit_size
        if value:
            result[key] = value
    return result


@cache
def jack_p(partition: tuple[int, ...]) -> dict[tuple[int, ...], Fraction]:
    """Monic Jack P_partition^(2), obtained from the triangular eigenproblem."""
    size = sum(partition)
    lower = [p for p in partitions(size) if dominates(partition, p)]
    # Strict dominance strictly decreases the sum of squares.
    lower.sort(key=lambda p: (sum(x * x for x in p), p), reverse=True)
    assert lower[0] == partition

    coefficients: dict[tuple[int, ...], Fraction] = {partition: Fraction(1)}
    target_eigenvalue = jack_eigenvalue(partition)
    for column in lower:
        scalar = coefficients.get(column, Fraction(0))
        if not scalar:
            continue
        for row, entry in jack_operator_column(column).items():
            if row == column:
                continue
            assert dominates(column, row)
            gap = target_eigenvalue - jack_eigenvalue(row)
            assert gap != 0
            coefficients[row] = coefficients.get(row, Fraction(0)) + scalar * entry / gap
    return {key: value for key, value in coefficients.items() if value}


Poly = dict[tuple[int, int], Fraction]  # exponent of s=t^2, exponent of w


def poly_add(left: Poly, right: Poly, scale: Fraction = Fraction(1)) -> Poly:
    result = left.copy()
    for key, value in right.items():
        result[key] = result.get(key, Fraction(0)) + scale * value
        if not result[key]:
            del result[key]
    return result


def poly_mul(left: Poly, right: Poly) -> Poly:
    result: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            result[a + c, b + d] += x * y
    return {key: value for key, value in result.items() if value}


def poly_pow(value: Poly, exponent: int) -> Poly:
    result: Poly = {(0, 0): Fraction(1)}
    for _ in range(exponent):
        result = poly_mul(result, value)
    return result


def monomial_orbit_size(partition: tuple[int, ...]) -> int:
    return len(set(permutations(pad(partition))))


@cache
def jack_spherical_polynomial(partition: tuple[int, ...]) -> Poly:
    """P^(2)(1,1,x,y)/P^(2)(1,1,1,1), in s,w coordinates.

    The two nontrivial eigenvalues have sum 2(s+w) and product (s-w)^2.
    """
    s: Poly = {(1, 0): Fraction(1)}
    w: Poly = {(0, 1): Fraction(1)}
    eigen_sum = {(1, 0): Fraction(2), (0, 1): Fraction(2)}
    eigen_product = poly_pow(poly_add(s, w, Fraction(-1)), 2)
    power_sums = [{(0, 0): Fraction(2)}, eigen_sum]
    for _ in range(2, max(2, sum(partition) + 1)):
        power_sums.append(
            poly_add(poly_mul(eigen_sum, power_sums[-1]),
                     poly_mul(eigen_product, power_sums[-2]), Fraction(-1))
        )

    def xy_symmetric(a: int, b: int) -> Poly:
        if a < b:
            a, b = b, a
        return poly_mul(poly_pow(eigen_product, b), power_sums[a - b])

    value: Poly = {}
    normalization = Fraction(0)
    for exponent_partition, coefficient in jack_p(partition).items():
        orbit: defaultdict[tuple[int, int], int] = defaultdict(int)
        for alpha in set(permutations(pad(exponent_partition))):
            orbit[alpha[2], alpha[3]] += 1
        seen: set[tuple[int, int]] = set()
        evaluation: Poly = {}
        for (a, b), multiplicity in orbit.items():
            if (a, b) in seen:
                continue
            if a == b:
                term = poly_pow(eigen_product, a)
                scalar = multiplicity
                seen.add((a, b))
            else:
                assert multiplicity == orbit[b, a]
                term = xy_symmetric(a, b)
                scalar = multiplicity
                seen.update(((a, b), (b, a)))
            evaluation = poly_add(
                evaluation,
                {key: scalar * item for key, item in term.items()},
            )
        value = poly_add(value, evaluation, coefficient)
        normalization += coefficient * monomial_orbit_size(exponent_partition)

    return {key: coefficient / normalization for key, coefficient in value.items()
            if coefficient}


def hook_dimension(partition: tuple[int, ...]) -> int:
    hooks = 1
    for i, row in enumerate(partition):
        for j in range(row):
            below = sum(1 for k in range(i + 1, len(partition))
                        if partition[k] > j)
            hooks *= row - j + below
    return factorial(sum(partition)) // hooks


@cache
def spherical_blocks(semilength: int):
    blocks = []
    for nu in partitions(semilength):
        mu = tuple(2 * value for value in pad(nu))
        mu = tuple(value for value in mu if value)
        blocks.append((mu, hook_dimension(mu), jack_spherical_polynomial(nu)))
    if semilength >= 2:
        for nu in partitions(semilength - 2):
            mu = tuple(2 * value + 1 for value in pad(nu))
            blocks.append((mu, hook_dimension(mu),
                           poly_mul(jack_spherical_polynomial(nu),
                                    {(1, 0): Fraction(1), (0, 1): Fraction(-1)})))
    return blocks


def rank_pairs(semilength: int):
    blocks = spherical_blocks(semilength)
    result = []
    for a in range(semilength + 1):
        for b in range(semilength - a + 1):
            c = semilength - a - b
            row = sum(dim for _, dim, poly in blocks if poly.get((b, c), 0))
            col = sum(dim for _, dim, poly in blocks if poly.get((b, a), 0))
            result.append((a, b, c, row, col))
    return result


def decimal_cost(pairs: list[tuple[int, int, int, int, int]]) -> str:
    with localcontext() as context:
        context.prec = 40
        total = sum((Decimal(row) * Decimal(col)).sqrt()
                    for _, _, _, row, col in pairs)
        return f"{total:.12f}"


def dyck_paths(length: int):
    result = []

    def walk(prefix, height):
        if len(prefix) == length:
            if height == 0:
                result.append(tuple(prefix))
            return
        walk(prefix + [1], height + 1)
        if height:
            walk(prefix + [-1], height - 1)

    walk([], 0)
    return result


def recursive_vector(path: tuple[int, ...], prime: int) -> dict[int, int]:
    vector, height = {0: 1}, 0
    for position, step in enumerate(path):
        if step == -1:
            following: defaultdict[int, int] = defaultdict(int)
            for word, value in vector.items():
                following[word | (1 << position)] += height * value
                for j in range(position):
                    if not (word >> j) & 1:
                        following[word | (1 << j)] -= value
            vector = {word: value % prime for word, value in following.items()
                      if value % prime}
        height += step
    return vector


def path_norm(path: tuple[int, ...], prime: int) -> int:
    height, value = 0, 1
    for step in path:
        if step == -1:
            value = value * height * (height + 1) % prime
        height += step
    return value


def rank_mod(matrix, prime: int) -> int:
    import numpy as np

    matrix = np.asarray(matrix, dtype=np.int64).copy() % prime
    pivot_row = 0
    for column in range(matrix.shape[1]):
        hits = np.flatnonzero(matrix[pivot_row:, column])
        if not len(hits):
            continue
        hit = pivot_row + int(hits[0])
        matrix[[pivot_row, hit]] = matrix[[hit, pivot_row]]
        inverse = pow(int(matrix[pivot_row, column]), -1, prime)
        matrix[pivot_row, column:] = matrix[pivot_row, column:] * inverse % prime
        factors = matrix[pivot_row + 1:, column].copy()
        matrix[pivot_row + 1:, column:] = (
            matrix[pivot_row + 1:, column:]
            - factors[:, None] * matrix[pivot_row, column:]
        ) % prime
        pivot_row += 1
        if pivot_row == matrix.shape[0]:
            break
    return pivot_row


def direct_eight_site_ranks(prime: int = 101):
    """Independent word-basis coefficient matrices and finite-field ranks."""
    import numpy as np

    length, semilength = 8, 4
    paths = dyck_paths(length)
    words = [word for word in range(1 << length)
             if word.bit_count() == semilength]
    word_index = {word: i for i, word in enumerate(words)}
    frame = np.zeros((len(paths), len(words)), dtype=np.int64)
    for i, path in enumerate(paths):
        vector = recursive_vector(path, prime)
        for word, value in vector.items():
            if word in word_index:
                frame[i, word_index[word]] = value
    norms = [path_norm(path, prime) for path in paths]
    assert all(norms)

    groups: defaultdict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for x, y in product(words, repeat=2):
        groups[x & y, x ^ y].append((x, y))
    monomials = [(a, b, semilength - a - b)
                 for a in range(semilength + 1)
                 for b in range(semilength - a + 1)]
    slots = {monomial: i for i, monomial in enumerate(monomials)}
    q = len(paths)
    tensors = np.zeros((len(monomials), q, q, q, q), dtype=np.int64)

    for (_, difference), pairs in groups.items():
        s = difference.bit_count()
        assert s % 2 == 0
        encoded = np.stack([
            np.outer(frame[:, word_index[x]], frame[:, word_index[y]]).ravel()
            for x, y in pairs
        ], axis=1) % prime
        local = np.zeros((s // 2 + 1, len(pairs), len(pairs)), dtype=np.int64)
        for i, (x, _) in enumerate(pairs):
            for j, (_, u) in enumerate(pairs):
                h = (x ^ u).bit_count()
                assert h % 2 == 0
                local[h // 2, i, j] = 1
        for k in range(s // 2 + 1):
            monomial = (semilength - s // 2, s // 2 - k, k)
            block = (encoded @ local[k] % prime) @ encoded.T % prime
            slot = slots[monomial]
            tensors[slot] = (
                tensors[slot].reshape(q * q, q * q) + block
            ).reshape(q, q, q, q) % prime

    ranks = []
    for monomial, tensor in zip(monomials, tensors):
        row = rank_mod(tensor.reshape(q * q, q * q), prime)
        column = rank_mod(
            tensor.transpose(0, 3, 2, 1).reshape(q * q, q * q), prime
        )
        ranks.append((*monomial, row, column))
    return ranks


def run() -> dict:
    summaries = []
    full_pairs = {}
    for semilength in range(1, 13):
        pairs = rank_pairs(semilength)
        cat = comb(2 * semilength, semilength) // (semilength + 1)
        block_dimension = sum(dim for _, dim, _ in spherical_blocks(semilength))
        assert block_dimension == cat * cat
        summaries.append({
            "semilength": semilength,
            "catalan_labels": cat,
            "block_dimension": block_dimension,
            "coefficient_count": len(pairs),
            "natural_piece_rank_cost": decimal_cost(pairs),
        })
        if semilength == 4:
            full_pairs = {f"{a},{b},{c}": [row, col]
                          for a, b, c, row, col in pairs}

    direct = direct_eight_site_ranks()
    expected_direct = [
        (int(a), int(b), int(c), int(row), int(col))
        for a, b, c, row, col in rank_pairs(4)
    ]
    assert direct == expected_direct, (direct, expected_direct)
    assert decimal_cost(rank_pairs(4)) == "1608.377200168389"
    return {
        "scope": "full Catalan paths, characteristic-zero zonal rank formula; direct F101 control at eight sites",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "research_note_sha256": sha256(
            (ROOT / "research/catalan-spherical-ranks.md").read_bytes()
        ).hexdigest(),
        "jack_parameter": 2,
        "spherical_eigenvalues": ["1", "1", "(t+sqrt(w))^2", "(t-sqrt(w))^2"],
        "summaries": summaries,
        "semilength_four_rank_pairs": full_pairs,
        "f101_direct_eight_site_rank_pairs": [list(row) for row in direct],
        "f101_direct_rank_cost": decimal_cost([
            (a, b, c, row, col) for a, b, c, row, col in direct
        ]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
