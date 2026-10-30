#!/usr/bin/env python3
"""Independent Boolean/global ranks and exact intervals for the reciprocal pencil."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from check_hamming_relabeling import kron, xor
from check_modular_crossing import rank_two
from hamming_relabeling import normalized, seed_data
from sl2_boolean_rank import (
    boolean_rank,
    reciprocal_character,
    reciprocal_optimum,
    reciprocal_rank,
    reciprocal_tiltings,
    subset_weights,
    tensor_tiltings,
    tilting_weights,
)

ROOT = Path(__file__).resolve().parents[1]


def boolean_checks():
    cases = 0
    for m in range(10):
        size = 1 << m
        elementary = [[0] * size for _ in range(m + 1)]
        for row in range(size):
            free = (size - 1) ^ row
            subset = free
            while True:
                elementary[subset.bit_count()][row] |= 1 << (row | subset)
                if not subset:
                    break
                subset = (subset - 1) & free
        rows = [0] * size
        assert boolean_rank(m, 0) == 0
        cases += 1
        for i in range(1, 1 << (m + 1)):
            bit = (i & -i).bit_length() - 1
            rows = xor(rows, elementary[bit])
            code = i ^ (i >> 1)
            assert boolean_rank(m, code) == rank_two(rows)
            cases += 1
    # Weyl factors Sym^2(V) and determinant do not form a rank decomposition.
    assert boolean_rank(2, 2) == 2
    assert rank_two([2, 0, 0]) + rank_two([0]) == 1
    return cases


def pencil_checks():
    _, _, pieces, _ = seed_data()
    r, c, x = pieces[0], pieces[1], pieces[5]
    a = pieces[2] ^ pieces[3] ^ pieces[4]
    local = [normalized(r ^ x), normalized(a ^ x), normalized(c ^ x)]
    coefficients = [[1]]
    cases = 0
    for t in (1, 2):
        following = [[0] * (16**t) for _ in range(2 * t + 1)]
        for i, matrix in enumerate(coefficients):
            for j, seed in enumerate(local):
                following[i + j] = xor(following[i + j], kron(matrix, seed))
        coefficients = following
        for mask in range(1 << (2 * t + 1)):
            rows = [0] * (16**t)
            for j, matrix in enumerate(coefficients):
                if mask >> j & 1:
                    rows = xor(rows, matrix)
            assert reciprocal_rank(t, mask) == rank_two(rows)
            cases += 1
    return cases


def series_bounds(value, terms=48):
    z = (value - 1) / (value + 1)
    lower = 2 * sum(z ** (2 * j + 1) / Fraction(2 * j + 1) for j in range(terms))
    error = 2 * z ** (2 * terms + 1) / (Fraction(2 * terms + 1) * (1 - z * z))
    return lower, lower + error


def log_bounds(value):
    value = Fraction(value)
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    while value < 2**exponent:
        exponent -= 1
    while value >= 2 ** (exponent + 1):
        exponent += 1
    l2, u2 = series_bounds(Fraction(2))
    lo, hi = series_bounds(value / Fraction(2**exponent))
    assert exponent >= 0
    return lo + exponent * l2, hi + exponent * u2


def decimal_interval(lo, hi):
    scale = 10**12
    a = lo.numerator * scale // lo.denominator
    b = -((-hi.numerator * scale) // hi.denominator)
    return {
        "lower": f"{a//scale}.{a%scale:012d}",
        "upper": f"{b//scale}.{b%scale:012d}",
    }


def saddle_bounds():
    lo, hi = Fraction(1), Fraction(2)

    def equation(x):
        return 3 * x**4 - 8 * x * x - 2 * x - 9

    for _ in range(80):
        mid = (lo + hi) / 2
        if equation(mid) < 0:
            lo = mid
        else:
            hi = mid
    assert equation(lo) < 0 < equation(hi)
    lower = 1 + 8 / hi + 1 / hi**2 + 3 * lo + 3 / hi**3
    upper = 1 + 8 / lo + 1 / lo**2 + 3 * hi + 3 / lo**3
    ll, _ = log_bounds(lower)
    _, ul = log_bounds(upper)
    l2, u2 = log_bounds(2)
    rlo, rhi = ll / (4 * u2), ul / (4 * l2)
    assert 2 * rlo * rlo > 1 and rhi < 1
    return {
        "minimizer": decimal_interval(lo, hi),
        "rank_growth_base": decimal_interval(lower, upper),
        "optimal_leading_rho": decimal_interval(rlo, rhi),
    }


def character_checks():
    cases = 0
    for n in range(129):
        assert (1 << len(tilting_weights(n))) <= (n + 1) ** 2
        assert set(subset_weights(n)) == set(range(n + 1))
        cases += 1
    for m in range(129):
        assert (
            sum(mu * (1 << len(tilting_weights(n))) for n, mu in tensor_tiltings(m))
            == 1 << m
        )
        cases += 1
    for t in range(1, 65):
        char = reciprocal_character(t)
        decomposition = reciprocal_tiltings(t)
        assert sum(char.values()) == 16**t
        assert all(mu <= char[n] for n, mu in decomposition)
        multiplicity = sum(mu for n, mu in decomposition if n >= t)
        assert multiplicity * (2 * t + 1) ** 2 >= char[t]
        assert char[t] * (t + 1) >= 96 ** (t // 2)
        cases += 1
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--full", action="store_true", help="extend exact span optimization through t=8"
    )
    args = parser.parse_args()
    boolean = boolean_checks()
    pencils = pencil_checks()
    characters = character_checks()
    optima = [reciprocal_optimum(t) for t in range(1, 9 if args.full else 5)]
    assert [case["exact_full_span_cutoff"] for case in optima[:4]] == [
        49,
        6241,
        524142,
        69872881,
    ]
    # Cauchy bound (4) evaluated at the rational test point x=2.
    for case in optima:
        t = case["word_length"]
        assert case["exact_full_span_cutoff"] <= 4 * (2 * t + 1) ** 4 * Fraction(
            93, 8
        ) ** (2 * t)
        assert (
            case["exact_full_span_cutoff"] * (2 * t + 1) ** 4
            >= reciprocal_character(t)[t] ** 2
        )
    paths = [Path(__file__), ROOT / "tools/sl2_boolean_rank.py"]
    result = {
        "status": "exact rank engine and solved reciprocal-pencil rate; no factoring improvement",
        "full": args.full,
        "exhaustive_direct_boolean_rank_cases": boolean,
        "independent_global_pencil_rank_cases": pencils,
        "character_and_block_checks": characters,
        "saddle_intervals": saddle_bounds(),
        "exact_finite_span_optima": optima,
        "primary_sources": [
            {
                "title": "Schur-Weyl duality in positive characteristic",
                "url": "https://arxiv.org/pdf/math/0610591v3",
                "used": "Lemma 13",
                "sha256": "187a84c35624b1713db2f1b4e581ac4938a5bb0e12f2d99f26f4da834014184a",
            },
            {
                "title": "SL2 tilting modules in the mixed case",
                "url": "https://arxiv.org/pdf/2105.07724v2",
                "used": "Proposition 3.4 and Corollary 3.6, algebraic characteristic two",
                "sha256": "bc649dfaa1e76ec4682ea6882578517d48dd1d307336c855fec8337cf3253623",
            },
        ],
        "source_sha256": {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in result.items() if k != "exact_finite_span_optima"})
    )


if __name__ == "__main__":
    main()
