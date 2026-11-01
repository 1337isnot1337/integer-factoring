#!/usr/bin/env python3
"""Exact parity-minor and entropy certificate for flat cube tensors.

Uses only the Python standard library.  Finite polynomial examples replay
the indexing; the accompanying note proves the all-degree statement.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from random import Random


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / "evidence" / "cube-flat-parity-floor.json"


def antipode(word: int) -> int:
    return sum(1 << (v ^ 7) for v in range(8) if word >> v & 1)


def endpoint(word: int) -> int:
    return int(word == antipode(word))


def polynomial_product(left: set[tuple[int, int, int]],
                       right: set[tuple[int, int, int]]) -> set[tuple[int, int, int]]:
    result: set[tuple[int, int, int]] = set()
    for a in left:
        for b in right:
            monomial = tuple(a[i] + b[i] for i in range(3))
            result.symmetric_difference_update((monomial,))
    return result


def minor_layouts() -> list[dict]:
    result = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        vertices = [v for v in range(8) if not (v >> a & 1 or v >> b & 1)]
        assert len(vertices) == 2
        v, w = vertices
        for parity in (0, 1):
            rows = [0, (1 << v) | (1 << w)] if parity == 0 else [1 << v, 1 << w]
            columns = [antipode(word) for word in rows]
            entries = [[r | c for c in columns] for r in rows]
            assert [[endpoint(word) for word in row] for row in entries] == [[1, 0], [0, 1]]
            for axis in (a, b):
                assert all(all(not (v >> axis & 1) for v in range(8) if word >> v & 1)
                           for word in rows)
                assert all(all(v >> axis & 1 for v in range(8) if word >> v & 1)
                           for word in columns)
                assert all(word.bit_count() % 2 == parity for word in rows + columns)
            result.append({"axes": [a, b], "parity": parity, "row_words": rows,
                           "column_words": columns, "physical_entries": entries})
    assert len(result) == 6
    return result


def polynomial_replay(layouts: list[dict]) -> dict:
    rng = Random(20261101)
    cases = []
    checked = 0
    for degree in range(1, 7):
        monomials = [(i, j, degree - i - j) for i in range(degree + 1)
                     for j in range(degree - i + 1)]
        correction = monomials[0]
        for repetition in range(2):
            tensor = []
            for word in range(256):
                poly = set() if word.bit_count() % 2 else {
                    exponent for exponent in monomials if rng.randrange(2)}
                if len(poly) % 2 != endpoint(word):
                    poly.symmetric_difference_update((correction,))
                assert len(poly) % 2 == endpoint(word)
                assert word.bit_count() % 2 == 0 or not poly
                tensor.append(poly)
            minors = []
            for layout in layouts:
                ((u00, u01), (u10, u11)) = layout["physical_entries"]
                determinant = polynomial_product(tensor[u00], tensor[u11])
                determinant.symmetric_difference_update(polynomial_product(tensor[u01], tensor[u10]))
                assert determinant and len(determinant) % 2 == 1
                assert all(sum(exponent) == 2 * degree for exponent in determinant)
                valuations = [min(exponent[i] for exponent in determinant) for i in range(3)]
                a, b = layout["axes"]
                assert valuations[a] + valuations[b] <= 2 * degree
                assert sum(valuations) <= 2 * degree
                minors.append({"axes": [a, b], "parity": layout["parity"],
                               "monomials": len(determinant), "valuations": valuations})
                checked += 1
            cases.append({"degree": degree, "repetition": repetition, "minors": minors})
    return {"field": "F2", "tensors": len(cases), "endpoint_entries": 256 * len(cases),
            "shared_minor_determinants": checked, "cases": cases}


def integer_product(numerators: list[int], multiplicities: list[int]) -> int:
    return math.prod(pow(n, n * multiplicity) for n, multiplicity in zip(numerators, multiplicities))


def entropy_certificate(name: str, left: int, right: int) -> dict:
    assert left > right
    return {"name": name, "strict_integer_comparison": True,
            "left_bits": left.bit_length(), "right_bits": right.bit_length(),
            "positive_difference_bits": (left - right).bit_length(),
            "difference_sha256": hashlib.sha256((left - right).to_bytes(
                ((left - right).bit_length() + 7) // 8, "big")).hexdigest()}


def run() -> dict:
    f = Fraction
    same = [86, 4, 18, 2]
    same_mult = [1, 1, 2, 12]
    assert sum(n * m for n, m in zip(same, same_mult)) == 150
    a, b, c = f(4, 150), f(18, 150), f(2, 150)
    assert a >= c and b >= c and a + b + 14 * c == f(1, 3)

    A, B = [285, 23, 57, 6, 7], [315, 42, 8, 11, 7]
    A_mult, B_mult = [1, 1, 2, 6, 6], [1, 2, 6, 1, 6]
    assert sum(n * m for n, m in zip(A, A_mult)) == 500
    assert sum(n * m for n, m in zip(B, B_mult)) == 500
    aa = [f(n, 500) for n in A]
    bb = [f(n, 500) for n in B]
    k = 2 * (aa[1] - aa[3])
    ell = bb[3] - bb[4]
    mu = 2 * (aa[2] - aa[4]) - ell
    assert k == bb[1] - bb[2]
    assert min(k, ell, mu) >= 0
    constant = 16 * (aa[3] + aa[4]) + 8 * (bb[2] + bb[4])
    assert constant == f(328, 500)
    assert constant + 2 * k + mu + 2 * ell == 1
    assert [k, mu, ell] == [f(34, 500), f(96, 500), f(4, 500)]

    certificates = [
        entropy_certificate("same_parity", pow(150, 750),
                            pow(2, 1734) * pow(integer_product(same, same_mult), 5)),
        entropy_certificate("mixed_parity", pow(500, 1500),
                            pow(2, 3468) * pow(integer_product(A, A_mult), 2)
                            * integer_product(B, B_mult)),
    ]
    assert 3 * 289**2 > 500**2
    layouts = minor_layouts()
    return {
        "scope": "All even-occupation homogeneous flat binary cube tensors with antipodal equality endpoint; the common-degree pressure only.",
        "proved_floor": "rho > 289/500 > 1/sqrt(3)",
        "same_parity_weights": {"denominator": 150, "numerators": same, "multiplicities": same_mult},
        "mixed_parity_weights": {"denominator": 500, "A_numerators": A, "A_multiplicities": A_mult,
                                 "B_numerators": B, "B_multiplicities": B_mult},
        "mixed_nonnegative_combination": {"constant": str(constant), "pair_even": str(k),
                                          "pair_odd_repeated": str(mu), "pair_odd_mixed": str(ell)},
        "entropy_certificates": certificates,
        "shared_minor_layouts": layouts,
        "polynomial_replay": polynomial_replay(layouts),
        "floating_point_used_for_certificate": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-evidence", action="store_true")
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    result = run()
    if args.write_evidence:
        args.evidence.parent.mkdir(parents=True, exist_ok=True)
        args.evidence.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif args.evidence.exists():
        assert json.loads(args.evidence.read_text()) == result, "Stored evidence differs from exact replay"
    print("PASS: exact flat parity floor >289/500>1/sqrt(3); 72 shared polynomial minors; 3072 endpoints")


if __name__ == "__main__":
    main()
