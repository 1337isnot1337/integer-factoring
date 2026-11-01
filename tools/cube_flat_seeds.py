#!/usr/bin/env python3
"""Deterministic degree-six binary cube seeds, with exact polynomial checks.

Public constructors ``rank_two_seed()`` and ``odd_edge_face_seed()`` return
256 dictionaries mapping exponent triples to nonzero F2 coefficients.
``validate(tensor)`` checks every degree and the literal antipodal endpoint.
No optional packages, finite-field specialization, or temporary files are
needed.  Vertex v has coordinate a equal to (v >> a) & 1.

The rank-two seed is an even Gaussian quotient with common physical
profile (0,2,2,4,4,4,6,6,6,6,8,8,8,10,10,12).  The odd edge-face seed has
64 odd coefficients and rank-one coordinate faces.  Its known generic
profile is (0,2,2,2,2,2,2,6,6,10,10,10,10,10,10,12).  The present module
checks polynomial identities, not those generic Smith certificates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path


Exponent = tuple[int, int, int]
Polynomial = dict[Exponent, int]
Tensor = list[Polynomial]
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / "evidence" / "cube-flat-seeds.json"


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            value = result.get(exponent, 0) ^ coefficient
            if value:
                result[exponent] = value
            else:
                result.pop(exponent, None)
    return result


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for a, ca in left.items():
        for b, cb in right.items():
            assert ca == cb == 1
            exponent = tuple(a[i] + b[i] for i in range(3))
            if exponent in result:
                del result[exponent]
            else:
                result[exponent] = 1
    return result


def antipode(word: int) -> int:
    return sum(1 << (v ^ 7) for v in range(8) if word >> v & 1)


def endpoint(word: int) -> int:
    return int(word == antipode(word))


def validate(tensor: Tensor, degree: int = 6) -> None:
    """Raise ValueError unless this is a literal homogeneous F2 cube tensor."""
    if len(tensor) != 256:
        raise ValueError("Expected exactly 256 physical coefficients.")
    for word, polynomial in enumerate(tensor):
        if not isinstance(polynomial, dict):
            raise ValueError("Each coefficient must be an exponent dictionary.")
        for exponent, coefficient in polynomial.items():
            if (not isinstance(exponent, tuple) or len(exponent) != 3
                    or any(not isinstance(e, int) or e < 0 for e in exponent)
                    or sum(exponent) != degree):
                raise ValueError("Every monomial must have the declared total degree.")
            if coefficient != 1:
                raise ValueError("F2 dictionaries retain only nonzero coefficient one.")
        if len(polynomial) % 2 != endpoint(word):
            raise ValueError(f"Antipodal equality endpoint failed at word {word}.")


def rank_two_seed() -> Tensor:
    """Return the even rank-two-face Gaussian quotient, without rational division.

    Put P={0,1,2,4}, J the antipodal matching, u=1_P, and
    C=J+u wedge Ju.  With q=(xyz)^2 and
    a=x^3(e0+e1)+y^3(e0+e2)+z^3(e0+e4), the kernel is
    K=q*C+a wedge Ja.  The returned coefficient is
    q*Pf(C_U)+sum_{i<j in U}(a_i*a_bar(j)+a_j*a_bar(i))*Pf(C_(U-ij)).
    It equals q^(1-|U|/2)*Pf(K_U) on even U and vanishes on odd U.
    """
    positive = {0, 1, 2, 4}
    C = [[int(w == (v ^ 7)) ^ int((v in positive) != (w in positive))
          for w in range(8)] for v in range(8)]
    pfaffians = [0] * 256
    pfaffians[0] = 1
    for word in range(1, 256):
        if word.bit_count() % 2:
            continue
        first = word & -word
        v = first.bit_length() - 1
        rest = word ^ first
        for w in range(8):
            if rest >> w & 1:
                pfaffians[word] ^= C[v][w] * pfaffians[rest ^ (1 << w)]
    a: list[Polynomial] = [{} for _ in range(8)]
    for axis in range(3):
        exponent = tuple(3 if j == axis else 0 for j in range(3))
        a[1 << axis] = {exponent: 1}
        a[0] = add(a[0], a[1 << axis])
    q = {(2, 2, 2): 1}
    result = []
    for word in range(256):
        polynomial = dict(q) if pfaffians[word] else {}
        vertices = [v for v in range(8) if word >> v & 1]
        for v, w in combinations(vertices, 2):
            if pfaffians[word ^ (1 << v) ^ (1 << w)]:
                polynomial = add(polynomial, multiply(a[v], a[w ^ 7]),
                                 multiply(a[w], a[v ^ 7]))
        result.append(polynomial)
    validate(result)
    return result


SWAP = ((0, 1), (1, 0))
LOWER = ((1, 0), (1, 1))
ODD_EDGE_MATRICES = (SWAP,) * 8 + (LOWER,) * 4
EDGES = tuple(tuple((v, v ^ (1 << axis)) for v in range(8)
                    if not (v >> axis & 1)) for axis in range(3))


def _edge_matrices(matrices) -> tuple:
    matrices = tuple(tuple(tuple(row) for row in matrix) for matrix in matrices)
    if len(matrices) != 12:
        raise ValueError("Expected four edge matrices per axis.")
    for matrix in matrices:
        if (len(matrix) != 2 or any(len(row) != 2 for row in matrix)
                or any(c not in (0, 1) for row in matrix for c in row)
                or not ((matrix[0][0] * matrix[1][1])
                        ^ (matrix[0][1] * matrix[1][0]))):
            raise ValueError("Every edge matrix must be invertible over F2.")
    return matrices


def _half_value(matrices: tuple, edge_axis: int, face_axis: int,
                side: int, word: int) -> int:
    result = 1
    for number, (v, w) in enumerate(EDGES[edge_axis]):
        if (v >> face_axis & 1) == side:
            result &= matrices[4 * edge_axis + number][word >> v & 1][word >> w & 1]
    return result


def edge_face_seed(matrices=ODD_EDGE_MATRICES) -> Tensor:
    """Construct the equal-degree rank-one-face family over F2.

    For axis a, A_a is the product of its four edge matrices.  For a,b
    and the remaining axis c, B_ab is the sum of the two cross products
    of the A_a and A_b half-face factors at cut c.  Set
      boundary=sum x_a^6*A_a + sum x_a^3*x_b^3*B_ab,
      S=boundary+(xyz)^2*(EQ-boundary(1,1,1)).
    This is homogeneous of degree six and has the literal endpoint.
    """
    matrices = _edge_matrices(matrices)
    result = []
    for word in range(256):
        polynomial: Polynomial = {}
        for axis in range(3):
            value = 1
            for number, (v, w) in enumerate(EDGES[axis]):
                value &= matrices[4 * axis + number][word >> v & 1][word >> w & 1]
            if value:
                polynomial[tuple(6 if j == axis else 0 for j in range(3))] = 1
        for a, b in combinations(range(3), 2):
            c = 3 - a - b
            value = (_half_value(matrices, a, c, 0, word)
                     * _half_value(matrices, b, c, 1, word))
            value ^= (_half_value(matrices, b, c, 0, word)
                      * _half_value(matrices, a, c, 1, word))
            if value:
                polynomial[tuple(3 if j in (a, b) else 0 for j in range(3))] = 1
        if len(polynomial) % 2 != endpoint(word):
            polynomial[(2, 2, 2)] = 1
        result.append(polynomial)
    validate(result)
    return result


def odd_edge_face_seed() -> Tensor:
    """Return the fixed parity-breaking rank-one-face reference over F2."""
    return edge_face_seed(ODD_EDGE_MATRICES)


def verify_edge_faces(tensor: Tensor, matrices=ODD_EDGE_MATRICES) -> int:
    """Check all 768 face entries against their literal rank-one products."""
    matrices = _edge_matrices(matrices)
    validate(tensor)
    checked = 0
    for axis in range(3):
        other = [a for a in range(3) if a != axis]
        for word in range(256):
            factors = []
            for side in (0, 1):
                polynomial = {}
                for a in other:
                    if _half_value(matrices, a, axis, side, word):
                        polynomial[tuple(3 if j == a else 0 for j in range(3))] = 1
                factors.append(polynomial)
            face = {e: c for e, c in tensor[word].items() if e[axis] == 0}
            assert face == multiply(*factors), (axis, word)
            checked += 1
    return checked


def statistics(tensor: Tensor) -> dict:
    validate(tensor)
    serial = [[[list(e), c] for e, c in sorted(polynomial.items())]
              for polynomial in tensor]
    return {"degree": 6, "endpoint_entries_checked": 256,
            "nonzero_coefficients": sum(bool(p) for p in tensor),
            "nonzero_even_coefficients": sum(bool(p) for i, p in enumerate(tensor)
                                             if i.bit_count() % 2 == 0),
            "nonzero_odd_coefficients": sum(bool(p) for i, p in enumerate(tensor)
                                            if i.bit_count() % 2),
            "monomials": sum(map(len, tensor)),
            "polynomial_sha256": hashlib.sha256(json.dumps(serial,
                separators=(",", ":")).encode()).hexdigest()}


def check() -> dict:
    even, odd = rank_two_seed(), odd_edge_face_seed()
    assert all(not p for i, p in enumerate(even) if i.bit_count() % 2)
    assert sum(bool(p) for i, p in enumerate(odd) if i.bit_count() % 2) == 64
    assert even[0] == {(2, 2, 2): 1}
    report = {"rank_two_seed": statistics(even), "odd_edge_face_seed": statistics(odd),
              "odd_seed_face_entries_checked": verify_edge_faces(odd),
              "generic_smith_profiles_checked_by_this_module": False}
    # Deliberate endpoint and degree corruptions must be rejected.
    bad = [dict(p) for p in odd]
    bad[0] = add(bad[0], {(6, 0, 0): 1})
    try:
        validate(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("Endpoint corruption was accepted.")
    bad = [dict(p) for p in odd]
    bad[0][(7, 0, 0)] = 1
    try:
        validate(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("Degree corruption was accepted.")
    report["negative_validation_controls"] = 2
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-evidence", action="store_true")
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    report = check()
    if args.write_evidence:
        args.evidence.parent.mkdir(parents=True, exist_ok=True)
        args.evidence.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    elif args.evidence.exists():
        assert json.loads(args.evidence.read_text()) == report
    print("PASS: two deterministic flat cube seeds; 512 endpoints; 768 rank-one face entries")


if __name__ == "__main__":
    main()
