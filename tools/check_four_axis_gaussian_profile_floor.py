#!/usr/bin/env python3
"""Exact finite controls for the four-axis shared-minor profile floor.

This uses four-variable F2 polynomials, never a specialization of the other
coordinates to one when computing a coordinate valuation. The script checks
the proof's finite indexing and polynomial identities, not all kernels.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random


BITS = 6  # All checked determinant coordinate degrees are at most 48.
MASK = (1 << BITS) - 1
ZERO = frozenset()
ONE = frozenset({0})


def encode(exponents: tuple[int, ...]) -> int:
    assert len(exponents) == 4 and all(0 <= exponent <= MASK for exponent in exponents)
    return sum(exponent << (BITS * axis) for axis, exponent in enumerate(exponents))


def decode(monomial: int) -> tuple[int, ...]:
    return tuple((monomial >> (BITS * axis)) & MASK for axis in range(4))


def multiply(left: frozenset[int], right: frozenset[int]) -> frozenset[int]:
    result = set()
    for a in left:
        for b in right:
            monomial = a + b
            if monomial in result:
                result.remove(monomial)
            else:
                result.add(monomial)
    return frozenset(result)


def determinant(matrix: list[list[frozenset[int]]]) -> frozenset[int]:
    size = len(matrix)
    result = ZERO
    for permutation in permutations(range(size)):
        term = ONE
        for row, column in enumerate(permutation):
            term = multiply(term, matrix[row][column])
            if not term:
                break
        result = result.symmetric_difference(term)
    return frozenset(result)


def compositions(total: int, length: int = 4):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for tail in compositions(total - first, length - 1):
                yield (first,) + tail


def random_kernel(degree: int, seed: int, common_factor=(0, 0, 0, 0)):
    rng = random.Random(seed)
    residual_degree = degree - sum(common_factor)
    monomials = [encode(tuple(x + y for x, y in zip(exponents, common_factor)))
                 for exponents in compositions(residual_degree)]
    kernel = [[ZERO for _ in range(16)] for _ in range(16)]
    for u, v in combinations(range(16), 2):
        endpoint = int(v == (u ^ 15))
        allowed = [count for count in range(min(4, len(monomials)) + 1)
                   if count % 2 == endpoint]
        count = rng.choice(allowed)
        value = frozenset(rng.sample(monomials, count))
        kernel[u][v] = kernel[v][u] = value
    return kernel


def cayley_control():
    degree = 12
    coefficients = {}
    for subset in range(1, 16):
        exponent = degree // subset.bit_count()
        coefficients[subset] = frozenset({encode(tuple(
            exponent if subset >> axis & 1 else 0 for axis in range(4)
        ))})
    steps = [ZERO] * 16
    for step in range(1, 16):
        for subset, value in coefficients.items():
            if subset & step == step:
                steps[step] = steps[step].symmetric_difference(value)
    return degree, [[steps[u ^ v] for v in range(16)] for u in range(16)]


def is_cayley(kernel) -> bool:
    return all(kernel[u][v] == kernel[0][u ^ v] for u in range(16) for v in range(16))


def shared_minor_controls() -> dict:
    settings = [
        (1, 11, (0, 0, 0, 0)), (1, 17, (0, 0, 0, 0)),
        (2, 23, (0, 0, 0, 0)), (2, 29, (0, 0, 0, 0)),
        (3, 31, (0, 0, 0, 0)), (3, 37, (0, 0, 0, 0)),
        (4, 41, (0, 0, 0, 0)), (4, 43, (0, 0, 0, 0)),
        (5, 47, (1, 1, 1, 1)), (5, 53, (1, 1, 1, 1)),
        (8, 59, (1, 1, 1, 1)), (8, 61, (1, 1, 1, 1)),
    ]
    kernels = [(degree, f"seed-{seed}", random_kernel(degree, seed, common))
               for degree, seed, common in settings]
    degree, kernel = cayley_control()
    kernels.append((degree, "symmetric-cayley-control", kernel))
    records = []
    minors = non_cayley = 0
    for degree, name, kernel in kernels:
        for u, v in product(range(16), repeat=2):
            assert kernel[u][v] == kernel[v][u]
            assert bool(len(kernel[u][v]) % 2) == (v == (u ^ 15))
            assert all(sum(decode(monomial)) == degree for monomial in kernel[u][v])
            if u == v:
                assert not kernel[u][v]
        cayley = is_cayley(kernel)
        non_cayley += not cayley
        crossblocks = {}
        for axis in range(4):
            rows = [u for u in range(16) if not (u >> axis & 1)]
            columns = [v for v in range(16) if v >> axis & 1]
            matrix = [[kernel[u][v] for v in columns] for u in rows]
            crossblocks[axis] = rows, columns, matrix
        upper = [4 * degree] * 4
        pair_records = []
        for a, b in combinations(range(4), 2):
            rows = [u for u in range(16) if not (u >> a & 1) and not (u >> b & 1)]
            columns = [u ^ 15 for u in rows]
            matrix = [[kernel[u][v] for v in columns] for u in rows]
            endpoint = [[len(entry) % 2 for entry in row] for row in matrix]
            assert endpoint == [[int(i == j) for j in range(4)] for i in range(4)]
            for axis in (a, b):
                ambient_rows, ambient_columns, ambient = crossblocks[axis]
                extracted = [[ambient[ambient_rows.index(u)][ambient_columns.index(v)]
                              for v in columns] for u in rows]
                assert extracted == matrix
            polynomial = determinant(matrix)
            assert polynomial and len(polynomial) % 2 == 1
            exponents = [decode(monomial) for monomial in polynomial]
            assert all(sum(exponent) == 4 * degree for exponent in exponents)
            va, vb = min(exponent[a] for exponent in exponents), min(exponent[b] for exponent in exponents)
            assert va + vb <= 4 * degree
            upper[a] = min(upper[a], va)
            upper[b] = min(upper[b], vb)
            pair_records.append({"axes": [a, b], "rows": rows, "antipodal_columns": columns,
                                 "determinant_monomials": len(polynomial), "valuations": [va, vb]})
            minors += 1
        assert sum(upper) <= 8 * degree
        records.append({"name": name, "degree": degree, "cayley": cayley,
                        "shared_minor_upper_bounds_for_Delta4": upper,
                        "sum_upper_bounds": sum(upper), "budget": 8 * degree,
                        "minors": pair_records})
    assert non_cayley == len(settings)
    return {"kernels": len(kernels), "non_cayley_kernels": non_cayley,
            "exact_shared_determinants": minors, "cases": records}


def augmentation_minor_controls() -> dict:
    cases = 0
    for a in range(4):
        others = [axis for axis in range(4) if axis != a]
        subsets = [sum(((word >> j) & 1) << axis for j, axis in enumerate(others))
                   for word in range(8)]
        for b in others:
            rows = [subset for subset in subsets if subset >> b & 1]
            columns = [subset for subset in subsets if not (subset >> b & 1)]
            terms = set()
            for permutation in permutations(range(4)):
                monomial = []
                for row, column in enumerate(permutation):
                    u, v = rows[row], columns[column]
                    if u & v != v:
                        break
                    monomial.append((1 << a) | (u ^ v))
                else:
                    term = tuple(sorted(monomial))
                    if term in terms:
                        terms.remove(term)
                    else:
                        terms.add(term)
            assert terms == {((1 << a) | (1 << b),) * 4}
            cases += 1
    return {"symbolic_coefficient_minors": cases, "identity": "determinant = c_ab^4"}


def entropy_controls() -> dict:
    cases = equalities = 0
    for denominator in range(2, 17):
        for numerator in range(1, denominator):
            for s in range(2, 13):
                for r in range(1, s):
                    left = (r + s)**denominator * numerator**numerator * (
                        denominator - numerator
                    )**(denominator - numerator)
                    right = denominator**denominator * r**numerator * s**(denominator - numerator)
                    assert left >= right
                    equality = numerator * (r + s) == denominator * r
                    assert (left == right) == equality
                    equalities += equality
                    cases += 1
    return {"exact_rational_cases": cases, "equality_cases": equalities,
            "identity": "(r+s)^b a^a (b-a)^(b-a) >= b^b r^a s^(b-a)"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = {
        "status": "passed",
        "scope": "Exact shared-minor and entropy controls; the universal floor is the written proof.",
        "shared_minors": shared_minor_controls(),
        "augmentation_minors": augmentation_minor_controls(),
        "entropy": entropy_controls(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "kernels": result["shared_minors"]["kernels"],
        "non_cayley_kernels": result["shared_minors"]["non_cayley_kernels"],
        "shared_determinants": result["shared_minors"]["exact_shared_determinants"],
        "symbolic_augmentation_minors": result["augmentation_minors"]["symbolic_coefficient_minors"],
        "entropy_cases": result["entropy"]["exact_rational_cases"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
