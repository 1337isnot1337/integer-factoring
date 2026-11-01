#!/usr/bin/env python3
"""Exact one-defect profiles of the PR24 three-label crossing seed.

Only 15 equality patterns and 16 internal-color patterns are required.
The profiles are exact finite data, not an extrapolation of rank growth.
"""

from collections import Counter
from itertools import product
from pathlib import Path
import argparse
import json

from check_finite_seed_amplification import PIECES, SOURCE_COMMIT, crossing, rank


def canonical(values):
    names = {}
    return tuple(names.setdefault(value, len(names)) for value in values)


def equality_patterns():
    return tuple(sorted({canonical(values) for values in product(range(4), repeat=4)}))


def refines(first, second):
    return all(first[i] != first[j] or second[i] == second[j]
               for i in range(4) for j in range(4))


def local_polynomial(values):
    x, y, u, v = values
    result = set()
    for i, matrix in enumerate(PIECES):
        if (matrix[3 * x + y] >> (3 * u + v)) & 1:
            result.add(()) if i == 12 else result.add((i,))
    return result


def polynomial_multiply(first, second):
    answer = set()
    for left in first:
        for right in second:
            monomial = tuple(sorted(left + right))
            if monomial in answer:
                answer.remove(monomial)
            else:
                answer.add(monomial)
    return answer


def orbit_profiles():
    patterns = equality_patterns()
    assert len(patterns) == 15
    assert local_polynomial((2, 2, 2, 2)) == {()}
    profiles, vectors = [], {}
    for pattern_index, pattern in enumerate(patterns):
        for color_index, colors in enumerate(product(range(2), repeat=4)):
            polynomial = {()}
            for position in range(max(pattern) + 1):
                values = tuple(colors[j] if pattern[j] == position else 2 for j in range(4))
                polynomial = polynomial_multiply(polynomial, local_polynomial(values))
            assert all(len(term) <= 4 for term in polynomial)
            profiles.append(tuple(sorted(polynomial)))
            for term in polynomial:
                vectors[term] = vectors.get(term, 0) | (1 << (16 * pattern_index + color_index))
    return patterns, profiles, vectors


def diagram_coordinates(patterns, exact_vector):
    """Möbius inversion from exact equality orbits to equality diagrams."""
    coordinates = [0] * len(patterns)
    for i in sorted(range(len(patterns)), key=lambda j: -len(set(patterns[j]))):
        value = (exact_vector >> (16 * i)) & 65535
        for j, finer in enumerate(patterns):
            if j != i and refines(finer, patterns[i]):
                value ^= coordinates[j]
        coordinates[i] = value
    for i, pattern in enumerate(patterns):
        reconstructed = 0
        for j, finer in enumerate(patterns):
            if refines(finer, pattern):
                reconstructed ^= coordinates[j]
        assert reconstructed == ((exact_vector >> (16 * i)) & 65535)
    return tuple(coordinates)


def matrices(copies, patterns, profiles, terms):
    dimension = 2 * copies
    square = dimension ** 2
    rows = {term: [0] * square for term in terms}
    columns = {term: [0] * square for term in terms}
    lookup = {pattern: i for i, pattern in enumerate(patterns)}
    for x, y, u, v in product(range(dimension), repeat=4):
        values = (x, y, u, v)
        pattern = canonical(tuple(value // 2 for value in values))
        color = sum((value & 1) << (3 - j) for j, value in enumerate(values))
        for term in profiles[16 * lookup[pattern] + color]:
            rows[term][x * dimension + y] |= 1 << (u * dimension + v)
            columns[term][x * dimension + u] |= 1 << (y * dimension + v)
    total = [0] * square
    for values in rows.values():
        for j, value in enumerate(values):
            total[j] ^= value
    assert tuple(total) == crossing(dimension)
    return rows, columns


def verify(max_copies):
    patterns, profiles, vectors = orbit_profiles()
    terms = tuple(sorted(vectors, key=lambda term: (len(term), term)))
    diagrams = {term: diagram_coordinates(patterns, vector) for term, vector in vectors.items()}
    crossing_index = patterns.index((0, 1, 1, 0))
    witnesses = [term for term in terms if diagrams[term][crossing_index]]
    assert witnesses
    summed = [0] * 15
    for coordinates in diagrams.values():
        for j, value in enumerate(coordinates):
            summed[j] ^= value
    crossing_color = sum(1 << color for color, (x, y, u, v) in
                         enumerate(product(range(2), repeat=4)) if x == v and y == u)
    assert summed == [crossing_color if j == crossing_index else 0 for j in range(15)]
    # Exact all-copy witness: only the equality pattern xv|yu with its
    # two positions DISTINCT survives. The all-equal orbit is zero.
    assert vectors[(17, 19)] == crossing_color << (16 * crossing_index)
    assert diagrams[(17, 19)] == tuple(
        crossing_color if pattern in ((0, 0, 0, 0), (0, 1, 1, 0)) else 0
        for pattern in patterns)

    all_profiles = []
    for copies in range(1, max_copies + 1):
        rows, columns = matrices(copies, patterns, profiles, terms)
        records = []
        for term in terms:
            r, c = rank(rows[term]), rank(columns[term])
            if r or c:
                assert r and c
                records.append({"nonvacuum_type": term, "row_rank": r, "column_rank": c,
                                "crossing_diagram_color_mask": diagrams[term][crossing_index]})
        if copies >= 2:
            witness = next(item for item in records if item["nonvacuum_type"] == (17, 19))
            assert witness["row_rank"] == witness["column_rank"] == 4 * copies * (copies - 1)
        all_profiles.append({"copies": copies, "dimension": 2 * copies,
                             "nonzero_pieces": len(records),
                             "maximum_rank_product": max(item["row_rank"] * item["column_rank"]
                                                         for item in records),
                             "pieces": records})
    # The first tensor power is just the restriction of each seed matrix to U.
    one_rows, _ = matrices(1, patterns, profiles, terms)
    for index, matrix in enumerate(PIECES):
        term = () if index == 12 else (index,)
        if term not in one_rows:
            assert not any((matrix[3*x+y] >> (3*u+v)) & 1
                           for x, y, u, v in product(range(2), repeat=4))
        else:
            expected = [sum(((matrix[3*x+y] >> (3*u+v)) & 1) << (2*u+v)
                            for u, v in product(range(2), repeat=2))
                        for x, y in product(range(2), repeat=2)]
            assert one_rows[term] == expected
    return {
        "scope": "exact one-defect source profiles and equality-diagram coefficients; no empirical asymptotic rank claim",
        "source_commit": SOURCE_COMMIT,
        "equality_patterns": patterns,
        "internal_color_order": "x,y,u,v; lexicographic binary order",
        "orbit_entries": len(profiles),
        "type_count": len(terms),
        "type_span_dimension": rank(list(vectors.values())),
        "maximum_nonvacuum_degree": max(map(len, terms)),
        "exact_offdiagonal_crossing_witness": {
            "nonvacuum_type": (17, 19),
            "remaining_vacuum_power": "t-2",
            "both_ranks": "4*t*(t-1)",
            "identity": "[pos_x=pos_v != pos_y=pos_u] [color_x=color_v] [color_y=color_u]",
        },
        "crossing_diagram_witnesses": witnesses,
        "diagram_types": [{"nonvacuum_type": term, "exact_orbit_bits": str(vectors[term]),
                           "diagram_color_masks": diagrams[term]} for term in terms],
        "finite_profiles": all_profiles,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-copies", type=int, default=4)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if not 1 <= arguments.max_copies <= 8:
        raise ValueError("use 1 through 8 copies for this finite checker")
    result = verify(arguments.max_copies)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("diagram_types", "finite_profiles")}, indent=2, sort_keys=True))
    print(json.dumps([{"copies": item["copies"], "dimension": item["dimension"],
                       "nonzero_pieces": item["nonzero_pieces"],
                       "maximum_rank_product": item["maximum_rank_product"]}
                      for item in result["finite_profiles"]], indent=2))
