#!/usr/bin/env python3
"""Sparse P6 isotropic boundary sections and their exact active operators."""

from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import json
import math
import random
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local
from check_noncommutative_crossing_pfaffian_boundaries import omega_boundary


PAIRS = tuple(combinations(range(6), 2))
ALIGNED_FACES = (
    ((0, 1), 0), ((2, 3), 16), ((0, 5), 4), ((1, 3), 17),
    ((1, 5), 21), ((0, 3), 6), ((2, 5), 26), ((3, 4), 35),
    ((1, 2), 40), ((0, 4), 32), ((1, 4), 45), ((0, 2), 8),
    ((2, 4), 34),
)


def support_vector(terms):
    value = 0
    for term in terms:
        value ^= 1 << term
    return value


def monomials(vector):
    while vector:
        bit = vector & -vector
        yield bit.bit_length() - 1
        vector ^= bit


def full_pair(index):
    return 3 << (2 * index)


def single_word(label):
    return sum(1 << (2 * i + ((label >> i) & 1)) for i in range(6))


def boundary_basis():
    stars = []
    for sigma in PAIRS:
        first_even = min(set(range(6)) - set(sigma))
        for bits in range(4):
            source = full_pair(first_even)
            for offset, i in enumerate(sigma):
                source |= 1 << (2 * i + ((bits >> offset) & 1))
            stars.append(support_vector(omega_boundary(source, 6)))

    def b(i, j):
        return support_vector(omega_boundary(full_pair(i) | full_pair(j), 6))

    u1, v1 = b(0, 1), b(0, 2)
    u2 = b(0, 3) ^ u1 ^ v1
    v2 = b(0, 4) ^ u1 ^ v1
    return tuple(stars), (u1, v1, u2, v2)


def construction(aligned=False):
    stars, hyperbolic = boundary_basis()
    columns = [0] * 64
    assigned = set()
    star_assignments = {}
    if aligned:
        for sigma, base in ALIGNED_FACES:
            mask = sum(1 << j for j in sigma)
            assert not base & mask
            for state in range(4):
                label = base | sum(((state >> offset) & 1) << j
                                   for offset, j in enumerate(sigma))
                assert label not in assigned
                assigned.add(label)
                index = 4 * PAIRS.index(sigma) + state
                star_assignments[index] = label
                columns[label] = 1 << index
    free_stars = [i for i in range(60) if i not in star_assignments]
    free_labels = [i for i in range(64) if i not in assigned]
    first, second = free_stars[:2]
    ordinary = free_stars[2:]
    assert len(free_labels) == len(ordinary) + 6
    for index, label in zip(ordinary, free_labels):
        columns[label] = 1 << index
    tail = free_labels[len(ordinary):]
    for start, norm_one, left, right in ((0, first, 60, 61), (3, second, 62, 63)):
        columns[tail[start]] = (1 << norm_one) ^ (1 << left)
        columns[tail[start + 1]] = (1 << norm_one) ^ (1 << right)
        columns[tail[start + 2]] = (1 << norm_one) ^ (1 << left) ^ (1 << right)
    assert all(columns)
    active_rows = [sum(((column >> i) & 1) << label
                       for label, column in enumerate(columns)) for i in range(64)]
    actual = []
    ambient = stars + hyperbolic
    for column in columns:
        value = 0
        for i in monomials(column):
            value ^= ambient[i]
        actual.append(value)
    gamma = [active_rows[4 * i:4 * i + 4] for i in range(15)]
    even = [0] * 64
    for a, b in ((active_rows[60], active_rows[61]),
                 (active_rows[62], active_rows[63])):
        for i in range(64):
            if a >> i & 1:
                even[i] ^= b
            if b >> i & 1:
                even[i] ^= a
    return columns, active_rows, actual, gamma, even


def flatten_matrix(rows):
    return sum(row << (64 * i) for i, row in enumerate(rows))


def add_outer(rows, first, second):
    for i in monomials(first):
        rows[i] ^= second


def compressed_operators(gamma, even):
    off_diagonal = []
    diagonal_projection = []
    identity_coefficient = list(even)
    for site in range(6):
        raising, lowering, projection = [0] * 64, [0] * 64, [0] * 64
        for sigma, block in zip(PAIRS, gamma):
            if site not in sigma:
                continue
            position = sigma.index(site)
            for state in range(4):
                add_outer(projection, block[state], block[state])
                if state >> position & 1:
                    add_outer(identity_coefficient, block[state], block[state])
                else:
                    add_outer(raising, block[state], block[state | (1 << position)])
                    add_outer(lowering, block[state | (1 << position)], block[state])
        assert rank_two(raising) == rank_two(lowering) == 10
        assert rank_two(projection) == 20
        off_diagonal += [flatten_matrix(raising), flatten_matrix(lowering)]
        diagonal_projection.append(flatten_matrix(projection))
    coefficient_vectors = off_diagonal + diagonal_projection + [flatten_matrix(identity_coefficient)]
    assert rank_two(coefficient_vectors) == 18
    assert not __import__("functools").reduce(int.__xor__, diagonal_projection)
    return {"off_diagonal_input_rank_each": 10,
            "diagonal_input_projection_rank_each": 20,
            "column_rank_Z1": 18,
            "independent_diagonal_projection_count": rank_two(diagonal_projection)}


def local_two_site_span():
    operators = []
    for sigma in PAIRS:
        for a, b in product(range(4), repeat=2):
            value = 0
            for source in range(64):
                bits = sum(((source >> site) & 1) << offset
                           for offset, site in enumerate(sigma))
                if bits != b:
                    continue
                target = source
                for offset, site in enumerate(sigma):
                    target = (target & ~(1 << site)) | (((a >> offset) & 1) << site)
                value ^= 1 << (64 * target + source)
            operators.append(value)
    assert rank_two(operators) == 154
    return 154


def first_jet_rank():
    matrix = [0] * 4096
    for x, y in product(range(64), repeat=2):
        row = 0
        for j in range(6):
            different = ((x ^ y) >> j) & 1
            u, v = x ^ (different << j), y ^ (different << j)
            row ^= 1 << (64 * u + v)
        matrix[64 * x + y] = row
    assert rank_two(matrix) == 2016
    return 2016


def general_boundary_gram_ranks():
    records = []
    for even_pairs in range(2, 13, 2):
        half = even_pairs // 2
        degree = half - 1
        sources = [sum(1 << i for i in subset)
                   for subset in combinations(range(even_pairs), degree)]
        positions = {mask: i for i, mask in enumerate(sources)}
        rows = []
        for i, mask in enumerate(sources):
            row = (1 << i) if (even_pairs - degree) & 1 else 0
            for a in range(even_pairs):
                if mask >> a & 1:
                    for b in range(even_pairs):
                        if not mask >> b & 1:
                            row ^= 1 << positions[mask ^ (1 << a) ^ (1 << b)]
            rows.append(row)
        expected = math.comb(even_pairs - 2, half - 2) if half >= 2 else 0
        assert rank_two(rows) == expected
        records.append({"even_pairs": even_pairs, "source_dimension": len(sources),
                        "boundary_Gram_rank": expected})
    return records


def sparse_all_size_section(pairs):
    assert pairs >= 7
    sectors = list(combinations(range(pairs), 4))[:16]
    support = set()
    for label in range(1 << pairs):
        even = sectors[label & 15]
        tail = label >> 4
        source = 3 << (2 * even[0])
        odd = sorted(set(range(pairs)) - set(even))
        for offset, site in enumerate(odd):
            source |= 1 << (2 * site + ((tail >> offset) & 1))
        boundary = omega_boundary(source, pairs)
        assert len(boundary) == 3
        assert all(term.bit_count() == pairs for term in boundary)
        assert not support.intersection(boundary)
        support.update(boundary)
    assert len(support) == 3 * (1 << pairs)
    return {"pairs": pairs, "labels": 1 << pairs, "chosen_even_sectors": 16,
            "boundary_terms_per_label": 3, "pairwise_disjoint_boundary_supports": True}


def literal_polynomial_checks(actual, gamma, even, count=512):
    p, q, r = 2, 1 << 7, 1 << 49
    local = evaluate_local((p, q, r, 1, 1, 1))
    a, b = r ^ q, p ^ q
    A = multiply(p, p ^ q ^ r)
    A2 = multiply(A, A)
    a2, ab, b2 = multiply(a, a), multiply(a, b), multiply(b, b)
    single_masks = [single_word(i) for i in range(64)]
    lifts = [actual[i] ^ (1 << single_masks[i]) for i in range(64)]
    rng = random.Random(20261030)

    @lru_cache(None)
    def literal(x, u, v, y):
        value = 1
        for j in range(6):
            xx, uu = 1 + ((x >> j) & 1), 1 + ((u >> j) & 1)
            vv, yy = (v >> (2 * j)) & 3, (y >> (2 * j)) & 3
            value = multiply(value, local[((xx * 4 + yy) * 4 + uu) * 4 + vv])
            if not value:
                break
        return value

    nonzero = 0
    for trial in range(count):
        x, y, u, v = (rng.randrange(64) for _ in range(4))
        if trial % 4 == 0:
            u, v = x, y
        elif trial % 4 == 1:
            sector = rng.randrange(15)
            sigma, block = PAIRS[sector], gamma[sector]
            state_v, state_y = rng.randrange(4), rng.randrange(4)
            v = rng.choice(tuple(monomials(block[state_v])))
            y = rng.choice(tuple(monomials(block[state_y])))
            for offset, j in enumerate(sigma):
                x = (x & ~(1 << j)) | (((state_v >> offset) & 1) << j)
            u = x
            for offset, j in enumerate(sigma):
                u = (u & ~(1 << j)) | (((state_y >> offset) & 1) << j)
        elif trial % 4 == 2:
            u = x
        total = 0
        for vv, yy in product(monomials(lifts[v]), monomials(lifts[y])):
            total ^= literal(x, u, vv, yy)
        baseline = 1
        for j in range(6):
            xx, yy, uu, vv = ((word >> j) & 1 for word in (x, y, u, v))
            factor = ((a if xx == vv and yy == uu else 0)
                      ^ (b if xx == uu and yy == vv else 0))
            baseline = multiply(baseline, factor)
        z0, z1 = 0, int(x == u and even[y] >> v & 1)
        for sigma, block in zip(PAIRS, gamma):
            selected = sum(1 << j for j in sigma)
            if (x ^ u) & ~selected:
                continue
            xx = sum(((x >> j) & 1) << i for i, j in enumerate(sigma))
            uu = sum(((u >> j) & 1) << i for i, j in enumerate(sigma))
            z0 ^= ((block[xx] >> v) & 1) & ((block[uu] >> y) & 1)
            for position in range(2):
                other = position ^ 1
                if (xx ^ uu) >> other & 1:
                    continue
                for c in range(2):
                    aa = (xx & (1 << position)) | (c << other)
                    bb = (uu & (1 << position)) | (c << other)
                    z1 ^= ((block[aa] >> v) & 1) & ((block[bb] >> y) & 1)
        z2 = int(x == u and y == v)
        correction = (a2 if z0 else 0) ^ (ab if z1 else 0) ^ (b2 if z2 else 0)
        expected = baseline ^ multiply(A2, correction)
        assert total == expected, (trial, x, y, u, v)
        nonzero += bool(total)
    assert nonzero >= count // 8
    return {"checked": count, "nonzero_polynomial_entries": nonzero}


def audit(aligned=False):
    columns, active, actual, gamma, even = construction(aligned)
    assert rank_two(active) == 64
    assert rank_two(even) == 4
    for i, first in enumerate(actual):
        boundary = 0
        for term in monomials(first):
            assert term.bit_count() == 6
            boundary ^= support_vector(omega_boundary(term, 6))
        assert boundary == 0
        assert not any(first >> single_word(j) & 1 for j in range(64))
        for j, second in enumerate(actual):
            assert (first & second).bit_count() % 2 == int(i == j)
    gram = list(even)
    for block in gamma:
        for vector in block:
            add_outer(gram, vector, vector)
    assert gram == [1 << i for i in range(64)]
    lifts = [actual[i] ^ (1 << single_word(i)) for i in range(64)]
    assert all(not ((a & b).bit_count() & 1) for a, b in product(lifts, repeat=2))
    operators = compressed_operators(gamma, even)
    checks = literal_polynomial_checks(actual, gamma, even)
    return {"name": "thirteen_aligned_faces" if aligned else "canonical_sparse",
            "labels": 64, "boundary_gram_is_identity": True,
            "lift_gram_is_zero": True, "active_coordinate_matrix_rank": 64,
            "even_Gram_rank": 4, "gamma_maps": gamma, "even_Gram": even,
            "active_columns": columns,
            "maximum_boundary_support": max(v.bit_count() for v in actual),
            "literal_polynomial_entries": checks,
            **operators}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [audit(False), audit(True)]
    result = {
        "status": "Exact isotropic sections and compressed Pfaffian identities.",
        "scope": "A positive construction; no successful factoring seed asserted.",
        "cases": cases,
        "aligned_faces": [[list(sigma), base] for sigma, base in ALIGNED_FACES],
        "two_site_operator_span_dimension": local_two_site_span(),
        "first_binary_jet_row_rank": first_jet_rank(),
        "isometry_invariant_correction_ranks": {
            "Z0": [3840, 154], "Z1": [1696, 18],
            "Z0_plus_Z1": [4096, 154], "Z2": [4096, 1]},
        "first_jet_row_rank_lower": 320,
        "off_diagonal_first_jet_rank_lower": 22,
        "all_size_boundary_Gram_rank_checks": general_boundary_gram_ranks(),
        "all_size_sparse_section_checks": [sparse_all_size_section(p) for p in range(7, 11)],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
