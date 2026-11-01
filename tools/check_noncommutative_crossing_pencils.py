#!/usr/bin/env python3
"""Exact Hopf profiles for arbitrary alternating two-step commutator maps."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess

from check_noncommutative_crossing_hopf import (
    add, coefficients, column_ranks, compose, evaluate_one, multiply,
    parameter_translate, scale, variable,
)
from check_noncommutative_crossing_orbit import binary_rank
from check_noncommutative_crossing_ut import RestrictedUpperTriangular


MASK64 = (1 << 64) - 1


def mix64(value):
    value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
    return value ^ (value >> 31)


def matrix_digest(matrices):
    xor = total = count = 0
    for source in range(len(matrices[0])):
        entries = {}
        for degree, matrix in enumerate(matrices):
            bits = matrix[source]
            while bits:
                bit = bits & -bits
                target = bit.bit_length() - 1
                entries[target] = entries.get(target, 0) ^ (1 << degree)
                bits ^= bit
        for target, polynomial in entries.items():
            token = mix64(mix64(source ^ (target << 32)) ^ polynomial)
            xor ^= token
            total = (total + token) & MASK64
            count += 1
    return {"matrix_digest_xor": str(xor), "matrix_digest_sum": str(total),
            "nonzero_polynomial_entries": count}


class LiteralTwoStep(RestrictedUpperTriangular):
    """Independent PBW rewriting and literal Hopf operations for low weights."""
    def __init__(self, horizontal, center, edges):
        self.width = horizontal + center
        self.d = 1 << self.width
        self.weights = (1,) * horizontal + (2,) * center
        self.degree = sum(self.weights)
        self.basis_weights = tuple(sum(w for i, w in enumerate(self.weights) if mask >> i & 1)
                                   for mask in range(self.d))
        self.words = tuple(tuple(i for i in range(self.width) if mask >> i & 1)
                           for mask in range(self.d))
        self.edge_map = {(a, b): tuple(horizontal + j for j in range(center) if vector >> j & 1)
                         for a, b, vector in edges}

    def normal(self, word):
        # The bounded audit uses direct recursion; it shares neither the
        # quadratic substitution nor the native polynomial builder.
        for i in range(len(word) - 1):
            a, b = word[i:i + 2]
            if a == b:
                return ()
            if a > b:
                out = set(self.normal(word[:i] + (b, a) + word[i + 2:]))
                for central in self.edge_map.get((b, a), ()):
                    out.symmetric_difference_update(self.normal(word[:i] + (central,) + word[i + 2:]))
                return tuple(sorted(out))
        return (sum(1 << i for i in word),)


def generator_maps(horizontal, center, edges):
    width = horizontal + center
    identity = [variable(i) for i in range(2 * width)]
    R, M, K = [], [], []
    translations = []
    for side in range(2):
        for j in range(width):
            own, other = variable(side * width + j), variable((1 - side) * width + j)
            if j < horizontal:
                image = add(scale(own, 2), scale(other, 3))
                R.append(image)
                M.append(image)
                K.append(own)
                continue
            central = j - horizontal
            linear = add(scale(own, 4), scale(other, 5))
            rquad, kquad = {}, {}
            for a, b, vector in edges:
                assert 0 <= a < b < horizontal and 0 < vector < 1 << center
                if not vector >> central & 1:
                    continue
                xL, yL = variable(a), variable(b)
                xR, yR = variable(width + a), variable(width + b)
                terms = [multiply(xL, yL), multiply(xL, yR),
                         multiply(yL, xR), multiply(xR, yR)]
                rcoef = (0, 2, 4, 6) if side == 0 else (2, 0, 6, 4)
                kcoef = (10, 8, 14, 12) if side == 0 else (8, 10, 12, 14)
                rquad = add(rquad, *(scale(term, c) for term, c in zip(terms, rcoef)))
                kquad = add(kquad, *(scale(term, c) for term, c in zip(terms, kcoef)))
            R.append(add(linear, rquad))
            M.append(linear)
            K.append(add(own, kquad))
            translations.append((side * width + j, kquad))
    assert compose(M, K) == compose(K, M) == R
    return R, translations


def compressed_translation(horizontal, center, translations, diagnostics=None):
    width, d = horizontal + center, 1 << (horizontal + center)
    polynomial = {0: 1}
    for derivative, translated in translations:
        # Central variables here label derivative operators, not coordinates.
        factor = add({0: 1}, multiply(variable(derivative), translated))
        polynomial = multiply(polynomial, factor)
    shifted = parameter_translate(polynomial)
    maximum = max(p.bit_length() - 1 for p in shifted.values())
    assert maximum <= 3 * center
    row_labels = sorted({monomial & (d - 1) for monomial in shifted})
    column_labels = sorted({monomial >> width for monomial in shifted})
    row_positions = {label: i for i, label in enumerate(row_labels)}
    column_positions = {label: i for i, label in enumerate(column_labels)}
    matrices = [[0] * len(row_labels) for _ in range(maximum + 1)]
    for monomial, p in shifted.items():
        row = row_positions[monomial & (d - 1)]
        column = column_positions[monomial >> width]
        while p:
            bit = p & -p
            matrices[bit.bit_length() - 1][row] ^= 1 << column
            p ^= bit
    if diagnostics is not None:
        diagnostics.update({"row_operator_monomials": len(row_labels),
                            "column_operator_monomials": len(column_labels),
                            "polynomial_entries": len(shifted),
                            "translation_parameter_degree": maximum})
    return [binary_rank(matrix) for matrix in matrices]


def common_coordinate_abelian_size(horizontal, edges):
    best = 0
    for mask in range(1 << horizontal):
        if all(not (mask >> a & 1 and mask >> b & 1) for a, b, _ in edges):
            best = max(best, mask.bit_count())
    return best


def case(name, horizontal, center, edges, full, native=None):
    R, translations = generator_maps(horizontal, center, edges)
    diagnostics = {}
    compressed = compressed_translation(horizontal, center, translations, diagnostics)
    result = {"name": name, "horizontal_generators": horizontal, "central_generators": center,
              "edges": [list(e) for e in edges], "label_bits": horizontal + center,
              "labels": 1 << (horizontal + center), "degree": horizontal + 2 * center,
              "constant_column_rank": compressed[0],
              "translation_coefficient_column_ranks": compressed,
              "translation_compression": diagnostics,
              "largest_coordinate_abelian_horizontal_subspace": common_coordinate_abelian_size(horizontal, edges)}
    if full:
        literal = LiteralTwoStep(horizontal, center, edges)
        assert literal.generator_images() == R
        matrices = coefficients([parameter_translate(image) for image in R])
        assert len(matrices) == horizontal + 2 * center + 1
        rows = [binary_rank(matrix) for matrix in matrices]
        cols = column_ranks(matrices, horizontal + center)
        assert cols[0] == compressed[0]
        d = 1 << (horizontal + center)
        for source in range(d * d):
            total = 0
            for matrix in matrices:
                total ^= matrix[source]
            left, right = source & (d - 1), source >> (horizontal + center)
            assert total == 1 << (right | (left << (horizontal + center)))
        result.update({"row_ranks": rows, "column_ranks": cols,
                       "rank_products": [r * c for r, c in zip(rows, cols)],
                       "maximum_rank_product": max(r * c for r, c in zip(rows, cols))})
        digest = matrix_digest(matrices)
        result.update(digest)
        result["literal_Hopf_generator_audit"] = True
        if native is not None and horizontal % 2 == 0:
            events = [json.loads(line) for line in subprocess.check_output(
                [str(native), "--variables", str(horizontal), "--center", str(center),
                 "--edges", ",".join(":".join(map(str, edge)) for edge in edges)], text=True).splitlines()]
            profile = [event for event in events if event["event"] == "coefficient" and event["row_rank"]]
            assert [event["row_rank"] for event in profile] == rows
            assert [event["column_rank"] for event in profile] == cols
            assert all(events[-1][key] == value for key, value in digest.items())
            result["native_full_profile_and_matrix_digest_match"] = True
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--local-ring-max", type=int, default=0)
    args = parser.parse_args()
    examples = [
        ("split_F2_pairs", 4, 2, ((0, 1, 1), (2, 3, 2)), True),
        ("F4_field_pencil", 4, 2, ((0, 2, 1), (0, 3, 2), (1, 2, 2), (1, 3, 3)), True),
        ("nilpotent_pencil", 4, 2, ((0, 2, 1), (1, 3, 1), (0, 3, 2)), True),
        ("free_bipartite_1_by_2", 3, 2, ((0, 1, 1), (0, 2, 2)), True),
        ("F8_projected_pencil", 6, 2,
         ((0, 3, 1), (0, 4, 2), (1, 3, 2), (1, 5, 3), (2, 4, 3), (2, 5, 2)), False),
        ("path_four_unique_centers", 4, 3, ((0, 1, 1), (1, 2, 2), (2, 3, 4)), False),
    ]
    for k in range(1, args.local_ring_max + 1):
        edges = tuple((i, k + j, 1 << (i + j))
                      for i in range(k) for j in range(k) if i + j < k)
        examples.append((f"F2_t_mod_t{k}", 2 * k, k, edges, k <= 2))
    result = {"status": "exact alternating-commutator constructions; no stronger factoring exponent",
              "cases": [case(*example, native=args.native) for example in examples],
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
