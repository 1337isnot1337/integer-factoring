#!/usr/bin/env python3
"""Exploratory graded Hopf crossings for restricted upper-triangular algebras."""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

from check_noncommutative_crossing_orbit import binary_rank
from check_noncommutative_crossing_hopf import coefficients, column_ranks, parameter_translate


def subsets(value):
    sub = value
    while True:
        yield sub
        if not sub:
            break
        sub = (sub - 1) & value


class RestrictedUpperTriangular:
    def __init__(self, order):
        self.order = order
        self.generators = tuple((i, i + gap)
                                for gap in range(1, order) for i in range(order - gap))
        self.positions = {pair: i for i, pair in enumerate(self.generators)}
        self.width = len(self.generators)
        self.d = 1 << self.width
        self.weights = tuple(j - i for i, j in self.generators)
        self.degree = sum(self.weights)
        self.basis_weights = tuple(sum(w for i, w in enumerate(self.weights) if mask >> i & 1)
                                   for mask in range(self.d))
        self.words = tuple(tuple(i for i in range(self.width) if mask >> i & 1)
                           for mask in range(self.d))

    def bracket(self, left, right):
        i, j = self.generators[left]
        k, l = self.generators[right]
        if j == k:
            return self.positions[i, l]
        if l == i:
            return self.positions[k, j]
        return None

    @lru_cache(None)
    def normal(self, word):
        for i in range(len(word) - 1):
            left, right = word[i:i + 2]
            if left == right:
                return ()
            if left > right:
                result = set(self.normal(word[:i] + (right, left) + word[i + 2:]))
                commutator = self.bracket(left, right)
                if commutator is not None:
                    result.symmetric_difference_update(self.normal(
                        word[:i] + (commutator,) + word[i + 2:]))
                return tuple(sorted(result))
        return (sum(1 << i for i in word),)

    @lru_cache(None)
    def multiply(self, left, right):
        result = self.normal(self.words[left] + self.words[right])
        assert all(self.basis_weights[x] == self.basis_weights[left] + self.basis_weights[right]
                   for x in result)
        return result

    def multiply_support(self, left, right):
        result = set()
        for a in left:
            for b in right:
                result.symmetric_difference_update(self.multiply(a, b))
        return result

    @lru_cache(None)
    def antipode(self, value):
        return self.normal(self.words[value][::-1])

    def validate(self):
        for a, b, c in itertools.product(range(self.d), repeat=3):
            assert self.multiply_support(self.multiply(a, b), (c,)) == \
                self.multiply_support((a,), self.multiply(b, c))
        for a in range(self.d):
            convolution = set()
            for left in subsets(a):
                convolution.symmetric_difference_update(self.multiply_support(
                    (left,), self.antipode(a ^ left)))
            assert convolution == ({0} if a == 0 else set())

    def literal_image(self, u, v):
        result = {}
        for u2 in subsets(u):
            for v2 in subsets(v):
                degree = self.basis_weights[u2] + self.basis_weights[v2]
                for a in self.multiply_support((u2,), self.antipode(v2)):
                    for a1 in subsets(a):
                        for x in self.multiply(v ^ v2, a1):
                            for y in self.multiply(u ^ u2, a ^ a1):
                                key = x | (y << self.width)
                                result[key] = result.get(key, 0) ^ (1 << degree)
        return {key: p for key, p in result.items() if p}

    def generator_images(self):
        result = [{} for _ in range(2 * self.width)]
        for u in range(self.d):
            for v in range(self.d):
                weight = self.basis_weights[u] + self.basis_weights[v]
                if weight not in self.weights:
                    continue
                image = self.literal_image(u, v)
                for side in range(2):
                    for j, degree in enumerate(self.weights):
                        if degree != weight:
                            continue
                        source = 1 << (j + side * self.width)
                        coefficient = image.get(source, 0)
                        if coefficient:
                            result[j + side * self.width][u | (v << self.width)] = coefficient
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, choices=(3, 4), default=4)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    algebra = RestrictedUpperTriangular(args.order)
    algebra.validate()
    images = algebra.generator_images()
    result = {"status": "exact restricted-algebra and Hopf-generator construction; no stronger exponent",
              "matrix_order": args.order, "generators": algebra.generators,
              "weights": algebra.weights, "label_bits": algebra.width,
              "labels": algebra.d, "degree": algebra.degree,
              "generator_substitution": [
                  {"variable": i, "terms": [[mask, p] for mask, p in sorted(image.items())]}
                  for i, image in enumerate(images)]}
    if args.profile:
        matrices = coefficients([parameter_translate(image) for image in images])
        assert len(matrices) == algebra.degree + 1
        for source, _ in enumerate(matrices[0]):
            summed = 0
            for matrix in matrices:
                summed ^= matrix[source]
            left, right = source & (algebra.d - 1), source >> algebra.width
            assert summed == 1 << (right | (left << algebra.width))
        rows = [binary_rank(matrix) for matrix in matrices]
        columns = column_ranks(matrices, algebra.width)
        result["row_ranks"] = rows
        result["column_ranks"] = columns
        result["rank_products"] = [a * b for a, b in zip(rows, columns)]
        result["maximum_rank_product"] = max(result["rank_products"])
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
