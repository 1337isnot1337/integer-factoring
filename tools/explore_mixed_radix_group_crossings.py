#!/usr/bin/env python3
"""Exact mixed-radix group-word crossings in characteristic two.

This includes groups whose order is not a power of the coefficient
characteristic. No augmentation-unit or binary-conservation bound is assumed.
"""

import argparse
import hashlib
import json
from math import isqrt, prod, sqrt
from pathlib import Path
import random

from check_modular_crossing import rank_two
from explore_complement_edge_pencils import minimum_decomposition


class Group:
    def __init__(self, name):
        self.name = name
        self.radices = [3, 2] if name == "S3" else [int(s[1:]) for s in name.split("x")]
        self.d = prod(self.radices)
        self.generators = []
        scale = 1
        for radix in self.radices:
            self.generators.append(scale)
            scale *= radix
        self.table = [[self.multiply(x, y) for y in range(self.d)] for x in range(self.d)]
        self.inverse = [next(y for y in range(self.d) if self.table[x][y] == 0)
                        for x in range(self.d)]
        for x in range(self.d):
            assert self.table[self.inverse[x]][x] == 0

    def digits(self, x):
        result = []
        for radix in self.radices:
            result.append(x % radix)
            x //= radix
        return result

    def encode(self, digits):
        result, scale = 0, 1
        for digit, radix in zip(digits, self.radices):
            result += scale * (digit % radix)
            scale *= radix
        return result

    def multiply(self, x, y):
        a, b = self.digits(x), self.digits(y)
        if self.name == "S3":
            return self.encode(((a[0] + (-1 if a[1] else 1) * b[0]) % 3, (a[1] + b[1]) % 2))
        return self.encode([u + v for u, v in zip(a, b)])


def word_kernel(group, word, weights):
    values = {0: 1}
    for generator, digit, weight in zip(group.generators, group.digits(word), weights):
        for _ in range(digit):
            following = {}
            for value, polynomial in values.items():
                following[value] = following.get(value, 0) ^ (polynomial << weight)
                other = group.table[value][generator]
                following[other] = following.get(other, 0) ^ polynomial ^ (polynomial << weight)
            values = {a: polynomial for a, polynomial in following.items() if polynomial}
    assert all((polynomial.bit_count() & 1) == (a == 0) for a, polynomial in values.items())
    assert 0 in values and values[0].bit_count() & 1
    return values


def build(group, weights, seed):
    d = group.d
    controls = list(range(d))
    if seed:
        random.Random(seed).shuffle(controls)
    kernels = [word_kernel(group, control, weights) for control in controls]
    degree = sum((radix - 1) * weight for radix, weight in zip(group.radices, weights))
    coefficients = [0] * (degree + 1)
    for S, kernel in enumerate(kernels):
        for a, polynomial in kernel.items():
            while polynomial:
                bit = polynomial & -polynomial
                coefficients[bit.bit_length() - 1] ^= 1 << (S * d + a)
                polynomial ^= bit
    pivots, selected = {}, []
    for k, value in enumerate(coefficients):
        coordinates = 0
        while value:
            lead = value.bit_length() - 1
            if lead in pivots:
                row, old = pivots[lead]
                value ^= row
                coordinates ^= old
            else:
                coordinates ^= 1 << len(selected)
                pivots[lead] = value, coordinates
                selected.append(k)
                break
    endpoint = 0
    for value in coefficients:
        endpoint ^= value
    assert endpoint == sum(1 << (S * d) for S in range(d))
    target = 0
    while endpoint:
        value, coordinate = pivots[endpoint.bit_length() - 1]
        endpoint ^= value
        target ^= coordinate
    order = d * d
    rows = [[0] * order for _ in selected]
    columns = [[0] * order for _ in selected]
    entry_checks = 0
    for S, kernel in enumerate(kernels):
        for y in range(d):
            x = group.table[S][y]
            for a, polynomial in kernel.items():
                u, v = group.table[y][a], group.table[x][a]
                assert group.table[x][group.inverse[y]] == group.table[v][group.inverse[u]]
                for j, degree_index in enumerate(selected):
                    if polynomial >> degree_index & 1:
                        rows[j][x * d + y] ^= 1 << (v * d + u)
                        columns[j][x * d + u] ^= 1 << (v * d + y)
                entry_checks += 1
    for cut in (rows, columns):
        for i in range(order):
            value = 0
            for j, matrix in enumerate(cut):
                if target >> j & 1:
                    value ^= matrix[i]
            assert value == 1 << i
    name = f"mixed_radix_{group.name}_w{'_'.join(map(str, weights))}_seed{seed}"
    return {"name": name, "labels": d, "group": group.name, "radices": group.radices,
            "weights": weights, "control_word_permutation": controls,
            "basis_dimension": len(selected), "basis_coefficient_indices": selected,
            "endpoint_coordinates": [(target >> j) & 1 for j in range(len(selected))],
            "block_sizes": [order],
            "basis_row_blocks": [[[hex(row) for row in matrix]] for matrix in rows],
            "basis_column_blocks": [[[hex(row) for row in matrix]] for matrix in columns],
            "basis_ranks": [[rank_two(a), rank_two(b)] for a, b in zip(rows, columns)],
            "orbit_entry_checks": entry_checks,
            "group_table_sha256": hashlib.sha256(json.dumps(group.table).encode()).hexdigest(),
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, rows, columns


def binary_cost(data, rows, columns):
    dimension = len(rows)
    if dimension > 12:
        raise ValueError("use native screen for dimension above12")
    order = data["labels"] ** 2
    rwork, cwork = [0] * order, [0] * order
    ranks = [[0, 0] for _ in range(1 << dimension)]
    previous = 0
    for step in range(1, 1 << dimension):
        mask = step ^ (step >> 1)
        change = (mask ^ previous).bit_length() - 1
        for i in range(order):
            rwork[i] ^= rows[change][i]
            cwork[i] ^= columns[change][i]
        ranks[mask] = [rank_two(rwork), rank_two(cwork)]
        previous = mask
    target = sum(value << j for j, value in enumerate(data["endpoint_coordinates"]))
    cost, witness = minimum_decomposition([sqrt(r * c) for r, c in ranks], target=target)
    lower, _ = minimum_decomposition([isqrt(r * c) for r, c in ranks], target=target)
    return {"all_binary_ranks": ranks, "found_cost": cost, "witness": witness,
            "integer_shortest_lower": lower,
            "binary_exclusion_by_integer_comparison": lower ** 70 > data["labels"] ** 99}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", required=True)
    parser.add_argument("--weights")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--binary-search", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    group = Group(args.group)
    weights = list(map(int, args.weights.split(","))) if args.weights else [1] * len(group.radices)
    if len(weights) != len(group.radices) or min(weights) < 1:
        raise ValueError("one positive weight per generator is required")
    data, rows, columns = build(group, weights, args.seed)
    if args.binary_search:
        data["binary_search"] = binary_cost(data, rows, columns)
    args.output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in data.items()
                      if k not in ("basis_row_blocks", "basis_column_blocks", "binary_search")}, indent=2))
    if args.binary_search:
        print(json.dumps({k: v for k, v in data["binary_search"].items() if k != "all_binary_ranks"}))


if __name__ == "__main__":
    main()
