#!/usr/bin/env python3
"""Exact checks for the digit reduction from Schur kernels to a q=4 seed.

The output certifies finite identities and Smith profiles. It does not
establish a stronger factoring exponent.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb
from pathlib import Path

from schur_fermion import BinaryField, crossing_value_binary


ROOT = Path(__file__).resolve().parent.parent
SEED_KERNEL = [[0, 3, 3, 1], [3, 0, 4, 3],
               [3, 4, 0, 12], [1, 3, 12, 0]]


def polynomial_multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def polynomial_determinant(matrix):
    result = 0
    for order in permutations(range(len(matrix))):
        term = 1
        for i, j in enumerate(order):
            term = polynomial_multiply(term, matrix[i][j])
        result ^= term
    return result


def flip(state):
    return ((state & 1) << 1) | ((state & 2) >> 1)


def seed_block_entry(states):
    rows = [i for i, state in enumerate(states) if state & 1]
    columns = [i for i, state in enumerate(states) if state & 2]
    if len(rows) != len(columns):
        return 0
    return polynomial_determinant([[SEED_KERNEL[i][j] for j in columns] for i in rows])


def seed_data():
    local = []
    for x, y, u, v in product(range(4), repeat=4):
        value = seed_block_entry((flip(u), flip(x), v, y))
        assert value.bit_count() % 2 == int(x == v and y == u)
        local.append({"indices": [x, y, u, v], "matrix": [[value]]})
    return {"physical_alphabet": 4, "auxiliary_dimension": 1,
            "field": 2, "degree": 8,
            "label_encoding": "0=empty, 1=plus, 2=minus, 3=plus_minus",
            "kernel_polynomial_bits": SEED_KERNEL, "local": local}


def evaluate(bits, point, field):
    value = 0
    for i in range(bits.bit_length() - 1, -1, -1):
        value = field.multiply(value, point) ^ ((bits >> i) & 1)
    return value


def matrix_multiply(a, b, field):
    result = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(b[k]):
                    result[i][j] ^= field.multiply(value, other)
    return result


def scale(a, value, field):
    return [[field.multiply(x, value) for x in row] for row in a]


def transpose(a):
    return list(map(list, zip(*a)))


def kernel_data(m, z, field):
    q = z ^ 1
    ident = [[int(i == j) for j in range(m)] for i in range(m)]
    complement = [[int(i + j == m - 1) for j in range(m)] for i in range(m)]
    pascal = [[int(j & ~i == 0) for j in range(m)] for i in range(m)]
    disjoint = [[int(i & j == 0) for j in range(m)] for i in range(m)]
    disjoint_inverse = [[int(i | j == m - 1) for j in range(m)] for i in range(m)]
    zeta = [[field.power(z, i & ~j) for j in range(m)] for i in range(m)]
    right = [[field.power(q, -j - 1) if i & ~j == 0 else 0
              for j in range(m)] for i in range(m)]
    south = [[field.power(z, -1 - (i | j)) for j in range(m)] for i in range(m)]
    kernel = {(0, 1): zeta, (0, 2): pascal, (0, 3): ident,
              (1, 2): right, (1, 3): south, (2, 3): disjoint}
    kernel.update({(j, i): transpose(a) for (i, j), a in list(kernel.items())})
    return ident, complement, pascal, disjoint, disjoint_inverse, zeta, right, south, kernel


def check_gauges(field):
    checks = 0
    for m in (1, 2, 4, 8, 16):
        for w in (2, 3, 5):
            z = field.power(w, 2)
            q = z ^ 1
            ident, j, p, d, di, zz, r, ss, kernel = kernel_data(m, z, field)
            zm = field.power(z, m)
            small_w = field.power(w, m)
            root_qm = small_w ^ 1
            zi = scale(zz, field.power(q, 1 - m), field)
            u = [ident,
                 scale(matrix_multiply(zi, j, field), field.multiply(root_qm, zm), field),
                 scale(di, field.multiply(root_qm, zm), field),
                 scale(j, zm, field)]
            v = [scale(j, field.inverse[zm], field),
                 scale(zi, root_qm, field), scale(p, root_qm, field), ident]
            assert matrix_multiply(zz, zz, field) == scale(ident, field.power(q, m - 1), field)
            assert matrix_multiply(transpose(zz), d, field) == scale(r, field.power(q, m), field)
            assert matrix_multiply(ss, p, field) == scale(r, field.power(field.multiply(q, field.inverse[z]), m), field)
            checks += 3 * m * m
            for (a, b), block in kernel.items():
                transformed = matrix_multiply(matrix_multiply(u[a], block, field), v[b], field)
                expected = scale(ident, evaluate(SEED_KERNEL[a][b], small_w, field), field)
                assert transformed == expected, (m, w, a, b)
                checks += m * m
    return checks


def minor(matrix, rows, columns, field):
    return field.determinant([[matrix[i][j] for j in columns] for i in rows])


def lift(left, right, labels, field):
    # A left row change and right column change act on the two exterior powers.
    return [[field.multiply(minor(left, a, c, field), minor(right, d, b, field))
             for c, d in labels] for a, b in labels]


def apply_axis(tensor, matrix, axis, dimension, field):
    result = [0] * len(tensor)
    stride = dimension ** (3 - axis)
    block = stride * dimension
    sparse = [[(j, value) for j, value in enumerate(row) if value] for row in matrix]
    for start in range(0, len(tensor), block):
        for inside in range(stride):
            for i, row in enumerate(sparse):
                value = 0
                for j, entry in row:
                    value ^= field.multiply(entry, tensor[start + j * stride + inside])
                result[start + i * stride + inside] = value
    return result


def partition_from_frobenius(arms, legs):
    arms = sorted(arms, reverse=True)
    legs = sorted(legs, reverse=True)
    r = len(arms)
    parts = [a + i + 1 for i, a in enumerate(arms)]
    for row in range(r + 1, max((b + i + 1 for i, b in enumerate(legs)), default=0) + 1):
        parts.append(sum(row <= b + i + 1 for i, b in enumerate(legs)))
    assert parts == sorted(parts, reverse=True)
    return tuple(reversed(parts))


def check_projected_identity(m, r, w, field):
    # Full tensor verification, including every local exterior map.
    h = m * r
    z = field.power(w, 2)
    q = z ^ 1
    ident, j, p, d, _, zz, _, _, _ = kernel_data(m, z, field)
    subsets = tuple(combinations(range(m), r))
    labels = tuple(product(subsets, repeat=2))
    dimension = len(labels)
    maps = [lift(ident, j, labels, field),
            lift(matrix_multiply(j, zz, field), zz, labels, field),
            lift(d, p, labels, field), lift(j, ident, labels, field)]
    point = field.power(w, m)
    site = {states: evaluate(seed_block_entry(states), point, field)
            for states in product(range(4), repeat=4)}
    tensor = []
    for positions in product(range(dimension), repeat=4):
        value = 1
        for index in range(m):
            states = tuple(int(index in labels[k][0]) + 2 * int(index in labels[k][1])
                           for k in positions)
            value = field.multiply(value, site[states])
        tensor.append(value)
    for axis, matrix in enumerate(maps):
        tensor = apply_axis(tensor, matrix, axis, dimension, field)
    lookup = {label: i for i, label in enumerate(labels)}
    scalar = field.power(q, -r * (m - 1))
    count = 0
    physical = [(a, b, partition_from_frobenius(a, b)) for a, b in labels
                if sum(a) + sum(b) + r == h]
    for x, y, u, v in product(physical, repeat=4):
        positions = [lookup[(u[1], u[0])], lookup[(x[1], x[0])],
                     lookup[(v[0], v[1])], lookup[(y[0], y[1])]]
        index = 0
        for position in positions:
            index = index * dimension + position
        actual = field.multiply(tensor[index], scalar)
        expected = crossing_value_binary(x[2], y[2], u[2], v[2], z, field)
        assert actual == expected, (m, w, x, y, u, v)
        count += 1
    return count


def rank_binary(rows):
    basis = {}
    for row in rows:
        while row:
            lead = row.bit_length() - 1
            if lead not in basis:
                basis[lead] = row
                break
            row ^= basis[lead]
    return len(basis)


def cut_matrix(seed, cut):
    result = [[0] * 16 for _ in range(16)]
    for term in seed["local"]:
        x, y, u, v = term["indices"]
        row, column = (4 * x + y, 4 * v + u) if cut == 0 else (4 * x + u, 4 * v + y)
        result[row][column] = term["matrix"][0][0]
    return result


def reverse_matrix(matrix, degree):
    return [[sum(((p >> j) & 1) << (degree - j) for j in range(degree + 1))
             for p in row] for row in matrix]


def smith_orders_by_toeplitz(matrix, max_precision=20):
    n = len(matrix)
    degree = max(x.bit_length() - 1 for row in matrix for x in row)
    coefficients = [[sum(((matrix[i][j] >> power) & 1) << j for j in range(n))
                     for i in range(n)] for power in range(degree + 1)]
    previous_rank = previous_difference = 0
    orders = []
    records = []
    for precision in range(1, max_precision + 1):
        rows = []
        for a in range(precision):
            for i in range(n):
                rows.append(sum(coefficients[a - b][i] << (b * n)
                                for b in range(max(0, a - degree), a + 1)))
        rank = rank_binary(rows)
        difference = rank - previous_rank
        assert difference >= previous_difference
        orders += [precision - 1] * (difference - previous_difference)
        records.append({"precision": precision, "rank": rank})
        previous_rank, previous_difference = rank, difference
        if difference == n:
            break
    assert len(orders) == n
    return orders, records


def check_profiles(seed):
    expected_zero = [[0] * 4 + [4] * 8 + [8] * 4, [0] * 16]
    expected_degrees = [[8, 7, 7, 6, 5, 5, 4, 4, 4, 4, 3, 3, 2, 1, 1, 0],
                        [8, 6, 6, 4, 2, 2, 0, 0, 0, 0, -2, -2, -4, -6, -6, -8]]
    result = []
    for cut in range(2):
        matrix = cut_matrix(seed, cut)
        zero, zero_records = smith_orders_by_toeplitz(matrix)
        infinity, infinity_records = smith_orders_by_toeplitz(reverse_matrix(matrix, 8))
        degrees = [8 - order for order in infinity]
        assert zero == expected_zero[cut]
        assert degrees == expected_degrees[cut]
        result.append({"cut": cut, "zero_orders": zero, "infinity_degrees": degrees,
                       "zero_toeplitz_checks": zero_records,
                       "infinity_toeplitz_checks": infinity_records})
    return result


def check_dimensions():
    result = []
    for m in (2, 4, 8, 16):
        r = m // 2
        counts = {}
        for subset in combinations(range(m), r):
            counts[sum(subset)] = counts.get(sum(subset), 0) + 1
        dimension = sum(value * value for value in counts.values())
        full = comb(m, r) ** 2
        width = r * (m - r) + 1
        assert dimension * width >= full
        result.append({"M": m, "r": r, "h": m * r,
                       "occupation_dimension": full, "energy_dimension": dimension,
                       "possible_subset_sums": width})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed-output", type=Path)
    args = parser.parse_args()
    seed = seed_data()
    field = BinaryField(7)
    result = {
        "status": "exact Gaussian reduction and finite profiles; no stronger factoring bound",
        "seed_endpoint_entries": 256,
        "seed_nonzero_polynomial_entries": sum(bool(t["matrix"][0][0]) for t in seed["local"]),
        "one_particle_field_entry_checks": check_gauges(field),
        "field_modulus_bits": field.modulus,
        "projected_schur_identity_checks": sum(check_projected_identity(m, r, w, field)
                                               for m, r, w in ((2, 1, 2), (2, 1, 3),
                                                               (4, 1, 2), (4, 3, 2))),
        "non_prime_power_weight_check": {"M": 4, "r": 3, "h": 12},
        "dimension_checks": check_dimensions(),
        "smith_profiles": check_profiles(seed),
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("tools/schur_fermion.py",
                                       "tools/check_schur_gaussian_reduction.py")},
    }
    for path, data in ((args.output, result), (args.seed_output, seed)):
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    if args.output is None:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
