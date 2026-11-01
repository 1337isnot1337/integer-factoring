#!/usr/bin/env python3
"""Replay the product-16 relabeling witness and constant four-piece blocks."""

import argparse
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import random

from check_hamming_pfaffian import independent_ranks
from check_modular_crossing import rank_two
from check_modular_recombination import (
    crossing,
    decode,
    encode,
    insert,
    permute,
    ranks,
    represent,
)
from hamming_relabeling import action, block_basis, normalized, seed_data


def kron(left, right):
    result, size = [], len(right)
    for row in left:
        bits = [j for j in range(row.bit_length()) if row >> j & 1]
        result.extend(sum(v << (j * size) for j in bits) for v in right)
    return result


def xor(left, right):
    return [a ^ b for a, b in zip(left, right)]


def block_generators(kind):
    if kind == "scalar":
        return [[0], [1], [0], [0]]
    if kind == "special":
        return [[2, 0], [1, 2], [0, 0], [2, 0]]
    identity, x, xy = [1, 2, 4, 8], [4, 8, 0, 0], [8, 0, 0, 0]
    return [[0] * 4, xor(identity, x), xor(x, xy), xy]


def block_checks(generators, basis):
    lookup = {}
    for j, vector in enumerate(basis):
        assert insert(lookup, vector, 1 << j)
    checked = 0
    kinds = ("special", "jordan", "jordan", "jordan", "scalar", "scalar")
    for column in (False, True):
        for index in range(4):
            rows = generators[(1 - index) if column and index < 2 else index]
            columns = [represent(lookup, action(rows, vector)) for vector in basis]
            actual = [
                sum(((columns[j] >> i) & 1) << j for j in range(16)) for i in range(16)
            ]
            expected, start = [], 0
            letter = (1 - index) if column and index < 2 else index
            for kind in kinds:
                block = block_generators(kind)[letter]
                expected.extend(v << start for v in block)
                start += len(block)
            assert actual == expected
            checked += 16 * 16
    return checked


def type_matrices(left, right):
    return [
        (
            kron(left[i], right[j])
            if i == j
            else xor(kron(left[i], right[j]), kron(left[j], right[i]))
        )
        for i, j in combinations_with_replacement(range(4), 2)
    ]


def combination(matrices, mask):
    rows = [0] * len(matrices[0])
    for j, matrix in enumerate(matrices):
        if mask >> j & 1:
            rows = xor(rows, matrix)
    return rows


def tensor_checks(generators):
    rng = random.Random(160025)
    masks = [0, 1023] + [1 << j for j in range(10)]
    masks += [rng.randrange(1024) for _ in range(40)]
    checked = 0
    for column in (False, True):
        g = generators if not column else [generators[i] for i in (1, 0, 2, 3)]
        direct = type_matrices(g, g)
        data = []
        kinds = [("special", 1), ("jordan", 3), ("scalar", 2)]
        for i, (left, left_count) in enumerate(kinds):
            for j in range(i, len(kinds)):
                right, right_count = kinds[j]
                a, b = block_generators(left), block_generators(right)
                if column:
                    a = [a[k] for k in (1, 0, 2, 3)]
                    b = [b[k] for k in (1, 0, 2, 3)]
                count = left_count * right_count * (1 if i == j else 2)
                data.append((type_matrices(a, b), count))
        for mask in masks:
            expected = sum(
                rank_two(combination(matrices, mask)) * count
                for matrices, count in data
            )
            assert rank_two(combination(direct, mask)) == expected
            checked += 1
        assert rank_two(combination(direct, 1023)) == 256
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    original, images, pieces, masks = seed_data()
    total, witness = 0, []
    for vector, mask in zip(pieces, masks):
        replay = 0
        terms = []
        for j, (image, (permutation, source_index)) in enumerate(images):
            if mask >> j & 1:
                replay ^= image
                terms.append(
                    {
                        "label_permutation": list(permutation),
                        "joint_source_basis_index": source_index,
                    }
                )
        assert vector == replay
        actual = ranks(decode(vector, 4), 4)
        assert actual == independent_ranks(vector, 4)
        total ^= vector
        witness.append(
            {
                "row_rank": actual[0],
                "column_rank": actual[1],
                "tensor_hex": hex(vector),
                "source_terms": terms,
            }
        )
    assert total == encode(crossing(4), 4)
    assert [(p["row_rank"], p["column_rank"]) for p in witness] == [
        (1, 16),
        (16, 1),
    ] + [(4, 4)] * 4
    r, c, x = pieces[0], pieces[1], pieces[5]
    e = encode(crossing(4), 4)
    distinct = e ^ r ^ c ^ x
    direct = [0] * 16
    from itertools import product

    for a, b, u, v in product(range(4), repeat=4):
        direct[4 * a + b] |= int(len({a, b, u, v}) == 4) << (4 * u + v)
    assert encode(direct, 4) == distinct == pieces[2] ^ pieces[3] ^ pieces[4]
    generators = [normalized(v) for v in (r, c, distinct, x)]
    basis = block_basis()
    local = block_checks(generators, basis)
    tensors = tensor_checks(generators)
    # Independently check that the two cuts exchange only R,C in the aggregate.
    for i, vector in enumerate((r, c, distinct, x)):
        other = normalized(encode(permute(decode(vector, 4), 4, (0, 2, 1, 3)), 4))
        assert other == generators[1 - i if i < 2 else i]
    paths = [Path(__file__), Path(__file__).with_name("hamming_relabeling.py")]
    result = {
        "status": "finite product-16 witness and all-size block model; no improved factoring exponent",
        "binary_width": 2,
        "joint_source_dimension": len(original),
        "label_orbit_images": len(images),
        "label_orbit_span_dimension": 102,
        "maximum_rank_product": 16,
        "maximum_binary_width_sum": 4,
        "witness_piece_count": len(witness),
        "crossing_entries": 256,
        "constant_basis": basis,
        "local_block_entries_checked": local,
        "independent_tensor_power_rank_checks": tensors,
        "witness": witness,
        "source_sha256": {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "witness"}))


if __name__ == "__main__":
    main()
