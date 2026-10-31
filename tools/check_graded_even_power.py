#!/usr/bin/env python3
"""Exact signs for graded powers of even four-boundary polynomial tensors.

No rank-pressure improvement is claimed. The negative controls separate
ordinary local flattenings from the cut order used for graded powers.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_all_field_pascal_crossings import add, multiply, scale
from check_signed_residual_codes import build_residual, flatten, prime_rank


# Original physical order is (V,Y,U,X). Each tuple gives contiguous cut order.
CUTS = ((3, 1, 0, 2), (3, 2, 0, 1))


def inversion_sign(original, reordered):
    """Koszul sign of an explicit permutation of occupied variable identifiers."""
    positions = {value: i for i, value in enumerate(reordered)}
    permutation = [positions[value] for value in original]
    return (-1) ** sum(a > b for i, a in enumerate(permutation)
                      for b in permutation[i + 1:])


def boundary_sign(labels, parity, order):
    occupied = [a for a in range(4) if parity[labels[a]]]
    return inversion_sign(occupied, [a for a in order if a in occupied])


def local_cut(tensor, parity, order, prime):
    d = len(parity)
    result = [[()] * (d * d) for _ in range(d * d)]
    for labels, value in tensor.items():
        row = d * labels[order[0]] + labels[order[1]]
        column = d * labels[order[2]] + labels[order[3]]
        result[row][column] = scale(value, boundary_sign(labels, parity, order), prime)
    return result


def literal_global_entry(tensor, words, parity, prime, order=(0, 1, 2, 3)):
    """Multiply copy tensors, then literally reorder occupied odd letters."""
    m = len(words[0])
    value = (1,)
    for j in range(m):
        value = multiply(value, tensor.get(tuple(word[j] for word in words), ()), prime)
    original = [(j, a) for j in range(m) for a in range(4)
                if parity[words[a][j]]]
    reordered = [(j, a) for a in order for j in range(m)
                 if parity[words[a][j]]]
    return scale(value, inversion_sign(original, reordered), prime)


def cut_factorized_entry(tensor, words, parity, prime, order):
    """Row/column diagonal formula for an even seed, in contiguous cut order."""
    m = len(words[0])
    value = (1,)
    left_parities = []
    row_sign = column_sign = 1
    for j in range(m):
        labels = tuple(word[j] for word in words)
        value = multiply(value, scale(tensor.get(labels, ()),
                                     boundary_sign(labels, parity, order), prime), prime)
        left_parities.append(sum(parity[labels[a]] for a in order[:2]) % 2)
    row_sign *= (-1) ** sum(a * b for i, a in enumerate(left_parities)
                           for b in left_parities[i + 1:])
    # Regroup copies inside each cut side. These signs depend on that side only.
    for side, destination in ((order[:2], "row"), (order[2:], "column")):
        original = [(j, a) for j in range(m) for a in side if parity[words[a][j]]]
        reordered = [(j, a) for a in side for j in range(m) if parity[words[a][j]]]
        sign = inversion_sign(original, reordered)
        if destination == "row":
            row_sign *= sign
        else:
            column_sign *= sign
    return scale(value, row_sign * column_sign, prime)


def sample_even_tensor(d, parity, prime):
    tensor = {}
    for labels in product(range(d), repeat=4):
        if sum(parity[a] for a in labels) % 2:
            continue
        number = sum((a + 1) * (i + 2) for i, a in enumerate(labels))
        tensor[labels] = tuple((number * (j + 1) + j * j + 1) % prime for j in range(3))
        while tensor[labels] and not tensor[labels][-1]:
            tensor[labels] = tensor[labels][:-1]
    return tensor


def check_tensor_power(prime, parity, copies):
    d = len(parity)
    tensor = sample_even_tensor(d, parity, prime)
    words = tuple(product(range(d), repeat=copies))
    checked = 0
    for order in CUTS:
        for labels in product(words, repeat=4):
            direct = literal_global_entry(tensor, labels, parity, prime, order)
            predicted = cut_factorized_entry(tensor, labels, parity, prime, order)
            assert direct == predicted
            checked += 1
    return checked


def evaluate(polynomial, point, prime):
    result = 0
    for coefficient in reversed(polynomial):
        result = (result * point + coefficient) % prime
    return result


def check_typed_endpoints(prime, parity, copies):
    d = len(parity)
    # Endpoint ordinary crossing, with arbitrary even perturbation vanishing at 1.
    perturbation = sample_even_tensor(d, parity, prime)
    tensor = {}
    for labels in product(range(d), repeat=4):
        v, y, u, x = labels
        base = (int(v == x and y == u),)
        tensor[labels] = add(base, multiply((prime - 1, 1),
                                           perturbation.get(labels, ()), prime), prime)
    words = tuple(product(range(d), repeat=copies))
    types = {}
    for word in words:
        kind = tuple(word.count(a) for a in range(d))
        types.setdefault(kind, []).append(word)
    checked = 0
    for word_set in types.values():
        for labels in product(word_set, repeat=4):
            v, y, u, x = labels
            value = evaluate(literal_global_entry(tensor, labels, parity, prime), 1, prime)
            assert value == int(v == x and y == u)
            # Global boundary reorder is constant on a fixed physical type.
            totals = [sum(parity[a] for a in word) % 2 for word in labels]
            for order in CUTS:
                sign = boundary_sign(tuple(totals), (0, 1), order)
                cut_value = literal_global_entry(tensor, labels, parity, prime, order)
                assert cut_value == scale(literal_global_entry(tensor, labels, parity, prime),
                                          sign, prime)
            checked += 1
    return checked


def odd_negative_control():
    # Its correctly cut-ordered local matrix is all ones, hence has rank one.
    # The seed has both odd and even entries, so the even-power lemma cannot apply.
    tensor = {(v, 0, 0, x): ((-1) ** (v * x) % 3,)
              for v, x in product(range(2), repeat=2)}
    words = tuple(product(range(2), repeat=2))
    direct = [[evaluate(literal_global_entry(tensor, (v, (0, 0), (0, 0), x),
                                            (0, 1), 3, CUTS[0]), 0, 3)
               for v in words] for x in words]
    assert prime_rank([[1] * 4 for _ in range(4)], 3) == 1
    assert prime_rank(direct, 3) == 2
    return {"ordinary_power_rank": 1, "graded_power_rank": 2, "field": 3}


def residual_sign_control():
    residual, _, _ = build_residual()
    ordinary = flatten(residual, 0)
    corrected = [row[:] for row in ordinary]
    for x, y, u, v in product(range(4), repeat=4):
        corrected[4 * x + y][4 * u + v] *= (-1) ** (v.bit_count() * y.bit_count())
    assert prime_rank(ordinary, 3) == 15
    assert prime_rank(corrected, 3) == 16
    return {"field": 3, "ordinary_residual_row_rank": 15,
            "contiguous_Grassmann_row_rank": 16,
            "essential_sign": "(-1)^(parity(v)*parity(y))"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for prime in (2, 3, 5):
        for parity, copies in (((0, 1), 2), ((0, 1), 3), ((0, 1, 1, 0), 2)):
            records.append({"prime": prime, "parity": parity, "copies": copies,
                            "polynomial_entry_equalities": check_tensor_power(prime, parity, copies),
                            "typed_endpoint_entries": check_typed_endpoints(prime, parity, copies)})
    result = {"scope": "Exact finite controls for graded even tensor powers and their actual cut matrices.",
              "cases": records, "odd_seed_negative_control": odd_negative_control(),
              "signed_residual_control": residual_sign_control(),
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
