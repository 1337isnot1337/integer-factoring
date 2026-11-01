#!/usr/bin/env python3
"""Exact checks for finite-seed amplification and the PR24 common flags."""

import argparse
from collections import Counter
from functools import reduce
import json
from math import log, sqrt
from pathlib import Path


SOURCE_COMMIT = "e96ee50560284301428d72421330d8773155d25b"
# PR24 modular-higher-witnesses.json, first case, in witness order.
# Rows are (x,y); columns are (u,v), each pair in lexicographic order.
PIECES = (
    (273,273,0,273,273,0,0,0,0),
    (9,18,36,9,18,36,0,0,0),
    (3,3,0,24,24,0,192,192,0),
    (27,0,0,0,27,0,0,0,27),
    (0,0,0,0,0,0,24,3,0),
    (32,4,0,32,4,0,0,0,0),
    (192,192,27,192,192,27,0,0,0),
    (195,195,0,216,216,0,195,216,0),
    (36,36,219,36,36,219,0,0,0),
    (0,0,219,0,0,219,63,63,0),
    (36,36,256,36,36,256,292,292,0),
    (0,0,0,0,0,0,63,63,228),
    (465,465,465,465,465,465,465,465,465),
    (10,10,0,10,10,0,0,0,0),
    (24,24,0,3,3,0,0,0,0),
    (238,0,238,0,238,238,238,238,0),
    (1,2,0,8,16,0,0,0,0),
    (0,192,64,192,0,128,0,0,0),
    (46,46,0,46,46,0,0,0,46),
    (32,32,0,4,4,0,4,32,0),
)
FLAGS = {
    "row": (27,36,192,10,256,3,4,64,1),
    "column": (27,192,36,10,3,256,64,4,1),
}


def rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def crossing(d):
    return tuple(1 << (y*d+x) for x in range(d) for y in range(d))


def column(rows, d):
    answer = [0] * (d*d)
    for x in range(d):
        for y in range(d):
            bits = rows[x*d+y]
            while bits:
                bit = bits & -bits
                bits ^= bit
                u, v = divmod(bit.bit_length()-1, d)
                answer[x*d+u] |= 1 << (y*d+v)
    return tuple(answer)


def add(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def multiply(left, right):
    return tuple(reduce(int.__xor__,
                        (right[k] for k in range(len(right)) if (row >> k) & 1),
                        0) for row in left)


def apply(matrix, vector):
    return sum(((row & vector).bit_count() % 2) << j
               for j, row in enumerate(matrix))


def conjugate(matrix, basis):
    assert rank(basis) == len(basis)
    coordinates = {}
    for j, vector in enumerate(basis):
        mask = 1 << j
        while vector:
            pivot = vector.bit_length()-1
            if pivot not in coordinates:
                coordinates[pivot] = (vector, mask)
                break
            other, change = coordinates[pivot]
            vector ^= other
            mask ^= change
    columns = []
    for vector in basis:
        value, result = apply(matrix, vector), 0
        while value:
            other, change = coordinates[value.bit_length()-1]
            value ^= other
            result ^= change
        columns.append(result)
    return tuple(sum(((value >> i) & 1) << j
                     for j, value in enumerate(columns))
                 for i in range(len(basis)))


def algebra_dimension(generators):
    n = len(generators[0])
    pivots, matrices = {}, []

    def insert(matrix):
        value = sum(row << (n*i) for i, row in enumerate(matrix))
        while value:
            pivot = value.bit_length()-1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                matrices.append(matrix)
                return

    insert(tuple(1 << i for i in range(n)))
    for matrix in generators:
        insert(matrix)
    for matrix in matrices:
        for generator in generators:
            insert(multiply(matrix, generator))
    return len(pivots)


def tensor_two(left, right, d=3):
    """Tensor product with the four word-label coordinates grouped correctly."""
    D = d*d
    result = [0] * (D*D)
    for x1 in range(d):
        for y1 in range(d):
            for x2 in range(d):
                for y2 in range(d):
                    output_row = (x1*d+x2)*D + y1*d+y2
                    a = left[x1*d+y1]
                    while a:
                        bit = a & -a
                        a ^= bit
                        u1, v1 = divmod(bit.bit_length()-1, d)
                        b = right[x2*d+y2]
                        while b:
                            bit = b & -b
                            b ^= bit
                            u2, v2 = divmod(bit.bit_length()-1, d)
                            result[output_row] ^= 1 << ((u1*d+u2)*D+v1*d+v2)
    return tuple(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert reduce(add, PIECES, (0,)*9) == crossing(3)
    pairs = [(rank(A), rank(column(A, 3))) for A in PIECES]
    assert Counter(r*c for r, c in pairs) == {3: 5, 2: 2, 4: 13}
    flags = {}
    for side, basis in FLAGS.items():
        raw = PIECES if side == "row" else tuple(column(A, 3) for A in PIECES)
        operators = tuple(multiply(A, crossing(3)) for A in raw)
        changed = tuple(conjugate(A, basis) for A in operators)
        assert all(not (row & ((1 << i)-1))
                   for A in changed for i, row in enumerate(A))
        characters = [sum(((A[i] >> i) & 1) << j
                          for j, A in enumerate(changed)) for i in range(9)]
        expected = {1 << (16 if side == "row" else 15): 4,
                    1 << 17: 2, 1 << 19: 2, 1 << 12: 1}
        assert Counter(characters) == expected
        patterns = []
        for j, vector in enumerate(basis):
            support_patterns = {
                (x == 2, y == 2)
                for x in range(3) for y in range(3)
                if (vector >> (3*x+y)) & 1
            }
            assert len(support_patterns) == 1
            pattern = support_patterns.pop()
            patterns.append(pattern)
            common = characters[j] in (1 << 12, 1 << 17, 1 << 19)
            assert common == any(pattern)
        dimension = algebra_dimension(operators)
        assert dimension == (30 if side == "row" else 31)
        flags[side] = {"basis": basis, "characters": characters,
                       "individual_w_indicators": patterns,
                       "algebra_dimension": dimension,
                       "triangular_operators": changed}
    total = (0,)*81
    profile = []
    for i, A in enumerate(PIECES):
        for j in range(i+1):
            B = tensor_two(A, PIECES[j])
            if i != j:
                B = add(B, tensor_two(PIECES[j], A))
            total = add(total, B)
            r, c = rank(B), rank(column(B, 9))
            multiplicity = 1 if i == j else 2
            rb = multiplicity*pairs[i][0]*pairs[j][0]
            cb = multiplicity*pairs[i][1]*pairs[j][1]
            assert r <= rb and c <= cb
            profile.append({"indices": (j, i), "row_rank": r,
                            "column_rank": c, "row_bound": rb,
                            "column_bound": cb})
    assert total == crossing(9)
    assert add(total, tensor_two(PIECES[0], PIECES[0])) != crossing(9)
    for d in range(2, 6):
        assert (4*(d*d-1))**70 > d**198
    assert 99**2 > 2*70**2
    assert (4*(6*6-1))**5 < 6**14  # 2sqrt(35)<6^(7/5)<6^sqrt(2).
    for d in (2, 3, 5, 6):
        for length in range(1, 65):
            size, copies = 1, 0
            while size < (1 << length):
                size *= d
                copies += 1
            assert (1 << length) <= size < d*(1 << length)
    N = 9
    # b(3)>71/100>1/sqrt(2), using only integer inequalities.
    assert (2*N)**(2*N*100) > (N+1)**((N+1)*100)*3**(4*71*N)
    assert 2*71**2 > 100**2
    result = {
        "status": "exact finite amplification checks; no improved factoring exponent",
        "seed_source_commit": SOURCE_COMMIT,
        "seed_dimension": 3, "seed_piece_count": len(PIECES),
        "seed_rank_pairs": pairs,
        "rank_sum_expression": "26+5*sqrt(3)+2*sqrt(2)",
        "rank_sum_decimal": 26+5*sqrt(3)+2*sqrt(2),
        "small_seed_exclusions": [2, 3, 4, 5],
        "smith_three_label_lower_bound": (
            2*log(2*N)-(N+1)/N*log(N+1))/(4*log(3)),
        "unrestricted_type_rank_lower_coefficient": log(5)/(2*log(3)),
        "flags": flags,
        "two_copy_piece_count": len(profile),
        "two_copy_maximum_rank_product": max(
            x["row_rank"]*x["column_rank"] for x in profile),
        "two_copy_profile": profile,
        "negative_control": "omitting one nonzero tensor-square changes crossing",
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("flags", "two_copy_profile")}, indent=2))


if __name__ == "__main__":
    main()
