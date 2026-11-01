#!/usr/bin/env python3
"""Exact cyclic orbit-sum crossing and alignment checks for the PR24 seed."""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from check_finite_seed_amplification import PIECES, SOURCE_COMMIT, column, crossing, rank


def orbit_basis(copies):
    assert copies > 0 and copies % 2
    words = tuple(product(range(3), repeat=copies))
    representatives, sizes, index = [], [], {}
    for word in words:
        if word in index:
            continue
        rotations = {word[j:]+word[:j] for j in range(copies)}
        assert len(rotations) % 2
        number = len(representatives)
        representatives.append(min(rotations))
        sizes.append(len(rotations))
        for rotation in rotations:
            index[rotation] = number
    return words, tuple(representatives), tuple(sizes), index


def coefficient_profile(copies):
    words, representatives, sizes, orbit_index = orbit_basis(copies)
    dimension = len(representatives)
    entries = {
        q: tuple(i for i, A in enumerate(PIECES)
                 if (A[3*q[0]+q[1]] >> (3*q[2]+q[3])) & 1)
        for q in product(range(3), repeat=4)
    }

    @lru_cache(None)
    def polynomial(key):
        terms = set()
        for indices in product(*key):
            term = tuple(sorted(indices))
            if term in terms:
                terms.remove(term)
            else:
                terms.add(term)
        return tuple(terms)

    matrices = {}
    # Simultaneous cyclic invariance fixes the first representative. The
    # omitted first-orbit multiplicity is odd, and is therefore one in F_2.
    for xi, x in enumerate(representatives):
        for y in words:
            yi = orbit_index[y]
            row = xi*dimension+yi
            for u in words:
                ui = orbit_index[u]
                for v in words:
                    vi = orbit_index[v]
                    key = tuple(sorted(entries[q] for q in zip(x, y, u, v)))
                    if not all(key):
                        continue
                    bit = 1 << (ui*dimension+vi)
                    for term in polynomial(key):
                        matrices.setdefault(term, [0]*(dimension*dimension))[row] ^= bit

    total = [0]*(dimension*dimension)
    profile = []
    for term, A in sorted(matrices.items()):
        if not any(A):
            continue
        total = [a ^ b for a, b in zip(A, total)]
        r, c = rank(A), rank(column(A, dimension))
        profile.append({"indices": term, "row_rank": r, "column_rank": c})
    assert tuple(total) == crossing(dimension)
    maximum = max(item["row_rank"]*item["column_rank"] for item in profile)
    return {"copies": copies, "labels": dimension,
            "orbit_sizes": dict(Counter(sizes)), "nonzero_types": len(profile),
            "maximum_rank_product": maximum, "profile": profile}


def rotate(bits, shift, copies):
    full = (1 << copies)-1
    return ((bits << shift) | (bits >> (copies-shift))) & full


def good_partition(copies, common, first, second, third):
    q, r = divmod(copies, 3)
    values = (common, first, second, third)
    assert [x.bit_count() for x in values] == [r, q, q, q]
    assert sum(values) == (1 << copies)-1
    assert all(not (values[i] & values[j]) for i in range(4) for j in range(i))
    x, y, v = common | first, common | third, common | second
    assert (x & y) == (v & y) == common
    for shift in range(1, copies):
        shifted = rotate(y, shift, copies)
        if (x & shifted).bit_count() <= r or (v & shifted).bit_count() <= r:
            return False
    return True


def alignment_check(copies, common, first, second, third):
    assert good_partition(copies, common, first, second, third)
    _, r = divmod(copies, 3)
    full = (1 << copies)-1
    x, y0, v0 = common | first, common | third, common | second
    y_rotations = [rotate(y0, s, copies) for s in range(copies)]
    v_rotations = [rotate(v0, s, copies) for s in range(copies)]
    found = []
    # Fix the first leg's rotation. All t remaining common rotations have
    # identical coefficients and contribute one because t is odd.
    for ys, y in enumerate(y_rotations):
        for us, u in enumerate(y_rotations):
            if y != u:
                continue
            for vs, v in enumerate(v_rotations):
                intersection = x & y
                if (intersection.bit_count() == r
                        and (x & v) == intersection and (y & v) == intersection
                        and (x | y | v) == full):
                    found.append((ys, us, vs))
    assert found == [(0, 0, 0)]
    return {"copies": copies, "partition": (common, first, second, third),
            "valid_alignments_with_first_leg_fixed": found,
            "full_valid_alignments": copies}


def verify():
    local_checks = 0
    for x, y, u, v in product(range(3), repeat=4):
        swap = lambda a, b: a < 2 and b < 2 and a != b
        expected = {
            4: x == 2 and v < 2 and swap(y, u),
            5: x < 2 and v == 2 and swap(y, u),
            17: (u == 2) and ((y == 2 and x == v and x < 2)
                              ^ (swap(x, y) and v < 2)),
        }
        for piece, value in expected.items():
            assert ((PIECES[piece][3*x+y] >> (3*u+v)) & 1) == value
            local_checks += 1
    assert (PIECES[12][8] >> 8) & 1

    small = coefficient_profile(1)
    assert small["labels"] == 3 and small["maximum_rank_product"] == 4
    profile = coefficient_profile(3)
    assert (profile["labels"], profile["nonzero_types"],
            profile["maximum_rank_product"]) == (11, 1522, 210)
    alignments = [alignment_check(*case) for case in (
        (27, 0, 68722889, 4344118, 61150720),
        (31, 32, 286369042, 1658864264, 202250309),
        (39, 0, 215565227025, 281572149450, 52618437412),
        (41, 264, 1327688390660, 479372912, 870855491715),
    )]

    # Exact eventual counting threshold: epsilon_t < 1/2 for q >= 512.
    assert 4*1538*1539**4*8**511 < 9**511
    assert 8*1539**5 < 9*1536**5
    return {"source_commit": SOURCE_COMMIT,
            "local_support_checks": local_checks,
            "single_copy": small, "three_copy": profile,
            "alignment_examples": alignments,
            "counting_threshold": {
                "q_min": 512,
                "base_certificate": "4*1538*1539^4*8^511 < 9^511",
                "step_certificate": "8*1539^5 < 9*1536^5",
            }}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass",
                      "local_support_checks": result["local_support_checks"],
                      "three_copy_labels": result["three_copy"]["labels"],
                      "three_copy_types": result["three_copy"]["nonzero_types"],
                      "three_copy_maximum": result["three_copy"]["maximum_rank_product"],
                      "alignment_examples": len(result["alignment_examples"])}))


if __name__ == "__main__":
    main()
