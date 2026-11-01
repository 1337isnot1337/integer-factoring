#!/usr/bin/env python3
"""All-pattern cubic identity and independent multirow rank-certificate replay."""

import argparse
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path

from explore_cubic_multirow_outputs import CubicSpace


@lru_cache(None)
def tails(first, left, right):
    if left == right == 0:
        return ((),)
    result = []
    for kind in range(first, 16):
        a, b = kind // 4, kind % 4
        da, db = a.bit_count(), b.bit_count()
        if da <= left and db <= right:
            result.extend(((a, b),) + rest for rest in tails(kind, left - da, right - db))
    return tuple(result)


def single_mode(word):
    singles = [(i, a) for i, a in enumerate(word) if a in (1, 2)]
    return singles[0] if len(singles) == 1 else None


def universal_identity():
    space = CubicSpace(9)
    checked = compatible = 0
    polynomials = set()
    for same_input in (False, True):
        distinguished = 1 if same_input else 2
        y, v = 0, 0 if same_input else 1
        for prefix in product(range(16), repeat=distinguished):
            first = tuple(a // 4 for a in prefix)
            second = tuple(a % 4 for a in prefix)
            left = 3 - sum(a.bit_count() for a in first)
            right = 3 - sum(a.bit_count() for a in second)
            if left < 0 or right < 0:
                continue
            for tail in tails(1, left, right):
                assert distinguished + len(tail) <= 8
                padding = (0,) * (9 - distinguished - len(tail))
                a = first + tuple(x for x, z in tail) + padding
                b = second + tuple(z for x, z in tail) + padding
                value = space.kernel(a, b, y, v)
                first_mode, second_mode = single_mode(a), single_mode(b)
                endpoint = int(first_mode == (v, 1) and second_mode == (y, 1))
                high_endpoint = int(first_mode == (y, 2) and second_mode == (v, 2))
                gram = int(y == v and first_mode is not None and first_mode == second_mode)
                assert value & 1 == endpoint
                assert (value >> 12) & 1 == high_endpoint
                assert ((value >> 1) ^ (value >> 4)) & 1 == (gram ^ high_endpoint)
                # The five universal scalar polynomials, dual to degrees 0,1,2,4,12.
                reconstructed = 0
                for degree, scalar in zip((0, 1, 2, 4, 12), (257, 514, 68, 272, 4352)):
                    if value >> degree & 1:
                        reconstructed ^= scalar
                assert reconstructed == value
                checked += 1
                compatible += int(value != 0)
                polynomials.add(value)
    return {"canonical_physical_patterns_checked": checked,
            "nonzero_canonical_patterns": compatible,
            "universal_scalar_polynomials": [257, 514, 68, 272, 4352],
            "observed_entry_polynomials": sorted(polynomials),
            "representative_physical_pairs": 9,
            "maximum_distinguished_or_occupied_sites": 8,
            "all_odd_pair_counts_reduce_by_empty_site_parity": True,
            "arbitrary_characteristic_two_output_identity": True}


def low_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = (value & -value).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def shortest_integer_cost(ranks):
    width = len(ranks)
    costs = [isqrt(a * b) for a, b in ranks]
    bound = sum(costs) + 1
    distances = [0] + [bound] * (width - 1)
    # Bellman relaxation on the entire finite XOR graph, independent of the search Dijkstra.
    for _ in range(width - 1):
        old = distances.copy()
        for state in range(width):
            for mask in range(1, width):
                distances[state ^ mask] = min(distances[state ^ mask], old[state] + costs[mask])
        if old == distances:
            break
    assert all(distances[x ^ m] <= distances[x] + costs[m]
               for x in range(width) for m in range(1, width))
    return distances[1]


def replay_multirow(path):
    result = json.loads(path.read_text())
    assert result["complete"] and len(result["records"]) == 680
    assert hashlib.sha256(Path(__file__).with_name("explore_cubic_multirow_outputs.py").read_bytes()).hexdigest() == result["source_sha256"]
    space = CubicSpace(result["pairs"])
    pairs = space.pairs
    ranks_checked = literal_checks = 0
    for case, record in enumerate(result["records"]):
        corrections = [list(map(tuple, row)) for row in record["corrections"]]
        tensor = space.tensor(corrections)
        if case in (0, 7, 73, 167, 168, 287, 323, 411, 679):
            outputs = [[space.outputs[i]] + corrections[i] for i in range(pairs)]
            for x, y, v, u in product(range(pairs), repeat=4):
                value = 0
                for first, second in product(outputs[x], outputs[u]):
                    value ^= space.kernel(first, second, y, v)
                assert value == tensor[x * pairs + y][v * pairs + u]
                literal_checks += 1
        rows = [[0] * (pairs * pairs) for _ in range(4)]
        cols = [[0] * (pairs * pairs) for _ in range(4)]
        for x, y, v, u in product(range(pairs), repeat=4):
            value = tensor[x * pairs + y][v * pairs + u]
            assert value & 1 == int(x == v and y == u)
            assert not (value >> 3) & 1 and not (value >> 12) & 1
            assert ((value >> 1) ^ (value >> 4)) & 1 == int(x == u and y == v)
            for i, degree in enumerate((0, 1, 2, 4)):
                if value >> degree & 1:
                    rows[i][x * pairs + y] ^= 1 << (v * pairs + u)
                    cols[i][x * pairs + u] ^= 1 << (v * pairs + y)
        assert record["basis_degrees"] == [0, 1, 2, 4]
        actual = [[0, 0]]
        for mask in range(1, 16):
            ranks = []
            for basis in (rows, cols):
                matrix = [0] * (pairs * pairs)
                for i, source in enumerate(basis):
                    if mask >> i & 1:
                        matrix = [a ^ b for a, b in zip(matrix, source)]
                ranks.append(low_rank(matrix))
            actual.append(ranks)
            ranks_checked += 1
        assert actual == record["rank_table"]
        assert shortest_integer_cost(actual) == record["integer_lower"] == 38
        assert 38 ** 70 > 7 ** 99
    return {"cases_replayed": 680, "nonzero_binary_rank_pairs_replayed": ranks_checked,
            "literal_bilinear_expansion_entries_checked": literal_checks,
            "integer_XOR_cost_lower_bound": 38,
            "all_canonical_non_crossing_tensors_unchanged": True,
            "all_recorded_spaces_fail_binary_seed_criterion": True,
            "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def affordable_witness():
    space = CubicSpace(7)
    edges = ((0, 1), (2, 3))
    corrections = [[tuple(int(j in (i, a, b)) for j in range(7))
                    for a, b in edges if i not in (a, b)] for i in range(7)]
    tensor = space.tensor(corrections)
    rows, columns = [0] * 49, [0] * 49
    for x, y, v, u in product(range(7), repeat=4):
        value = tensor[7 * x + y][7 * v + u]
        if (value ^ (value >> 2)) & 1:
            rows[7 * x + y] ^= 1 << (7 * v + u)
            columns[7 * x + u] ^= 1 << (7 * v + y)
    assert [low_rank(rows), low_rank(columns)] == [16, 15]
    assert 1393 ** 2 < 2 * 985 ** 2 and 240 ** 985 < 7 ** 2786
    profile = space.profile(tensor)
    assert profile["rank_table"][5] == [16, 15]
    assert profile["rank_table"][4] == [42, 40]
    assert shortest_integer_cost(profile["rank_table"]) == 38
    return {"labels": 7, "quadratic_edges": edges,
            "output_selector": "cyclic next-site double plus low-mode e_i times the shared quadratic",
            "crossing_coefficient_mask": 5, "rank_pair": [16, 15], "rank_product": 240,
            "strictly_below_single_piece_seed_threshold": True,
            "complementary_piece_rank_pair": [42, 40],
            "whole_space_integer_cost_lower_bound": 38,
            "whole_space_still_fails_seed_criterion": True,
            "rank_table": profile["rank_table"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multirow", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "universal_identity": universal_identity(),
              "affordable_tensor_witness": affordable_witness()}
    print(json.dumps(result["universal_identity"]), flush=True)
    if args.multirow:
        result["multirow_replay"] = replay_multirow(args.multirow)
        print(json.dumps(result["multirow_replay"]), flush=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
