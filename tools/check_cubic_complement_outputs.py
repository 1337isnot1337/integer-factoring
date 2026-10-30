#!/usr/bin/env python3
"""Independent audits of binary and F4 complementary-output searches."""

import argparse
import hashlib
import heapq
from itertools import product
import json
from math import isqrt
from pathlib import Path
import random

from check_cubic_output_identity import low_rank, shortest_integer_cost
from check_pfaffian_span_f4 import independent_rank, scale_vector, span_basis, remainder
from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_fields import matrices, MULTIPLY


def literal_comparison(family, complement):
    d = family.pairs
    low = family.space.outputs
    high = [tuple(2 if a == 1 else a for a in word) for word in low]
    outputs = [[(low[i], 1)] + [(high[j], value) for j, value in enumerate(complement[i]) if value]
               for i in range(d)]
    rows, columns = matrices(family, complement)
    checked = 0
    for x, y, v, u in product(range(d), repeat=4):
        polynomial = 0
        for (first, a), (second, b) in product(outputs[x], outputs[u]):
            scalar = MULTIPLY[a][b]
            binary = family.space.kernel(first, second, y, v)
            spread = sum(((binary >> degree) & 1) << (2 * degree) for degree in range(13))
            polynomial ^= spread * scalar
        for k, degree in enumerate(family.degrees):
            value = (polynomial >> (2 * degree)) & 3
            assert ((rows[k][d * x + y] >> (2 * (d * v + u))) & 3) == value
            assert ((columns[k][d * x + u] >> (2 * (d * v + y))) & 3) == value
            checked += 1
    return checked


def audit_binary(path, family):
    data = json.loads(path.read_text())
    assert data["complete"] and len(data["cases"]) == 1195
    assert hashlib.sha256(Path(__file__).with_name("explore_cubic_complement_outputs.py").read_bytes()).hexdigest() == data["source_sha256"]
    checked = literal = 0
    for case, record in enumerate(data["cases"]):
        raw, raw_column = family.matrices(record["complement_rows"])
        indices = [family.degrees.index(degree) for degree in record["basis_degrees"]]
        rows, columns = [[cut[i] for i in indices] for cut in (raw, raw_column)]
        ranks = [[0, 0]]
        for mask in range(1, 1 << len(indices)):
            pair = []
            for cut in (rows, columns):
                matrix = [0] * 49
                for i, source in enumerate(cut):
                    if mask >> i & 1:
                        matrix = [a ^ b for a, b in zip(matrix, source)]
                pair.append(low_rank(matrix))
            ranks.append(pair)
            checked += 1
        assert ranks == record["rank_table"]
        lower = shortest_integer_cost(ranks)
        assert lower == record["integer_shortest_lower_bound"]
        assert lower ** 70 > 7 ** 99
        if case in (0, 1, 117, 555, 1194):
            complement = [[(row >> j) & 1 for j in range(7)] for row in record["complement_rows"]]
            literal += literal_comparison(family, complement)
    return {"cases": len(data["cases"]), "independent_binary_rank_pairs": checked,
            "direct_physical_coefficient_entries_checked": literal,
            "minimum_integer_cost_lower_bound": min(r["integer_shortest_lower_bound"] for r in data["cases"]),
            "all_recorded_binary_spaces_fail": True,
            "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def integer_distances(ranks):
    edges = [(isqrt(a * b), i) for i, (a, b) in enumerate(ranks) if i]
    distances = [10 ** 6] * len(ranks)
    distances[0] = 0
    queue = [(0, 0)]
    while queue:
        distance, state = heapq.heappop(queue)
        if distance != distances[state]:
            continue
        for cost, vector in edges:
            candidate, successor = distance + cost, state ^ vector
            if candidate < distances[successor]:
                distances[successor] = candidate
                heapq.heappush(queue, (candidate, successor))
    assert all(distances[x ^ v] <= distances[x] + cost
               for x in range(len(ranks)) for cost, v in edges)
    return distances


def audit_field4(path, family):
    data = json.loads(path.read_text())
    assert data["complete"] and len(data["cases"]) == 56
    assert hashlib.sha256(Path(__file__).with_name("explore_cubic_complement_fields.py").read_bytes()).hexdigest() == data["source_sha256"]
    projective = checked = literal = 0
    for case, record in enumerate(data["cases"]):
        raw, raw_columns = matrices(family, record["complement_matrix"])
        indices = [family.degrees.index(degree) for degree in record["basis_degrees"]]
        rows, columns = [[cut[i] for i in indices] for cut in (raw, raw_columns)]
        dimension = len(indices)
        assert len(record["rank_table"]) == 4 ** dimension
        count = (4 ** dimension - 1) // 3
        assert record["projective_points_checked"] == count
        projective += count
        affordable = [code for code, (a, b) in enumerate(record["rank_table"])
                      if code and a * b < 256]
        residual = remainder(1, span_basis(affordable, dimension), dimension)
        assert residual, "The affordable span needs a cost audit."
        # Check every affordable projective point, the displayed witness, and controls.
        selected = set(record["witness"])
        for code in affordable:
            pivot = (code.bit_length() - 1) // 2
            if (code >> (2 * pivot)) & 3 == 1:
                selected.add(code)
        rng = random.Random(20261030 + case)
        selected.update(rng.randrange(1, 4 ** dimension) for _ in range(2))
        for code in selected:
            pair = []
            for cut in (rows, columns):
                matrix = [0] * 49
                for i, source in enumerate(cut):
                    scalar = (code >> (2 * i)) & 3
                    matrix = [a ^ scale_vector(b, scalar, 49) for a, b in zip(matrix, source)]
                pair.append(independent_rank(matrix, 49))
            assert pair == record["rank_table"][code]
            checked += 1
        if case in (0, 2, 3, 34, 55):
            literal += literal_comparison(family, record["complement_matrix"])
        if case in (2, 3):
            distances = integer_distances(record["rank_table"])
            assert distances[1] == record["integer_shortest_lower_bound"] == 36
            assert sorted(record["rank_table"][i] for i in record["witness"]) == [[1, 49], [22, 22], [49, 1]]
    return {"cases": len(data["cases"]), "complete_projective_points": projective,
            "scalar_F4_elimination_controls": checked,
            "direct_physical_coefficient_entries_checked": literal,
            "all_affordable_spans_miss_endpoint": True,
            "scalar_omega_identity_exact_cost": 36,
            "exact_integer_shortest_path_potentials_checked": 2,
            "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--field4", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(7)
    result = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "binary": audit_binary(args.binary, family)}
    print(json.dumps(result["binary"]), flush=True)
    result["field4"] = audit_field4(args.field4, family)
    print(json.dumps(result["field4"]), flush=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
