#!/usr/bin/env python3
"""Independent F4 audit and finite rank-cost certificate for Pfaffian spans."""

import argparse
import hashlib
import heapq
import json
from math import isqrt
from pathlib import Path
import random

MULTIPLY = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
INVERSE = (0, 1, 3, 2)


def scale_vector(vector, scalar, width):
    return sum(MULTIPLY[(vector >> (2 * j)) & 3][scalar] << (2 * j)
               for j in range(width))


def span_basis(vectors, width):
    basis = {}
    for vector in vectors:
        while vector:
            pivot = (vector.bit_length() - 1) // 2
            scalar = (vector >> (2 * pivot)) & 3
            if pivot not in basis:
                basis[pivot] = scale_vector(vector, INVERSE[scalar], width)
                break
            vector ^= scale_vector(basis[pivot], scalar, width)
    return basis


def remainder(vector, basis, width):
    while vector:
        pivot = (vector.bit_length() - 1) // 2
        if pivot not in basis:
            break
        vector ^= scale_vector(basis[pivot], (vector >> (2 * pivot)) & 3, width)
    return vector


def independent_rank(rows, width):
    # Scalar F4 elimination with the opposite pivot direction from the native
    # implementation. Each entry is a two-bit field element, not a binary rank.
    basis = {}
    for row in rows:
        while row:
            pivot = ((row & -row).bit_length() - 1) // 2
            scalar = (row >> (2 * pivot)) & 3
            if pivot not in basis:
                basis[pivot] = scale_vector(row, INVERSE[scalar], width)
                break
            row ^= scale_vector(basis[pivot], scalar, width)
    return len(basis)


def prepare_basis(data):
    result = []
    for name in ("basis_row_blocks", "basis_column_blocks"):
        cut = []
        for generator in data[name]:
            cut.append([[sum(((int(row, 16) >> j) & 1) << (2 * j)
                              for j in range(n)) for row in block]
                        for n, block in zip(data["block_sizes"], generator)])
        result.append(cut)
    return result


def exact_ranks(data, prepared, code):
    result = []
    dimension = data["basis_dimension"]
    for cut in prepared:
        total = 0
        for b, n in enumerate(data["block_sizes"]):
            rows = [0] * n
            for i in range(dimension):
                scalar = (code >> (2 * i)) & 3
                if not scalar:
                    continue
                for j, row in enumerate(cut[i][b]):
                    # Every supplied generator entry is binary, so scaling its
                    # packed row is ordinary multiplication by 1,2,or3.
                    rows[j] ^= row * scalar
            total += independent_rank(rows, n)
        result.append(total)
    return result


def shortest_lower_bound(records, dimension, target, scale):
    weights = {}
    for record in records:
        weight = isqrt(record["product"] * scale * scale)
        for scalar in (1, 2, 3):
            weights[scale_vector(record["code"], scalar, dimension)] = weight
    queue = [(0, 0)]
    distance = {0: 0}
    while queue:
        cost, state = heapq.heappop(queue)
        if distance[state] != cost:
            continue
        if state == target:
            return cost, len(distance)
        for vector, weight in weights.items():
            new, successor = cost + weight, state ^ vector
            if new < distance.get(successor, new + 1):
                distance[successor] = new
                heapq.heappush(queue, (new, successor))
    return None, len(distance)


def verify(basis_path, screen_path, samples):
    data = json.loads(basis_path.read_text())
    events = [json.loads(line) for line in screen_path.read_text().splitlines()]
    start, complete = events[0], events[-1]
    dimension, labels = data["basis_dimension"], data["labels"]
    assert start["event"] == "start" and complete["event"] == "complete"
    assert start["labels"] == labels and start["dimension"] == dimension
    assert start["code_start"] == 1 and start["code_count"] == 4 ** dimension - 1
    assert complete["projective_tested"] == (4 ** dimension - 1) // 3
    gate = start["integer_cost_gate"]
    # sqrt(2)<99/70 and gate^70>d^99 give an entirely integer upper
    # certificate for the irrational finite-seed threshold.
    assert 99 ** 2 > 2 * 70 ** 2 and gate ** 70 > labels ** 99
    assert start["product_cap"] == gate * gate - 1
    eligible = [x for x in events if x["event"] == "eligible"]
    assert len(eligible) == complete["eligible"]
    records = {x["code"]: x for x in eligible}
    assert len(records) == len(eligible)
    for record in eligible:
        assert record["product"] == record["row_rank"] * record["column_rank"] < gate * gate
        assert (record["code"] >> (2 * ((record["code"].bit_length() - 1) // 2))) == 1
    rng = random.Random(8675309)
    selected = set(records)
    while len(selected) < len(records) + samples:
        code = rng.randrange(1, 4 ** dimension)
        pivot = (code.bit_length() - 1) // 2
        selected.add(scale_vector(code, INVERSE[(code >> (2 * pivot)) & 3], dimension))
    prepared = prepare_basis(data)
    for code in sorted(selected):
        r, c = exact_ranks(data, prepared, code)
        assert (code in records) == (r * c < gate * gate)
        if code in records:
            assert (r, c) == (records[code]["row_rank"], records[code]["column_rank"])
    target = sum(c << (2 * i) for i, c in enumerate(data["endpoint_coordinates"]))
    assert exact_ranks(data, prepared, target) == [labels * labels] * 2
    basis = span_basis(records, dimension)
    residue = remainder(target, basis, dimension)
    scale = 10 ** 12
    lower, visited = (None, 0) if residue else shortest_lower_bound(eligible, dimension, target, scale)
    excluded = bool(residue) or lower >= gate * scale
    assert excluded, "no exclusion certificate: investigate this candidate"
    return {"name": data["name"], "labels": labels, "coefficient_dimension": dimension,
            "field": "F2[t]/(t^2+t+1)", "integer_cost_gate": gate,
            "projective_vectors_tested": complete["projective_tested"],
            "eligible_count": len(eligible), "eligible_vectors": eligible,
            "independently_checked_vectors": len(selected),
            "eligible_span_basis": {str(k): v for k, v in sorted(basis.items())},
            "eligible_span_dimension": len(basis), "endpoint_remainder": residue,
            "restricted_shortest_cost_lower_numerator": lower,
            "cost_denominator": scale, "shortest_path_states": visited,
            "all_F4_decompositions_fail_seed_criterion": excluded,
            "basis_sha256": hashlib.sha256(basis_path.read_bytes()).hexdigest(),
            "screen_sha256": hashlib.sha256(screen_path.read_bytes()).hexdigest(),
            "scope": "the supplied finite coefficient space over F4; other fields and label retracts are not excluded"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.basis, args.screen, args.samples)
    result["source_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "screen_pfaffian_span_f4.cpp")}
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "eligible_vectors"}, indent=2))


if __name__ == "__main__":
    main()
