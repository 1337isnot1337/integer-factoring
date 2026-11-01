#!/usr/bin/env python3
"""Decide the finite cost target from a complete capped F2 coefficient screen."""

import argparse
import hashlib
import heapq
import json
from math import isqrt
from pathlib import Path


def unpack_binary(code, dimension):
    assert all(((code >> (2 * i)) & 3) in (0, 1) for i in range(dimension))
    return sum(((code >> (2 * i)) & 1) << i for i in range(dimension))


def analyze(basis_file, screen_file):
    data = json.loads(basis_file.read_text())
    events = [json.loads(line) for line in screen_file.read_text().splitlines()]
    start, complete = events[0], events[-1]
    dimension, labels = data["basis_dimension"], data["labels"]
    assert start["event"] == "start" and complete["event"] == "complete"
    assert start["field_size"] == 2 and start["dimension"] == dimension
    assert start["labels"] == labels
    assert start["code_start"] == 1 and start["code_count"] == (1 << dimension) - 1
    assert complete["projective_tested"] == (1 << dimension) - 1
    gate = start["integer_cost_gate"]
    assert gate ** 70 > labels ** 99
    assert start["product_cap"] == gate * gate - 1
    records = {}
    for record in events:
        if record["event"] == "eligible":
            mask = unpack_binary(record["code"], dimension)
            assert mask and mask not in records
            assert record["product"] == record["row_rank"] * record["column_rank"] < gate * gate
            records[mask] = record
    assert len(records) == complete["eligible"]
    target = sum(value << i for i, value in enumerate(data["endpoint_coordinates"]))
    span = {}
    for value in records:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in span:
                span[pivot] = value
                break
            value ^= span[pivot]
    residual = target
    while residual:
        pivot = residual.bit_length() - 1
        if pivot not in span:
            break
        residual ^= span[pivot]
    result = {"labels": labels, "dimension": dimension,
              "integer_cost_gate": gate, "tested_vectors": complete["projective_tested"],
              "eligible_vectors": len(records), "eligible_span_dimension": len(span),
              "endpoint_remainder": residual,
              "basis_sha256": hashlib.sha256(basis_file.read_bytes()).hexdigest(),
              "screen_sha256": hashlib.sha256(screen_file.read_bytes()).hexdigest()}
    if residual:
        result["status"] = "excluded: affordable vectors do not span crossing"
        return result
    scale = 10 ** 12
    weights = {mask: isqrt(record["product"] * scale * scale)
               for mask, record in records.items()}
    distances, previous = {0: 0}, {}
    queue = [(0, 0)]
    bound = gate * scale
    while queue:
        cost, state = heapq.heappop(queue)
        if distances[state] != cost:
            continue
        if cost >= bound:
            result["status"] = "excluded: optimistic shortest cost exceeds integer gate"
            result["cost_lower_numerator"] = bound
            break
        if state == target:
            witness = []
            while state:
                parent, mask = previous[state]
                witness.append(mask)
                state = parent
            upper = sum(weights[m] + (weights[m] ** 2 != records[m]["product"] * scale * scale)
                        for m in witness)
            result.update(cost_lower_numerator=cost, witness_upper_numerator=upper,
                          witness=[{"mask": m, **records[m]} for m in witness])
            # 1393/985 < sqrt(2) < 577/408. Integer powers avoid floating
            # comparisons with the irrational finite-seed threshold.
            if upper ** 985 < labels ** 1393 * scale ** 985:
                result["status"] = "native seed candidate: independent tensor and rank audit required"
            elif cost ** 408 > labels ** 577 * scale ** 408:
                result["status"] = "excluded: exact lower cost exceeds seed threshold"
            else:
                result["status"] = "unresolved threshold interval: refine the radical and exponent bounds"
            break
        for mask, weight in weights.items():
            successor, value = state ^ mask, cost + weight
            if value < distances.get(successor, value + 1):
                distances[successor] = value
                previous[successor] = state, mask
                heapq.heappush(queue, (value, successor))
    result["cost_denominator"] = scale
    result["shortest_path_states"] = len(distances)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.basis, args.screen)
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(result, sort_keys=True, indent=2) + "\n"
    args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
