#!/usr/bin/env python3
"""Certify optimal binary recombinations of a completed native coefficient span."""

import argparse
from collections import Counter
import hashlib
import heapq
import json
from math import isqrt
from pathlib import Path


def shortest_xor(weights, target):
    count = target + 1
    infinity = sum(weights.values()) + 1
    distance = [infinity] * count
    previous = [None] * count
    distance[0] = 0
    queue = [(0, 0)]
    while queue:
        cost, state = heapq.heappop(queue)
        if distance[state] != cost:
            continue
        if state == target:
            witness = []
            while state:
                parent, mask = previous[state]
                witness.append(mask)
                state = parent
            return cost, witness
        for mask, weight in weights.items():
            successor = state ^ mask
            new = cost + weight
            if new < distance[successor]:
                distance[successor] = new
                previous[successor] = state, mask
                heapq.heappush(queue, (new, successor))
    raise AssertionError("crossing is not generated")


def bottleneck(records, target):
    basis = {}
    for mask, record in sorted(records.items(), key=lambda pair: (pair[1]["product"], pair[0])):
        value, witness = mask, 1 << (mask - 1)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value, witness
                break
            other, other_witness = basis[pivot]
            value ^= other
            witness ^= other_witness
        value, witness = target, 0
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                break
            other, other_witness = basis[pivot]
            value ^= other
            witness ^= other_witness
        if not value:
            masks = [i + 1 for i in range(target) if witness >> i & 1]
            assert sum_xor(masks) == target
            return record["product"], masks
    raise AssertionError("incomplete span")


def sum_xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def summarize(paths):
    records = {}
    reference = None
    intervals = []
    for path in paths:
        events = [json.loads(line) for line in path.read_text().splitlines()]
        assert events[0]["event"] == "start" and events[-1]["event"] == "complete"
        start, complete = events[0], events[-1]
        coefficients = [x for x in events if x["event"] == "coefficient"]
        signature = (start["labels"], start["label_bits"],
                     complete["matrix_digest_xor"], complete["matrix_digest_sum"],
                     complete["nonzero_polynomial_entries"], coefficients)
        if reference is None:
            reference = signature
        assert reference == signature, "coefficient chunks are from different tensors"
        chunk = [x for x in events if x["event"] == "recombination"]
        assert len(chunk) == complete["span_count"]
        assert [x["mask"] for x in chunk] == list(range(complete["span_start"],
                                                      complete["span_start"] + len(chunk)))
        for record in chunk:
            mask = record["mask"]
            assert mask not in records, "overlapping chunks"
            assert record["product"] == record["row_rank"] * record["column_rank"]
            records[mask] = record
        intervals.append({"start": complete["span_start"], "count": complete["span_count"]})
    dimension = len(reference[-1])
    target = (1 << dimension) - 1
    assert set(records) == set(range(1, target + 1)), "coefficient span is not complete"
    best_product, best_masks = bottleneck(records, target)
    scale = 10 ** 12
    lower, upper = {}, {}
    for mask, record in records.items():
        radicand = record["product"] * scale * scale
        value = isqrt(radicand)
        lower[mask] = value
        upper[mask] = value + (value * value != radicand)
    lo, lo_masks = shortest_xor(lower, target)
    hi, hi_masks = shortest_xor(upper, target)
    assert sum_xor(lo_masks) == sum_xor(hi_masks) == target
    assert lo <= hi
    labels = reference[0]
    return {"status": "complete binary coefficient-span certificate",
            "labels": labels, "label_bits": reference[1],
            "coefficient_dimension": dimension, "target_mask": target,
            "matrix_digest_xor": reference[2], "matrix_digest_sum": reference[3],
            "nonzero_polynomial_entries": reference[4],
            "coefficient_profile": reference[-1],
            "span_intervals": intervals, "nonzero_combinations": len(records),
            "minimum_bottleneck_product": best_product, "bottleneck_witness_masks": best_masks,
            "minimum_rank_sum_interval": {"lower_numerator": lo, "upper_numerator": hi,
                                          "denominator": scale},
            "rank_sum_upper_witness_masks": hi_masks,
            "rank_sum_lower_exceeds_d_to_three_halves": lo * lo > labels ** 3 * scale ** 2,
            "rank_product_histogram": {str(k): v for k, v in sorted(Counter(
                r["product"] for r in records.values()).items())},
            "records": [records[m] for m in range(1, target + 1)],
            "scope": "all binary recombinations of this finite coefficient span; other fields and crossing tensors are not excluded"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = summarize(args.inputs)
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["input_sha256"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in args.inputs}
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ("records", "rank_product_histogram", "coefficient_profile")}, indent=2))


if __name__ == "__main__":
    main()
