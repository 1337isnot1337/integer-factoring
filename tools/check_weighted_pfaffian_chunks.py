#!/usr/bin/env python3
"""Independent binary audit of a complete chunked weighted coefficient screen."""

import argparse
import hashlib
import json
from pathlib import Path
import random

from check_symplectic_seed_tensor_fields import independent_field_rank
from check_weighted_pfaffian_screen import read_native
from symplectic_seed_weighted_modules import load_certificate, verify_certificate


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def flatten(matrix, size):
    return sum(row << (size * i) for i, row in enumerate(matrix))


def make_ranker(cuts):
    prepared = {}
    for name, blocks in cuts.items():
        prepared[name] = [(block["size"], block["multiplicity"],
                           [flatten(matrix, block["size"]) for matrix in block["generator_rows"]], {})
                          for block in blocks]

    def exact(mask):
        coordinates = list(bits(mask))
        result = []
        for name in ("row", "column"):
            total = 0
            for size, multiplicity, matrices, cache in prepared[name]:
                value = 0
                for i in coordinates:
                    value ^= matrices[i]
                rank = cache.get(value)
                if rank is None:
                    bitmask = (1 << size) - 1
                    rows = [(value >> (size * i)) & bitmask for i in range(size)]
                    rank = independent_field_rank(rows, [0] * size)
                    if len(cache) >= 256:
                        cache.clear()
                    cache[value] = rank
                total += multiplicity * rank
            result.append(total)
        return result

    return exact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--export-metadata", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    labels, dimension, gate, cuts = read_native(args.input)
    assert 99 ** 2 > 2 * 70 ** 2 and gate ** 70 >= labels ** 99
    metadata = json.loads(args.export_metadata.read_text())
    if "native_export" in metadata:
        metadata = metadata["native_export"]
    assert metadata["input_sha256"] == hashlib.sha256(args.input.read_bytes()).hexdigest()
    assert (metadata["labels"], metadata["dimension"], metadata["gate"]) == (labels, dimension, gate)
    verification = verify_certificate(load_certificate())

    records, chunks = {}, []
    next_start, current, eligible_in_chunk = 1, None, 0
    for line in args.screen.read_text().splitlines():
        event = json.loads(line)
        if event["event"] == "start":
            assert current is None
            assert event["code_start"] == next_start and event["field_size"] == 2
            assert not event["all_ranks"]
            assert (event["labels"], event["dimension"], event["integer_cost_gate"]) == (labels, dimension, gate)
            current, eligible_in_chunk = event, 0
        elif event["event"] == "eligible":
            assert current is not None
            index = event["index"]
            assert next_start <= index < next_start + current["code_count"]
            assert index not in records
            assert event["code"] == sum(1 << (2 * bit) for bit in bits(index))
            assert event["product"] == event["row_rank"] * event["column_rank"] < gate * gate
            records[index] = event
            eligible_in_chunk += 1
        elif event["event"] == "complete":
            assert current is not None
            assert event["projective_tested"] == current["code_count"]
            assert event["eligible"] == eligible_in_chunk
            chunks.append({"start": next_start, "count": current["code_count"],
                           "eligible": eligible_in_chunk, "milliseconds": event["milliseconds"]})
            next_start += current["code_count"]
            current = None
        else:
            raise AssertionError("unexpected native event")
    assert current is None and next_start == 1 << dimension

    selected = set(records)
    rng = random.Random(287720)
    while len(selected) < len(records) + args.samples:
        selected.add(rng.randrange(1, 1 << dimension))
    exact = make_ranker(cuts)
    for count, index in enumerate(sorted(selected), 1):
        row, column = exact(index)
        assert (index in records) == (row * column < gate * gate)
        if index in records:
            assert (row, column) == (records[index]["row_rank"], records[index]["column_rank"])
        if count % 250 == 0:
            print(json.dumps({"independent_vectors_checked": count, "total": len(selected)}), flush=True)
    pivots = {}
    for index in records:
        value = index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    remainder = 1
    while remainder and remainder.bit_length() - 1 in pivots:
        remainder ^= pivots[remainder.bit_length() - 1]
    result = {
        "labels": labels, "coefficient_dimension": dimension, "field_size": 2,
        "projective_vectors_tested": (1 << dimension) - 1, "integer_cost_gate": gate,
        "rank_bound": metadata["rank_bound"], "maximum_block_size": metadata["maximum_block_size"],
        "omitted_weighted_dimensions": {cut: record["omitted_weighted_dimension"]
                                         for cut, record in metadata["cuts"].items()},
        "eligible_vectors": list(records.values()), "independently_checked_vectors": len(selected),
        "eligible_span_dimension": len(pivots), "eligible_span_basis": list(pivots.values()),
        "endpoint_remainder": remainder, "excluded_by_eligible_span": bool(remainder),
        "chunks": chunks, "certificate_verification": verification,
        "native_input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "screen_sha256": hashlib.sha256(args.screen.read_bytes()).hexdigest(),
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "symplectic_seed_weighted_modules.py",
                                       "check_symplectic_seed_tensor_fields.py", "check_weighted_pfaffian_screen.py",
                                       "screen_weighted_pfaffian_wide.cpp", "screen_pfaffian_span_f4.cpp")},
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("eligible_vectors", "chunks", "certificate_verification")}, indent=2))


if __name__ == "__main__":
    main()
