#!/usr/bin/env python3
"""Independent scalar audit of a complete weighted F2/F4 coefficient screen."""

import argparse
import hashlib
import json
from pathlib import Path
import random

from check_pfaffian_span_f4 import INVERSE, remainder, scale_vector, span_basis
from check_symplectic_seed_tensor_fields import exact_ranks


def read_native(path):
    values = iter(path.read_text().split())
    labels, dimension, gate = (int(next(values)) for _ in range(3))
    cuts = {}
    for cut in ("row", "column"):
        blocks = []
        for _ in range(int(next(values))):
            size, multiplicity = int(next(values)), int(next(values))
            assert size > 0 and multiplicity > 0
            matrices = []
            for _ in range(dimension):
                rows = []
                for _ in range(size):
                    value = sum(int(next(values), 16) << shift for shift in range(0, size, 64))
                    assert value >> size == 0
                    rows.append(value)
                matrices.append(rows)
            blocks.append({"size": size, "multiplicity": multiplicity,
                           "generator_rows": matrices})
        assert sum(block["size"] * block["multiplicity"] for block in blocks) == labels * labels
        cuts[cut] = blocks
    assert next(values, None) is None
    return labels, dimension, gate, cuts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--target-code", type=int, default=1)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert args.samples >= 0
    labels, dimension, gate, cuts = read_native(args.input)
    events = [json.loads(line) for line in args.screen.read_text().splitlines()]
    start, complete = events[0], events[-1]
    assert start["event"] == "start" and complete["event"] == "complete"
    assert (start["labels"], start["dimension"], start["integer_cost_gate"]) == (labels, dimension, gate)
    field = start["field_size"]
    assert field in (2, 4)
    assert start["code_start"] == 1 and start["code_count"] == field ** dimension - 1
    expected = (field ** dimension - 1) // (field - 1)
    assert complete["projective_tested"] == expected
    assert 99 ** 2 > 2 * 70 ** 2 and gate ** 70 > labels ** 99
    eligible = [record for record in events if record["event"] in ("eligible", "rank")
                and record["product"] < gate * gate]
    assert len(eligible) == complete["eligible"]
    records = {record["code"]: record for record in eligible}
    assert len(records) == len(eligible)
    selected = set(records)
    rng = random.Random(285050)
    target_count = min(expected, len(selected) + args.samples)
    while len(selected) < target_count:
        if field == 2:
            mask = rng.randrange(1, 1 << dimension)
            code = sum(((mask >> i) & 1) << (2 * i) for i in range(dimension))
        else:
            code = rng.randrange(1, 4 ** dimension)
            pivot = (code.bit_length() - 1) // 2
            code = scale_vector(code, INVERSE[(code >> (2 * pivot)) & 3], dimension)
        selected.add(code)
    for code in sorted(selected):
        row, column = exact_ranks(cuts, code)
        assert (code in records) == (row * column < gate * gate)
        if code in records:
            assert (row, column) == (records[code]["row_rank"], records[code]["column_rank"])
            assert row * column == records[code]["product"]
    basis = span_basis(records, dimension)
    residue = remainder(args.target_code, basis, dimension)
    target_ranks = exact_ranks(cuts, args.target_code)
    result = {
        "labels": labels, "coefficient_dimension": dimension, "field_size": field,
        "projective_vectors_tested": expected, "integer_cost_gate": gate,
        "eligible_vectors": eligible, "independently_checked_vectors": len(selected),
        "eligible_span_dimension": len(basis), "endpoint_remainder": residue,
        "target_ranks": target_ranks, "excluded_by_eligible_span": bool(residue),
        "native_input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "screen_sha256": hashlib.sha256(args.screen.read_bytes()).hexdigest(),
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "check_symplectic_seed_tensor_fields.py",
                                       "check_pfaffian_span_f4.py", "screen_weighted_pfaffian_span.cpp",
                                       "screen_pfaffian_span_f4.cpp")},
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "eligible_vectors"}, indent=2))


if __name__ == "__main__":
    main()
