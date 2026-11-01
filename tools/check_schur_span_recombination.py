#!/usr/bin/env python3
"""Replay a complete finite F2 span-recombination screen.

The native search uses exact matrix ranks and integer-rounded-down edge
costs. Its complete Dijkstra distance is a lower bound for every
decomposition in the supplied coefficient span, with arbitrary length.
Python independently prepares the span and rechecks each returned witness.
"""

import argparse
from hashlib import sha256
import json
from math import isqrt, sqrt
from pathlib import Path
import subprocess
import tempfile


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            k = row.bit_length() - 1
            if k not in pivots:
                pivots[k] = row
                break
            row ^= pivots[k]
    return len(pivots)


def prepare(data):
    q = data["q"]
    n = q*q
    pivots = {}
    basis, original_masks = [], []
    endpoint = 0
    for index, rows in enumerate(data["pieces"]):
        assert len(rows) == n
        assert all(0 <= row < 1 << n for row in rows)
        value = sum(row << (n*j) for j, row in enumerate(rows))
        endpoint ^= value
        original = 1 << index
        while value:
            key = value.bit_length() - 1
            if key in pivots:
                old, old_original, _ = pivots[key]
                value ^= old
                original ^= old_original
            else:
                pivots[key] = value, original, len(basis)
                basis.append(value)
                original_masks.append(original)
                break
    target = 0
    value = endpoint
    while value:
        old, _, index = pivots[value.bit_length() - 1]
        value ^= old
        target ^= 1 << index
    assert endpoint == sum(1 << (n*(x*q+y)+y*q+x)
                           for x in range(q) for y in range(q))
    return basis, original_masks, target


def matrices(tensor, q):
    n = q*q
    rows = [(tensor >> (n*i)) & ((1 << n)-1) for i in range(n)]
    columns = [0] * n
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    if (rows[x*q+y] >> (u*q+v)) & 1:
                        columns[x*q+u] |= 1 << (y*q+v)
    return rows, columns


def run(path, native):
    raw = path.read_bytes()
    data = json.loads(raw)
    assert data["field"] == 2
    basis, original_masks, target = prepare(data)
    q, n = data["q"], data["q"]**2
    with tempfile.TemporaryDirectory(prefix="schur_span_") as temporary:
        input_path = Path(temporary) / "input.txt"
        lines = [f"{q} {len(basis)} {target}"]
        for tensor in basis:
            rows, _ = matrices(tensor, q)
            lines.append(" ".join(map(str, rows)))
        input_path.write_text("\n".join(lines) + "\n")
        completed = subprocess.run([str(native), str(input_path)], check=True,
                                   text=True, capture_output=True)
    result = json.loads(completed.stdout)
    assert result["q"] == q and result["target"] == target
    assert result["span_dimension"] == len(basis)
    assert result["combinations"] == 1 << len(basis)
    scale = result["integer_scale"]
    upper = 0
    total = 0
    for piece in result["pieces"]:
        mask = piece["mask"]
        tensor = original = 0
        for i, vector in enumerate(basis):
            if mask >> i & 1:
                tensor ^= vector
                original ^= original_masks[i]
        rows, columns = matrices(tensor, q)
        assert rank(rows) == piece["row_rank"]
        assert rank(columns) == piece["column_rank"]
        square = piece["row_rank"]*piece["column_rank"]*scale*scale
        lower = isqrt(square)
        upper += lower + (lower*lower != square)
        piece["original_shift_mask"] = original
        total ^= mask
    assert total == target
    assert upper == result["selected_cost_upper_integer"]
    assert result["minimum_cost_lower_integer"] <= upper
    result.update({"input_sha256": sha256(raw).hexdigest(),
                   "input_name": path.name,
                   "weight": data.get("h"),
                   "basis_original_shift_masks": original_masks,
                   "lower_bound_exceeds_q_to_three_halves":
                       result["minimum_cost_lower_integer"]**2 > q**3*scale**2,
                   "half_threshold_approximation": q**sqrt(2),
                   "status": "pass"})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--pieces", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    native = args.native.resolve()
    result = {"status": "pass",
              "scope": "Complete finite coefficient spans only; no general factoring lower bound.",
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "native_source_sha256": sha256(Path(__file__).with_name(
                  "schur_span_recombination.cpp").read_bytes()).hexdigest(),
              "records": [run(path, native) for path in args.pieces]}
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
