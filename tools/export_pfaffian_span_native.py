#!/usr/bin/env python3
"""Encode an exact binary tensor basis for the F2/F4 native span screen."""

import argparse
import hashlib
import json
from pathlib import Path


def encode(data, gate=None):
    labels = data["labels"]
    if gate is None:
        gate = 1
        while gate ** 70 <= labels ** 99:
            gate += 1
    assert gate >= 2
    dimension = data["basis_dimension"]
    sizes = data["block_sizes"]
    assert sum(sizes) == labels * labels and all(n > 0 for n in sizes)
    lines = [f"{labels} {dimension} {len(sizes)} {gate}", " ".join(map(str, sizes))]
    for name in ("basis_row_blocks", "basis_column_blocks"):
        assert len(data[name]) == dimension
        for basis in data[name]:
            assert len(basis) == len(sizes)
            for n, block in zip(sizes, basis):
                assert len(block) == n
                for row in block:
                    value = int(row, 16)
                    assert value >= 0 and value.bit_length() <= n
                    lines.append(" ".join(format((value >> (64 * j)) & ((1 << 64) - 1), "x")
                                          for j in range((n + 63) // 64)))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("basis", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cost-gate", type=int)
    args = parser.parse_args()
    data = json.loads(args.basis.read_text())
    text = encode(data, args.cost_gate)
    args.output.write_text(text)
    print(json.dumps({"labels": data["labels"], "dimension": data["basis_dimension"],
                      "basis_sha256": hashlib.sha256(args.basis.read_bytes()).hexdigest(),
                      "native_input_sha256": hashlib.sha256(text.encode()).hexdigest()}))


if __name__ == "__main__":
    main()
