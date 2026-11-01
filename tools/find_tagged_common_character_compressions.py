#!/usr/bin/env python3
"""Fast coordinate-subspace discovery using a tagged Gaussian basis.

This is separate from the earlier finder so its fixed certificates remain
exactly reproducible.  The independent checker constructs the left maps.
"""

import argparse
import json
from pathlib import Path
import random

from assemble_two_even_site_cores import assemble


def columns(rows):
    result = [0] * len(rows)
    for i, row in enumerate(rows):
        while row:
            bit = row & -row
            result[bit.bit_length() - 1] |= 1 << i
            row ^= bit
    return result


def insert_tagged(pivots, vector, label):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in pivots:
            pivots[pivot] = vector, label
            return True
        old, tag = pivots[pivot]
        vector ^= old
        label ^= tag
    return label == 0


def search_block(matrices, character, seed=0):
    size = len(matrices[0])
    assert character[0] == 1 and matrices[0] == [1 << i for i in range(size)]
    differences = [columns([row ^ ((1 << i) if scalar else 0)
                            for i, row in enumerate(matrix)])
                   for matrix, scalar in zip(matrices[1:], character[1:])]
    order = list(range(size))
    random.Random(seed).shuffle(order)
    pivots, selected = {}, []
    for index in order:
        trial = pivots.copy()
        if not all(insert_tagged(trial, matrix[index], 0) for matrix in differences):
            continue
        if not insert_tagged(trial, 1 << index, 1 << len(selected)):
            continue
        pivots = trial
        selected.append(index)
    return {"size": size, "compression": len(selected),
            "annihilated_span": len(pivots) - len(selected),
            "coordinate_indices": selected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--pairs", type=int)
    source.add_argument("--basis", type=Path)
    parser.add_argument("--character", required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    character = [int(bit) for bit in args.character]
    assert set(character) <= {0, 1} and character[0] == 1
    if args.basis:
        data = json.loads(args.basis.read_text())
        cuts = [[[[int(row, 16) for row in block] for block in matrix]
                 for matrix in data[name]] for name in ("basis_row_blocks", "basis_column_blocks")]
        pairs = data.get("pairs", data.get("code", {}).get("pairs"))
    else:
        data, cuts = assemble(args.pairs)
        pairs = args.pairs
    assert len(character) == data["basis_dimension"]
    result = {"labels": data["labels"], "pairs": pairs, "character": character,
              "discovery": {"algorithm": "tagged coordinate-subspace search", "seed": args.seed},
              "cuts": {}}
    for name, matrices in zip(("row", "column"), cuts):
        records = []
        for block in range(len(matrices[0])):
            witness = search_block([matrix[block] for matrix in matrices], character, args.seed)
            records.append(witness)
            print(name, block, witness["compression"], flush=True)
        result["cuts"][name] = records
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
