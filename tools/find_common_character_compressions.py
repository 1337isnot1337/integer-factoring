#!/usr/bin/env python3
"""Find directly verifiable scalar compressions without closing a matrix algebra.

The independent checker constructs the left maps and checks all matrix
entries.  This discovery tool only proposes right-vector witnesses.
"""

import argparse
import json
from pathlib import Path
import random

from check_common_character_compressions import coefficient_digest


def insert(pivots, vector):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in pivots:
            pivots[pivot] = vector
            return
        vector ^= pivots[pivot]


def columns(rows):
    return [sum(((row >> col) & 1) << i for i, row in enumerate(rows))
            for col in range(len(rows))]


def apply(matrix_columns, vector):
    result = 0
    while vector:
        bit = vector & -vector
        result ^= matrix_columns[bit.bit_length() - 1]
        vector ^= bit
    return result


def search_block(matrices, character, seed, trials):
    ambient = len(matrices[0])
    assert matrices[0] == [1 << i for i in range(ambient)] and character[0] == 1
    differences = [columns([row ^ ((1 << i) if scalar else 0)
                            for i, row in enumerate(matrix)])
                   for matrix, scalar in zip(matrices[1:], character[1:])]
    rng = random.Random(seed)
    order = list(range(ambient))
    rng.shuffle(order)
    candidates = [1 << i for i in order] + [rng.getrandbits(ambient) for _ in range(trials)]
    vectors, zero_space = [], {}
    for vector in candidates:
        enlarged = zero_space.copy()
        images = [apply(matrix, vector) for matrix in differences]
        for image in images:
            insert(enlarged, image)
        difference_dimension = len(enlarged)
        for old in vectors:
            insert(enlarged, old)
        insert(enlarged, vector)
        if len(enlarged) != difference_dimension + len(vectors) + 1:
            continue
        for image in images:
            insert(zero_space, image)
        vectors.append(vector)
    return {"size": ambient, "compression": len(vectors),
            "annihilated_span": len(zero_space), "vectors": [hex(vector) for vector in vectors]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", type=Path, required=True)
    parser.add_argument("--character", required=True, help="Binary functional, e.g.11110000")
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--random-trials", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    character = [int(bit) for bit in args.character]
    assert set(character) <= {0, 1} and character[0] == 1
    data = json.loads(args.basis.read_text())
    assert len(character) == data["basis_dimension"]
    rows, columns_ = [[[[int(row, 16) for row in block] for block in matrix]
                       for matrix in data[name]]
                      for name in ("basis_row_blocks", "basis_column_blocks")]
    cuts = {}
    for name, matrices in (("row", rows), ("column", columns_)):
        records = []
        for block in range(len(data["block_sizes"])):
            local = [matrix[block] for matrix in matrices]
            choices = [search_block(local, character, seed, args.random_trials)
                       for seed in range(args.restarts)]
            best = max(choices, key=lambda item: item["compression"])
            records.append(best)
            print(name, block, best["compression"], flush=True)
        cuts[name] = records
    result = {"labels": data["labels"], "pairs": data["code"]["pairs"],
              "coefficient_space_dimension": len(character), "character": character,
              "coefficient_matrices_sha256": coefficient_digest(rows, columns_), "cuts": cuts}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print({cut: sum(item["compression"] for item in records) for cut, records in cuts.items()})


if __name__ == "__main__":
    main()
