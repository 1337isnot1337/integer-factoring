#!/usr/bin/env python3
"""A single-coordinate certificate for mixing the recorded cubic pieces."""

import argparse
import hashlib
import json
from pathlib import Path

from explore_cubic_multirow_outputs import CubicSpace


def span_and_residual(vectors, target):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    while target:
        pivot = target.bit_length() - 1
        if pivot not in pivots:
            break
        target ^= pivots[pivot]
    return len(pivots), target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("evidence/cubic-multirow-outputs.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    assert data["complete"] and len(data["records"]) == 680
    assert hashlib.sha256(Path(__file__).with_name("explore_cubic_multirow_outputs.py").read_bytes()).hexdigest() == data["source_sha256"]
    family = CubicSpace(7)
    endpoint = sum(1 << (49 * i + i) for i in range(49))
    pool, variable = {}, []
    for case, record in enumerate(data["records"]):
        tensor = family.tensor([list(map(tuple, row)) for row in record["corrections"]])
        basis = [sum(((tensor[i][j] >> degree) & 1) << (49 * i + j)
                     for i in range(49) for j in range(49)) for degree in (0, 1, 2, 4)]
        assert basis[0] == endpoint
        variable.extend((basis[1], basis[2], basis[1] ^ basis[3]))
        for mask, (row_rank, column_rank) in enumerate(record["rank_table"]):
            if not mask or row_rank * column_rank >= 256:
                continue
            value = 0
            for i, matrix in enumerate(basis):
                if mask >> i & 1:
                    value ^= matrix
            assert not (value >> (49 * 3 + 3)) & 1
            if value in pool:
                assert pool[value]["rank_pair"] == [row_rank, column_rank]
            else:
                pool[value] = {"case": case, "mask": mask,
                               "rank_pair": [row_rank, column_rank],
                               "tensor_sha256": hashlib.sha256(value.to_bytes(301, "little")).hexdigest()}
    dimension, residual = span_and_residual(pool, endpoint)
    whole_dimension, whole_residual = span_and_residual(variable, endpoint)
    assert len(pool) == 659 and dimension == 341 and residual
    assert whole_dimension == 819 and whole_residual == 0
    assert 16 ** 70 > 7 ** 99 and 99 ** 2 > 2 * 70 ** 2
    result = {"input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "source_spaces": 680, "unique_recorded_affordable_tensors": len(pool),
              "integer_cost_gate": 16, "affordable_span_dimension": dimension,
              "separating_logical_entry_x_y_v_u": [0, 3, 0, 3],
              "every_pool_tensor_entry": 0, "crossing_entry": 1,
              "all_mixtures_of_recorded_binary_affordable_pieces_excluded": True,
              "full_variable_span_dimension": whole_dimension,
              "crossing_in_full_variable_span": True,
              "pieces": list(pool.values()),
              "scope": "The recorded binary coefficient tensors, with arbitrary characteristic-two scalars; new cross-space recombinations and untested field coefficients remain open."}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "pieces"}), flush=True)


if __name__ == "__main__":
    main()
