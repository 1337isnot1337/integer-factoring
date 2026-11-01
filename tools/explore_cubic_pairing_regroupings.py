#!/usr/bin/env python3
"""Bounded exact F2/F4 regroupings inside the 177-dimensional cubic space."""

import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_fields import MULTIPLY, rank4
from check_cubic_varying_image_escape import xor_matrices


def pool():
    source = json.loads(Path("evidence/cubic-varying-image-escape.json").read_text())
    family = ComplementOutputs(7)
    candidates, names, seen, q3 = [], [], {}, []

    def insert(tensor, name):
        key = tuple(tensor[0])
        if key not in seen:
            seen[key] = len(candidates)
            candidates.append(tensor)
            names.append(name)

    insert(([0] * 49, [0] * 49), ["zero"])
    for i, record in enumerate(source["rank_one_output_maps"]):
        rows, columns = family.matrices(record["complement_rows"])
        q3.append((rows[4], columns[4]))
        for k in range(5):
            insert((rows[k], columns[k]), ["raw", i, k])
    diagonal = sum(1 << (8 * i) for i in range(7))
    identity = [1 << i for i in range(49)]
    rank_one = [diagonal if i // 7 == i % 7 else 0 for i in range(49)]
    swap = [1 << (7 * (i % 7) + i // 7) for i in range(49)]
    endpoint, h, n = (identity, identity), (rank_one, swap), (swap, rank_one)
    balanced = tuple(xor_matrices(e, first, second) for e, first, second in zip(endpoint, h, n))
    for name, tensor in (("H", h), ("N", n), ("K", balanced)):
        insert(tensor, [name])
    for record in source["affordable_edge_pieces"]:
        tensor = tuple(xor_matrices(*(q3[i][cut] for i in record["q3_map_indices"])) for cut in range(2))
        insert(tensor, ["edge"] + record["labels"])
    targets = [endpoint, tuple(xor_matrices(e, value) for e, value in zip(endpoint, h)),
               tuple(xor_matrices(e, value) for e, value in zip(endpoint, n)), balanced]
    assert len(candidates) == 538
    return candidates, names, targets


def payload(candidates, targets):
    lines = [str(len(candidates))]
    for tensor in targets + candidates:
        for matrix in tensor:
            lines.extend(format(value, "x") for value in matrix)
    return "\n".join(lines) + "\n"


def field_combination(first, second, a, b, target=None):
    result = []
    for cut in range(2):
        matrix = []
        for i in range(49):
            row = 0
            for j in range(49):
                value = (a if (first[cut][i] >> j) & 1 else 0) ^ (b if (second[cut][i] >> j) & 1 else 0)
                if target is not None:
                    value ^= (target[cut][i] >> j) & 1
                row |= value << (2 * j)
            matrix.append(row)
        result.append(matrix)
    return result


def improves(first, second, offset):
    difference = (36 - offset) ** 2 - first - second
    return difference > 0 and 4 * first * second < difference ** 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.native and not args.replay:
        parser.error("--native or --replay is required")
    candidates, names, targets = pool()
    text = payload(candidates, targets)
    hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in (Path(__file__).name, "search_cubic_pairing_regroupings.cpp",
                           "check_cubic_varying_image_escape.py", "explore_cubic_complement_outputs.py",
                           "explore_cubic_complement_fields.py")}
    saved = json.loads(args.replay.read_text()) if args.replay else None
    if saved:
        assert saved["source_sha256"] == hashes
    offsets = (0, 7, 7, 14)
    results, controls = [], []
    representatives = [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3)]
    covered = {(a, b) for a, b in representatives}
    covered.update((MULTIPLY[a][a], MULTIPLY[b][b]) for a, b in representatives)
    assert covered == {(a, b) for a in range(1, 4) for b in range(1, 4)}
    for field, expected_cases, expected_comparisons in ((2, 144453, 533823), (4, 722265, 2495874)):
        if args.native:
            screen = subprocess.run([str(args.native), str(field)], input=text, text=True,
                                    stdout=subprocess.PIPE, check=True).stdout
        else:
            screen = saved["fields"][len(results)]["native_output"]
        events = [json.loads(line) for line in screen.splitlines()]
        assert len(events) == 1
        assert events[0] == {"event": "complete", "field": field, "candidates": 538,
                             "coefficient_cases": expected_cases,
                             "exact_target_comparisons": expected_comparisons,
                             "improvements_below_36": 0}
        results.append({"field": field, "native_output": screen,
                        "native_output_sha256": hashlib.sha256(screen.encode()).hexdigest()})
        rng = random.Random(20261030 + field)
        for _ in range(32):
            i, j = rng.sample(range(len(candidates)), 2)
            a, b = (1, 1) if field == 2 else representatives[rng.randrange(5)]
            target = rng.randrange(4)
            first = field_combination(candidates[i], candidates[j], a, b)
            second = field_combination(candidates[i], candidates[j], a, b, targets[target])
            ranks = [[rank4(matrix, 49, low_pivot=True) for matrix in tensor] for tensor in (first, second)]
            assert ranks == [[rank4(matrix, 49) for matrix in tensor] for tensor in (first, second)]
            assert not improves(ranks[0][0] * ranks[0][1], ranks[1][0] * ranks[1][1], offsets[target])
            controls.append({"field": field, "indices": [i, j], "scalars": [a, b],
                             "target": target, "rank_pairs": ranks})
    result = {"source_sha256": hashes, "native_input_sha256": hashlib.sha256(text.encode()).hexdigest(),
              "candidate_count": len(candidates), "candidate_definitions": names,
              "targets": ["E", "E+H", "E+N", "E+H+N"], "fixed_costs": list(offsets),
              "fields": results, "independent_rank_controls": controls,
              "all_recorded_regroupings_fail_to_improve_36": True,
              "scope": "At most two scalar-weighted candidates, over F2 or F4, against the four displayed targets with their stated fixed costs. This is not optimization over the full 177-dimensional space."}
    if saved:
        assert result == saved
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("candidate_definitions", "independent_rank_controls")}), flush=True)


if __name__ == "__main__":
    main()
