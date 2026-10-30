#!/usr/bin/env python3
"""Physical expansion and archive audits for the bounded cubic F8 probe."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import TABLE, matrices, native_input, rank8


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("evidence/cubic-complement-F8.json"))
    parser.add_argument("--native", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    assert data["complete"] and data["all_recorded_spaces_excluded_over_F8"]
    assert len(data["cases"]) == 428 and data["projective_vectors_tested"] == 2003468
    for name, digest in data["source_sha256"].items():
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest
    for a, b, c in product(range(8), repeat=3):
        assert TABLE[a][TABLE[b][c]] == TABLE[TABLE[a][b]][c]
        assert TABLE[a][b ^ c] == TABLE[a][b] ^ TABLE[a][c]
    total = affordable = 0
    for record in data["cases"]:
        events = [json.loads(line) for line in record["native_screen"].splitlines()]
        assert hashlib.sha256(record["native_screen"].encode()).hexdigest() == record["native_screen_sha256"]
        assert events[0]["dimension"] == 5 and events[0]["field_size"] == 8
        assert events[0]["product_cap"] == 255 and events[-1]["projective_tested"] == 4681
        points = [event for event in events if event["event"] == "eligible"]
        assert len(points) == events[-1]["eligible"] == record["affordable_count"]
        codes = [point["code"] for point in points]
        assert len(codes) == len(set(codes))
        dimension = rank8(codes, 5)
        assert dimension == record["affordable_span_dimension"]
        assert rank8(codes + [1], 5) == dimension + 1
        total += 4681
        affordable += len(points)
    assert total == 2003468 and affordable == 4376
    assert 16 ** 70 > 7 ** 99 and 99 ** 2 > 2 * 70 ** 2

    family = ComplementOutputs(7)
    low = family.space.outputs
    high = [tuple(2 if a == 1 else a for a in word) for word in low]
    literal = native_replays = 0
    for index in (0, 3, 131, 260, 391, 427):
        record = data["cases"][index]
        complement = record["complement_matrix"]
        rows, columns = matrices(family, complement)
        payload = native_input(rows, columns)
        assert hashlib.sha256(payload.encode()).hexdigest() == record["native_input_sha256"]
        if args.native:
            screen = subprocess.run([str(args.native)], input=payload, text=True,
                                    stdout=subprocess.PIPE, check=True).stdout
            assert screen == record["native_screen"]
            native_replays += 1
        outputs = [[(low[i], 1)] + [(high[j], value) for j, value in enumerate(complement[i]) if value]
                   for i in range(7)]
        for x, y, v, u in product(range(7), repeat=4):
            polynomial = 0
            for (first, a), (second, b) in product(outputs[x], outputs[u]):
                binary = family.space.kernel(first, second, y, v)
                spread = sum(((binary >> degree) & 1) << (3 * degree) for degree in range(13))
                polynomial ^= spread * TABLE[a][b]
            for k, degree in enumerate(family.degrees):
                value = (polynomial >> (3 * degree)) & 7
                assert ((rows[k][7 * x + y] >> (3 * (7 * v + u))) & 7) == value
                assert ((columns[k][7 * x + u] >> (3 * (7 * v + y))) & 7) == value
                literal += 1
    result = {"source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "complement_maps": 428, "complete_formal_projective_vectors": total,
              "affordable_projective_vectors": affordable,
              "independent_affordable_span_checks": 428,
              "source_independent_rank_vectors": data["independently_checked_vectors"],
              "literal_physical_coefficient_entries_checked": literal,
              "deterministic_native_control_replays": native_replays,
              "all_recorded_F8_spaces_fail_seed_criterion": True}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
