#!/usr/bin/env python3
"""Replay the six complete F4 symplectic-code coefficient-space certificates."""

import argparse
import hashlib
import json
from pathlib import Path

from check_pfaffian_span_f4 import verify


CASES = (
    "omega_radical_P5_L3",
    "omega_boundary_P4_L3",
    "omega_radical_P4_L4_half",
    "omega_radical_P4_L4",
    "omega_radical_P6_L4",
    "omega_radical_P5_L5_half",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis-dir", type=Path, required=True,
                        help="export directory from explore_symplectic_boundary_codes.py")
    parser.add_argument("--screen-dir", type=Path,
                        default=Path("evidence/symplectic-code-screens"))
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {"cases": [], "source_sha256": {
        name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "check_pfaffian_span_f4.py",
                     "screen_pfaffian_span_f4.cpp", "explore_symplectic_boundary_codes.py")}}
    for name in CASES:
        record = verify(args.basis_dir / (name + ".json"),
                        args.screen_dir / (name + "-F4.jsonl"), args.samples)
        result["cases"].append(record)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({key: record[key] for key in (
            "name", "labels", "coefficient_dimension", "projective_vectors_tested",
            "eligible_count", "eligible_span_dimension", "endpoint_remainder")}), flush=True)


if __name__ == "__main__":
    main()
