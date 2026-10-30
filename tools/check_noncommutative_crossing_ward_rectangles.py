#!/usr/bin/env python3
"""Audit stable rectangular bounds against independently literal Pfaffian ranks."""

from pathlib import Path
import argparse
import hashlib
import json

from check_noncommutative_crossing_ward_rank_descent import (
    rectangle_profile, actual_code_profile, partial_form_code,
)
from explore_symplectic_boundary_codes import boundary_code, half_orientation_code
from ward_rectangular_bounds import rectangular_bounds, ward_score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    specifications = set()
    for pairs, length in ((2, 2), (3, 3), (3, 2)):
        for j in range(2):
            for h in range(2 - j):
                ra, ry = length - 1 - j, length - 3 + j + 2 * h
                specifications.add((pairs, ra, ry, length, length))
                specifications.add((pairs, ra, length, length, ry))
    reports = []
    for spec in sorted(specifications):
        literal = rectangle_profile(*spec)
        upper, largest = rectangular_bounds(*spec)
        assert all(rank <= upper.get(key, 0) for key, rank in literal.items()), spec
        reports.append({"grades_P_RA_RY_CV_CU": list(spec),
                        "literal_profile": [[a, b, rank] for (a, b), rank in sorted(literal.items())],
                        "stable_upper_profile": [[a, b, rank] for (a, b), rank in sorted(upper.items()) if rank],
                        "largest_stable_block": largest})
    full = boundary_code(3, 3, True)
    cases = ((full, "global_P3_L3"),
             (half_orientation_code(full), "global_P3_L3_half"))
    complete = []
    for code, name in cases:
        exact = actual_code_profile(code)
        result = ward_score(code["pairs"], code["length"], 1, code["labels"])
        bound = {(a, b): (r, c) for a, b, r, c in result["profile"]}
        assert all(r <= bound.get(key, (0, 0))[0] and c <= bound.get(key, (0, 0))[1]
                   for key, (r, c) in exact.items())
        complete.append({"name": name, "labels": code["labels"],
                         "certified_integer_rank_sum_upper": result["integer_rank_sum_upper_bound"],
                         "exact_jet_ranks": [[a, b, r, c] for (a, b), (r, c) in sorted(exact.items())]})
    target = Path(__file__).with_name("ward_rectangular_bounds.py")
    result = {
        "status": "Every literal rectangle and complete restricted-code control satisfies the stable Ward bounds.",
        "rectangle_controls": reports,
        "complete_code_controls": complete,
        "rectangular_engine_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "literal_engine_sha256": hashlib.sha256(Path(__file__).with_name(
            "check_noncommutative_crossing_ward_rank_descent.py").read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()

