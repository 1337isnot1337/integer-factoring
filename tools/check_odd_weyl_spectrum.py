#!/usr/bin/env python3
"""Exact characteristic-three component and constant-column certificate."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic

from odd_weyl_projectors import make_algebra


def labels_at(n, weights):
    if len(weights) == 1:
        return [(n,)]
    return [a+(k,) for k in range(n//weights[-1]+1)
            for a in labels_at(n-k*weights[-1], weights[:-1])]


def digit_disjoint(a, b, p):
    while a or b:
        if a % p and b % p:
            return False
        a //= p
        b //= p
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy", type=int, choices=(27, 81), default=81)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = monotonic()
    here = Path(__file__).resolve().parent
    algebra = make_algebra(3, args.energy)
    labels = labels_at(algebra["reduced_energy"], algebra["weights"])
    witnesses = {}
    eligible = 0
    for a in labels:
        for b in labels:
            if all(digit_disjoint(x, y, 3) for x, y in zip(a, b)):
                witnesses.setdefault(tuple(x+y for x, y in zip(a, b)), (a, b))
                eligible += 1
    expected = {27: (23, 399, 93, 88), 81: (239, 33933, 1632, 1617)}
    with tempfile.TemporaryDirectory(prefix="odd-weyl-") as directory:
        temporary = Path(directory)
        binary = temporary / "rank"
        command = ["g++", "-O3", "-std=c++17", str(here / "explore_odd_weyl.cpp"), "-o", str(binary)]
        subprocess.run(command, check=True, capture_output=True, text=True)
        mask_file = temporary / "mask.txt"
        raw0 = next(value for name, value in algebra["masks"] if name == "raw_0")
        mask_file.write_text("1\nraw_0 "+" ".join(map(str, raw0))+"\n")
        completed = subprocess.run([str(binary), "3", str(args.energy), str(mask_file)],
                                   check=True, capture_output=True, text=True)
        ranks = json.loads(completed.stdout)
    row, col = ranks["pieces"][0]["row"], ranks["pieces"][0]["col"]
    assert (len(labels), row, col, len(witnesses)) == expected[args.energy]
    product = row*len(witnesses)
    assert product > len(labels)**3
    result = {"characteristic": 3, "energy": args.energy,
              "reduced_energy": algebra["reduced_energy"], "weights": algebra["weights"],
              "labels": len(labels), "abstract_dimension": algebra["abstract_dimension"],
              "basis_associativity_checks": algebra["basis_associativity_checks"],
              "spectral_projectors": algebra["projector_count"],
              "frobenius_exponent": algebra["frobenius_exponent"],
              "zero_character_projector": algebra["projectors"][0],
              "local_unit_and_nilpotence_verified": True,
              "constant_operator_ranks": {"row": row, "column": col},
              "constant_one_column_blocks": len(witnesses), "eligible_digit_disjoint_pairs": eligible,
              "forced_rank_product": product, "strict_target_inequality": "product > d^3 > d^(2 sqrt(2))",
              "native_rank_evidence": ranks,
              "source_sha256": {name: hashlib.sha256((here/name).read_bytes()).hexdigest()
                                 for name in (Path(__file__).name, "odd_weyl_projectors.py", "explore_odd_weyl.cpp")},
              "scope": "All decompositions from this original coefficient span over every characteristic-three extension field; opposite symmetries and mixed products are not covered.",
              "elapsed_seconds": round(monotonic()-started, 3)}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
