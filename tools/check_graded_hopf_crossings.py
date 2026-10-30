#!/usr/bin/env python3
"""Independently audit native graded-Hopf profiles against literal Hopf maps."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from explore_graded_hopf_crossings import RestrictedHeisenberg


MASK64 = (1 << 64) - 1


def mix64(value):
    value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
    return value ^ (value >> 31)


def matrix_digest(algebra, matrices):
    xor, total, count = 0, 0, 0
    d = algebra.d
    for row in range(d * d):
        entries = {}
        for degree, matrix in enumerate(matrices):
            bits = matrix[row]
            while bits:
                bit = bits & -bits
                column = bit.bit_length() - 1
                entries[column] = entries.get(column, 0) ^ (1 << degree)
                bits ^= bit
        x, y = divmod(row, d)
        output = x | (y << algebra.width)
        for column, polynomial in entries.items():
            u, v = divmod(column, d)
            input = u | (v << algebra.width)
            token = mix64(mix64(output ^ (input << 32)) ^ polynomial)
            xor ^= token
            total = (total + token) & MASK64
            count += 1
    return {"matrix_digest_xor": str(xor), "matrix_digest_sum": str(total),
            "nonzero_polynomial_entries": count}


def check_case(binary, columns, center):
    algebra = RestrictedHeisenberg(len(columns), center, tuple(columns))
    algebra.check_algebra()
    literal = algebra.coefficient_matrices()
    substitution = algebra.substitution_matrices()
    assert literal == substitution
    profile = algebra.rank_profile(literal)
    command = [str(binary.resolve()), "--horizontal", str(len(columns)), "--center", str(center),
               "--columns", ",".join(map(str, columns))]
    events = [json.loads(line) for line in subprocess.check_output(command, text=True).splitlines()]
    assert events[-1]["event"] == "complete"
    native = [[r["degree"], r["row_rank"], r["column_rank"], r["product"]]
              for r in events if r["event"] == "coefficient" and r["row_rank"]]
    assert native == profile["coefficient_profiles"]
    digest = matrix_digest(algebra, literal)
    assert all(events[-1][name] == value for name, value in digest.items())
    return {**profile, **digest,
            "literal_and_substitution_matrices_equal": True,
            "native_complete_profile_and_digest_equal": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [check_case(args.binary, (1,), 1), check_case(args.binary, (1, 1), 1),
             check_case(args.binary, (1, 2), 2), check_case(args.binary, (1, 3), 2)]
    # A central basis change is a graded Hopf isomorphism. It preserves both
    # ranks of every coefficient, although the individual matrix entries differ.
    assert cases[2]["coefficient_profiles"] == cases[3]["coefficient_profiles"]
    root = Path(__file__).parent
    result = {"status": "pass", "scope": "literal Hopf formula, polynomial substitution, both native cuts, and central basis naturality; no stronger factoring exponent",
              "cases": cases,
              "source_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "explore_graded_hopf_crossings.py",
                                             "explore_graded_hopf_crossings.cpp")}}
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
