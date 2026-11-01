#!/usr/bin/env python3
"""Reproduce bounded dense-grade and independent-energy Gram screens."""

import argparse
import hashlib
import json
from pathlib import Path

from folded_gram_aggregate import gaussian_distributions, shell_dimension
from unequal_grade_gram import score_case


def first_state_energy(modulus, states, wanted):
    """Lowest homogeneous energy whose shell has at least wanted labels."""
    distribution = gaussian_distributions(modulus // 2)[states]
    for index, count in enumerate(distribution):
        if count >= wanted:
            return states ** 2 + 2 * index
    raise ValueError("the state grade has no sufficiently large shell")


def central_cases():
    for modulus, lengths in ((64, (8, 12, 16)), (128, (8, 16, 24, 32)),
                             (256, (16, 32, 48, 64))):
        modes = modulus // 2
        for length in lengths:
            energy = length * modulus // 2
            dimension = shell_dimension(modulus, length, energy)
            for states in sorted({length + (modes - 2 * length) * i // 4
                                  for i in range(5)}):
                yield modulus, length, states, energy, first_state_energy(modulus, states, dimension)


def joint_energy_cases():
    for modulus, length, grades in ((64, 8, (8, 10, 12, 14, 16)),
                                    (128, 16, (16, 20, 24, 28, 32)),
                                    (256, 32, (32, 40, 48, 56, 64))):
        for fraction in range(1, 5):
            energy = length ** 2 + (length * modulus // 2 - length ** 2) * fraction // 4
            energy -= (energy - length) % 2
            dimension = shell_dimension(modulus, length, energy)
            for states in grades:
                first = first_state_energy(modulus, states, dimension)
                central = states * modulus // 2
                central -= (central - states) % 2
                for state_energy in sorted({first, first + (central - first) // 4 * 2, central}):
                    yield modulus, length, states, energy, state_energy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for name, cases in (("central_operator", central_cases()),
                        ("joint_energy", joint_energy_cases())):
        group = []
        for case in cases:
            record = score_case(*case)
            record.pop("complete_uncapped_channels")
            group.append(record)
        records.append({"grid": name, "cases": group,
                        "minimum_rho_case": min(group, key=lambda row: row["rho_of_upper_bound"])})
        print(json.dumps({"grid": name, "case_count": len(group),
                          "minimum_rho_case": records[-1]["minimum_rho_case"]}), flush=True)
    report = {
        "status": "pass",
        "scope": "bounded unequal-grade ambient Gram upper-bound screen; no improved factoring exponent",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "module_sha256": hashlib.sha256(Path(__file__).with_name("unequal_grade_gram.py").read_bytes()).hexdigest(),
        "grids": records,
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
