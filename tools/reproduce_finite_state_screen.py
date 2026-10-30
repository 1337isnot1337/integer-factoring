#!/usr/bin/env python3
"""Rebuild the finite-state screen and every required larger profile outside the repo."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True,
                        help="Directory for temporary native binaries and complete raw profiles.")
    parser.add_argument("--output", type=Path,
                        help="Compact evidence destination; defaults inside --directory.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    small, compressed = directory / "markov_search", directory / "two_state_profile"
    for source, executable in (("markov_pairing_search.cpp", small),
                               ("two_state_run_profile.cpp", compressed)):
        subprocess.run(["g++", "-O3", "-std=c++17", str(root / source), "-o", str(executable)],
                       check=True)

    screens = [
        ("screen", 5, "markov-cost-screen-t5.jsonl", [2]),
        ("screen", 6, "markov-cost-screen-t6.jsonl", [5]),
        ("polyscreen", 6, "markov-polynomial-screen-012-t6.jsonl", [0, 1, 2]),
        ("polyscreen", 6, "markov-polynomial-screen-001-t6.jsonl", [0, 0, 1]),
        *[("auxscreen", 6, f"auxiliary-matrix-screen-{a}-b{b}-t6.jsonl",
           [*map(int, a), b]) for a in ("001", "012") for b in (0, 1)],
    ]
    for mode, t, name, parameters in screens:
        print(f"Rebuilding {name}", flush=True)
        subprocess.run([str(small), mode, str(t), str(directory / name),
                        *map(str, parameters)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # The compressed engine reproduces exactly the same integer ranks as
    # the original XOR-circulant engine and avoids its much larger storage.
    profiles = [
        (8, "markov-run-cost-t8.json", [0, 2, 0, 0, 0, 1]),
        (10, "markov-run-cost-t10.json", [0, 2, 0, 0, 0, 1]),
        (8, "markov-domainwall-cost-t8.json", [0, 5, 0, 1, 1, 0]),
        (10, "markov-domainwall-cost-t10.json", [0, 5, 0, 1, 1, 0]),
        (12, "markov-domainwall-cost-t12.json", [0, 5, 0, 1, 1, 0]),
        (10, "markov-run-count-bigrading-t10.json", [0, 6, 0, 1, 0, 0]),
        (12, "markov-run-count-bigrading-t12.json", [0, 7, 0, 1, 0, 0]),
        (8, "markov-two-pair-baseline-t8.json", [0, 1, 0, 0, 0, 0]),
    ]
    for t, name, parameters in profiles:
        print(f"Rebuilding {name}", flush=True)
        subprocess.run([str(compressed), "profile", str(t), str(directory / name),
                        *map(str, parameters)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    command = [sys.executable, str(root / "audit_finite_state_screen.py"),
               "--directory", str(directory)]
    if args.output:
        command += ["--output", str(args.output.resolve())]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
