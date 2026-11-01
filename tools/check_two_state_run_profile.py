#!/usr/bin/env python3
"""Check the run-compressed engine against independent tensor and kernel ranks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess
from tempfile import TemporaryDirectory

from check_markov_pairing_transfer import direct_tensors, ranks_of_coefficients


def lift(vertex, edge):
    a, b = vertex
    ss, sc, cs, cc = edge
    return [a, a, b], [ss, ss, sc, ss, ss, sc, cs, cs, cc]


def invoke(executable, directory, t, vertex, edge):
    path = directory / "profile.json"
    subprocess.run([str(executable), "profile", str(t), str(path),
                    *map(str, vertex + edge)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    rng = random.Random(20261030)
    cases = [([0, 5], [0, 1, 1, 0]), ([0, 2], [0, 0, 0, 1])]
    cases += [([rng.randrange(3), rng.randrange(3)],
               [rng.randrange(3) for _ in range(4)]) for _ in range(10)]
    profile_path = root.parent / "evidence" / "finite-state-pairing-screen.json"
    if not profile_path.exists():
        profile_path = root / "finite-state-pairing-screen.json"
    profiles = json.loads(profile_path.read_text())["selected_larger_profiles"]
    with TemporaryDirectory(prefix="two-state-ranks-") as temporary:
        directory = Path(temporary)
        small, compressed = directory / "small", directory / "compressed"
        for source, executable in (("markov_pairing_search.cpp", small),
                                   ("two_state_run_profile.cpp", compressed)):
            subprocess.run(["g++", "-O2", "-std=c++17", str(root / source), "-o", str(executable)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        result = subprocess.run([str(compressed), "selftest", "4"], check=True,
                                text=True, stdout=subprocess.PIPE)
        exhaustive = json.loads(result.stdout)["polynomials"]
        for vertex, edge in cases:
            triple, matrix = lift(vertex, edge)
            result = invoke(compressed, directory, 3, vertex, edge)
            tensors = direct_tensors(3, triple, matrix, "profile")
            rows, columns = ranks_of_coefficients(3, tensors, len(result["row_ranks"]))
            assert rows == result["row_ranks"] and columns == result["column_ranks"]
            result = invoke(compressed, directory, 6, vertex, edge)
            independent = invoke(small, directory, 6, triple, matrix)
            assert independent["row_ranks"] == result["row_ranks"]
            assert independent["column_ranks"] == result["column_ranks"]
        larger = [
            (8, [0, 5], [0, 1, 1, 0], "markov-domainwall-cost-t8.json"),
            (12, [0, 5], [0, 1, 1, 0], "markov-domainwall-cost-t12.json"),
            (12, [0, 7], [0, 1, 0, 0], "markov-run-count-bigrading-t12.json"),
        ]
        for t, vertex, edge, name in larger:
            result = invoke(compressed, directory, t, vertex, edge)
            independent = profiles[name]
            assert independent["row_ranks"] == result["row_ranks"]
            assert independent["column_ranks"] == result["column_ranks"]
    evidence = {
        "status": "passed", "exhaustive_four_variable_polynomials": exhaustive,
        "independent_full_tensor_cases": len(cases),
        "independent_six_position_kernel_cases": len(cases),
        "stored_independent_larger_profiles": len(larger),
        "maximum_independent_comparison_length": 12,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((root / "two_state_run_profile.cpp").read_bytes()).hexdigest(),
        "scope": "Exact compression and rank checks; no asymptotic contraction claim.",
    }
    encoded = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
