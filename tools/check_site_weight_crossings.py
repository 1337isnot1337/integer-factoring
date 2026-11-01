#!/usr/bin/env python3
"""Independent endpoint and full-cut rank checks for bivariate site weights."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import random
import subprocess
from tempfile import TemporaryDirectory

from check_markov_pairing_transfer import ranks_of_coefficients, xor_all


def direct_tensors(weights):
    t = len(weights)
    dimension = 1 << t
    support_c, support_s = [], []
    for q in product(range(2), repeat=4):
        x, y, u, v = q
        c = x == u and y == v
        e = x == v and y == u
        if c:
            support_c.append(q)
        if c != e:
            support_s.append(q)
    tensors = {}
    for word in range(1 << t):
        key = (word.bit_count(), sum(weights[i] for i in range(t) if word >> i & 1))
        tensor = 0
        choices = [support_c if word >> i & 1 else support_s for i in range(t)]
        for quads in product(*choices):
            x, y, u, v = (sum(quads[i][leg] << i for i in range(t)) for leg in range(4))
            index = ((x * dimension + y) * dimension + u) * dimension + v
            tensor ^= 1 << index
        tensors[key] = tensors.get(key, 0) ^ tensor
    crossing = sum(1 << (((x * dimension + y) * dimension + y) * dimension + x)
                   for x, y in product(range(dimension), repeat=2))
    assert xor_all(tensors.values()) == crossing
    return tensors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    rng = random.Random(20261030)
    cases = [[0, 0, 0], [1, 1, 1], [1, 4, 9], [0, 1, 3]]
    cases += [[rng.randrange(8) for _ in range(3)] for _ in range(20)]
    with TemporaryDirectory(prefix="site-weight-ranks-") as temporary:
        directory = Path(temporary)
        binary = directory / "site_weights"
        subprocess.run(["g++", "-O2", "-std=c++17", str(root / "site_weight_crossings.cpp"),
                        "-o", str(binary)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        for weights in cases:
            path = directory / "current.json"
            subprocess.run([str(binary), "profile", str(len(weights)), str(path), *map(str, weights)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            native = json.loads(path.read_text())
            tensors = direct_tensors(weights)
            keys = [(k, w) for k, w, _, _ in native["pieces"]]
            assert keys == sorted(tensors)
            rows, columns = ranks_of_coefficients(3, {i: tensors[key] for i, key in enumerate(keys)}, len(keys))
            assert rows == [entry[2] for entry in native["pieces"]]
            assert columns == [entry[3] for entry in native["pieces"]]
        path = directory / "uniform.json"
        subprocess.run([str(binary), "profile", "8", str(path), *(["1"] * 8)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        uniform = json.loads(path.read_text())
        evidence_path = root.parent / "evidence" / "finite-state-pairing-screen.json"
        if not evidence_path.exists():
            evidence_path = root / "finite-state-pairing-screen.json"
        baseline = json.loads(evidence_path.read_text())["selected_larger_profiles"]["markov-two-pair-baseline-t8.json"]
        assert [x[2] for x in uniform["pieces"]] == baseline["row_ranks"]
        assert [x[3] for x in uniform["pieces"]] == baseline["column_ranks"]
    # Central symmetry is needed for the full-span exclusion. Test arbitrary
    # coefficient recombinations, not just the individual graded pieces.
    projection_cases = 0
    for weights in ([0, 1, 2, 3], [0, 1, 4, 5], [0, 0, 3, 3]):
        tensors = direct_tensors(weights)
        for _ in range(8):
            chosen = [key for key in tensors if rng.randrange(2)]
            tensor = xor_all(tensors[key] for key in chosen)
            central = xor_all(tensors[key] for key in chosen if key[0] == 2)
            rows, columns = ranks_of_coefficients(4, {0: tensor, 1: central}, 2)
            assert rows[1] == columns[1]
            assert rows[1] <= 5 * rows[0] and columns[1] <= 5 * columns[0]
            projection_cases += 1
    asymmetric = direct_tensors([0, 1, 1, 1])
    central = {i: tensor for i, (key, tensor) in enumerate(sorted(asymmetric.items())) if key[0] == 2}
    # Use a dense index list for the rank helper.
    central = dict(enumerate(central.values()))
    rows, columns = ranks_of_coefficients(4, central, len(central))
    assert any(r != c for r, c in zip(rows, columns))
    evidence = {
        "status": "passed", "field": "F_2",
        "independent_full_tensor_cases": len(cases),
        "independent_eight_position_uniform_profile": 1,
        "symmetric_weight_recombination_checks": projection_cases,
        "asymmetric_weight_negative_control": [0, 1, 1, 1],
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((root / "site_weight_crossings.cpp").read_bytes()).hexdigest(),
        "scope": "Exact finite coefficient identities and ranks; no asymptotic contraction claim.",
    }
    encoded = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
