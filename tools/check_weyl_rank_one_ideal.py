#!/usr/bin/env python3
"""Regenerate the Weyl curves and verify their finite full-algebra orbit ideal."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command):
    return subprocess.run([str(x) for x in command], check=True,
                          capture_output=True, text=True)


def labels_at(h, variables):
    if variables == 0:
        return [()] if h == 0 else []
    weight = 1 << (variables-1)
    return [low+(a,) for a in range(h//weight+1)
            for low in labels_at(h-a*weight, variables-1)]


def parity_minor(h):
    labels = labels_at(h, h.bit_length())
    selected = [b for b in labels if any(x & 1 for x in b) and b[-1] != 1]
    witnesses = []
    for beta in selected:
        sigma = tuple(2*x for x in beta)
        bits = sorted(((1 << (j+k), j, 1 << k)
                       for j, x in enumerate(sigma) for k in range(x.bit_length())
                       if x >> k & 1), reverse=True)
        left = [0]*len(beta)
        remaining = h
        for weight, j, value in bits:
            if weight <= remaining:
                left[j] += value
                remaining -= weight
        assert remaining == 0
        right = tuple(s-a for s, a in zip(sigma, left))
        assert not any(a & u for a, u in zip(left, right))
        assert not any((a | u) & 1 for a, u in zip(left, right))
        witnesses.append((tuple(left), right))
    for i, (alpha, uu) in enumerate(witnesses):
        for j, beta in enumerate(selected):
            conservation = all(a+u == 2*b for a, u, b in zip(alpha, uu, beta))
            assert conservation == (i == j)
            assert any((a+b) & 1 for a, b in zip(alpha, beta))
    d = len(labels)
    half = len(labels_at(h//2, h.bit_length()))
    assert len(selected) == d-half-1
    n0 = sum(len(labels_at(k, h.bit_length()))**2 for k in range(h//2+1))
    assert n0 <= d*half
    product = (d*d-n0)*len(selected)
    if h == 32:
        assert 17**2 > 8*6**2 and product**6 > d**17
    return {"protected_minor": len(selected), "verified_minor_entries": len(selected)**2,
            "row_rank_bound": d*d-n0, "rank_product_bound": product}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy", type=int, choices=(2, 4, 8, 16, 32), default=32)
    parser.add_argument("--generator-source", type=Path)
    parser.add_argument("--dump-prefix", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = monotonic()
    here = Path(__file__).resolve().parent
    generator = (args.generator_source or here / "explore_weighted_binomial_crossing.cpp").resolve()
    orbit_source = here / "explore_weyl_rank_one_orbit.cpp"
    symmetry_source = here / "explore_weyl_symmetry_orbit.cpp"
    expected = {2: (2, 2, 1, 1, 0), 4: (4, 6, 2, 1, 0),
                8: (10, 26, 6, 6, 0), 16: (36, 214, 74, 77, 0),
                32: (202, 4154, 1912, 1894, 0)}
    symmetry_expected = {2: 1, 4: 2, 8: 10, 16: 99, 32: 2051}
    with tempfile.TemporaryDirectory(prefix="weyl-orbit-") as directory:
        temporary = Path(directory)
        orbit_binary = temporary / "orbit"
        symmetry_binary = temporary / "symmetry-orbit"
        run(["g++", "-O3", "-std=c++17", orbit_source, "-o", orbit_binary])
        run(["g++", "-O3", "-std=c++17", symmetry_source, "-o", symmetry_binary])
        prefix = args.dump_prefix.resolve() if args.dump_prefix else temporary / "weyl"
        generated = None
        if args.dump_prefix is None:
            generator_binary = temporary / "generator"
            run(["g++", "-O3", "-std=c++17", generator, "-o", generator_binary])
            generated = json.loads(run([generator_binary, args.energy, "--dump-prefix", prefix]).stdout)
            assert generated["labels"] == expected[args.energy][0]
        result = json.loads(run([orbit_binary, prefix]).stdout)
        actual = (result["labels"], result["parity_zero_dimension"],
                  result["left_orbit"], result["right_orbit"], result["pairing_rank"])
        assert result["energy"] == args.energy and actual == expected[args.energy]
        result["ideal_dimension"] = result["left_orbit"]*result["right_orbit"]
        result["maximum_ideal_row_rank"] = min(result["left_orbit"], result["right_orbit"])
        result["middle_quotient_dimension"] = (result["parity_zero_dimension"]
                                               -result["left_orbit"]-result["right_orbit"])
        result["square_zero_ideal_verified"] = True
        result["parity_zero_affine_minor"] = parity_minor(args.energy)
        symmetry = json.loads(run([symmetry_binary, prefix]).stdout)
        assert symmetry["energy"] == args.energy
        assert symmetry["transpose_closed_orbit"] == symmetry_expected[args.energy]
        assert symmetry["pairing_rank"] == 0
        symmetry["ideal_dimension"] = symmetry["transpose_closed_orbit"]**2
        symmetry["middle_quotient_dimension"] = (symmetry["parity_zero_dimension"]
                                                 -2*symmetry["transpose_closed_orbit"])
        symmetry["square_zero_ideal_verified"] = True
        result["full_crossing_symmetry_algebra"] = symmetry
        result["source_sha256"] = {"checker": digest(__file__), "generator": digest(generator),
                                    "orbit_engine": digest(orbit_source),
                                    "symmetry_engine": digest(symmetry_source)}
        result["generator_dump_sha256"] = {cut: digest(str(prefix)+"."+cut+".bin")
                                           for cut in ("row", "col")}
        result["fresh_generator_replay"] = generated is not None
        result["scope"] = ("Exact finite orbit closure for P,Q and for P,P^T,Q,Q^T. "
                           "The latter generates all eight crossing-preserving external images. "
                           "The rank-one generated ideals are square-zero for arbitrary words and "
                           "extension-field coefficients at this energy; no all-energy radical theorem "
                           "or score gain is claimed.")
        result["total_seconds"] = round(monotonic()-started, 3)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
