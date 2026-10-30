#!/usr/bin/env python3
"""Audit per-leg Hopf gauges by generic exterior-algebra composition."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from check_graded_hopf_crossings import matrix_digest
from check_noncommutative_crossing_orbit import polynomial_mul
from explore_graded_hopf_crossings import RestrictedHeisenberg


def multiply(left, right):
    result = {}
    for a, p in left.items():
        for b, q in right.items():
            if not a & b:
                monomial = a | b
                result[monomial] = result.get(monomial, 0) ^ polynomial_mul(p, q)
    return {m: p for m, p in result.items() if p}


def add_scaled(target, source, coefficient):
    for monomial, p in source.items():
        target[monomial] = target.get(monomial, 0) ^ polynomial_mul(p, coefficient)
    return {m: p for m, p in target.items() if p}


def gauged_images(algebra, code):
    choices = (0, 3, 5, 6)  # 0,t,t²,t+t² after t=1+z.
    gauges = [choices[(code >> (2 * j)) & 3] for j in range(4)]
    width = algebra.width
    original = algebra.generator_images()
    output_images = [dict(image) for image in original]
    input_images = [{1 << i: 1} for i in range(2 * width)]
    for side in range(2):
        for center in range(algebra.t):
            index = side * width + 2 * algebra.m + center
            for pair, column in enumerate(algebra.columns):
                if not (column >> center) & 1:
                    continue
                x = side * width + pair
                y = side * width + algebra.m + pair
                output_images[index] = add_scaled(output_images[index],
                    multiply(original[x], original[y]), gauges[side])
                input_images[index] = add_scaled(input_images[index],
                    {(1 << x) | (1 << y): 1}, gauges[side + 2])
    result = []
    for polynomial in output_images:
        composed = {}
        for monomial, coefficient in polynomial.items():
            product = {0: coefficient}
            while monomial:
                bit = monomial & -monomial
                product = multiply(product, input_images[bit.bit_length() - 1])
                monomial ^= bit
            composed = add_scaled(composed, product, 1)
        result.append(composed)
    return result


def verify_case(binary, columns, center, codes):
    algebra = RestrictedHeisenberg(len(columns), center, tuple(columns))
    all_events = [json.loads(line) for line in subprocess.check_output(
        [str(binary.resolve()), "--horizontal", str(len(columns)), "--center", str(center),
         "--columns", ",".join(map(str, columns)), "--sweep"], text=True).splitlines()]
    native = {r["gauge_code"]: r for r in all_events if r["event"] == "gauge_profile"}
    assert len(native) == 256
    checked = []
    baseline_constant = native[0]["profile"][0][1:]
    for record in native.values():
        assert record["profile"][0][1:] == baseline_constant
    for code in codes:
        matrices = algebra.substitution_matrices(gauged_images(algebra, code))
        for row in range(algebra.d * algebra.d):
            x, y = divmod(row, algebra.d)
            total = 0
            for matrix in matrices:
                total ^= matrix[row]
            assert total == 1 << (y * algebra.d + x)
        expected = algebra.rank_profile(matrices)
        actual = native[code]
        assert [[i, r, c] for i, r, c, _ in expected["coefficient_profiles"]] == actual["profile"]
        digest = matrix_digest(algebra, matrices)
        assert all(actual[key] == value for key, value in digest.items())
        checked.append({"gauge_code": code, "profile": actual["profile"], **digest})
    best_product = min(native.values(), key=lambda x: (x["maximum_rank_product"], x["gauge_code"]))
    best_sum = min(native.values(), key=lambda x: (x["rank_sum_approx"], x["gauge_code"]))
    return {"columns": list(columns), "center_bits": center, "labels": algebra.d,
            "native_gauges": 256, "independent_generic_compositions": checked,
            "best_product_gauge": best_product["gauge_code"],
            "best_product": best_product["maximum_rank_product"],
            "best_rank_sum_gauge": best_sum["gauge_code"],
            "rank_sum_comparison_is_numerical": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--all-small", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    selected = (0, 1, 4, 16, 64, 85, 93, 170, 255)
    cases = [verify_case(args.binary, (1,), 1, range(256) if args.all_small else selected),
             verify_case(args.binary, (1, 1), 1, selected)]
    root = Path(__file__).parent
    result = {"status": "pass", "scope": "endpoint-preserving per-leg gauges; independent complete tensors and exact native cuts; no exponent improvement",
              "cases": cases,
              "source_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "explore_graded_hopf_crossings.py",
                                             "explore_graded_hopf_crossings.cpp", "explore_graded_hopf_gauges.cpp")}}
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
