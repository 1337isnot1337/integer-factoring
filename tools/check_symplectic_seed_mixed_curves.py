#!/usr/bin/env python3
"""Complete binary recombinations for two different endpoint-preserving curves."""

import argparse
import hashlib
from itertools import combinations, product
import json
from math import log2
from pathlib import Path

from check_modular_crossing import rank_two
from check_symplectic_seed_tensor_modules import add, polynomial_tensor, shortest_cost
from symplectic_seed_weighted_modules import (
    CERTIFICATE, load_certificate, specialize_types, verify_certificate,
)


def mixed_rank_table(first, second):
    coordinates = len(first[0]["generator_rows"]) + len(second[0]["generator_rows"]) - 1
    ranks = [0] * (1 << coordinates)
    for left in first:
        for right in second:
            coefficients = polynomial_tensor(left["generator_rows"], right["generator_rows"])
            multiplicity = left["multiplicity"] * right["multiplicity"]
            rows, previous = [0] * len(coefficients[0]), 0
            for step in range(1, 1 << coordinates):
                mask = step ^ (step >> 1)
                changed = (mask ^ previous).bit_length() - 1
                rows = add(rows, coefficients[changed])
                previous = mask
                ranks[mask] += multiplicity * rank_two(rows)
    return ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    certificate = load_certificate()
    verification = verify_certificate(certificate)
    parameters = list(product((0, 1), repeat=3))
    specialized = {
        cut: [specialize_types(types, 1 | (a << 1), 2, (b << 1) | (c << 2))
              for a, b, c in parameters]
        for cut, types in certificate["types"].items()}
    result = {
        "curve": "K(1+a*t,t,b*t+c*t^2)", "labels": 676,
        "source_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in (Path(__file__).name, "symplectic_seed_weighted_modules.py",
                         "check_symplectic_seed_tensor_modules.py", "check_modular_crossing.py")},
        "certificate_sha256": hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest(),
        "certificate_verification": verification, "cases": []}
    for first, second in combinations(range(len(parameters)), 2):
        tables = {cut: mixed_rank_table(types[first], types[second])
                  for cut, types in specialized.items()}
        row, column = tables["row"], tables["column"]
        assert row[1] == column[1] == 676 ** 2
        lower = shortest_cost(row, column, "floor")
        upper = shortest_cost(row, column, "ceiling")
        numerical = shortest_cost(row, column, "float")
        assert lower["cost"] ** 2 > 676 ** 3
        record = {"parameters": [parameters[first], parameters[second]],
                  "rank_tables": tables, "integer_lower": lower, "integer_upper": upper,
                  "numerical_optimum": numerical,
                  "numerical_rho": log2(numerical["cost"]) / (4 * log2(26)),
                  "integer_lower_squared_exceeds_labels_cubed": True}
        result["cases"].append(record)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"parameters": record["parameters"],
                          "integer_interval": [lower["cost"], upper["cost"]],
                          "rho": record["numerical_rho"]}), flush=True)


if __name__ == "__main__":
    main()
