#!/usr/bin/env python3
"""Finite exact profiles of ternary Bell-diagonal linear pencils.

Exhausts the specified small alphabets and/or runs a seeded local search.
Every retained tensor has s(z)=1+(z-1)v and the exact crossing at z=1.
Search output is finite evidence, never an all-length factoring rate.
"""

import argparse
from collections import Counter
from itertools import product
import json
from math import log2
import random
import time

from bell_crossing import linear_profile_ternary, sign, transform


def counts(values: tuple[int, ...]) -> tuple[int, int, int]:
    histogram = Counter(values)
    return tuple(histogram[i] for i in range(3))


def profile_record(values: tuple[int, ...], bits: int, lengths: list[int]) -> dict:
    other = transform(values, bits, 3)
    row_counts, col_counts = counts(values), counts(other)
    profiles = []
    for length in lengths:
        row = linear_profile_ternary(row_counts, length)
        col = linear_profile_ternary(col_counts, length)
        maximum = max(r * c for r, c in zip(row, col))
        profiles.append({"length": length, "row_ranks": row, "column_ranks": col,
                         "maximum_rank_product": maximum,
                         "finite_rho_without_slack": log2(maximum) / (4 * bits * length)})
    return {"values": values, "transformed_values": other,
            "row_value_counts": row_counts, "column_value_counts": col_counts,
            "profiles": profiles}


def run(bits: int, length: int, exhaustive: str | None,
        steps: int, seed: int) -> dict:
    start = time.monotonic()
    d, q = 1 << bits, 1 << (2 * bits)
    rng = random.Random(seed)
    signs = [[pow(d, -1, 3) * sign((a ^ c) & (b ^ e)) % 3
              for a in range(d) for b in range(d)]
             for c in range(d) for e in range(d)]
    catalogue = {}
    best_key = None
    best_values = None
    evaluations = 0

    def inspect(values: tuple[int, ...], other: tuple[int, ...]) -> tuple[int, int]:
        nonlocal evaluations, best_key, best_values
        evaluations += 1
        rc, cc = counts(values), counts(other)
        key = tuple(sorted((rc, cc)))
        if key not in catalogue:
            row = linear_profile_ternary(rc, length)
            col = linear_profile_ternary(cc, length)
            catalogue[key] = max(r * c for r, c in zip(row, col))
        # Prefer lower endpoint pressure when finite products tie.
        objective = (catalogue[key], (q - rc[1]) * (q - cc[1])
                     + (q - rc[0]) * (q - cc[0]))
        if best_key is None or objective < best_key:
            best_key, best_values = objective, values
        return objective

    if exhaustive:
        alphabet = range(2 if exhaustive == "binary" else 3)
        if q > 16 or exhaustive == "ternary" and q > 4:
            raise ValueError("exhaustive schedule is limited to 2^16 or 3^4")
        for values in product(alphabet, repeat=q):
            inspect(values, transform(values, bits, 3))

    # Sparse, complementary, affine and quadratic initial configurations.
    starts = [tuple(0 for _ in range(q)), tuple(1 for _ in range(q))]
    for scalar in (1, 2):
        for index in range(q):
            starts.append(tuple(scalar if i == index else 0 for i in range(q)))
            starts.append(tuple((1 - scalar) % 3 if i == index else 1 for i in range(q)))
    starts += [tuple(int(a == b) for a in range(d) for b in range(d)),
               tuple(int(a == 0) for a in range(d) for _ in range(d)),
               tuple(int(sign(a & b) == -1) for a in range(d) for b in range(d))]
    for values in ([] if exhaustive and not steps else starts):
        inspect(values, transform(values, bits, 3))

    values = best_values
    other = transform(values, bits, 3)
    current = inspect(values, other)
    accepted = 0
    for step in range(steps):
        if step % 512 == 0:
            if step % 2048 == 0:
                values = tuple(rng.randrange(3) for _ in range(q))
            else:
                values = best_values
            other = transform(values, bits, 3)
            current = inspect(values, other)
        index = rng.randrange(q)
        delta = rng.choice((1, 2))
        proposal = list(values)
        proposal[index] = (proposal[index] + delta) % 3
        proposal = tuple(proposal)
        next_other = tuple((x + delta * row[index]) % 3 for x, row in zip(other, signs))
        objective = inspect(proposal, next_other)
        temperature = 0.12 * (1 - (step % 512) / 512) + 0.002
        worsening = log2(objective[0]) - log2(current[0])
        if objective <= current or rng.random() < 2 ** (-worsening / temperature):
            values, other, current = proposal, next_other, objective
            accepted += 1
    assert transform(best_values, bits, 3) == tuple(
        sum(row[j] * best_values[j] for j in range(q)) % 3 for row in signs)
    lengths = sorted(set((1, 2, 4, 8, 12, 16, 24, 32, length)))
    return {"status": "finite search completed", "bits_per_local_label": bits,
            "local_dimension": d, "coefficient_prime": 3, "objective_length": length,
            "exhaustive_alphabet": exhaustive, "seed": seed, "local_steps": steps,
            "candidate_evaluations": evaluations, "distinct_histogram_pairs": len(catalogue),
            "accepted_local_steps": accepted, "best": profile_record(best_values, bits, lengths),
            "elapsed_seconds": time.monotonic() - start,
            "scope": "Specified finite schedules only; no asymptotic contraction proved"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=3)
    parser.add_argument("--length", type=int, default=12)
    parser.add_argument("--exhaustive", choices=("binary", "ternary"))
    parser.add_argument("--steps", type=int, default=8192)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = run(args.bits, args.length, args.exhaustive, args.steps, args.seed)
    with open(args.output, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in (
        "status", "local_dimension", "candidate_evaluations", "distinct_histogram_pairs",
        "elapsed_seconds")}))
    print(json.dumps({"row_value_counts": result["best"]["row_value_counts"],
                      "column_value_counts": result["best"]["column_value_counts"],
                      "profiles": [{"length": x["length"],
                                    "rho": x["finite_rho_without_slack"]}
                                   for x in result["best"]["profiles"]]}))


if __name__ == "__main__":
    main()
