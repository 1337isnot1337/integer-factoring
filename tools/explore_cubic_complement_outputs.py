#!/usr/bin/env python3
"""Cubic dual outputs with an arbitrary nonzero complementary-mode image."""

import argparse
import hashlib
from itertools import product
import json
from math import isqrt, sqrt
from pathlib import Path
import random

from check_double_endpoint_retracts import bit_indices, independent_coefficients
from check_modular_crossing import rank_two
from explore_complement_edge_pencils import minimum_decomposition
from explore_cubic_multirow_outputs import CubicSpace


class ComplementOutputs:
    def __init__(self, pairs):
        self.pairs = pairs
        self.space = CubicSpace(pairs)
        low = self.space.outputs
        high = [tuple(2 if x == 1 else x for x in word) for word in low]
        self.degrees = (0, 1, 2, 4, 12)
        self.parent = [[[[0] * (2 * pairs) for _ in range(pairs)]
                        for _ in range(pairs)] for _ in self.degrees]
        for x, u, y, v in product(range(2 * pairs), range(2 * pairs), range(pairs), range(pairs)):
            value = self.space.kernel((low + high)[x], (low + high)[u], y, v)
            for k, degree in enumerate(self.degrees):
                if value >> degree & 1:
                    self.parent[k][y][v][x] ^= 1 << u

    def matrices(self, complement):
        d = self.pairs
        output = [(1 << i) | (complement[i] << d) for i in range(d)]
        supports = list(map(bit_indices, output))
        rows = [[0] * (d * d) for _ in self.degrees]
        columns = [[0] * (d * d) for _ in self.degrees]
        for k, y, v in product(range(5), range(d), range(d)):
            parent = self.parent[k][y][v]
            for x in range(d):
                value = 0
                for index in supports[x]:
                    value ^= parent[index]
                for u in range(d):
                    if (value & output[u]).bit_count() & 1:
                        rows[k][x * d + y] ^= 1 << (v * d + u)
                        columns[k][x * d + u] ^= 1 << (v * d + y)
        identity = [1 << i for i in range(d * d)]
        assert rows[0] == columns[0] == identity
        # Full-entry checks for the two fixed physical identities after applying [I C].
        gram = [sum((((complement[x] & complement[u]).bit_count() & 1) ^ int(x == u)) << u
                    for u in range(d)) for x in range(d)]
        for x, y, v, u in product(range(d), repeat=4):
            high = ((complement[x] >> y) & 1) & ((complement[u] >> v) & 1)
            assert ((rows[4][x * d + y] >> (v * d + u)) & 1) == high
            wanted = (int(y == v) & ((gram[x] >> u) & 1)) ^ high
            assert (((rows[1][x * d + y] ^ rows[3][x * d + y]) >> (v * d + u)) & 1) == wanted
        return rows, columns

    def profile(self, complement):
        rows, columns = self.matrices(complement)
        chosen = independent_coefficients([[matrix] for matrix in rows])
        rows, columns = [[cut[i] for i in chosen] for cut in (rows, columns)]
        work = [[0] * (self.pairs * self.pairs) for _ in range(2)]
        ranks = [[0, 0] for _ in range(1 << len(chosen))]
        previous = 0
        for step in range(1, 1 << len(chosen)):
            mask = step ^ (step >> 1)
            changed = (mask ^ previous).bit_length() - 1
            for target, source in zip(work, (rows[changed], columns[changed])):
                for i, value in enumerate(source):
                    target[i] ^= value
            ranks[mask] = [rank_two(matrix) for matrix in work]
            previous = mask
        upper, witness = minimum_decomposition([sqrt(a * b) for a, b in ranks], target=1)
        lower, _ = minimum_decomposition([isqrt(a * b) for a, b in ranks], target=1)
        return {"complement_rows": complement, "complement_rank": rank_two(complement),
                "basis_degrees": [self.degrees[i] for i in chosen], "rank_table": ranks,
                "found_decomposition_cost": upper, "witness": witness,
                "integer_shortest_lower_bound": lower,
                "binary_seed_criterion_excluded": lower ** 70 > self.pairs ** 99}


def candidates(d, seed, count):
    for mask in range(1 << d):
        yield "diagonal", [((mask >> i) & 1) << i for i in range(d)]
        yield "circulant", [sum(((mask >> ((j - i) % d)) & 1) << j for j in range(d))
                            for i in range(d)]
        if mask:
            yield "symmetric_rank_one", [mask if mask >> i & 1 else 0 for i in range(d)]
    for mask in range(1 << (d - 1)):
        yield "nilpotent_toeplitz", [sum(((mask >> (j - i - 1)) & 1) << j for j in range(i + 1, d))
                                    for i in range(d)]
    rng = random.Random(seed)
    for case in range(count):
        if case % 3 == 0:
            u, v = rng.randrange(1, 1 << d), rng.randrange(1, 1 << d)
            yield "random_rank_one", [v if u >> i & 1 else 0 for i in range(d)]
        elif case % 3 == 1:
            yield "random_dense", [rng.randrange(1 << d) for _ in range(d)]
        else:
            yield "random_strict_upper", [rng.randrange(1 << d) & ~((1 << (i + 1)) - 1)
                                          for i in range(d)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=7)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--random-cases", type=int, default=768)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(args.pairs)
    seen = set()
    result = {"pairs": args.pairs, "labels": args.pairs, "seed": args.seed,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "dependency_sha256": hashlib.sha256(Path(__file__).with_name("explore_cubic_multirow_outputs.py").read_bytes()).hexdigest(),
              "cases": [], "complete": False}
    best = None
    for description, complement in candidates(args.pairs, args.seed, args.random_cases):
        key = tuple(complement)
        if key in seen:
            continue
        seen.add(key)
        record = dict(family.profile(complement), family=description)
        result["cases"].append(record)
        if best is None or record["found_decomposition_cost"] < best["found_decomposition_cost"]:
            best = record
            print(json.dumps({"new_best_cost": best["found_decomposition_cost"],
                              "dimension": len(best["basis_degrees"]), "family": description,
                              "complement_rows": complement}), flush=True)
        if len(result["cases"]) % 50 == 0:
            print(json.dumps({"cases": len(result["cases"]),
                              "best_cost": best["found_decomposition_cost"]}), flush=True)
            args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    result["best"] = best
    result["complete"] = True
    result["all_recorded_spaces_excluded_over_F2"] = all(r["binary_seed_criterion_excluded"] for r in result["cases"])
    result["scope"] = "The recorded binary complement matrices; arbitrary complements and extension fields remain open."
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"complete": True, "cases": len(result["cases"]),
                      "best_cost": best["found_decomposition_cost"],
                      "all_recorded_spaces_excluded_over_F2": result["all_recorded_spaces_excluded_over_F2"]}), flush=True)


if __name__ == "__main__":
    main()
