#!/usr/bin/env python3
"""Exact binary searches for simultaneous off-syndrome cubic output changes."""

import argparse
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from math import isqrt, sqrt
from pathlib import Path
import random

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local, minimum_decomposition
from check_double_endpoint_retracts import independent_coefficients


class CubicSpace:
    def __init__(self, pairs):
        assert pairs >= 3 and pairs % 2
        self.pairs = pairs
        self.width = pairs * pairs
        self.degrees = [p + 4 * q for q in range(4) for p in range(4 - q)]
        self.local = evaluate_local((2, 16, 1, 1, 1, 1))
        def monomial(i, j):
            return tuple(1 if h == i else 3 if h == j else 0 for h in range(pairs))
        self.inputs = [[monomial(i, j) for j in range(pairs) if i != j] for i in range(pairs)]
        self.outputs = [monomial(i, (i + 1) % pairs) for i in range(pairs)]
        self.triples = [tuple(1 + ((orientation >> triple.index(j)) & 1) if j in triple else 0
                             for j in range(pairs))
                        for triple, orientation in product(combinations(range(pairs), 3), range(8))]
        self.baseline = [[self.kernel(self.outputs[x], self.outputs[u], y, v)
                          for v in range(pairs) for u in range(pairs)]
                         for x in range(pairs) for y in range(pairs)]

    @lru_cache(None)
    def pattern_value(self, pattern):
        value = 1
        for index in pattern:
            value = multiply(value, self.local[index])
            if not value:
                break
        return value

    @lru_cache(None)
    def kernel(self, first, second, y, v):
        syndrome = sum(((a.bit_count() + b.bit_count()) % 2) << i
                       for i, (a, b) in enumerate(zip(first, second)))
        if syndrome ^ (1 << y) ^ (1 << v):
            return 0
        value = 0
        for yy, vv in product(self.inputs[y], self.inputs[v]):
            pattern = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + d
                                   for a, b, c, d in zip(first, yy, second, vv)))
            value ^= self.pattern_value(pattern)
        return value

    def tensor(self, corrections):
        pairs = self.pairs
        tensor = [row.copy() for row in self.baseline]
        changed = [i for i, words in enumerate(corrections) if words]
        for i in changed:
            for extra in corrections[i]:
                assert sum(word.bit_count() for word in extra) == 3
                assert sum(word in (1, 2) for word in extra) == 3
                for y, v, u in product(range(pairs), repeat=3):
                    tensor[i * pairs + y][v * pairs + u] ^= self.kernel(extra, self.outputs[u], y, v)
                for x, y, v in product(range(pairs), repeat=3):
                    tensor[x * pairs + y][v * pairs + i] ^= self.kernel(self.outputs[x], extra, y, v)
        for i, j in product(changed, repeat=2):
            for first, second in product(corrections[i], corrections[j]):
                for y, v in product(range(pairs), repeat=2):
                    tensor[i * pairs + y][v * pairs + j] ^= self.kernel(first, second, y, v)
        return tensor

    def profile(self, tensor):
        pairs, width, degrees = self.pairs, self.width, self.degrees
        rows = [[0] * width for _ in degrees]
        columns = [[0] * width for _ in degrees]
        positions = {degree: i for i, degree in enumerate(degrees)}
        for x, y, v, u in product(range(pairs), repeat=4):
            value = tensor[x * pairs + y][v * pairs + u]
            assert (value & 1) == int(x == v and y == u)
            assert not (value >> 3) & 1 and not (value >> 12) & 1
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                i = positions[degree]
                rows[i][x * pairs + y] ^= 1 << (v * pairs + u)
                columns[i][x * pairs + u] ^= 1 << (v * pairs + y)
                value ^= bit
        chosen = independent_coefficients([[matrix] for matrix in rows])
        assert chosen[0] == 0
        row, column = [rows[i] for i in chosen], [columns[i] for i in chosen]
        ranks = [[0, 0] for _ in range(1 << len(chosen))]
        work = [[0] * width for _ in range(2)]
        previous = 0
        for step in range(1, 1 << len(chosen)):
            mask = step ^ (step >> 1)
            index = (mask ^ previous).bit_length() - 1
            for target, basis in zip(work, (row, column)):
                for i, value in enumerate(basis[index]):
                    target[i] ^= value
            ranks[mask] = [rank_two(matrix) for matrix in work]
            previous = mask
        cost, witness = minimum_decomposition([sqrt(r * c) for r, c in ranks], target=1)
        lower, _ = minimum_decomposition([isqrt(r * c) for r, c in ranks], target=1)
        return {"basis_degrees": [degrees[i] for i in chosen], "rank_table": ranks,
                "cost": cost, "integer_lower": lower, "witness": witness,
                "exclusion": lower ** 70 > pairs ** 99,
                "rank_table_sha256": hashlib.sha256(json.dumps(ranks, separators=(",", ":")).encode()).hexdigest()}


def candidates(space, seed, random_cases):
    d = space.pairs
    # A fixed pair of sites creates coherent triple corrections on d-2 rows.
    for a, b in combinations(range(d), 2):
        for orientation in range(8):
            corrections = [[] for _ in range(d)]
            for i in range(d):
                if i not in (a, b):
                    corrections[i] = [tuple((1 + (orientation & 1)) if j == i else
                                             (1 + ((orientation >> 1) & 1)) if j == a else
                                             (1 + ((orientation >> 2) & 1)) if j == b else 0
                                             for j in range(d))]
            yield {"family": "coherent_pair", "sites": [a, b], "orientation": orientation}, corrections
    rng = random.Random(seed)
    for case in range(random_cases):
        corrections = [[] for _ in range(d)]
        if case % 4 == 3:
            # One low-mode alternating quadratic form shared by every row.
            edges = [edge for edge in combinations(range(d), 2) if rng.randrange(4) == 0]
            for i in range(d):
                corrections[i] = [tuple(int(j in (i, a, b)) for j in range(d))
                                  for a, b in edges if i not in (a, b)]
            descriptor = {"family": "shared_quadratic", "case": case, "edges": edges}
        else:
            number = (2, 3, d)[case % 4]
            indices = sorted(rng.sample(range(d), number))
            for i in indices:
                corrections[i] = [space.triples[rng.randrange(len(space.triples))]]
            descriptor = {"family": "independent_triples", "case": case, "changed_rows": indices}
        yield descriptor, corrections


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=7)
    parser.add_argument("--random-cases", type=int, default=512)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    space = CubicSpace(args.pairs)
    baseline = space.profile(space.baseline)
    result = {"pairs": args.pairs, "labels": args.pairs, "length": 3, "seed": args.seed,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "baseline": baseline, "best": baseline, "records": [], "complete": False}
    for descriptor, corrections in candidates(space, args.seed, args.random_cases):
        tensor = space.tensor(corrections)
        record = dict(space.profile(tensor), **descriptor,
                      corrections=corrections,
                      changed_output_rows=sum(bool(words) for words in corrections),
                      changed_polynomial_entries=sum(a != b for row, old in zip(tensor, space.baseline)
                                                     for a, b in zip(row, old)))
        if record["cost"] < result["best"]["cost"]:
            result["best"] = record
        result["records"].append(record)
        if len(result["records"]) % 25 == 0:
            args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
            print(json.dumps({"cases": len(result["records"]), "best_cost": result["best"]["cost"],
                              "latest_cost": record["cost"], "latest_dimension": len(record["basis_degrees"])}), flush=True)
    result["complete"] = True
    result["all_recorded_spaces_excluded_over_F2"] = all(r["exclusion"] for r in result["records"])
    result["distinct_rank_tables"] = len({r["rank_table_sha256"] for r in result["records"]})
    result["endpoint_entries_checked"] = (len(result["records"]) + 1) * args.pairs ** 4
    result["distinct_local_patterns"] = space.pattern_value.cache_info().currsize
    result["scope"] = "The explicitly recorded multirow output maps over F2; no exhaustive claim about all multirow maps."
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("records", "baseline", "best")}), flush=True)


if __name__ == "__main__":
    main()
