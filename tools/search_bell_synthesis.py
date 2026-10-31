#!/usr/bin/env python3
"""Synthesize Bell polynomial paths by exact two-cut filtration constraints.

GF(3) rows use disjoint bit sets for coefficients one and two. Feasibility
is exact. The greedy order of tightening constraints is a finite heuristic.
"""

import argparse
import json
from math import log
import random
import time

from bell_crossing import sign, transform
from search_bell_profiles import optimum, profiles, record, tensor_start


class AffineTernary:
    def __init__(self, variables: int):
        self.variables = variables
        self.variable_mask = (1 << variables) - 1
        self.mask = (1 << (variables + 1)) - 1
        self.basis = {}

    def add(self, one: int, two: int) -> bool:
        """Commit an equation if consistent; leave the basis intact otherwise."""
        while (one | two) & self.variable_mask:
            pivot_bit = ((one | two) & self.variable_mask) & -((one | two) & self.variable_mask)
            pivot = pivot_bit.bit_length() - 1
            if pivot not in self.basis:
                if two & pivot_bit:
                    one, two = two, one
                self.basis[pivot] = one, two
                return True
            b1, b2 = self.basis[pivot]
            if one & pivot_bit:
                b1, b2 = b2, b1  # subtract the stored row
            a0, b0 = self.mask ^ one ^ two, self.mask ^ b1 ^ b2
            one, two = ((a0 & b1) | (one & b0) | (two & b2),
                        (a0 & b2) | (two & b0) | (one & b1))
        return not (one | two)

    def solution(self) -> tuple[int, ...]:
        values = [0] * self.variables
        value_one = value_two = 0
        for pivot in sorted(self.basis, reverse=True):
            one, two = self.basis[pivot]
            rhs = ((one >> self.variables) & 1) + 2 * ((two >> self.variables) & 1)
            dot = ((one & value_one).bit_count() + (two & value_two).bit_count()
                   + 2 * ((one & value_two).bit_count() + (two & value_one).bit_count()))
            value = (rhs - dot) % 3
            values[pivot] = value
            if value == 1:
                value_one |= 1 << pivot
            elif value == 2:
                value_two |= 1 << pivot
        return tuple(values)


def synthesis(bits: int, degree: int, scale: float, seed: int,
              jitter: float = 0.25, tensor_initial: bool = False,
              relax_probability: float = 0.15,
              initial_coefficients=None) -> dict:
    d, q = 1 << bits, 1 << (2 * bits)
    rng = random.Random(seed)
    rows = []
    for c in range(d):
        for e in range(d):
            one = two = 0
            for a in range(d):
                for b in range(d):
                    value = pow(d, -1, 3) * sign((a ^ c) & (b ^ e)) % 3
                    if value == 1:
                        one |= 1 << (a * d + b)
                    else:
                        two |= 1 << (a * d + b)
            rows.append((one, two))
    system = AffineTernary((degree + 1) * q)
    rhs = 1 << system.variables
    for j in range(q):
        assert system.add(sum(1 << (i * q + j) for i in range(degree + 1)) | rhs, 0)
    lows, highs = [0] * q, [degree] * q
    if tensor_initial:
        initial = initial_coefficients or tensor_start(bits, degree)
        initial_dual = tuple(transform(row, bits, 3) for row in initial)
        for j in range(q):
            lows[j] = next(i for i in range(degree + 1) if initial[i][j])
            highs[j] = next(i for i in range(degree, -1, -1) if initial_dual[i][j])
            if rng.random() < relax_probability:
                lows[j] = rng.randrange(lows[j] + 1)
            if rng.random() < relax_probability:
                highs[j] = rng.randrange(highs[j], degree + 1)
            for i in range(lows[j]):
                assert system.add(1 << (i * q + j), 0)
            for i in range(highs[j] + 1, degree + 1):
                one, two = rows[j]
                assert system.add(one << (i * q), two << (i * q))
    blocked_low, blocked_high = set(), set()
    proposals = accepted = 0
    while True:
        a = sum(scale ** (-x) for x in lows)
        b = sum(scale ** x for x in highs)
        choices = []
        for j in range(q):
            if lows[j] < degree and j not in blocked_low:
                gain = scale ** (-lows[j]) / a
                choices.append((gain * (1 + rng.uniform(-jitter, jitter)), False, j))
            if highs[j] and j not in blocked_high:
                gain = scale ** highs[j] / b
                choices.append((gain * (1 + rng.uniform(-jitter, jitter)), True, j))
        if not choices:
            break
        _, column, j = max(choices)
        proposals += 1
        if column:
            one, two = rows[j]
            shift = highs[j] * q
            good = system.add(one << shift, two << shift)
            if good:
                highs[j] -= 1
            else:
                blocked_high.add(j)
        else:
            good = system.add(1 << (lows[j] * q + j), 0)
            if good:
                lows[j] += 1
            else:
                blocked_low.add(j)
        accepted += good
    solution = system.solution()
    coefficients = tuple(tuple(solution[i * q:(i + 1) * q]) for i in range(degree + 1))
    other = tuple(transform(row, bits, 3) for row in coefficients)
    assert all(sum(col) % 3 == 1 for col in zip(*coefficients))
    assert all(coefficients[i][j] == 0 for j in range(q) for i in range(lows[j]))
    assert all(other[i][j] == 0 for j in range(q) for i in range(highs[j] + 1, degree + 1))
    rl, rh = profiles(coefficients)
    cl, ch = profiles(other)
    value = min(optimum(rl, ch)[0], optimum(cl, rh)[0])
    return {"coefficients": coefficients, "value": value,
            "proposals": proposals, "accepted": accepted,
            "affine_rank": len(system.basis)}


def run(bits: int, degree: int, trials: int, seed: int,
        tensor_initial: bool = False, adaptive: bool = False) -> dict:
    start = time.monotonic()
    best = None
    improvements = []
    profiles_seen = set()
    total_proposals = 0
    scales = (1.25, 1.5, 2., 3., 4., 8., 16., 32.)
    for trial in range(trials):
        scale = scales[trial % len(scales)]
        initial = None
        if adaptive and best is not None and trial % 16:
            initial = best["coefficients"]
            dual = tuple(transform(row, bits, 3) for row in initial)
            rl, rh = profiles(initial)
            cl, ch = profiles(dual)
            if optimum(cl, rh)[0] < optimum(rl, ch)[0]:
                initial = dual
        result = synthesis(bits, degree, scale, seed + trial,
                           jitter=(0., .15, .4, .8)[trial // len(scales) % 4],
                           tensor_initial=tensor_initial,
                           relax_probability=(.02, .05, .1, .2, .4, .6)[trial // 32 % 6],
                           initial_coefficients=initial)
        total_proposals += result["proposals"]
        coefficients = result["coefficients"]
        other = tuple(transform(row, bits, 3) for row in coefficients)
        profiles_seen.add((profiles(coefficients), profiles(other)))
        if best is None or result["value"] < best["value"] - 1e-12:
            best = result
            improvements.append({"trial": trial, "seed": seed + trial,
                                 "scale": scale,
                                 "rho_approx": result["value"] / (4 * bits * log(2))})
    return {"status": "exact synthesis schedule completed", "bits": bits,
            "local_dimension": 1 << bits, "degree": degree, "trials": trials,
            "tensor_initial": tensor_initial,
            "adaptive": adaptive,
            "seed": seed, "distinct_profile_pairs": len(profiles_seen),
            "total_exact_constraint_proposals": total_proposals,
            "improvements": improvements, "best": record(best["coefficients"], bits),
            "elapsed_seconds": time.monotonic() - start,
            "scope": "Exact affine feasibility in a specified heuristic target order"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=3)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--tensor-initial", action="store_true")
    parser.add_argument("--adaptive", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = run(args.bits, args.degree, args.trials, args.seed, args.tensor_initial, args.adaptive)
    with open(args.output, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({key: result[key] for key in ("status", "local_dimension", "degree",
          "trials", "distinct_profile_pairs", "total_exact_constraint_proposals", "elapsed_seconds")}))
    print(json.dumps({"bounds": result["best"]["bounds"], "improvements": result["improvements"]}))


if __name__ == "__main__":
    main()
