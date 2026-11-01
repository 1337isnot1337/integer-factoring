#!/usr/bin/env python3
"""Exact geometry and cost controls for the fixed-axis cube transfer.

The channel assignments below are abstract widths, not constructed tensors.
Only the geometric transfer and schedule arithmetic are being tested.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path
import random


def geometry(d: int, g: int) -> dict:
    vertices = list(product(range(g), repeat=d))
    cubes: dict[tuple, list] = {}
    other = [Counter() for _ in range(d)]
    generic = 0
    for first, second in combinations(vertices, 2):
        if all(x != y for x, y in zip(first, second)):
            cube = tuple(tuple(sorted((x, y))) for x, y in zip(first, second))
            cubes.setdefault(cube, []).append((first, second))
            generic += 1
        else:
            for axis in range(d):
                if first[axis] != second[axis]:
                    other[axis][tuple(sorted((first[axis], second[axis])))] += 1
    pairs = list(combinations(range(g), 2))
    expected_cubes = len(pairs) ** d
    assert len(cubes) == expected_cubes
    assert generic == (1 << (d - 1)) * expected_cubes
    assert all(len(edges) == 1 << (d - 1) for edges in cubes.values())
    for cube, edges in cubes.items():
        seen = set()
        for first, second in edges:
            assert all({first[a], second[a]} == set(cube[a]) for a in range(d))
            seen.add(first)
            seen.add(second)
        assert len(seen) == 1 << d
    plane_cubes = [Counter(cube[a] for cube in cubes) for a in range(d)]
    expected_meet = len(pairs) ** (d - 1)
    expected_other = g ** (2 * d - 2) - (g * (g - 1)) ** (d - 1)
    for axis in range(d):
        for pair in pairs:
            assert plane_cubes[axis][pair] == expected_meet
            assert other[axis][pair] == expected_other

    lam = 7
    cap = (1 << (d - 1)) * lam
    volume = g ** (2 * d) * lam
    assert cap * expected_meet + expected_other * lam == volume // (g * g)
    rng = random.Random(319018 + 101 * d + g)
    assignments = 0
    for mode in range(5):
        totals = [0] * d
        widths = [Counter() for _ in range(d)]
        max_sum = 0
        for number, cube in enumerate(sorted(cubes)):
            if mode == 0:
                local = [0] * d
            elif mode == 1:
                local = [cap] * d
            elif mode == 2:
                local = [cap if a == number % d else 0 for a in range(d)]
            elif mode == 3:
                local = [(number * (a + 1) + 3 * a) % (cap + 1) for a in range(d)]
            else:
                local = [rng.randrange(cap + 1) for _ in range(d)]
            max_sum = max(max_sum, sum(local))
            for axis, width in enumerate(local):
                totals[axis] += width
                widths[axis][cube[axis]] += width
        rho = Fraction(max_sum, d * cap)
        axis = min(range(d), key=totals.__getitem__)
        actual = 2 * sum(widths[axis][pair] + other[axis][pair] * lam for pair in pairs)
        formula = (Fraction(g - 1, g) - (1 - rho) * Fraction(g - 1, g) ** d) * volume
        assert actual <= formula
        assert formula <= (rho + Fraction(d - 1, g)) * volume
        assert actual <= Fraction(g - 1, g) * volume
        assert all(widths[axis][pair] + other[axis][pair] * lam <= volume // (g * g)
                   for pair in pairs)
        assignments += 1
    pure_assignment_exponent = cap * expected_cubes
    assert pure_assignment_exponent == Fraction(1, 2) * Fraction(g - 1, g) ** d * volume
    return {"axes": d, "side": g, "vertices": len(vertices), "cubes": expected_cubes,
            "generic_edges": generic,
            "other_edges": len(vertices) * (len(vertices) - 1) // 2 - generic,
            "plane_pair_other_edges": expected_other, "width_assignments": assignments,
            "pure_piece_exponent": pure_assignment_exponent}


def schedules() -> list[dict]:
    out = []
    for d in range(2, 17):
        a = (4 * d).bit_length()
        assert 1 << a > 4 * d
        x0 = max(a, 32)
        logs = [x0]
        for j in range(9):
            old = logs[-1]
            increment = a * (j + 1) + x0
            new = d * (old + increment)
            assert increment <= 2 * old
            assert new <= 3 * d * old
            exact = d ** (j + 1) * (
                Fraction(x0) + sum(Fraction(d * (a * (i + 1) + x0), d ** (i + 1))
                                   for i in range(j + 1)))
            assert new == exact
            # ln(K) < log2(K), so this stronger integer test implies (7).
            assert (1 << increment) >= 16 * (1 << (2 * j + 16)) * (2 * old + 1)
            logs.append(new)
        out.append({"axes": d, "a": a, "x0": x0, "log2_K": logs})
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = ([(d, 2) for d in range(2, 10)] + [(d, 3) for d in range(2, 6)]
             + [(d, 4) for d in range(2, 4)])
    data = {"status": "exact finite geometry and cost controls; local tensor family remains conditional",
            "complete_factoring_score": "1/2 unchanged",
            "geometry": [geometry(d, g) for d, g in cases], "schedules": schedules(),
            "transfer_exponent": "log_d(1/rho)"}
    rendered = json.dumps(data, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(json.dumps({"geometry_cases": len(cases), "schedule_dimensions": len(data["schedules"]),
                      "width_assignments": sum(row["width_assignments"] for row in data["geometry"]),
                      "status": "passed; no stronger factoring exponent claimed"}, sort_keys=True))


if __name__ == "__main__":
    main()
