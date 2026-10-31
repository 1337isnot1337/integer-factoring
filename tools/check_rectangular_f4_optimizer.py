#!/usr/bin/env python3
"""Independent controls for paired Gibbs derivatives and minor separation.

This optional check needs NumPy, SciPy, and the native C++17 F4 oracle.
Floating checks test the implementation; they do not certify any search rate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random

import numpy as np

from check_rectangular_f4_repair import build
from native_rectangular_smith_f4 import RectangularOracle
from search_rectangular_f4_convex import RectangularConvexSearch
from search_rectangular_f4_pressure import paired_value_gradient


def check() -> dict:
    rng = random.Random(8041)
    seed = build((0, 1), 2, check_types=False)
    derivative_errors = []
    root_errors = []
    minor_slacks = []
    determinant_errors = []
    constraints_per_cut = []
    heldout_comparisons = 0

    for cut in (0, 1):
        search = RectangularConvexSearch(
            seed, cut=cut, height_bound=20.0, normalized_bound=40.0
        )
        oracle = search.oracle
        earlier_points = []
        for _ in range(4):
            parameters = np.array(
                [rng.uniform(-0.7, 0.7) for _ in range(oracle.parameters - 1)]
                + [rng.uniform(0.1, 0.3)]
            )
            rho, derivative, _, _ = paired_value_gradient(oracle, parameters, cut)
            reference_rho, *_ = oracle.types(parameters, cut)
            root_errors.append(abs(rho - reference_rho))

            for coordinate in range(oracle.parameters):
                delta = 1e-7
                plus = parameters.copy()
                minus = parameters.copy()
                plus[coordinate] += delta
                minus[coordinate] -= delta
                difference = (
                    paired_value_gradient(oracle, plus, cut)[0]
                    - paired_value_gradient(oracle, minus, cut)[0]
                ) / (2 * delta)
                derivative_errors.append(abs(float(derivative[coordinate]) - difference))

            profile, _, _ = search.separate(parameters)
            earlier_points.append(np.r_[parameters, profile["orders"]])
            inequalities = np.array(list(search.constraints.values()))
            equalities = np.array(list(search.equalities.values()))
            # Minors found at another gauge remain valid at every earlier gauge.
            for point in earlier_points:
                minor_slacks.append(float(min(inequalities @ point)))
                determinant_errors.append(float(max(abs(equalities @ point))))
                heldout_comparisons += len(inequalities)
        constraints_per_cut.append(len(search.constraints))

    assert max(root_errors) < 2e-10
    assert max(derivative_errors) < 2e-6
    assert min(minor_slacks) > -1e-8
    assert max(determinant_errors) < 1e-8

    # A separate unequal-alphabet seed has six disconnected support components.
    # Each component consists of one nonzero entry, so its complete spectra and
    # global prefix determinant orders follow directly from polynomial degrees.
    alphabet_x, alphabet_y = 3, 2
    values = [(0, 0)] * (alphabet_x**2 * alphabet_y**2)
    terms = []
    for x in range(alphabet_x):
        for y in range(alphabet_y):
            exponent = 2 * x + y + 1
            index = x + alphabet_x * (
                y + alphabet_y * (y + alphabet_y * x)
            )
            values[index] = (1, (1 << exponent) ^ (1 << (exponent + 1)))
            terms.append((x, y, exponent + 1))
    oracle = RectangularOracle(values, alphabet_x, alphabet_y)
    parameters = [rng.uniform(-1, 1) for _ in range(oracle.parameters - 1)] + [0.17]
    component_orders = 0
    component_prefixes = 0
    for cut in (0, 1):
        prefixes = oracle.prefixes(parameters, cut)
        for side in (0, 1):
            expected = []
            for x, y, degree in terms:
                labels = (x, y, y, x)
                height = sum(
                    parameters[side * oracle.stride + oracle.offsets[leg] + label - 1]
                    for leg, label in enumerate(labels)
                    if label
                )
                base = -degree if side == 0 else 0
                expected.append((height + base * parameters[-1], base))
            expected.sort()
            cumulative = 0
            for index, (order, base) in enumerate(expected):
                cumulative += base
                position = side * oracle.n + index
                assert abs(order - prefixes["orders"][position]) < 1e-10
                assert cumulative == prefixes["base_prefix"][position]
                component_orders += 1
                component_prefixes += 1

    source_names = (
        "native_rectangular_smith_f4.cpp",
        "native_rectangular_smith_f4.py",
        "search_rectangular_f4_pressure.py",
        "search_rectangular_f4_convex.py",
        Path(__file__).name,
    )
    source_directory = Path(__file__).parent
    return {
        "status": "pass",
        "random_seed": 8041,
        "actual_mixed_seed_points": len(root_errors),
        "paired_gradient_coordinates": len(derivative_errors),
        "paired_gradient_maximum_error": max(derivative_errors),
        "gibbs_root_maximum_error": max(root_errors),
        "heldout_minor_comparisons": heldout_comparisons,
        "heldout_minor_minimum_slack": min(minor_slacks),
        "determinant_equality_maximum_error": max(determinant_errors),
        "constraints_per_cut": constraints_per_cut,
        "disconnected_component_orders": component_orders,
        "disconnected_component_prefixes": component_prefixes,
        "source_sha256": {
            name: hashlib.sha256((source_directory / name).read_bytes()).hexdigest()
            for name in source_names
        },
        "scope": (
            "Finite implementation checks for both cuts, independent Gibbs roots, "
            "the envelope derivative, attaining-minor inequalities across gauges, "
            "and disconnected component assembly. Floating comparisons do not "
            "certify an exponent, a global optimum, or arbitrary future profiles."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps({"status": "fail", "error": repr(error)}, indent=2) + "\n"
            )
        raise
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
