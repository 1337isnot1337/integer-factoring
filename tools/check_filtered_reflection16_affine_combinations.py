#!/usr/bin/env python3
"""Replay a completed bounded cross-direction q16 affine search."""
from __future__ import annotations

import argparse
from datetime import datetime
import gzip
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import numpy as np

from dense_coupled_reflection_frames import digest, rank_one_terms, transformed_from_terms
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from reflection_laurent_physical_frames import add, bell_normalize, scale
from search_dense_coupled_reflection_frames import profile_rho
from search_filtered_reflection16_affine import scaled_terms
from search_filtered_reflection16_affine_combinations import component_catalogue, combined_seed
from search_physical_code_gauges import start_params


def field_power(value, exponent):
    result = 1
    for _ in range(exponent):
        result = F4_PRODUCT[result][value]
    return result


def scalar_inverse_row(component):
    """Invert the three-scalar evaluation map independently by enumeration."""
    wanted = [int(index == component) for index in range(3)]
    solutions = []
    for coefficients in product(range(4), repeat=3):
        values = []
        for exponent in (1, 2, 3):
            value = 0
            for coefficient, point in zip(coefficients, (1, 2, 3)):
                value ^= F4_PRODUCT[coefficient][field_power(point, exponent)]
            values.append(value)
        if values == wanted:
            solutions.append(coefficients)
    assert len(solutions) == 1
    return solutions[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(
        "evidence/filtered-reflection16-affine-combinations-search.json.gz"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    compressed = args.archive.read_bytes()
    raw = gzip.decompress(compressed) if args.archive.suffix == ".gz" else compressed
    archive = json.loads(raw)
    assert archive["completed"]
    assert archive["stopping_reason"] == "Declared time or trial limit reached."
    elapsed = (datetime.fromisoformat(archive["finished_utc"])
               - datetime.fromisoformat(archive["started_utc"])).total_seconds()
    assert abs(elapsed - archive["elapsed_seconds"]) < .01
    assert (archive["elapsed_seconds"] >= archive["requested_seconds"]
            or archive["counters"]["proposed"] == archive["requested_trials"])
    for name, expected in archive["source_sha256"].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name

    base = tensor(8)
    metadata, components = component_catalogue(base, archive["requested_poles"])
    assert metadata == archive["components"]
    assert len(components) == 52
    assert (len(archive["records"]) + archive["counters"]["duplicate_tensors"]
            == archive["counters"]["proposed"])
    hashes, trials = {digest(base["values"])}, set()
    for profile in archive["profiles"]:
        assert len(profile["orders"]) == 512
        assert all(isinstance(value, int) for value in profile["orders"])
    for row in archive["records"]:
        assert row["trial"] not in trials
        trials.add(row["trial"])
        assert row["word"] == sorted(row["word"])
        assert len({index for index, _ in row["word"]}) == len(row["word"])
        assert all(0 <= index < len(components) and coefficient in (1, 2, 3)
                   for index, coefficient in row["word"])
        changed = combined_seed(base, components, row["word"])
        actual = digest(changed["values"])
        assert actual == row["seed_sha256"] and actual not in hashes
        hashes.add(actual)
        assert changed["degree"] == row["degree"] == 96
        assert changed["monomial_content"] == row["monomial_content"] == 0
        if "profile_index" in row:
            assert row["natural_rho"] == archive["profiles"][row["profile_index"]]["rho"]
    unresolved = sum("profile_index" not in row for row in archive["records"])
    assert unresolved == archive["counters"]["unresolved_profiles"]
    best = min((row for row in archive["records"] if "natural_rho" in row),
               key=lambda row: row["natural_rho"])
    assert best["natural_rho"] == archive["best_natural_rho"]

    # Reconstruct the winner's components from the three complete scalar
    # frames. This does not use the affine_terms component formula.
    _, factors = bell_normalize(base)
    frame_cache = {}
    independent_entries = 0
    for index, _ in best["word"]:
        record = metadata[index]
        key = record["direction"], record["pole"]
        if key not in frame_cache:
            terms = rank_one_terms(base, record["u"], record["v"], factors)
            frame_cache[key] = [
                [add(value, original) for value, original in zip(
                    transformed_from_terms(base, scaled_terms(terms, scalar),
                                           record["pole"])["values"], base["values"])]
                for scalar in (1, 2, 3)
            ]
        coefficients = scalar_inverse_row(record["component"])
        for position in range(16**4):
            value = (0, 0)
            for coefficient, frame in zip(coefficients, frame_cache[key]):
                value = add(value, scale(frame[position], coefficient))
            value = scale(value, record["normalizing_scalar"])
            assert value == components[index].get(position, (0, 0))
            independent_entries += 1

    weights = physical_bit_energies(8)
    energy = np.asarray([sum(value for bit, value in enumerate(weights)
                             if state >> bit & 1) for state in range(16)], float)
    natural = start_params(energy, 1)
    fresh = []
    for label, changed, profile in (
        ("original", base, archive["profiles"][0]),
        ("best_retained", combined_seed(base, components, best["word"]),
         archive["profiles"][best["profile_index"]]),
    ):
        raw_orders = RectangularOracle(changed["values"], 16, 16).evaluate(
            natural, [1 / 16] * 16, [1 / 16] * 16, want_spectra=True)[2]
        orders = tuple(round(value) for value in raw_orders)
        assert max(abs(first - second) for first, second in zip(raw_orders, orders)) < 1e-8
        assert list(orders) == profile["orders"]
        rho = profile_rho(orders, energy, profile["scale"])
        assert abs(rho - profile["rho"]) < 1e-10
        fresh.append({"label": label, "seed_sha256": digest(changed["values"]),
                      "rho": rho, "scale": profile["scale"],
                      "complete_channel_count": len(orders)})
    report = {
        "status": "passed",
        "scope": "all retained complete tensor hashes, endpoints, degrees, source bindings; independent three-scalar winner reconstruction; fresh original and winner numerical profiles",
        "archive_sha256": sha256(compressed).hexdigest(),
        "uncompressed_archive_sha256": sha256(raw).hexdigest(),
        "recorded_search_source_hashes_checked": len(archive["source_sha256"]),
        "declared_seconds": archive["requested_seconds"],
        "actual_search_seconds": archive["elapsed_seconds"],
        "started_utc": archive["started_utc"],
        "finished_utc": archive["finished_utc"],
        "exact_components_reconstructed": len(components),
        "complete_retained_tensors_reconstructed": len(archive["records"]),
        "complete_polynomial_entries_reconstructed": len(archive["records"]) * 16**4,
        "independent_three_scalar_component_entries": independent_entries,
        "counters": archive["counters"],
        "fresh_profile_replays": fresh,
        "no_natural_pressure_improvement_over_original": best["natural_rho"] >= archive["baseline_rho"],
        "elapsed_seconds": time.monotonic() - started,
        "source_sha256": {Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()},
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "source_sha256"}))


if __name__ == "__main__":
    main()
