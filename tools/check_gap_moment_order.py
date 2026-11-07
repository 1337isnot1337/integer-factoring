#!/usr/bin/env python3
"""Exact checks for the gap/order tradeoff in the cofactor moment criterion.

The tool records two finite facts on already used arithmetic orbits:

1. The maximum inadmissible gap on every saved complete orbit is below
   ``2 * ceil(1/delta) * bit_length(period)``.  A clustered abstract orbit
   violates the same bound, so density alone does not imply it.
2. At the period-scale moment order ``q = input_bits + 1`` the sufficient
   window is much shorter than at the previously tested small orders, and
   the exact lower-tail bound holds with a small constant on the saved
   generous-predicate cycles.  On the strict cycle the raw process fits
   inside the period while the cluster-entry process still wraps.

This is a finite diagnostic.  It proves no uniform gap, moment, or window
theorem and does not improve the complete factoring exponent.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from check_cofactor_orbit_moments import (
    actual_orbit,
    maximum_zero_run,
    moment_record,
    sufficient_moment_window,
)
from prime_form_factor_base import build_prime_form_factor_base


GENEROUS = (1000036000099, 19, 10000)
STRICT = (1000036000099, 19, 1)
TEST_PRIMES = (5, 7)
STRICT_TEST_PRIMES = (7,)


def prime_form(integer: int, cutoff: int, prime: int) -> tuple[int, int, int]:
    delta = -integer if integer % 4 == 3 else -4 * integer
    base = build_prime_form_factor_base(delta, cutoff)
    for item in base["records"]:
        if item.prime == prime:
            return item.positive
    raise ValueError(f"Prime {prime} is not in the split base at cutoff {cutoff}")


_ORBIT_CACHE: dict[tuple[int, int, int, int], tuple[dict, bytes]] = {}


def orbit(integer: int, cutoff: int, cap: int, prime: int) -> tuple[dict, bytes]:
    key = (integer, cutoff, cap, prime)
    if key not in _ORBIT_CACHE:
        generator = {"kind": "prime", "norm": prime,
                     "form": prime_form(integer, cutoff, prime)}
        record, flags = actual_orbit(integer, cutoff, cap, generator, 1_000_000)
        if flags is None:
            raise RuntimeError(f"Orbit for N={integer}, prime={prime} did not close")
        _ORBIT_CACHE[key] = (record, flags)
    return _ORBIT_CACHE[key]


def entry_flags(flags: bytes) -> bytes:
    period = len(flags)
    return bytes(flags[index] and not flags[index - 1] for index in range(period))


def panel_rows() -> list[dict]:
    rows: dict[tuple, dict] = {}
    specifications = (
        ("evidence/cofactor-orbit-moments.json", "cases"),
        ("evidence/entry-gap-cluster-loss.json", "complete_cases"),
    )
    for name, key in specifications:
        data = json.loads(Path(name).read_text())
        for case in data[key]:
            if not case.get("complete_period"):
                continue
            generator = case["generator"]
            period = case["period"]
            good = case["good_classes"]
            gap = case["maximum_inadmissible_gap"]
            identity = (case["integer"], generator["kind"], generator["norm"],
                        period, good, gap)
            denominator = (period + good - 1) // good
            rows[identity] = {
                "integer": case["integer"],
                "kind": generator["kind"],
                "norm": generator["norm"],
                "period": period,
                "good_classes": good,
                "density_denominator": denominator,
                "maximum_inadmissible_gap": gap,
                "period_bit_length": period.bit_length(),
                "integer_gap_bound": 2 * denominator * period.bit_length(),
                "gap_ratio_over_d_log2_h": gap / (denominator * math.log2(period)),
                "saved_integer_bound_holds": gap < 2 * denominator * period.bit_length(),
            }
    return sorted(rows.values(), key=lambda row: (
        row["integer"], row["kind"], row["norm"], row["period"]))


def clustered_control() -> dict:
    period, good = 4096, 2048
    flags = bytes([1] * good + [0] * (period - good))
    denominator = (period + good - 1) // good
    gap = maximum_zero_run(flags)
    return {
        "period": period,
        "good_classes": good,
        "density": good / period,
        "density_denominator": denominator,
        "maximum_inadmissible_gap": gap,
        "integer_gap_bound": 2 * denominator * period.bit_length(),
        "gap_ratio_over_d_log2_h": gap / (denominator * math.log2(period)),
        "saved_integer_bound_holds": gap < 2 * denominator * period.bit_length(),
        "claim": "An interval-clustered admissible set satisfies the density and violates the recorded gap bound.",
    }


def scan(flags: bytes, period_bound: int, orders: tuple[int, ...],
         constants: tuple[int, ...]) -> list[dict]:
    period, good, gap = len(flags), sum(flags), maximum_zero_run(flags)
    denominator = (period + good - 1) // good
    records = []
    for order in orders:
        for constant in constants:
            window = sufficient_moment_window(period_bound, denominator, order, constant)
            check = moment_record(flags, order, window, constant)
            records.append({
                "order": order,
                "constant": constant,
                "window": window,
                "period": period,
                "wrapping": window >= period,
                "minimum_hits": check["minimum_hits"],
                "maximum_hits": check["maximum_hits"],
                "lower_tail_least_positive_integer_constant":
                    check["lower_tail_least_positive_integer_constant"],
                "lower_tail_moment_bound_holds":
                    check["lower_tail_moment_bound_holds"],
                "maximum_gap": gap,
                "gap_below_window": gap < window,
            })
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = panel_rows()
    for row in rows:
        assert row["saved_integer_bound_holds"], row
        assert row["gap_ratio_over_d_log2_h"] < 2.0, row
    control = clustered_control()
    assert not control["saved_integer_bound_holds"]
    assert control["gap_ratio_over_d_log2_h"] > 2.0

    integer, cutoff, generous_cap = GENEROUS
    generous = {}
    for prime in TEST_PRIMES:
        record, flags = orbit(integer, cutoff, generous_cap, prime)
        generous[str(prime)] = {
            "generator": record["generator"],
            "period": record["period"],
            "good_classes": record["good_classes"],
            "density_denominator":
                (record["period"] + record["good_classes"] - 1) // record["good_classes"],
            "maximum_inadmissible_gap": record["maximum_inadmissible_gap"],
            "scans": scan(flags, 2 * integer, (6, 20, 40, 60), (1, 2)),
        }
    # The small order fails at C=1 for generator five; the period-scale
    # order passes at C=1 with a much shorter window.
    first = generous["5"]["scans"]
    small = [item for item in first if item["order"] == 6 and item["constant"] == 1][0]
    large = [item for item in first if item["order"] == 40 and item["constant"] == 1][0]
    assert not small["lower_tail_moment_bound_holds"] and not small["wrapping"]
    assert large["lower_tail_moment_bound_holds"] and large["gap_below_window"]
    assert large["window"] < small["window"]
    for prime in TEST_PRIMES:
        final = [item for item in generous[str(prime)]["scans"]
                 if item["order"] == 40 and item["constant"] == 1][0]
        assert final["lower_tail_moment_bound_holds"]
        assert final["gap_below_window"]

    strict_integer, strict_cutoff, strict_cap = STRICT
    strict = {}
    for prime in STRICT_TEST_PRIMES:
        record, flags = orbit(strict_integer, strict_cutoff, strict_cap, prime)
        entries = entry_flags(flags)
        entry_count = sum(entries)
        entry_denominator = (len(entries) + entry_count - 1) // entry_count
        raw_scan = scan(flags, 2 * strict_integer, (20, 40), (1,))[1]
        entry_window = sufficient_moment_window(
            2 * strict_integer, entry_denominator, 40, 1)
        strict[str(prime)] = {
            "generator": record["generator"],
            "period": record["period"],
            "good_classes": record["good_classes"],
            "density_denominator":
                (record["period"] + record["good_classes"] - 1) // record["good_classes"],
            "maximum_inadmissible_gap": record["maximum_inadmissible_gap"],
            "raw_order_40": raw_scan,
            "entry_count": entry_count,
            "entry_density_denominator": entry_denominator,
            "entry_maximum_gap": maximum_zero_run(entries),
            "entry_order_40_window": entry_window,
            "entry_wraps_period": entry_window >= len(entries),
        }
        assert raw_scan["lower_tail_moment_bound_holds"]
        assert not raw_scan["wrapping"]
        assert strict[str(prime)]["entry_wraps_period"]

    # The prescribed window, not the gap itself, is what makes an empty
    # window contradict the moment bound.  The clustered control has a gap
    # above its own prescribed window and fails; the saved strict cycle has
    # a gap below the prescribed window and passes, even at W = gap, where
    # the mean count is still below the moment order.
    clustered = clustered_control()
    clustered_flags = bytes([1] * clustered["good_classes"]
                            + [0] * (clustered["period"] - clustered["good_classes"]))
    clustered_window = sufficient_moment_window(
        clustered["period"], clustered["density_denominator"], 40, 1)
    clustered_check = moment_record(clustered_flags, 40, clustered_window, 1)
    assert clustered["maximum_inadmissible_gap"] >= clustered_window
    assert clustered_check["minimum_hits"] == 0
    assert not clustered_check["lower_tail_moment_bound_holds"]

    _, strict_flags = orbit(*STRICT, 7)
    strict_gap = maximum_zero_run(strict_flags)
    strict_gap_check = moment_record(strict_flags, 40, strict_gap, 1)
    assert strict_gap_check["minimum_hits"] == 0
    assert strict_gap_check["lower_tail_moment_bound_holds"]

    empty_controls = {
        "abstract_clustered": {
            "order": 40,
            "constant": 1,
            "window": clustered_window,
            "maximum_gap": clustered["maximum_inadmissible_gap"],
            "minimum_hits": clustered_check["minimum_hits"],
            "lower_tail_moment_bound_holds":
                clustered_check["lower_tail_moment_bound_holds"],
            "note": "Gap exceeds the prescribed window, so an empty window contradicts the bound.",
        },
        "saved_strict_seven_at_gap": {
            "order": 40,
            "constant": 1,
            "window": strict_gap,
            "maximum_gap": strict_gap,
            "minimum_hits": strict_gap_check["minimum_hits"],
            "lower_tail_moment_bound_holds":
                strict_gap_check["lower_tail_moment_bound_holds"],
            "note": ("The empty window alone does not contradict the bound here because "
                     "delta*W is still below the moment order; the prescribed window is larger."),
        },
    }

    result = {
        "claim": "Finite gap/order tradeoff for the cofactor lower-tail moment criterion",
        "scope": ("Saved complete orbits and one abstract clustered control. "
                  "No uniform gap, density, moment, or window theorem is proved; "
                  "the complete factoring exponent remains 1/2."),
        "panel": {
            "complete_orbits": len(rows),
            "maximum_gap_ratio_over_d_log2_h":
                max(row["gap_ratio_over_d_log2_h"] for row in rows),
            "minimum_gap_ratio_over_d_log2_h":
                min(row["gap_ratio_over_d_log2_h"] for row in rows),
            "all_integer_bounds_hold": True,
            "rows": rows,
        },
        "abstract_clustered_control": control,
        "generous_cycle_scans": generous,
        "strict_cycle_scans": strict,
        "empty_window_controls": empty_controls,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
