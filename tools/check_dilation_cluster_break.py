#!/usr/bin/env python3
"""Finite checks for exponent dilation and cluster-entry density.

On the saved strict generator-seven cycle, the cluster-entry process at
the generator step wraps the period.  A small coprime exponent dilation
increases the number of entry clusters, lowers the entry density
denominator, and makes the same entry window fit.  The tool also records
the exact unit-dilation average for the sparse strict cycle and a
deterministic sample for the dense large-period cycles.

This is a finite diagnostic.  It proves no uniform cluster-breaking or
entry-gap theorem and does not improve the complete factoring exponent.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from check_cofactor_orbit_moments import (
    maximum_zero_run,
    moment_record,
    sufficient_moment_window,
)
from check_gap_moment_order import (
    GENEROUS,
    STRICT,
    entry_flags,
    orbit,
)


LARGE_N = 1000036000099
ORDER = 40
CONSTANT = 1
SAMPLE_SIZE = 32
SEED = 20261102


def dilate(flags: bytes, step: int) -> bytes:
    period = len(flags)
    out = bytearray(period)
    index = 0
    for position in range(period):
        out[position] = flags[index]
        index = (index + step) % period
    return bytes(out)


def runs(flags: bytes) -> int:
    return sum(entry_flags(flags))


def unit_list(period: int) -> list[int]:
    return [step for step in range(1, period) if math.gcd(step, period) == 1]


def exact_unit_average(flags: bytes) -> dict:
    period = len(flags)
    support = [index for index, bit in enumerate(flags) if bit]
    units = len(unit_list(period))
    unit_pairs = 0
    for left in support:
        for right in support:
            if math.gcd((left - right) % period, period) == 1:
                unit_pairs += 1
    average = len(support) - unit_pairs / units
    return {
        "method": "pair_enumeration",
        "period": period,
        "good_classes": len(support),
        "unit_count": units,
        "unit_difference_pairs": unit_pairs,
        "all_shift_average_runs": len(support) * (period - len(support)) / period,
        "unit_average_runs": average,
        "unit_average_over_good_classes": average / len(support),
    }


def exact_unit_average_bitset(flags: bytes) -> dict:
    period = len(flags)
    mask = (1 << period) - 1
    bits = 0
    for index, bit in enumerate(flags):
        if bit:
            bits |= 1 << index
    units = unit_list(period)
    unit_pairs = 0
    for step in units:
        rotated = ((bits << step) | (bits >> (period - step))) & mask
        unit_pairs += (bits & rotated).bit_count()
    support = sum(flags)
    average = support - unit_pairs / len(units)
    return {
        "method": "bitset_enumeration",
        "period": period,
        "good_classes": support,
        "unit_count": len(units),
        "unit_difference_pairs": unit_pairs,
        "all_shift_average_runs": support * (period - support) / period,
        "unit_average_runs": average,
        "unit_average_over_good_classes": average / support,
        "unit_pair_error_fraction":
            (unit_pairs - len(units) * support * support / period) / (support * support),
    }


def sampled_unit_average(flags: bytes, size: int) -> dict:
    period = len(flags)
    units = unit_list(period)
    sample = random.Random(SEED).sample(units, min(size, len(units)))
    values = [runs(dilate(flags, step)) for step in sample]
    support = sum(flags)
    return {
        "sample_size": len(sample),
        "sample_steps": sample,
        "minimum_runs": min(values),
        "average_runs": sum(values) / len(values),
        "maximum_runs": max(values),
        "all_shift_average_runs": support * (period - support) / period,
        "sample_average_over_good_classes": (sum(values) / len(values)) / support,
    }


def entry_stats(flags: bytes, step: int, integer: int, check_moment: bool) -> dict:
    period = len(flags)
    dilated = flags if step == 1 else dilate(flags, step)
    entries = entry_flags(dilated)
    entry_count = sum(entries)
    density_denominator = (period + entry_count - 1) // entry_count
    entry_gap = maximum_zero_run(entries)
    window = sufficient_moment_window(2 * integer, density_denominator, ORDER, CONSTANT)
    record = {
        "step": step,
        "raw_gap": maximum_zero_run(dilated),
        "runs": entry_count,
        "entry_density_denominator": density_denominator,
        "entry_gap": entry_gap,
        "moment_order": ORDER,
        "moment_constant": CONSTANT,
        "entry_window": window,
        "entry_window_fits_period": window < period,
    }
    if check_moment:
        check = moment_record(entries, ORDER, window, CONSTANT)
        record.update({
            "minimum_entry_hits": check["minimum_hits"],
            "lower_tail_moment_bound_holds":
                check["lower_tail_moment_bound_holds"],
            "lower_tail_least_positive_integer_constant":
                check["lower_tail_least_positive_integer_constant"],
        })
    return record


def strict_record() -> dict:
    integer, cutoff, cap = STRICT
    _, flags = orbit(integer, cutoff, cap, 7)
    period = len(flags)
    steps = (1, 3, 7, 13)
    entries = {
        str(step): entry_stats(flags, step, integer, check_moment=True)
        for step in steps
    }
    assert not entries["1"]["entry_window_fits_period"]
    assert entries["3"]["entry_window_fits_period"]
    assert entries["13"]["entry_window_fits_period"]
    assert entries["3"]["lower_tail_moment_bound_holds"]
    assert entries["13"]["lower_tail_moment_bound_holds"]
    assert entries["13"]["runs"] == sum(flags)
    assert entries["3"]["runs"] > entries["1"]["runs"]
    average = exact_unit_average(flags)
    assert average["unit_average_over_good_classes"] > 0.99
    return {
        "integer": integer,
        "cutoff": cutoff,
        "max_cofactor": cap,
        "generator_norm": 7,
        "period": period,
        "good_classes": sum(flags),
        "density_denominator": (period + sum(flags) - 1) // sum(flags),
        "unit_average": average,
        "steps": entries,
    }


def generous_record(prime: int) -> dict:
    integer, cutoff, cap = GENEROUS
    _, flags = orbit(integer, cutoff, cap, prime)
    period = len(flags)
    support = sum(flags)
    steps = {
        str(step): entry_stats(flags, step, integer, check_moment=False)
        for step in (1, 3, 7)
    }
    assert steps["1"]["entry_window_fits_period"]
    assert steps["3"]["entry_window_fits_period"]
    sample = sampled_unit_average(flags, SAMPLE_SIZE)
    assert sample["sample_average_over_good_classes"] > 0.7
    average = exact_unit_average_bitset(flags)
    assert abs(average["unit_pair_error_fraction"]) < 1e-3
    return {
        "integer": integer,
        "cutoff": cutoff,
        "max_cofactor": cap,
        "generator_norm": prime,
        "period": period,
        "good_classes": support,
        "density_denominator": (period + support - 1) // support,
        "unit_average": average,
        "sample_average": sample,
        "steps": steps,
    }


def smallest_odd_prime_not_dividing(period: int) -> int:
    candidate = 3
    while math.gcd(candidate, period) != 1:
        candidate += 2
    return candidate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    strict = strict_record()
    generous_five = generous_record(5)
    generous_seven = generous_record(7)
    for record in (strict, generous_five, generous_seven):
        assert record["period"] == 198850
    assert (generous_five["unit_average"]["unit_difference_pairs"]
            == generous_seven["unit_average"]["unit_difference_pairs"])
    assert smallest_odd_prime_not_dividing(198850) == 3
    output = {
        "claim": (
            "On the saved strict generator-seven cycle the entry window wraps "
            "at step one, while a small coprime exponent dilation increases "
            "the entry run count and makes the entry window fit.  Exact "
            "unit averages are near the all-shift baseline."
        ),
        "scope": (
            "Finite complete-orbit evidence only.  No uniform cluster-breaking "
            "lemma, entry-gap theorem, or complete bound improvement is proved."
        ),
        "rule": {
            "step": smallest_odd_prime_not_dividing(198850),
            "description": "smallest odd prime not dividing the period",
        },
        "strict": strict,
        "large_period_generous": {
            "generator_five": generous_five,
            "generator_seven": generous_seven,
        },
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n")
    else:
        print(text)


if __name__ == "__main__":
    main()
