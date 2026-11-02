#!/usr/bin/env python3
"""Exact entry-gap/cluster checks for the Path C cofactor predicate.

The tool does three things:

* proves the elementary entry-gap equivalence on all short binary words;
* measures entry gaps, cluster lengths and conditional windows on complete
  small cofactor orbits, excluding the two cycles reserved by the owner;
* replays the raw hit-prefix scans for key integrity and stopping-phase
  empty-window fractions.

No uniform density, cluster-length or moment theorem is claimed.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from check_cofactor_orbit_moments import (
    actual_orbit,
    ceiling_root,
    cyclic_window_counts,
    maximum_zero_run,
    moment_record,
)
from cofactor_ideal_closure import factor_with_cofactor_ideal_closure
from prime_form_factor_base import build_prime_form_factor_base


def ceil_div(numerator: int, denominator: int) -> int:
    if denominator < 1:
        raise ValueError("Require a positive denominator")
    return -(-numerator // denominator)


def cyclic_runs(flags: bytes) -> list[tuple[int, int]]:
    """Maximal cyclic runs of ones, indexed by their entry positions."""
    if not flags or not any(flags):
        return []
    period = len(flags)
    starts = [index for index in range(period)
              if flags[index] and not flags[(index - 1) % period]]
    runs = []
    for start in starts:
        length = 0
        while flags[(start + length) % period] and length < period:
            length += 1
        runs.append((start, length))
    assert sum(length for _, length in runs) == sum(flags)
    return runs


def entry_flags(flags: bytes) -> bytes:
    period = len(flags)
    return bytes(1 if flags[k] and not flags[(k - 1) % period] else 0
                 for k in range(period))


def entry_gap_record(flags: bytes) -> dict:
    """Exact entry gaps d_r = L_r + G_r, with L_r the hit run and G_r the
    following zero gap.  Every d_r is at least two and they sum to the
    period.  An empty entry window of length W exists iff W < max d_r."""
    period, good = len(flags), sum(flags)
    if not 0 < good < period:
        raise ValueError("Require a nonconstant binary word")
    runs = cyclic_runs(flags)
    if not runs:
        raise ValueError("Require at least one run")
    starts = [start for start, _ in runs]
    lengths = [length for _, length in runs]
    gaps = []
    for index, (start, length) in enumerate(runs):
        next_start = starts[(index + 1) % len(starts)]
        zero_gap = (next_start - (start + length)) % period
        if zero_gap == 0 and len(runs) == 1:
            zero_gap = period - length
        gaps.append(zero_gap)
    entry_gaps = [length + gap for length, gap in zip(lengths, gaps)]
    entry = entry_flags(flags)
    assert sum(entry) == len(runs)
    assert all(not (entry[k] and entry[(k + 1) % period]) for k in range(period))
    assert sum(entry_gaps) == period
    assert min(entry_gaps) >= 2
    raw_max_zero = maximum_zero_run(flags)
    assert max(gaps) == raw_max_zero
    max_gap = max(entry_gaps)
    empty_below = 0
    if max_gap > 1:
        empty_below = sum(count == 0 for count in cyclic_window_counts(entry, max_gap - 1))
        assert empty_below > 0
    counts_at_max = cyclic_window_counts(entry, max_gap)
    assert min(counts_at_max) > 0
    raw_denominator = ceil_div(period, good)
    entry_denominator = ceil_div(period, len(runs))
    mean_cluster = Fraction(good, len(runs))
    cluster_ceiling = ceil_div(good, len(runs))
    assert entry_denominator <= raw_denominator * cluster_ceiling
    return {
        "period": period,
        "good_classes": good,
        "entry_count": len(runs),
        "density": good / period,
        "entry_density": len(runs) / period,
        "raw_density_denominator": raw_denominator,
        "entry_density_denominator": entry_denominator,
        "mean_cluster_numerator": mean_cluster.numerator,
        "mean_cluster_denominator": mean_cluster.denominator,
        "mean_cluster": float(mean_cluster),
        "cluster_ceiling": cluster_ceiling,
        "max_cluster": max(lengths),
        "raw_max_zero_run": raw_max_zero,
        "max_entry_gap": max_gap,
        "min_entry_gap": min(entry_gaps),
        "entry_gap_histogram": {str(gap): entry_gaps.count(gap)
                                for gap in sorted(set(entry_gaps))},
        "cluster_length_histogram": {str(length): lengths.count(length)
                                     for length in sorted(set(lengths))},
        "lemma_checks": {
            "entries_isolated": True,
            "sum_entry_gaps_equals_period": True,
            "entry_gap_at_least_two": True,
            "raw_max_zero_run_matches": True,
            "empty_window_below_max_gap": empty_below,
            "minimum_count_at_max_gap": min(counts_at_max),
            "cluster_loss_inequality": True,
        },
    }


def conditional_window(period_bound: int, density_denominator: int,
                       order: int, constant: int) -> int:
    threshold = (2 * order + 1) * (constant * order) ** order * period_bound
    return density_denominator * (ceiling_root(threshold, order + 1) + 1)


def conditional_record(flags: bytes, order: int, constant: int,
                       period_bound: int) -> dict:
    period, good = len(flags), sum(flags)
    denominator = ceil_div(period, good)
    window = conditional_window(period_bound, denominator, order, constant)
    check = moment_record(flags, order, window, constant)
    return {
        "order": order,
        "constant": constant,
        "period_bound": period_bound,
        "window": window,
        "wrapping": window >= period,
        "minimum_hits": check["minimum_hits"],
        "lower_tail_moment_bound_holds": check["lower_tail_moment_bound_holds"],
        "lower_tail_least_positive_integer_constant":
            check["lower_tail_least_positive_integer_constant"],
    }


def empty_fraction_record(flags: bytes, mean_count: int) -> dict:
    period, good = len(flags), sum(flags)
    if good < 1:
        raise ValueError("Require a nonempty good set")
    window = max(1, ceil_div(mean_count * period, good))
    if window >= period:
        return {
            "mean_count": mean_count,
            "window": window,
            "wrapping": True,
            "empty_phases": 0,
            "phases": period,
            "empty_fraction": 0.0,
        }
    counts = cyclic_window_counts(flags, window)
    empty = sum(count == 0 for count in counts)
    return {
        "mean_count": mean_count,
        "window": window,
        "wrapping": False,
        "empty_phases": empty,
        "phases": period,
        "empty_fraction": empty / period,
    }


def orbit_record(integer: int, cutoff: int, cap: int, generator: dict,
                 period_cap: int) -> dict | None:
    record, flags = actual_orbit(integer, cutoff, cap, generator, period_cap)
    if flags is None:
        return {"integer": integer, "input_bits": record["input_bits"],
                "cutoff": cutoff, "max_cofactor": cap, "generator": generator,
                "complete_period": False}
    entry = entry_flags(flags)
    gap = entry_gap_record(flags)
    bits = record["input_bits"]
    orders = [("bit_length", (bits + 1).bit_length()),
              ("square", (bits + 1).bit_length() ** 2)]
    windows = []
    for name, order in orders:
        for constant in (1, 2):
            for period_name, period_bound in (("period", gap["period"]),
                                              ("2N", 2 * integer)):
                for kind, sequence in (("raw", flags), ("entry", entry)):
                    item = conditional_record(sequence, order, constant, period_bound)
                    item.update({"order_name": name, "kind": kind,
                                 "period_bound_name": period_name})
                    windows.append(item)
    fractions = []
    for mean_count in (1, 2, 4, 8):
        for kind, sequence in (("raw", flags), ("entry", entry)):
            item = empty_fraction_record(sequence, mean_count)
            item.update({"kind": kind})
            fractions.append(item)
    record.update({
        "complete_period": True,
        "entry": gap,
        "conditional_windows": windows,
        "empty_window_fractions": fractions,
        "first_hit_phase": 0,
        "entry_flags_sha256": hashlib.sha256(entry).hexdigest(),
        "flags_sha256": hashlib.sha256(flags).hexdigest(),
        "inverse_symmetry_checked": record.get("inverse_symmetry_checked"),
    })
    return record


def complete_orbit_cases(period_cap: int, include_large: bool) -> list[dict]:
    specifications = ((1009 * 1013, 5, 100), (10007 * 10009, 7, 1000),
                      (1000003 * 1000033, 19, 10000))
    cases = []
    for integer, cutoff, cap in specifications:
        if not include_large and integer > 10 ** 9:
            continue
        delta = -integer if integer % 4 == 3 else -4 * integer
        base = build_prime_form_factor_base(delta, cutoff)
        generators = [{"kind": "prime", "norm": item.prime, "form": item.positive}
                      for item in base["records"][:2]]
        closure = factor_with_cofactor_ideal_closure(
            integer, cutoff, 32, max_cofactor=cap, max_generators=256,
            max_queries=10000)
        residuals = [item for item in closure["generators"]
                     if item["kind"] == "residual"][:2]
        generators += [{key: item[key] for key in ("kind", "norm", "form", "key")}
                       for item in residuals]
        for generator in generators:
            if (integer == 1000036000099 and generator["kind"] == "prime"
                    and generator["norm"] == 5):
                continue
            item = orbit_record(integer, cutoff, cap, generator, period_cap)
            if item is not None:
                cases.append(item)
    return cases


def raw_prefix_records() -> tuple[list[dict], dict]:
    records = []
    validation_errors = []
    key_map = {}
    file_records = []
    for path in sorted(glob.glob("tools/window_lemma/data/scan_*.jsonl")):
        digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        count = 0
        with open(path) as handle:
            for line_number, line in enumerate(handle, 1):
                data = json.loads(line)
                integer, cutoff = int(data["N"]), int(data["cutoff"])
                for target in data["targets"]:
                    prime = int(target["p"])
                    hits = list(map(int, target["hits"]))
                    scanned = int(target["scanned"])
                    key = (integer, cutoff, prime)
                    if key in key_map:
                        if key_map[key] != hits:
                            validation_errors.append(("conflicting_prefix", key))
                    key_map[key] = hits
                    if hits != sorted(set(hits)):
                        validation_errors.append(("unsorted_or_duplicate_hits", key))
                    if any(hit < 0 or hit >= scanned for hit in hits):
                        validation_errors.append(("hit_outside_prefix", key))
                    if scanned < 1:
                        validation_errors.append(("empty_prefix", key))
                    count += 1
                    records.append({
                        "integer": integer,
                        "bits": int(data["bits"]),
                        "cutoff": cutoff,
                        "prime": prime,
                        "scanned": scanned,
                        "hits": hits,
                        "censored": bool(target.get("censored")),
                        "start_exponent": int(target["K_p"]),
                        "target_hits": None,
                    })
        file_records.append({"path": path, "sha256": digest, "records": count})
    return records, {"files": file_records,
                     "distinct_keys": len(key_map),
                     "validation_errors": validation_errors}


def raw_prefix_summary(records: list[dict], integrity: dict) -> dict:
    counts = []
    start_offsets = []
    scanned_values = []
    primes = set()
    n_values = set()
    for item in records:
        hits, scanned = item["hits"], item["scanned"]
        if not hits or scanned < 1:
            continue
        counts.append(len(hits))
        start_offsets.append(hits[0])
        scanned_values.append(scanned)
        primes.add(item["prime"])
        n_values.add(item["integer"])

    def median(values):
        ordered = sorted(values)
        if not ordered:
            return None
        middle = len(ordered) // 2
        if len(ordered) % 2:
            return ordered[middle]
        return (ordered[middle - 1] + ordered[middle]) / 2

    summary = {
        "records": len(records),
        "usable_records": len(counts),
        "distinct_integers": len(n_values),
        "distinct_primes": len(primes),
        "hits_total": sum(counts),
        "scanned_total": sum(scanned_values),
        "prefix_density": sum(counts) / sum(scanned_values),
        "hits_per_record_median": median(counts),
        "first_hit_offset_median": median(start_offsets),
        "first_hit_offset_max": max(start_offsets) if start_offsets else None,
        "scanned_median": median(scanned_values),
        "scanned_max": max(scanned_values) if scanned_values else None,
        "censored_records": sum(item["censored"] for item in records),
        "trailing_zero_gap_max": max(
            (item["scanned"] - 1 - item["hits"][-1]
             for item in records if item["hits"]), default=None),
        "start_exponent_median": median([item["start_exponent"] for item in records]),
        "start_exponent_max": max(item["start_exponent"] for item in records),
        "integrity": integrity,
    }
    mean_counts = []
    for mean_count in (1, 2, 4, 8):
        considered = skipped = windows = empty = records_with_empty = 0
        initial_empty = 0
        widths = []
        for item in records:
            hits, scanned = item["hits"], item["scanned"]
            if not hits or scanned < 1:
                continue
            width = max(1, ceil_div(mean_count * scanned, len(hits)))
            if width > scanned:
                skipped += 1
                continue
            considered += 1
            widths.append(width)
            # Zero runs of the prefix: before the first hit, between hits,
            # and after the last hit.  A zero run of length g contains
            # max(0, g-width+1) empty non-wrapping windows.
            zero_runs = [hits[0]]
            for previous, current in zip(hits, hits[1:]):
                zero_runs.append(current - previous - 1)
            zero_runs.append(scanned - 1 - hits[-1])
            empty_here = sum(max(0, gap - width + 1) for gap in zero_runs)
            total_here = scanned - width + 1
            windows += total_here
            empty += empty_here
            records_with_empty += empty_here > 0
            initial_empty += hits[0] >= width
        mean_counts.append({
            "mean_count": mean_count,
            "considered_records": considered,
            "skipped_wrapping_or_too_long": skipped,
            "windows": windows,
            "empty_windows": empty,
            "empty_fraction": empty / windows if windows else None,
            "records_with_any_empty_window": records_with_empty,
            "initial_K_p_window_empty": initial_empty,
            "initial_K_p_empty_fraction": initial_empty / considered if considered else None,
            "window_median": median(widths),
            "window_max": max(widths) if widths else None,
        })
    summary["mean_count_windows"] = mean_counts
    return summary


def synthetic_controls() -> dict:
    checked = 0
    clustered = 0
    maximum_period = 12
    for period in range(2, maximum_period + 1):
        for mask in range(1 << period):
            flags = bytes((mask >> index) & 1 for index in range(period))
            good = sum(flags)
            if not 0 < good < period:
                continue
            gap = entry_gap_record(flags)
            entry = entry_flags(flags)
            max_gap = gap["max_entry_gap"]
            for window in range(1, period + 1):
                empty = min(cyclic_window_counts(entry, window)) == 0
                assert empty == (window < max_gap)
            checked += 1
    for period in (17, 31, 32, 63, 64):
        for run_length in range(1, min(period // 2, 9) + 1):
            flags = bytes([1] * run_length + [0] * (period - run_length))
            gap = entry_gap_record(flags)
            entry = entry_flags(flags)
            for window in (1, max(1, gap["max_entry_gap"] - 1), gap["max_entry_gap"]):
                empty = min(cyclic_window_counts(entry, window)) == 0
                assert empty == (window < gap["max_entry_gap"])
            clustered += 1
    return {"short_words_tested": checked, "clustered_words_tested": clustered,
            "maximum_short_period": maximum_period,
            "entry_gap_equivalence_checked": True}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--period-cap", type=int, default=1_000_000)
    parser.add_argument("--skip-large", action="store_true")
    args = parser.parse_args()

    controls = synthetic_controls()
    cases = complete_orbit_cases(args.period_cap, not args.skip_large)
    raw_records, integrity = raw_prefix_records()
    raw_summary = raw_prefix_summary(raw_records, integrity)
    assert not integrity["validation_errors"], integrity["validation_errors"][:3]

    nonwrapping_entry_c1 = []
    nonwrapping_raw_c1 = []
    for case in cases:
        if not case["complete_period"]:
            continue
        for item in case["conditional_windows"]:
            if (item["order_name"] == "bit_length" and item["constant"] == 1
                    and item["period_bound_name"] == "2N" and not item["wrapping"]):
                if item["kind"] == "entry":
                    nonwrapping_entry_c1.append(item["lower_tail_moment_bound_holds"])
                else:
                    nonwrapping_raw_c1.append(item["lower_tail_moment_bound_holds"])
    complete = [case for case in cases if case["complete_period"]]
    summary = {
        "complete_orbits": len(complete),
        "incomplete_orbits": len(cases) - len(complete),
        "lemma_checks_all_pass": True,
        "entry_c1_nonwrapping_cases": len(nonwrapping_entry_c1),
        "entry_c1_nonwrapping_holds": sum(nonwrapping_entry_c1),
        "raw_c1_nonwrapping_cases": len(nonwrapping_raw_c1),
        "raw_c1_nonwrapping_holds": sum(nonwrapping_raw_c1),
        "max_entry_gap_over_bits": max(
            (case["entry"]["max_entry_gap"] / case["input_bits"] for case in complete),
            default=None),
        "max_mean_cluster": max(
            (case["entry"]["mean_cluster"] for case in complete), default=None),
        "max_cluster": max((case["entry"]["max_cluster"] for case in complete),
                           default=None),
        "max_cluster_loss_factor": max(
            (case["entry"]["entry_density_denominator"]
             / case["entry"]["raw_density_denominator"] for case in complete),
            default=None),
        "all_max_entry_gap_at_most_2b": all(
            case["entry"]["max_entry_gap"] <= 2 * case["input_bits"]
            for case in complete),
        "all_mean_cluster_at_most_b": all(
            case["entry"]["mean_cluster"] <= case["input_bits"] for case in complete),
        "reserved_cycle_excluded": "N=1000036000099,Y=19,g=5",
        "reserved_cycle_owner_check": "worker 5",
    }
    result = {
        "claim": "Exact entry-gap equivalence and finite entry/cluster diagnostics for Path C",
        "scope": ("Complete independent cofactor orbits exclude the two owner-named cycles; "
                  "raw prefixes are censored and no uniform theorem is claimed."),
        "synthetic_controls": controls,
        "summary": summary,
        "complete_cases": cases,
        "raw_prefix": raw_summary,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(json.dumps({"complete_orbits": len(complete),
                      "entry_c1_nonwrapping": [summary["entry_c1_nonwrapping_holds"],
                                               summary["entry_c1_nonwrapping_cases"]],
                      "raw_c1_nonwrapping": [summary["raw_c1_nonwrapping_holds"],
                                             summary["raw_c1_nonwrapping_cases"]],
                      "raw_records": raw_summary["records"],
                      "raw_empty_fraction_c4": next(
                          item["empty_fraction"] for item in raw_summary["mean_count_windows"]
                          if item["mean_count"] == 4)}, sort_keys=True))


if __name__ == "__main__":
    main()
