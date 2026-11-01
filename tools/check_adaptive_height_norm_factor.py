#!/usr/bin/env python3
"""Exact certificates and finite comparisons for adaptive target powers."""

import argparse
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path

from adaptive_height_norm_factor import factor_with_adaptive_height_norms
from arithmetic_two_group import power
from binary_quadratic_forms import reduce_form
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from smooth_norm_words import exponent_vector, norm_to_word
from weighted_norm_relations import factor_from_weighted_relations, least_power_above


ROOT = Path(__file__).resolve().parents[1]


def canonical_base(integer, cutoff):
    records = build_prime_form_factor_base(-4 * integer, cutoff)["records"]
    primes, generators = [], []
    for record in records:
        prime = record.prime
        root = min(record.root, 2 * prime - record.root)
        primes.append(prime)
        generators.append(reduce_form(
            (prime, root, (root * root + 4 * integer) // (4 * prime))))
    return primes, generators


def brute_acceptable_norm(form, bound, primes):
    """Independent square enumeration, used only on small failed targets.

    A reduced form satisfies Q(x,y)>=A*(x*x+y*y)/2. This supplies a
    complete square that is deliberately different from the ellipse scan.
    """
    a, b, c = form
    discriminant = b * b - 4 * a * c
    radius = isqrt((2 * bound) // a)
    for x in range(-radius, radius + 1):
        for y in range(-radius, radius + 1):
            if gcd(x, y) != 1:
                continue
            norm = a * x * x + b * x * y + c * y * y
            if not 0 < norm <= bound or gcd(norm, discriminant) != 1:
                continue
            remaining = norm
            for prime in primes:
                while remaining % prime == 0:
                    remaining //= prime
            if remaining == 1:
                return norm, (x, y)
    return None


def verify_record(record, *, exhaust_failed=False):
    integer, height = record["integer"], record["height_bound"]
    if record["method"] == "factor-base gcd":
        assert 1 < record["factor"] < integer and integer % record["factor"] == 0
        return {"verified_rows": 0, "lower_power_rows": 0,
                "failed_powers_exhausted": 0, "ordinary_dominance_failures": 0}
    primes, generators = canonical_base(integer, record["prime_cutoff"])
    assert primes == record["factor_base_primes"]
    identity, rows, pivots = (1, 0, integer), [], []
    lower_powers = ordinary_failures = exhausted = 0
    for witness in record["witnesses"]:
        index = witness["generator_index"]
        prime = primes[index]
        exponent = witness["target_exponent"]
        maximum = least_power_above(prime, height)
        assert prime == witness["generator_prime"]
        assert 1 <= exponent <= maximum == witness["maximum_exponent"]
        target = power(generators[index], exponent, identity, compose_forms_hnf)
        assert target == tuple(witness["target_form"])
        x, y = witness["primitive_vector"]
        norm = target[0] * x * x + target[1] * x * y + target[2] * y * y
        assert norm == witness["norm"] and 0 < norm <= height and norm < prime**exponent
        assert gcd(x, y) == 1 and gcd(norm, -4 * integer) == 1
        word = norm_to_word(target, x, y, record["prime_cutoff"])
        assert word is not None and word.primitive_norm == norm
        assert [(term.prime, term.exponent) for term in word.terms] == [
            tuple(term) for term in witness["signed_prime_word"]]
        encoded = exponent_vector(word, primes)
        row = tuple(exponent * (column == index) - entry
                    for column, entry in enumerate(encoded))
        rows.append(row)
        pivots.append(index)
        lower_powers += exponent < maximum
        ordinary_failures += abs(row[index]) <= sum(
            abs(entry) for column, entry in enumerate(row) if column != index)
    assert len(rows) == record["accepted_relations"] == record["stats"]["relations_verified"]
    assert pivots == record["relation_pivot_columns"] and len(set(pivots)) == len(pivots)
    # This invocation re-evaluates every class-group relation independently
    # of the production engine's accepted-row cache.
    replay = factor_from_weighted_relations(
        integer, generators, rows, primes, pivots,
        height_bound=height, verify_relations=True)
    if record["factor"] is not None:
        assert 1 < record["factor"] < integer and integer % record["factor"] == 0
        if "weighted relation halving" in record["method"]:
            assert replay["factor"] == record["factor"]
    else:
        assert replay["factor"] is None
    if exhaust_failed:
        for failed in record["failed_targets"]:
            index = failed["generator_index"]
            prime = primes[index]
            maximum = least_power_above(prime, height)
            assert maximum == failed["maximum_exponent"]
            for exponent in range(1, maximum + 1):
                target = power(generators[index], exponent, identity, compose_forms_hnf)
                bound = min(height, prime**exponent - 1)
                assert brute_acceptable_norm(target, bound, primes) is None
                exhausted += 1
    return {"verified_rows": len(rows), "lower_power_rows": lower_powers,
            "failed_powers_exhausted": exhausted,
            "ordinary_dominance_failures": ordinary_failures,
            "replayed_extraction": replay}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    controls = []
    for integer, cutoff, height, interval in (
        (1009, 13, 1, 0), (1009, 13, 64, 0), (1009, 13, 128, 0),
        (1013, 3, 64, 0), (1021, 7, 64, 0),
        (10807, 7, 64, 0), (13493, 3, 256, 0), (1022117, 3, 1024, 0),
        (10807, 31, 512, 0), (13493, 13, 512, 1),
        (1022117, 127, 2**16, 4), (15, 31, 128, 0), (3, 3, 1, 0),
    ):
        record = factor_with_adaptive_height_norms(
            integer, cutoff, height, extraction_interval=interval,
            keep_certificates=True)
        record["independent_verification"] = verify_record(record, exhaust_failed=True)
        controls.append(record)
    assert any(record["failed_targets"] for record in controls)
    assert any(record["factor"] is not None for record in controls)
    assert any(record["stats"]["leading_norm_successes"] for record in controls)
    rejections = 0
    for integer, cutoff, height, interval in (
        (2, 3, 1, 0), (1, 3, 1, 0), (15, 2, 1, 0),
        (15, 3, 0, 0), (15, 3, 1, -1),
    ):
        try:
            factor_with_adaptive_height_norms(
                integer, cutoff, height, extraction_interval=interval)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError("An invalid search request was accepted.")
    baseline_path = ROOT / "evidence/height-selected-norm-factoring.json"
    baseline = json.loads(baseline_path.read_text())
    records = []
    for previous in baseline["records"][:1 if args.small_only else None]:
        record = factor_with_adaptive_height_norms(
            previous["integer"], previous["prime_cutoff"], previous["height_bound"],
            keep_certificates=True)
        verification = verify_record(record)
        assert record["factor"] is not None
        assert record["accepted_relations"] == len(record["factor_base_primes"])
        assert verification["lower_power_rows"] > 0
        prior_points = previous["search"]["coordinate_points_visited"]
        assert record["stats"]["ellipse_coordinate_points"] < prior_points
        record["independent_verification"] = verification
        record["fixed_power_comparison"] = {
            "source": str(baseline_path.relative_to(ROOT)),
            "same_integer_cutoff_and_height": True,
            "coordinate_points_visited": prior_points,
            "search": previous["search"],
            "scope": "Adaptive powers plus leading tests and the new cap schedule; no isolated attribution or asymptotic claim",
        }
        records.append(record)
    result = {
        "status": "pass",
        "scope": "Adaptive positive target powers with exact weighted relations; all finite search work charged; uniform fast selection remains open",
        "malformed_request_rejections": rejections,
        "small_controls": controls,
        "records": records,
        "engine_sha256": hashlib.sha256((ROOT / "tools/adaptive_height_norm_factor.py").read_bytes()).hexdigest(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "baseline_sha256": hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
        print(json.dumps({"status": "pass", "records": len(records),
                          "verified_rows": sum(record["independent_verification"]["verified_rows"]
                                               for record in records + controls),
                          "failed_powers_exhausted": sum(record["independent_verification"]["failed_powers_exhausted"]
                                                          for record in controls)}))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
