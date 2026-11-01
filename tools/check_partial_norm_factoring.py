#!/usr/bin/env python3
"""Check partial norm search and replay its supplied relation certificates.

Small skipped-target controls independently enumerate every claimed empty
ellipse. Large fixtures are replayed from their integer representations;
their original exhaustive-search counters are not recomputed.
"""

import argparse
import json
from math import gcd
from pathlib import Path

from arithmetic_partial_norm_factor import factor_with_partial_short_norms
from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form
from check_smooth_norm_sieve import brute_witnesses
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from relation_two_saturation import factor_from_partial_relations
from smooth_norm_sieve import SmoothNormSieve
from smooth_norm_words import exponent_vector, norm_to_word


def replay(record, *, check_empty_ellipses=False):
    integer = record["integer"]
    discriminant, identity = -4 * integer, (1, 0, integer)
    bound, budget = record["smooth_bound"], record["word_budget"]
    diagonal = budget + 1
    assert record["discriminant"] == discriminant
    assert record["diagonal_parameter"] == diagonal
    assert not record["determinant_computed"]
    base = []
    for prime, root in SmoothNormSieve(discriminant, bound).prime_roots:
        root = root if root % 2 == 0 else root + prime
        root = min(root, 2 * prime - root)
        form = (prime, root, (root * root - discriminant) // (4 * prime))
        base.append((prime, reduce_with_transform(form)[0]))
    primes = [prime for prime, _ in base]
    generators = [form for _, form in base]
    assert primes == record["factor_base_primes"]
    rows, pivots = [], []
    for witness in record["witnesses"]:
        index = witness["generator_index"]
        assert index not in pivots
        assert witness["generator_prime"] == primes[index]
        target = power(generators[index], diagonal, identity, compose_forms_hnf)
        assert target == tuple(witness["target_form"])
        x, y = witness["primitive_vector"]
        a, b, c = target
        norm = a*x*x + b*x*y + c*y*y
        assert norm == witness["represented_norm"] > 0
        assert gcd(x, y) == 1 and gcd(norm, discriminant) == 1
        residual, multiplicity = norm, 0
        for prime in primes:
            while residual % prime == 0:
                residual //= prime
                multiplicity += 1
        assert residual == 1 and multiplicity == witness["multiplicity"] <= budget
        word = norm_to_word(target, x, y, bound, budget)
        assert word is not None and word.primitive_norm == norm
        assert [[term.prime, term.exponent] for term in word.terms] == [list(term) for term in witness["prime_word"]]
        represented_class = identity
        for term in word.terms:
            value = term.form if term.exponent >= 0 else inverse_form(term.form)
            represented_class = compose_forms_hnf(
                represented_class,
                power(value, abs(term.exponent), identity, compose_forms_hnf),
            )
        assert represented_class == target
        exponents = exponent_vector(word, primes)
        rows.append(tuple(diagonal * (column == index) - exponent
                          for column, exponent in enumerate(exponents)))
        pivots.append(index)
    assert pivots == record["relation_pivot_columns"]
    assert len(rows) == record["relation_rows_completed"]
    if "relation_matrix" in record:
        assert rows == [tuple(row) for row in record["relation_matrix"]]

    empty_controls = 0
    all_targets = record["witnesses"] + record["failed_targets"]
    assert sorted(target["generator_index"] for target in all_targets) == list(range(record["targets_visited"]))
    counters = {key: 0 for key in ("coordinate_points_visited", "modular_strikes", "prime_divisions")}
    for target_record in all_targets:
        index = target_record["generator_index"]
        form = power(generators[index], diagonal, identity, compose_forms_hnf)
        assert form == tuple(target_record["target_form"])
        cap = 4 * max(form[0], form[2])
        attempts = target_record["search_attempts"]
        for number, attempt in enumerate(attempts):
            assert attempt["maximum_norm"] == cap * 4**number
            for key in counters:
                counters[key] += attempt[key]
            if attempt["exhausted"]:
                assert attempt["vector"] is None and attempt["norm"] is None
                if check_empty_ellipses:
                    assert not any(value[3] <= budget for value in brute_witnesses(
                        form, attempt["maximum_norm"], bound))
                    empty_controls += 1
            else:
                assert number == len(attempts) - 1
                assert list(attempt["vector"]) == list(target_record["primitive_vector"])
                assert attempt["norm"] == target_record["represented_norm"]
        if index not in pivots:
            assert len(attempts) == record["maximum_expansions_per_target"] + 1
            assert all(attempt["exhausted"] for attempt in attempts)
    assert len(record["failed_targets"]) == record["failed_target_count"]
    assert all(counters[key] == record[key] for key in counters)

    # Every row identity was established above through its oriented word.
    # The independent extractor also rechecks the designated dominant minor.
    extraction = factor_from_partial_relations(
        integer, generators, rows, pivots, verify_relations=False)
    assert json.loads(json.dumps(extraction)) == json.loads(json.dumps(record["relation_extraction"]))
    factor = extraction["factor"]
    assert (factor is not None) == (record["status"] == "factor")
    if factor is not None:
        assert factor == record["factor"] and 1 < factor < integer
        assert integer % factor == 0
    return {
        "integer": integer, "input_bits": integer.bit_length(),
        "factor_base_size": len(base), "verified_rows": len(rows),
        "skipped_targets": record["failed_target_count"],
        "empty_ellipses_independently_checked": empty_controls,
        "factor": factor, "cofactor": None if factor is None else integer // factor,
        "original_coordinate_points": record["coordinate_points_visited"],
        "extraction_stats": extraction["stats"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificates", type=Path,
                        default=Path(__file__).resolve().parents[1] / "evidence" / "partial-norm-factoring-certificates.json")
    args = parser.parse_args()
    controls = []
    for integer in (1022117, 100160063, 10002200057):
        record = factor_with_partial_short_norms(
            integer, 127, extraction_interval=4, keep_relations=True)
        assert record["status"] == "factor"
        controls.append(replay(record))
    skipped_records = []
    for integer, bound in ((10807, 31), (13493, 13)):
        record = factor_with_partial_short_norms(
            integer, bound, maximum_expansions=0,
            extraction_interval=1, keep_relations=True)
        assert record["status"] == "factor" and record["failed_target_count"] > 0
        assert record["relation_rows_completed"] < len(record["factor_base_primes"])
        controls.append(replay(record, check_empty_ellipses=True))
        skipped_records.append(record)
    prime = factor_with_partial_short_norms(1009, 31, extraction_interval=1, keep_relations=True)
    assert prime["status"] == "unresolved"
    controls.append(replay(prime))
    direct = factor_with_partial_short_norms(15, 3)
    assert direct["status"] == "factor" and direct["factor"] == 3
    empty = factor_with_partial_short_norms(7, 3)
    assert empty["status"] == "unresolved" and empty["reason"] == "empty finite factor base"
    zero_budget = factor_with_partial_short_norms(
        10807, 31, word_budget=0, maximum_expansions=0,
        extraction_interval=1, keep_relations=True)
    assert zero_budget["status"] == "unresolved"
    controls.append(replay(zero_budget, check_empty_ellipses=True))
    certificates = []
    if args.certificates.exists():
        data = json.loads(args.certificates.read_text())
        certificates = [replay(record) for record in data["records"]]
    print(json.dumps({
        "scope": "Exact finite partial-relation factoring; large supplied norms are replayed without repeating their search",
        "complete_factoring_score_improved": False,
        "small_controls": controls, "skipped_target_records": skipped_records,
        "direct_gcd_control": direct, "empty_factor_base_control": empty,
        "supplied_certificates": certificates,
    }, indent=2))


if __name__ == "__main__":
    main()
