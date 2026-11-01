#!/usr/bin/env python3
"""Replay composite-ideal closure words and exact weighted dominance."""

from __future__ import annotations

import argparse
import json
from math import gcd, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form, reduce_form
from cofactor_collision_relations import RationalRowBasis, verify_partial_identity
from cofactor_ideal_closure import factor_with_cofactor_ideal_closure
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above


def normalize_record(result: dict) -> dict:
    result["prime_records"] = {int(prime): (record[0], tuple(record[1]))
                               for prime, record in result["prime_records"].items()}
    for generator in result["generators"]:
        generator["form"] = tuple(generator["form"])
        if generator["key"] is not None:
            generator["key"] = tuple(generator["key"])
        if generator["definition"] is not None:
            generator["definition"]["word"] = {int(prime): exponent for prime, exponent
                                                in generator["definition"]["word"].items()}
    for witness in result["witnesses"]:
        witness["target"] = tuple(witness["target"])
        witness["word"] = {int(prime): exponent for prime, exponent in witness["word"].items()}
        if witness["key"] is not None:
            witness["key"] = tuple(witness["key"])
            witness["residual_form"] = tuple(witness["residual_form"])
    return result


def verify_result(result: dict, *, large: bool = False) -> dict:
    n, delta = result["integer"], result["discriminant"]
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    generators = result["generators"]
    rows = [dict(row) for row in result["sparse_augmented_relation_matrix"]]
    basis = None if large else RationalRowBasis(len(generators))
    prime_indices = {prime: index for index, prime in enumerate(result["active_primes"])}
    words = []
    for index, generator in enumerate(generators):
        definition = generator["definition"]
        if definition is None:
            assert index < result["initial_base_size"] and generator["depth"] == 0
            prime = result["active_primes"][index]
            root, form = result["prime_records"][prime]
            assert generator["norm"] == prime and generator["form"] == form
            assert form == reduce_form((prime, root, (root * root - delta) // (4 * prime)))
            words.append([int(j == index) for j in range(result["initial_base_size"])])
        else:
            norm, root = generator["key"]
            assert norm == generator["norm"] and gcd(norm, 2 * delta) == 1
            assert norm > 1 and 0 <= root <= norm and (root * root - delta) % (4 * norm) == 0
            assert generator["form"] == reduce_form((norm, root, (root * root - delta) // (4 * norm)))
            parent = definition["source"]
            assert parent < index and generator["depth"] == generators[parent]["depth"] + 1
            sign = definition["orientation"]
            word = [sign * definition["exponent"] * value for value in words[parent]]
            for prime, exponent in definition["word"].items():
                word[prime_indices[prime]] -= sign * exponent
            words.append(word)
    fixed_by_source = dict(zip(result["fixed_relation_sources"], result["relation_matrix"]))
    self_references = negative_orientations = 0
    assert result["new_residual_generators"] <= result["max_cofactor"] * result["max_cofactor"].bit_length()
    for index, (row, witness) in enumerate(zip(rows, result["witnesses"])):
        assert witness["generator_index"] == index
        generator = generators[index]
        exponent = witness["exponent"]
        assert exponent >= least_power_above(generator["norm"], result["height_bound"])
        assert power(generator["form"], exponent, identity, compose_forms_hnf) == tuple(witness["target"])
        assert verify_partial_identity(n, witness["target"], witness, result["prime_records"])
        assert witness["smooth_part"] * witness["residual"] == witness["target"][0]
        assert prod(prime ** abs(value) for prime, value in witness["word"].items()) == witness["smooth_part"]
        expected = {index: exponent}
        for prime, value in witness["word"].items():
            column = prime_indices[prime]
            expected[column] = expected.get(column, 0) - value
        if witness["residual_index"] is not None:
            column = witness["residual_index"]
            expected[column] = expected.get(column, 0) - witness["orientation"]
        assert row == {column: value for column, value in expected.items() if value}
        other = prod(generators[column]["norm"] ** abs(value)
                     for column, value in row.items() if column != index)
        assert row.get(index, 0) >= least_power_above(generator["norm"], other)
        if not large:
            assert generator["norm"] ** row[index] > other
            value = identity
            for column, coefficient in row.items():
                form = generators[column]["form"]
                if coefficient < 0:
                    form, coefficient = inverse_form(form), -coefficient
                value = compose_forms_hnf(value, power(form, coefficient, identity, compose_forms_hnf))
            assert value == identity
            assert basis.add([row.get(column, 0) for column in range(len(generators))])
        eliminated = [sum(coefficient * words[column][j] for column, coefficient in row.items())
                      for j in range(result["initial_base_size"])]
        assert eliminated == fixed_by_source.get(index, [0] * result["initial_base_size"])
        if witness["residual_index"] is not None:
            residual = generators[witness["residual_index"]]
            assert residual["key"] == witness["key"] and residual["norm"] == witness["residual"]
            assert residual["form"] == witness["residual_form"]
            self_references += witness["residual_index"] == index
            negative_orientations += witness["orientation"] < 0
    assert result["augmented_relation_rank"] == len(rows)
    if not large:
        assert len(basis.rows) == len(rows)
    fixed_basis = None if large else RationalRowBasis(result["initial_base_size"])
    for row in result["relation_matrix"]:
        if not large:
            assert fixed_basis.add(row)
            value = identity
            for column, coefficient in enumerate(row):
                form = generators[column]["form"]
                if coefficient < 0:
                    form, coefficient = inverse_form(form), -coefficient
                value = compose_forms_hnf(value, power(form, coefficient, identity, compose_forms_hnf))
            assert value == identity
    assert len(result["relation_matrix"]) == result["relation_rank"] == len(rows) - result["new_residual_generators"]
    assert result["closed"] == (len(rows) == len(generators))
    factor = result["factor"]
    if factor is not None:
        assert 1 < factor < n and n % factor == 0
    for field in ("factoring_oracle_calls", "residual_primality_tests", "new_primes_introduced"):
        assert result[field] == 0
    return {"literal_relation_checks": 0 if large else len(rows),
            "partial_identity_checks": len(rows), "weighted_dominance_checks": len(rows),
            "fixed_base_relation_checks": len(result["relation_matrix"]),
            "integer_elimination_replays": len(rows),
            "self_references": self_references, "negative_orientations": negative_orientations}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--replay", type=Path)
    args = parser.parse_args()
    if args.replay:
        record = normalize_record(json.loads(args.replay.read_text()))
        result = {"claim": "Replay of stored composite-ideal closure certificates",
                  "integer": record["integer"], "factor": record["factor"],
                  "norm_queries": record["norm_queries"],
                  "verification": verify_result(record, large=True),
                  "scope": "Rank certified by exact weighted dominance and integral elimination; no dense augmented rank computation."}
        text = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.write_text(text)
        print(text, end="")
        return
    specifications = ((1009 * 1013, 5), (1019 * 1021, 5), (1031 * 1033, 5),
                      (10007 * 10009, 7), (10037 * 10039, 7), (1000003 * 1000033, 19))
    cases = []
    counts = {"literal_relation_checks": 0, "weighted_dominance_checks": 0,
              "partial_identity_checks": 0,
              "fixed_base_relation_checks": 0,
              "integer_elimination_replays": 0,
              "self_references": 0, "negative_orientations": 0}
    for n, cutoff in specifications:
        result = factor_with_cofactor_ideal_closure(
            n, cutoff, 32, max_cofactor=10000, max_generators=256,
            max_queries=10000, verify_candidates=True)
        assert result["factor"] is not None
        verification = verify_result(result)
        for field, value in verification.items():
            counts[field] += value
        cases.append(result)
    # The engine receives no factorization of these residual norms.
    residual_norms = {item["norm"] for result in cases for item in result["generators"]
                      if item["kind"] == "residual"}
    known_composites = {49: (7, 7), 413: (7, 59), 91: (7, 13), 121: (11, 11), 77: (7, 11)}
    assert all(norm in residual_norms and prod(factors) == norm for norm, factors in known_composites.items())
    self_reference = factor_with_cofactor_ideal_closure(
        1013, 5, 1, max_generators=80, max_queries=100, verify_candidates=True)
    control = verify_result(self_reference)
    assert self_reference["closed"] and self_reference["factor"] is None and control["self_references"] >= 1
    for field, value in control.items():
        counts[field] += value
    failures = []
    for options in ({"max_cofactor": 1, "window": 1},
                    {"max_cofactor": 10000, "window": 32, "max_queries": 0},
                    {"max_cofactor": 10000, "window": 32, "max_generators": 2}):
        options = dict(options)
        window = options.pop("window")
        result = factor_with_cofactor_ideal_closure(1009 * 1013, 5, window, **options)
        assert result["factor"] is None and not result["closed"]
        verify_result(result)
        failures.append(result)
    nongenerating = factor_with_cofactor_ideal_closure(31 * 37, 7, 16, max_cofactor=100)
    verify_result(nongenerating)
    assert nongenerating["closed"] and nongenerating["factor"] is None
    result = {"claim": "Dominant composite-ideal closure without residual factorization",
              "successful_factoring_controls": len(cases), "verification": counts,
              "known_composite_residual_controls": known_composites,
              "cases": cases, "self_reference_control": self_reference,
              "budget_failure_controls": failures, "small_base_nongeneration_control": nongenerating,
              "scope": "Uniform admissible-hit and closure-size bounds are not proved; finite small bases need not generate the class group."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
