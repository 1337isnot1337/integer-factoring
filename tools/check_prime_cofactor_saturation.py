#!/usr/bin/env python3
"""Verify prime residual certificates, dominance, and DAG elimination."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd, prod
from pathlib import Path
from time import perf_counter

from arithmetic_two_group import power
from binary_quadratic_forms import reduce_form
from cofactor_collision_relations import RationalRowBasis, partial_leading_word, verify_partial_identity
from cofactor_ideal_closure import factor_with_cofactor_ideal_closure
from prime_cofactor_saturation import factor_with_prime_cofactor_saturation, row_digest
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]


def normalize_record(record):
    for generator in record["generators"]:
        generator["form"] = tuple(generator["form"])
        definition = generator["definition"]
        if definition is not None:
            definition["word"] = {int(p): v for p, v in definition["word"].items()}
    for witness in record["witnesses"]:
        witness["target"] = tuple(witness["target"])
        witness["word"] = {int(p): v for p, v in witness["word"].items()}
        if witness["key"] is not None:
            witness["key"] = tuple(witness["key"])
            witness["residual_form"] = tuple(witness["residual_form"])
        witness["factors"] = [tuple(item) for item in witness["factors"]]
    return record


def verify_result(record, *, small=False):
    n, delta = record["integer"], record["discriminant"]
    cutoff, cap, m = record["cutoff"], record["max_cofactor"], record["initial_base_size"]
    promotion_cutoff = record.get("promotion_cutoff", cap)
    assert 1 <= cap <= cutoff * cutoff
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {item.prime: (item.root, item.positive) for item in base["records"]}
    if n % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (n + 1) // 2)))
    primes = sorted(records)
    assert primes == record["initial_primes"] and m == len(primes)
    indices = {p: i for i, p in enumerate(primes)}
    prime_product = initial_product = prod(primes)
    generators = record["generators"]
    definitions = [None] * m
    for index, p in enumerate(primes):
        item = generators[index]
        assert item["norm"] == p and item["root"] == records[p][0]
        assert item["form"] == records[p][1] and item["definition"] is None and item["depth"] == 0
    rows = [dict(row) for row in record["sparse_augmented_relation_matrix"]]
    rank = RationalRowBasis(len(generators)) if small else None
    fixed, fixed_sources = [], []
    counts = {"partial_identity_checks": 0, "prime_residual_certificates": 0,
              "weighted_dominance_checks": 0, "integer_dag_eliminations": 0,
              "negative_residual_orientations": 0, "multiple_noninitial_dependencies": 0,
              "promoted_prime_occurrences": 0}
    maximum_depth = 0
    for source, (row, witness) in enumerate(zip(rows, record["witnesses"])):
        assert source == witness["generator_index"]
        generator = generators[source]
        first = least_power_above(generator["norm"], record["height_bound"])
        exponent = witness["exponent"]
        assert first <= exponent < first + record["window"]
        target = power(generator["form"], exponent, identity, compose_forms_hnf)
        assert target == witness["target"]
        split = partial_leading_word(n, target, primes, records, prime_product)
        for name in ("smooth_part", "residual", "word", "key", "orientation", "residual_form"):
            assert witness[name] == split[name]
        assert verify_partial_identity(n, target, split, records)
        counts["partial_identity_checks"] += 1
        expected = {source: exponent}
        for prime, value in split["word"].items():
            column = indices[prime]
            expected[column] = expected.get(column, 0) - value
            counts["promoted_prime_occurrences"] += prime > cutoff
        residual, new = split["residual"], False
        assert residual <= cap
        if residual > 1:
            assert residual > cutoff and gcd(residual, initial_product * 2 * delta) == 1
            assert (split["key"][1] ** 2 - delta) % (4 * residual) == 0
            # These elementary conditions prove primality: a hypothetical
            # prime factor <=Y would be an initial split prime already removed.
            if residual not in indices:
                index = len(indices)
                assert witness["residual_index"] == index
                item = generators[index]
                assert item["norm"] == residual and item["root"] == split["key"][1]
                assert item["form"] == split["residual_form"]
                dependencies = {j: split["orientation"] * v for j, v in expected.items() if v}
                definition = item["definition"]
                assert dict(definition["dependencies"]) == dependencies
                assert definition["source"] == source and definition["exponent"] == exponent
                assert definition["word"] == split["word"] and definition["orientation"] == split["orientation"]
                assert all(j < index for j in dependencies)
                depth = 1 + max((generators[j]["depth"] for j in dependencies), default=0)
                assert item["depth"] == depth
                maximum_depth = max(maximum_depth, depth)
                definitions.append(dependencies)
                counts["multiple_noninitial_dependencies"] += sum(j >= m for j in dependencies) > 1
                indices[residual] = index
                records[residual] = (item["root"], item["form"])
                counts["prime_residual_certificates"] += 1
                if record["promote_residuals"] and residual <= promotion_cutoff:
                    primes = sorted(primes + [residual])
                    prime_product *= residual
                new = True
            else:
                assert not (record["promote_residuals"] and residual <= promotion_cutoff)
                index = indices[residual]
                assert records[residual] == (split["key"][1], split["residual_form"])
            assert witness["residual_index"] == index
            expected[index] = expected.get(index, 0) - split["orientation"]
            counts["negative_residual_orientations"] += split["orientation"] < 0
        else:
            assert witness["residual_index"] is None
        expected = {j: v for j, v in expected.items() if v}
        assert expected == row
        other = prod(generators[j]["norm"] ** abs(v) for j, v in row.items() if j != source)
        assert row.get(source, 0) >= least_power_above(generator["norm"], other)
        counts["weighted_dominance_checks"] += 1
        if small:
            assert rank.add([row.get(j, 0) for j in range(len(generators))])
        if not new:
            # Independent dense-index reverse substitution; the production
            # engine instead uses a sparse priority queue.
            eliminated = dict(row)
            for j in range(len(definitions) - 1, m - 1, -1):
                value = eliminated.pop(j, 0)
                if value:
                    for previous, coefficient in definitions[j].items():
                        eliminated[previous] = eliminated.get(previous, 0) + value * coefficient
            eliminated = {j: v for j, v in eliminated.items() if v}
            assert eliminated and all(j < m for j in eliminated)
            fixed.append(eliminated)
            fixed_sources.append(source)
            counts["integer_dag_eliminations"] += 1
    assert len(indices) == len(generators) == record["generator_count"]
    assert len(rows) == record["augmented_relation_rank"] == record["admissible_queries"]
    assert record["relation_rank"] == len(fixed) == len(rows) - len(generators) + m
    assert record["new_residual_generators"] == len(generators) - m
    assert record["fixed_relation_sources"] == fixed_sources
    assert record["fixed_relation_sha256"] == row_digest(fixed)
    assert record["maximum_definition_depth"] == maximum_depth
    assert record["maximum_saturation_product_bits"] == prime_product.bit_length()
    assert record["saturation_prime_count"] == len(primes)
    if "sparse_fixed_relation_matrix" in record:
        assert [dict(row) for row in record["sparse_fixed_relation_matrix"]] == fixed
    if small:
        fixed_rank = RationalRowBasis(m)
        for row in fixed:
            assert fixed_rank.add([row.get(j, 0) for j in range(m)])
    if record["factor"] is not None:
        assert 1 < record["factor"] < n and n % record["factor"] == 0
    assert record["factoring_oracle_calls"] == record["residual_primality_tests"] == 0
    return counts


def summary(record):
    fields = ("input_bits", "factor", "norm_queries", "generator_count", "initial_base_size",
              "relation_rank", "augmented_relation_rank", "maximum_definition_depth",
              "maximum_relation_entry_bits", "maximum_elimination_coefficient_bits",
              "maximum_saturation_product_bits", "fixed_relation_sha256")
    return {field: record[field] for field in fields}


def source_hashes():
    return {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
            for name in (Path(__file__).name, "prime_cofactor_saturation.py",
                         "cofactor_collision_relations.py", "leading_power_smooth_relations.py")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--benchmark", type=int, choices=(92, 124, 144))
    parser.add_argument("--promotion-cutoff", type=int)
    args = parser.parse_args()
    if args.replay:
        raw = json.loads(args.replay.read_text())
        record = normalize_record(raw.get("dynamic", raw))
        report = {"claim": "Prime cofactor saturation witness replay", "summary": summary(record),
                  "verification": verify_result(record)}
    elif args.benchmark:
        n, y, window, maxq, maxg = {
            92: (4951760154835678088235319297, 1021, 10000, 300000, 20000),
            124: (21267647932558653899591465697288388633, 4093, 100000, 3000000, 20000),
            144: (17355670471043633865640786569418441084204301, 8191, 200000, 12000000, 40000)
        }[args.benchmark]
        fixed_start = perf_counter()
        fixed = factor_with_prime_cofactor_saturation(n, y, window, max_queries=maxq,
            max_generators=maxg, promote_residuals=False, verify_candidates=True, include_relations=False)
        fixed_seconds = perf_counter() - fixed_start
        published_path = ROOT / "evidence" / f"cofactor-ideal-closure-{args.benchmark}.json"
        published = json.loads(published_path.read_text())
        for field in ("factor", "norm_queries", "generator_count", "relation_rank",
                      "augmented_relation_rank", "sparse_augmented_relation_matrix"):
            assert json.loads(json.dumps(fixed[field])) == published[field]
        assert fixed["fixed_relation_sha256"] == row_digest([
            {j: value for j, value in enumerate(row) if value} for row in published["relation_matrix"]])
        dynamic_start = perf_counter()
        dynamic = factor_with_prime_cofactor_saturation(n, y, window, max_queries=maxq,
            max_generators=maxg, promote_residuals=True, promotion_cutoff=args.promotion_cutoff,
            verify_candidates=True, include_relations=False)
        dynamic_seconds = perf_counter() - dynamic_start
        report = {"claim": "Matched fixed-versus-promoted prime cofactor closure",
                  "fixed": summary(fixed), "dynamic": dynamic,
                  "fixed_seconds": fixed_seconds, "dynamic_seconds": dynamic_seconds,
                  "fixed_matches_published_augmented_and_eliminated_rows": True,
                  "published_fixed_record": published_path.name,
                  "published_fixed_sha256": sha256(published_path.read_bytes()).hexdigest(),
                  "verification": verify_result(dynamic), "source_sha256": source_hashes(),
                  "scope": "Finite matched traversals; no uniform factoring exponent improvement."}
    else:
        cases = []
        for n, y in ((1009 * 1013, 5), (1019 * 1021, 5), (1031 * 1033, 5),
                     (10007 * 10009, 7), (10037 * 10039, 7), (1000003 * 1000033, 19)):
            old = factor_with_cofactor_ideal_closure(n, y, 1000, max_cofactor=y*y, verify_candidates=True)
            fixed = factor_with_prime_cofactor_saturation(n, y, 1000, promote_residuals=False, verify_candidates=True)
            dynamic = factor_with_prime_cofactor_saturation(n, y, 1000, verify_candidates=True)
            selective = factor_with_prime_cofactor_saturation(n, y, 1000,
                promotion_cutoff=2*y, verify_candidates=True)
            for field in ("factor", "norm_queries", "generator_count", "relation_rank",
                          "sparse_augmented_relation_matrix"):
                assert old[field] == fixed[field]
            assert row_digest([{j: v for j, v in enumerate(row) if v} for row in old["relation_matrix"]]) == fixed["fixed_relation_sha256"]
            cases.append({"fixed": fixed, "dynamic": dynamic,
                          "selective": selective,
                          "fixed_verification": verify_result(fixed, small=True),
                          "selective_verification": verify_result(selective, small=True),
                          "dynamic_verification": verify_result(dynamic, small=True)})
        failures = []
        for options in ({"max_queries": 0}, {"max_generators": 1}, {"max_cofactor": 1}):
            result = factor_with_prime_cofactor_saturation(1009*1013, 5, 1, **options)
            assert result["factor"] is None
            verify_result(result, small=True)
            failures.append({"options": options, "method": result["method"]})
        try:
            factor_with_prime_cofactor_saturation(1009*1013, 5, 1, max_cofactor=26)
        except ValueError:
            pass
        else:
            raise AssertionError("An unsupported composite residual cap was accepted")
        nonmonotone_fixed = factor_with_prime_cofactor_saturation(163*167, 5, 200,
            promote_residuals=False, verify_candidates=True)
        nonmonotone_dynamic = factor_with_prime_cofactor_saturation(163*167, 5, 200,
            promote_residuals=True, verify_candidates=True)
        assert nonmonotone_fixed["norm_queries"] == 8
        assert nonmonotone_dynamic["norm_queries"] == 9
        assert nonmonotone_fixed["factor"] and nonmonotone_dynamic["factor"]
        report = {"claim": "Exact prime residual promotion and DAG elimination controls",
                  "cases": cases, "budget_failures": failures, "oversized_cap_rejected": True,
                  "nonmonotone_query_control": {
                      "fixed": nonmonotone_fixed, "dynamic": nonmonotone_dynamic,
                      "fixed_verification": verify_result(nonmonotone_fixed, small=True),
                      "dynamic_verification": verify_result(nonmonotone_dynamic, small=True)},
                  "source_sha256": source_hashes(), "scope": "Finite controls; no uniform admissibility claim."}
    encoded = json.dumps(report, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        args.output.write_text(encoded)
    if "dynamic" in report:
        print(json.dumps({"fixed": report["fixed"], "dynamic": summary(report["dynamic"]),
                          "verification": report["verification"]}, indent=2))
    else:
        print(json.dumps({key: value for key, value in report.items() if key != "cases"}, indent=2))


if __name__ == "__main__":
    main()
