#!/usr/bin/env python3
"""Promote proven prime residuals into the cofactor saturation product.

At cap B <= Y^2, every nontrivial residual is prime: all its prime
factors are split, unramified, and above the initial cutoff Y.
No unknown factorization or residual primality test is performed.
"""

from __future__ import annotations

from hashlib import sha256
from heapq import heapify, heappop, heappush
import json
from math import gcd, isqrt, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from cofactor_collision_relations import partial_leading_word, verify_partial_identity
from leading_power_smooth_relations import compose_with_prime_form
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def row_digest(rows: list[dict[int, int]]) -> str:
    digest = sha256()
    for row in rows:
        digest.update(json.dumps(sorted(row.items()), separators=(",", ":")).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def eliminate_definitions(row: dict[int, int], generators: list[dict],
                          initial_size: int) -> tuple[dict[int, int], int, int]:
    """Eliminate a sparse row by descending indices in the definition DAG."""
    result = {index: value for index, value in row.items() if value}
    pending = [-index for index in result if index >= initial_size]
    heapify(pending)
    visits = 0
    maximum_bits = max((abs(value).bit_length() for value in result.values()), default=0)
    while pending:
        index = -heappop(pending)
        coefficient = result.pop(index, 0)
        if not coefficient:
            continue
        visits += 1
        definition = generators[index]["definition"]
        if definition is None:
            raise ArithmeticError("A new generator lacks its unit definition")
        for previous, multiplier in definition["dependencies"]:
            assert previous < index
            old = result.get(previous, 0)
            value = old + coefficient * multiplier
            if value:
                result[previous] = value
                maximum_bits = max(maximum_bits, abs(value).bit_length())
                if previous >= initial_size and not old:
                    heappush(pending, -previous)
            elif old:
                del result[previous]
    assert all(index < initial_size for index in result)
    return result, visits, maximum_bits


def factor_with_prime_cofactor_saturation(integer: int, cutoff: int, window: int, *,
                                          max_cofactor: int | None = None,
                                          max_generators: int | None = None,
                                          max_queries: int | None = None,
                                          promote_residuals: bool = True,
                                          promotion_cutoff: int | None = None,
                                          verify_candidates: bool = False,
                                          include_relations: bool = True,
                                          smooth_batch_size: int = 1,
                                          progress=None) -> dict:
    if not isinstance(integer, int) or integer < 3 or integer % 2 == 0:
        raise ValueError("Require an odd integer of at least three")
    if not isinstance(cutoff, int) or cutoff < 2 or not isinstance(window, int) or window < 1:
        raise ValueError("Require cutoff >= 2 and a positive window")
    if max_cofactor is None:
        max_cofactor = cutoff * cutoff
    if not isinstance(max_cofactor, int) or not 1 <= max_cofactor <= cutoff * cutoff:
        raise ValueError("The residual cap must lie between one and cutoff squared")
    if promotion_cutoff is None:
        promotion_cutoff = max_cofactor
    if not isinstance(promotion_cutoff, int) or promotion_cutoff < 0:
        raise ValueError("The promotion cutoff must be a nonnegative integer")
    if not isinstance(promote_residuals, bool) or not isinstance(include_relations, bool):
        raise ValueError("Promotion and relation-storage options must be boolean")
    if not isinstance(smooth_batch_size, int) or smooth_batch_size < 1:
        raise ValueError("The smoothness batch size must be a positive integer")
    for name, value in (("max_generators", max_generators), ("max_queries", max_queries)):
        if value is not None and (not isinstance(value, int) or value < 0):
            raise ValueError(f"{name} must be nonnegative")
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {item.prime: (item.root, item.positive) for item in base["records"]}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    initial_primes = sorted(records)
    initial_size = len(initial_primes)
    indices = {prime: index for index, prime in enumerate(initial_primes)}
    saturation_primes = list(initial_primes)
    prime_product = prod(saturation_primes)
    generators = [{"kind": "prime", "norm": prime, "root": records[prime][0],
                   "form": records[prime][1], "definition": None, "depth": 0}
                  for prime in initial_primes]
    rows, witnesses, fixed_rows, fixed_sources = [], [], [], []
    parity_pivots = {}
    last_extraction, last_extraction_rows = None, -1
    counters = {"norm_queries": 0, "cofactor_rejections": 0, "admissible_queries": 0,
                "new_residual_generators": 0, "trial_divisions": 0,
                "initial_power_computations": 0, "successive_power_compositions": 0,
                "extractions": 0, "maximum_exponent_tested": 0,
                "maximum_residual_bits": 0, "maximum_definition_depth": 0,
                "elimination_vertices": 0, "maximum_elimination_coefficient_bits": 0,
                "maximum_relation_entry_bits": 0, "partial_identities_verified": 0,
                "maximum_saturation_product_bits": prime_product.bit_length()}
    if smooth_batch_size > 1:
        from batch_smooth_parts import smooth_sequence_batches, smoothness_exponent
        batch_exponent = smoothness_exponent(height.bit_length())
        counters.update({"scanned_norm_queries": 0, "smooth_batches": 0})

    def candidates(item, first, target, limit):
        if smooth_batch_size == 1:
            for offset in range(limit):
                yield offset, target, None
                if offset + 1 < limit:
                    target = compose_with_prime_form(target, item["norm"],
                                                     item["root"], item["form"])
                    counters["successive_power_compositions"] += 1
            return
        successor = lambda form: compose_with_prime_form(
            form, item["norm"], item["root"], item["form"])
        for offset, block, parts in smooth_sequence_batches(
                target, successor, limit, prime_product, smooth_batch_size,
                norm=lambda form: form[0], exponent=batch_exponent):
            before = counters["norm_queries"]
            counters["norm_queries"] += len(block)
            counters["smooth_batches"] += 1
            counters["successive_power_compositions"] += len(block) - int(offset == 0)
            counters["maximum_exponent_tested"] = max(
                counters["maximum_exponent_tested"], first + offset + len(block) - 1)
            if progress is not None and counters["norm_queries"] // 1000 > before // 1000:
                progress({"generator_count": len(generators), "relation_rank": len(fixed_rows),
                          "augmented_relation_rank": len(rows), "current_generator": current, **counters})
            for index, (form, smooth_part) in enumerate(zip(block, parts)):
                yield offset + index, form, smooth_part

    def extract():
        nonlocal last_extraction, last_extraction_rows
        counters["extractions"] += 1
        dense = [[row.get(j, 0) for j in range(initial_size)] for row in fixed_rows]
        _, factor, last_extraction = first_useful_relation_top(
            [records[prime][1] for prime in initial_primes], dense, identity,
            compose_forms_hnf, inverse_form, lambda form: ambiguous_factor(form, integer),
            verify_relations=not verify_candidates)
        last_extraction_rows = len(fixed_rows)
        return factor

    def finish(factor, method, terminal=None):
        assert len(fixed_rows) == len(rows) - counters["new_residual_generators"]
        result = {"integer": integer, "input_bits": integer.bit_length(), "discriminant": delta,
                  "height_bound": height, "cutoff": cutoff, "window": window,
                  "max_cofactor": max_cofactor, "max_generators": max_generators,
                  "max_queries": max_queries, "promote_residuals": promote_residuals,
                  "promotion_cutoff": promotion_cutoff,
                  "factor": factor, "method": method, "initial_base_size": initial_size,
                  "generator_count": len(generators), "relation_rank": len(fixed_rows),
                  "augmented_relation_rank": len(rows), "closed": len(rows) == len(generators),
                  "generators": generators, "initial_primes": initial_primes,
                  "saturation_prime_count": len(saturation_primes),
                  "sparse_augmented_relation_matrix": [sorted(row.items()) for row in rows],
                  "fixed_relation_sources": fixed_sources,
                  "fixed_relation_sha256": row_digest(fixed_rows), "witnesses": witnesses,
                  "last_extraction": last_extraction, "terminal_witness": terminal,
                  "factoring_oracle_calls": 0, "residual_primality_tests": 0,
                  "norm_queries_policy": "FIFO generators; first admissible power from the norm threshold",
                  **counters}
        if include_relations:
            result["sparse_fixed_relation_matrix"] = [sorted(row.items()) for row in fixed_rows]
        if smooth_batch_size > 1:
            result.update({"smooth_batch_size": smooth_batch_size,
                           "smoothness_exponent": batch_exponent,
                           "modular_smoothness_tests": counters["norm_queries"],
                           "generated_ahead_norms": counters["norm_queries"] - counters["scanned_norm_queries"],
                           "norm_queries_policy": "FIFO generators; bounded geometric batches; first admissible power from the norm threshold"})
        return result

    def stop(method, terminal=None):
        if fixed_rows and last_extraction_rows != len(fixed_rows):
            factor = extract()
            if factor is not None:
                return finish(factor, "partial prime-cofactor relation halving", terminal)
        return finish(None, method, terminal)

    for prime in base["excluded_primes"]:
        factor = gcd(prime, integer)
        if 1 < factor < integer:
            return finish(factor, "factor-base gcd", {"prime": prime})
    if max_generators is not None and initial_size > max_generators:
        return finish(None, "initial generator cap exhausted")

    current = 0
    while current < len(generators):
        item = generators[current]
        first = least_power_above(item["norm"], height)
        target = power(item["form"], first, identity, compose_forms_hnf)
        counters["initial_power_computations"] += 1
        selected = False
        limit = window
        if smooth_batch_size > 1 and max_queries is not None:
            limit = min(window, max_queries - counters["norm_queries"])
        for offset, target, batched_smooth in candidates(item, first, target, limit):
            if (smooth_batch_size == 1 and max_queries is not None
                    and counters["norm_queries"] >= max_queries):
                return stop("query cap exhausted")
            exponent = first + offset
            if smooth_batch_size == 1:
                counters["norm_queries"] += 1
            else:
                counters["scanned_norm_queries"] += 1
            if (smooth_batch_size == 1 and progress is not None
                    and counters["norm_queries"] % 1000 == 0):
                progress({"generator_count": len(generators), "relation_rank": len(fixed_rows),
                          "augmented_relation_rank": len(rows), "current_generator": current, **counters})
            counters["maximum_exponent_tested"] = max(counters["maximum_exponent_tested"], exponent)
            terminal = {"generator_index": current, "exponent": exponent, "target": target}
            if target != identity and target == inverse_form(target):
                factor = ambiguous_factor(target, integer)
                if factor is not None:
                    return finish(factor, "target involution", terminal)
            factor = gcd(target[0], integer)
            if 1 < factor < integer:
                return finish(factor, "leading-norm gcd", terminal)
            if factor != 1 or not 1 <= target[0] <= height:
                raise ArithmeticError("The reduced leading-norm invariant failed")
            smooth = (batched_smooth if smooth_batch_size > 1 else
                      gcd(target[0], pow(prime_product, target[0].bit_length(), target[0])))
            residual = target[0] // smooth
            if residual > max_cofactor:
                counters["cofactor_rejections"] += 1
            else:
                split = partial_leading_word(integer, target, saturation_primes, records,
                                             prime_product, known_smooth_part=smooth)
                counters["trial_divisions"] += split["trial_divisions"]
                if verify_candidates:
                    if not verify_partial_identity(integer, target, split, records):
                        raise ArithmeticError("A partial ideal identity failed")
                    counters["partial_identities_verified"] += 1
                row = {current: exponent}
                for prime, count in split["word"].items():
                    column = indices[prime]
                    row[column] = row.get(column, 0) - count
                residual_index, new_definition = None, False
                if residual > 1:
                    # All prime factors of R exceed the initial cutoff; the
                    # residual cap therefore proves primality without a test.
                    assert residual > cutoff
                    counters["maximum_residual_bits"] = max(
                        counters["maximum_residual_bits"], residual.bit_length())
                    if residual not in indices:
                        if max_generators is not None and len(generators) >= max_generators:
                            return stop("residual generator cap exhausted", {**terminal, **split})
                        residual_index = len(generators)
                        indices[residual] = residual_index
                        dependencies = {column: split["orientation"] * value
                                        for column, value in row.items() if value}
                        assert all(column < residual_index for column in dependencies)
                        depth = 1 + max((generators[column]["depth"] for column in dependencies), default=0)
                        generators.append({"kind": "residual", "norm": residual,
                                           "root": split["key"][1], "form": split["residual_form"],
                                           "depth": depth,
                                           "definition": {"source": current, "exponent": exponent,
                                                          "word": split["word"],
                                                          "orientation": split["orientation"],
                                                          "dependencies": sorted(dependencies.items())}})
                        records[residual] = (split["key"][1], split["residual_form"])
                        counters["new_residual_generators"] += 1
                        counters["maximum_definition_depth"] = max(counters["maximum_definition_depth"], depth)
                        if promote_residuals and residual <= promotion_cutoff:
                            saturation_primes = sorted(saturation_primes + [residual])
                            prime_product *= residual
                            counters["maximum_saturation_product_bits"] = prime_product.bit_length()
                        new_definition = True
                    else:
                        assert not (promote_residuals and residual <= promotion_cutoff)
                        residual_index = indices[residual]
                        if records[residual] != (split["key"][1], split["residual_form"]):
                            raise ArithmeticError("A prime residual changed its canonical ideal key")
                    row[residual_index] = row.get(residual_index, 0) - split["orientation"]
                row = {column: value for column, value in row.items() if value}
                assert current == len(rows)
                rows.append(row)
                witnesses.append({**terminal, **split, "residual_index": residual_index})
                counters["admissible_queries"] += 1
                if not new_definition:
                    fixed, visits, coefficient_bits = eliminate_definitions(row, generators, initial_size)
                    assert fixed and len(fixed_rows) < initial_size
                    fixed_rows.append(fixed)
                    fixed_sources.append(current)
                    counters["elimination_vertices"] += visits
                    counters["maximum_elimination_coefficient_bits"] = max(
                        counters["maximum_elimination_coefficient_bits"], coefficient_bits)
                    counters["maximum_relation_entry_bits"] = max(
                        counters["maximum_relation_entry_bits"],
                        max(abs(value).bit_length() for value in fixed.values()))
                    bits = sum((value & 1) << column for column, value in fixed.items())
                    while bits:
                        pivot = bits.bit_length() - 1
                        if pivot not in parity_pivots:
                            parity_pivots[pivot] = bits
                            break
                        bits ^= parity_pivots[pivot]
                    if not bits or len(fixed_rows) % 8 == 0 or len(fixed_rows) == initial_size:
                        factor = extract()
                        if factor is not None:
                            return finish(factor, "prime-cofactor relation halving")
                selected = True
                break
        if not selected:
            if limit < window:
                return stop("query cap exhausted")
            return stop("power window exhausted", {"generator_index": current,
                                                    "start_exponent": first, "window": window})
        current += 1
    return stop("closed prime-cofactor presentation")

