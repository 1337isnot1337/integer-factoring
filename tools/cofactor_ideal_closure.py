#!/usr/bin/env python3
"""Close a dominant relation system with unfactored composite ideal generators.

One admissible leading-norm query is selected per generator. Residual
ideal keys, not new primes, extend the base. Every finite budget is
explicit; no uniform admissible-window or closure-size gain is claimed.
"""

from __future__ import annotations

from math import gcd, isqrt, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from cofactor_collision_relations import partial_leading_word, verify_partial_identity
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def factor_with_cofactor_ideal_closure(integer: int, cutoff: int, window: int, *,
                                     max_cofactor: int | None = None,
                                     max_generators: int | None = None,
                                     max_queries: int | None = None,
                                     verify_candidates: bool = False,
                                     progress=None) -> dict:
    if not isinstance(integer, int) or integer < 3 or integer % 2 == 0:
        raise ValueError("Require an odd integer of at least three")
    if not isinstance(cutoff, int) or cutoff < 2 or not isinstance(window, int) or window < 1:
        raise ValueError("Require cutoff >= 2 and window >= 1")
    for name, value in (("max_generators", max_generators), ("max_queries", max_queries)):
        if value is not None and (not isinstance(value, int) or value < 0):
            raise ValueError(f"{name} must be a nonnegative integer")
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    if max_cofactor is None:
        max_cofactor = height
    if not isinstance(max_cofactor, int) or max_cofactor < 1:
        raise ValueError("Require a positive residual cap")
    base = build_prime_form_factor_base(delta, cutoff)
    records = {item.prime: (item.root, item.positive) for item in base["records"]}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes = sorted(records)
    indices = {prime: index for index, prime in enumerate(primes)}
    prime_product = prod(primes)
    generators = [{"kind": "prime", "norm": prime, "form": records[prime][1], "key": None,
                   "definition": None, "depth": 0}
                  for prime in primes]
    residual_indices: dict[tuple[int, int], int] = {}
    rows: list[dict[int, int]] = []
    fixed_rows: list[dict[int, int]] = []
    fixed_sources = []
    word_cache = {index: {index: 1} for index in range(len(primes))}
    witnesses = []
    parity_pivots = {}
    counters = {"norm_queries": 0, "cofactor_rejections": 0, "admissible_queries": 0,
                "new_residual_generators": 0, "trial_divisions": 0,
                "initial_power_computations": 0, "successive_power_compositions": 0,
                "extractions": 0, "maximum_exponent_tested": 0,
                "maximum_residual_bits": 0, "word_expansions": 0,
                "definition_nodes_expanded": 0, "maximum_word_coefficient_bits": 1,
                "maximum_definition_depth": 0, "maximum_relation_entry_bits": 0,
                "partial_identities_verified": 0}
    last_extraction, last_extraction_rows = None, -1

    def expand_word(index):
        if index in word_cache:
            return word_cache[index]
        counters["word_expansions"] += 1
        trail, current_word = [], index
        while current_word not in word_cache:
            trail.append(current_word)
            current_word = generators[current_word]["definition"]["source"]
        word = dict(word_cache[current_word])
        for vertex in reversed(trail):
            definition = generators[vertex]["definition"]
            sign = definition["orientation"]
            multiplier = sign * definition["exponent"]
            word = {column: multiplier * value for column, value in word.items()}
            for prime, count in definition["word"].items():
                column = indices[prime]
                word[column] = word.get(column, 0) - sign * count
            word = {column: value for column, value in word.items() if value}
        counters["definition_nodes_expanded"] += len(trail)
        counters["maximum_word_coefficient_bits"] = max(
            counters["maximum_word_coefficient_bits"],
            max((abs(value).bit_length() for value in word.values()), default=0))
        word_cache[index] = word
        return word

    def extract():
        nonlocal last_extraction, last_extraction_rows
        counters["extractions"] += 1
        dense = [[row.get(j, 0) for j in range(len(primes))] for row in fixed_rows]
        _, factor, last_extraction = first_useful_relation_top(
            [records[prime][1] for prime in primes], dense, identity,
            compose_forms_hnf, inverse_form,
            lambda form: ambiguous_factor(form, integer),
            # Literal partial identities plus exact integral substitutions
            # establish every emitted row. Avoid repeating their long word
            # evaluation at each extraction after they have been established.
            verify_relations=not verify_candidates)
        last_extraction_rows = len(fixed_rows)
        return factor

    def finish(factor, method, terminal=None):
        assert len(fixed_rows) == len(rows) - len(residual_indices)
        return {"integer": integer, "input_bits": integer.bit_length(), "discriminant": delta,
                "height_bound": height, "cutoff": cutoff, "window": window,
                "max_cofactor": max_cofactor, "max_generators": max_generators,
                "max_queries": max_queries, "factor": factor, "method": method,
                "initial_base_size": len(primes), "generator_count": len(generators),
                "relation_rank": len(fixed_rows), "augmented_relation_rank": len(rows),
                "closed": len(rows) == len(generators),
                "generators": generators, "active_primes": primes, "prime_records": records,
                "sparse_augmented_relation_matrix": [sorted(row.items()) for row in rows],
                "relation_matrix": [[row.get(j, 0) for j in range(len(primes))] for row in fixed_rows],
                "fixed_relation_sources": fixed_sources, "cached_word_count": len(word_cache),
                "witnesses": witnesses, "last_extraction": last_extraction,
                "terminal_witness": terminal, "factoring_oracle_calls": 0,
                "residual_primality_tests": 0, "new_primes_introduced": 0, **counters}

    def stop(method, terminal=None):
        if fixed_rows and last_extraction_rows != len(fixed_rows):
            factor = extract()
            if factor is not None:
                return finish(factor, "partial ideal-closure relation halving", terminal)
        return finish(None, method, terminal)

    for prime in base["excluded_primes"]:
        factor = gcd(prime, integer)
        if 1 < factor < integer:
            return finish(factor, "factor-base gcd", {"prime": prime})
    if max_generators is not None and len(generators) > max_generators:
        return finish(None, "initial generator cap exhausted")

    current = 0
    while current < len(generators):
        item = generators[current]
        first = least_power_above(item["norm"], height)
        target = power(item["form"], first, identity, compose_forms_hnf)
        counters["initial_power_computations"] += 1
        selected = False
        for offset in range(window):
            if max_queries is not None and counters["norm_queries"] >= max_queries:
                return stop("query cap exhausted")
            exponent = first + offset
            counters["norm_queries"] += 1
            if progress is not None and counters["norm_queries"] % 1000 == 0:
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
            smooth = gcd(target[0], pow(prime_product, target[0].bit_length(), target[0]))
            residual = target[0] // smooth
            if residual > max_cofactor:
                counters["cofactor_rejections"] += 1
            else:
                split = partial_leading_word(integer, target, primes, records, prime_product,
                                             known_smooth_part=smooth)
                counters["trial_divisions"] += split["trial_divisions"]
                if verify_candidates:
                    if not verify_partial_identity(integer, target, split, records):
                        raise ArithmeticError("A partial ideal identity failed")
                    counters["partial_identities_verified"] += 1
                residual_index = None
                new_definition = False
                if split["key"] is not None:
                    counters["maximum_residual_bits"] = max(counters["maximum_residual_bits"], residual.bit_length())
                    if split["key"] not in residual_indices:
                        if max_generators is not None and len(generators) >= max_generators:
                            return stop("residual generator cap exhausted", {**terminal, **split})
                        residual_index = len(generators)
                        residual_indices[split["key"]] = residual_index
                        generators.append({"kind": "residual", "norm": residual,
                                           "form": split["residual_form"], "key": split["key"],
                                           "depth": item["depth"] + 1,
                                           "definition": {"source": current, "exponent": exponent,
                                                          "word": split["word"],
                                                          "orientation": split["orientation"]}})
                        counters["new_residual_generators"] += 1
                        counters["maximum_definition_depth"] = max(
                            counters["maximum_definition_depth"], item["depth"] + 1)
                        new_definition = True
                    else:
                        residual_index = residual_indices[split["key"]]
                row = {current: exponent}
                for prime, count in split["word"].items():
                    column = indices[prime]
                    row[column] = row.get(column, 0) - count
                if residual_index is not None:
                    row[residual_index] = row.get(residual_index, 0) - split["orientation"]
                row = {column: value for column, value in row.items() if value}
                assert current == len(rows)
                rows.append(row)
                witnesses.append({**terminal, **split, "residual_index": residual_index})
                counters["admissible_queries"] += 1
                if not new_definition:
                    fixed = {column: exponent * value for column, value in expand_word(current).items()}
                    for prime, count in split["word"].items():
                        column = indices[prime]
                        fixed[column] = fixed.get(column, 0) - count
                    if residual_index is not None:
                        for column, value in expand_word(residual_index).items():
                            fixed[column] = fixed.get(column, 0) - split["orientation"] * value
                    fixed = {column: value for column, value in fixed.items() if value}
                    assert fixed and len(fixed_rows) < len(primes)
                    fixed_rows.append(fixed)
                    counters["maximum_relation_entry_bits"] = max(
                        counters["maximum_relation_entry_bits"],
                        max(abs(value).bit_length() for value in fixed.values()))
                    fixed_sources.append(current)
                    bits = sum((value & 1) << column for column, value in fixed.items())
                    while bits:
                        pivot = bits.bit_length() - 1
                        if pivot not in parity_pivots:
                            parity_pivots[pivot] = bits
                            break
                        bits ^= parity_pivots[pivot]
                    if not bits or len(fixed_rows) % 8 == 0 or len(fixed_rows) == len(primes):
                        factor = extract()
                        if factor is not None:
                            return finish(factor, "ideal-closure relation halving")
                selected = True
                break
            if offset + 1 < window:
                target = compose_forms_hnf(target, item["form"])
                counters["successive_power_compositions"] += 1
        if not selected:
            return stop("admissible power window exhausted", {"generator_index": current})
        current += 1
    return stop("complete composite-ideal presentation")
