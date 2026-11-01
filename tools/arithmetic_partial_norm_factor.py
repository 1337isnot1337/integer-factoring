#!/usr/bin/env python3
"""Bounded norm search with partial relation-halving extraction.

Every successful row is checked as a class-group relation. Failed finite
norm targets are recorded and skipped. Periodic extraction uses only the
independent rows collected so far and never constructs a determinant.
This remains a finite research prototype, not an asymptotic search bound.
"""

import argparse
from dataclasses import asdict
import json
from math import gcd
import sys

from arithmetic_short_norm_factor import _vector_value
from arithmetic_two_group import power
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from relation_two_saturation import factor_from_partial_relations
from smooth_norm_sieve import SmoothNormSieve, primes_through
from smooth_norm_words import dominant_relation, exponent_vector, norm_to_word


def factor_with_partial_short_norms(integer, smooth_bound=127, word_budget=None,
                                    maximum_expansions=8, extraction_interval=16,
                                    keep_relations=False, progress=None):
    """Run a fixed finite sieve budget and inspect independent partial rows."""
    if integer < 3 or integer % 2 == 0:
        raise ValueError("an odd integer at least three is required")
    if smooth_bound < 3 or maximum_expansions < 0 or extraction_interval < 1:
        raise ValueError("invalid finite search or extraction budget")
    if word_budget is None:
        word_budget = integer.bit_length()
    if word_budget < 0:
        raise ValueError("the multiplicity budget must be nonnegative")
    record = {
        "integer": integer, "input_bits": integer.bit_length(),
        "discriminant": -4 * integer, "smooth_bound": smooth_bound,
        "word_budget": word_budget, "diagonal_parameter": word_budget + 1,
        "maximum_expansions_per_target": maximum_expansions,
        "extraction_interval": extraction_interval,
        "complete_class_group_enumerated": False,
        "generator_orders_searched": False,
        "finite_cutoff_certified_to_generate_full_group": False,
        "determinant_computed": False,
    }
    for prime in primes_through(smooth_bound):
        divisor = gcd(prime, integer)
        if 1 < divisor < integer:
            record.update(status="factor", factor=divisor, method="factor-base gcd")
            return record
    discriminant, identity = -4 * integer, (1, 0, integer)
    sieve = SmoothNormSieve(discriminant, smooth_bound)
    base = []
    for prime, root in sieve.prime_roots:
        root = root if root % 2 == 0 else root + prime
        root = min(root, 2 * prime - root)
        form = (prime, root, (root * root - discriminant) // (4 * prime))
        base.append((prime, reduce_with_transform(form)[0]))
    primes = [prime for prime, _ in base]
    generators = [form for _, form in base]
    record["factor_base_primes"] = primes
    if not base:
        record.update(status="unresolved", reason="empty finite factor base")
        return record

    rows, pivots, witnesses, failed, extractions = [], [], [], [], []
    total_points = total_strikes = total_divisions = 0
    extraction = None

    def try_extract():
        value = factor_from_partial_relations(
            integer, generators, rows, pivots, verify_relations=False,
        )
        extractions.append({"rows": len(rows), "result": value})
        return value

    for index, (prime, generator) in enumerate(base):
        target = power(generator, word_budget + 1, identity, compose_forms_hnf)
        maximum_norm = 4 * max(target[0], target[2])
        attempts = []
        for expansion in range(maximum_expansions + 1):
            found = sieve.search(target, maximum_norm, word_budget)
            total_points += found.coordinate_points_visited
            total_strikes += found.modular_strikes
            total_divisions += found.prime_divisions
            attempts.append(asdict(found))
            if found.vector is not None:
                break
            maximum_norm *= 4
        if found.vector is None:
            failed.append({
                "generator_index": index, "generator_prime": prime,
                "target_form": target, "search_attempts": attempts,
            })
        else:
            word = norm_to_word(target, *found.vector, smooth_bound, word_budget)
            if word is None:
                raise ArithmeticError("the exact norm converter rejected a sieve witness")
            encoded = exponent_vector(word, primes)
            row = dominant_relation(index, [0] * len(base), encoded, word_budget, 0)
            if _vector_value(row, base, identity) != identity:
                raise ArithmeticError("a constructed row is not a class-group relation")
            rows.append(row)
            pivots.append(index)
            witnesses.append({
                "generator_index": index, "generator_prime": prime,
                "target_form": target, "represented_norm": found.norm,
                "primitive_vector": found.vector, "multiplicity": found.multiplicity,
                "prime_word": [(term.prime, term.exponent) for term in word.terms],
                "search_attempts": attempts,
            })
            if len(rows) % extraction_interval == 0:
                extraction = try_extract()
        if progress is not None:
            progress({
                "generator_index": index, "generators": len(base),
                "relation_rows": len(rows), "failed_targets": len(failed),
                "coordinate_points": total_points,
                "factor": None if extraction is None else extraction["factor"],
            })
        if extraction is not None and extraction["factor"] is not None:
            break
    if not extractions or extractions[-1]["rows"] != len(rows):
        extraction = try_extract()
    factor = extraction["factor"]
    if factor is not None:
        assert 1 < factor < integer and integer % factor == 0
    record.update(
        status="factor" if factor is not None else "unresolved",
        method="bounded norm search and determinant-free partial relation extraction",
        relation_rows_completed=len(rows), relation_pivot_columns=pivots,
        failed_target_count=len(failed), failed_targets=failed,
        targets_visited=index + 1, coordinate_points_visited=total_points,
        modular_strikes=total_strikes, prime_divisions=total_divisions,
        maximum_witness_norm=max((witness["represented_norm"] for witness in witnesses), default=None),
        partial_extractions=extractions, relation_extraction=extraction,
    )
    if factor is not None:
        record["factor"] = factor
    else:
        record["reason"] = "no useful top in the torsion image of the collected finite relations"
    if keep_relations:
        record["witnesses"] = witnesses
        record["relation_matrix"] = rows
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("integer", type=int)
    parser.add_argument("--smooth-bound", type=int, default=127)
    parser.add_argument("--word-budget", type=int)
    parser.add_argument("--maximum-expansions", type=int, default=8)
    parser.add_argument("--extraction-interval", type=int, default=16)
    parser.add_argument("--keep-relations", action="store_true")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    progress = (lambda value: print(json.dumps(value), file=sys.stderr, flush=True)) if args.progress else None
    print(json.dumps(factor_with_partial_short_norms(
        args.integer, args.smooth_bound, args.word_budget,
        args.maximum_expansions, args.extraction_interval,
        args.keep_relations, progress,
    ), indent=2))


if __name__ == "__main__":
    main()
