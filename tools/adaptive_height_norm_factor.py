#!/usr/bin/env python3
"""Finite height-bounded norm selection across several target powers.

For each generator, the search accepts any primitive smooth norm n<p**k,
where 1<=k<=the least exponent with p**k>H. It tries leading norms first
and then interleaves complete ellipse caps. The counters include repeated
ellipse visits and the separate leading-coefficient tests.
No uniform fast-selector bound is claimed.
"""

import argparse
from copy import copy
import json
from math import gcd
import sys

from arithmetic_short_norm_factor import _vector_value
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from smooth_norm_sieve import SmoothNormSieve, primes_through
from smooth_norm_words import exponent_vector, norm_to_word
from weighted_norm_relations import least_power_above, factor_from_weighted_relations


def factor_with_adaptive_height_norms(integer, cutoff, height_bound, *,
                                      extraction_interval=0,
                                      keep_certificates=False, progress=None):
    """Run the complete finite power/cap schedule and try weighted extraction.

    Zero extraction_interval means extract once after the supplied finite
    search. A positive interval enables early partial extraction. Missing
    target witnesses are recorded and skipped. Full class-group generation
    is not certified at the chosen finite cutoff.
    """
    if not isinstance(integer, int) or integer < 3 or integer % 2 == 0:
        raise ValueError("An odd integer at least three is required.")
    if not isinstance(cutoff, int) or cutoff < 3:
        raise ValueError("The prime cutoff must be at least three.")
    if not isinstance(height_bound, int) or height_bound < 1:
        raise ValueError("The norm height must be positive.")
    if not isinstance(extraction_interval, int) or extraction_interval < 0:
        raise ValueError("The extraction interval must be nonnegative.")
    stats = {
        "power_group_multiplications": 0,
        "eligible_target_powers": 0,
        "leading_norm_tests": 0,
        "leading_division_tests": 0,
        "leading_prime_divisions": 0,
        "leading_norm_successes": 0,
        "ellipse_attempts": 0,
        "ellipse_rows_visited": 0,
        "ellipse_coordinate_points": 0,
        "ellipse_primitive_coprime_candidates": 0,
        "ellipse_modular_strikes": 0,
        "ellipse_prime_divisions": 0,
        "target_root_schedules_built": 0,
        "target_prime_root_schedules": 0,
        "relations_verified": 0,
        "extraction_calls": 0,
        "complete_class_group_enumerated": False,
        "finite_cutoff_certified_to_generate_full_group": False,
        "determinant_computed": False,
    }
    base, rows, pivots, witnesses, failed, extractions = [], [], [], [], [], []

    def finish(factor, method):
        if factor is not None:
            assert 1 < factor < integer and integer % factor == 0
        result = {
            "integer": integer, "input_bits": integer.bit_length(),
            "prime_cutoff": cutoff, "height_bound": height_bound,
            "extraction_interval": extraction_interval,
            "status": "factor" if factor is not None else "unresolved",
            "factor": factor, "method": method,
            "factor_base_primes": [prime for prime, _ in base],
            "relation_pivot_columns": list(pivots),
            "accepted_relations": len(rows), "failed_targets": failed,
            "extractions": extractions, "stats": stats,
        }
        if keep_certificates:
            result["witnesses"] = witnesses
        return result

    for prime in primes_through(cutoff):
        if 1 < gcd(prime, integer) < integer:
            return finish(prime, "factor-base gcd")
    discriminant, identity = -4 * integer, (1, 0, integer)
    sieve = SmoothNormSieve(discriminant, cutoff)
    for prime, root in sieve.prime_roots:
        root = root if root % 2 == 0 else root + prime
        root = min(root, 2 * prime - root)
        base.append((prime, reduce_with_transform(
            (prime, root, (root * root - discriminant) // (4 * prime)))[0]))
    primes = [prime for prime, _ in base]
    generators = [form for _, form in base]

    def extract():
        result = factor_from_weighted_relations(
            integer, generators, rows, primes, pivots,
            height_bound=height_bound, verify_relations=False,
        )
        stats["extraction_calls"] += 1
        extractions.append({"rows": len(rows), "result": result})
        return result["factor"]

    for index, (prime, generator) in enumerate(base):
        maximum_exponent = least_power_above(prime, height_bound)
        value, integer_power, powers = identity, 1, []
        for exponent in range(1, maximum_exponent + 1):
            value = compose_forms_hnf(value, generator)
            integer_power *= prime
            stats["power_group_multiplications"] += 1
            allowed = min(height_bound, integer_power - 1)
            # A reduced positive form has no nonzero norm below its A.
            if value[0] > allowed:
                continue
            divisor = gcd(value[0], integer)
            if 1 < divisor < integer:
                return finish(divisor, "target leading-norm gcd")
            powers.append({
                "exponent": exponent, "form": value, "allowed": allowed,
                "cap": min(allowed, 4 * value[0]),
                # Prime roots are immutable and shared; the form-specific
                # modular slopes are cached separately for each power.
                "sieve": copy(sieve),
            })
        stats["eligible_target_powers"] += len(powers)
        found = None
        for candidate in reversed(powers):
            norm = candidate["form"][0]
            if gcd(norm, discriminant) != 1:
                continue
            stats["leading_norm_tests"] += 1
            remaining = norm
            for divisor in primes:
                if remaining == 1:
                    break
                while True:
                    stats["leading_division_tests"] += 1
                    quotient, remainder = divmod(remaining, divisor)
                    if remainder:
                        break
                    remaining = quotient
                    stats["leading_prime_divisions"] += 1
            if remaining == 1:
                found = candidate, (1, 0), norm, "leading norm"
                stats["leading_norm_successes"] += 1
                break

        active = list(reversed(powers))
        while found is None and active:
            next_active = []
            for candidate in active:
                if candidate["sieve"]._scheduled_form != tuple(candidate["form"]):
                    stats["target_root_schedules_built"] += 1
                    stats["target_prime_root_schedules"] += len(sieve.prime_roots)
                searched = candidate["sieve"].search(
                    candidate["form"], candidate["cap"], height_bound.bit_length())
                stats["ellipse_attempts"] += 1
                stats["ellipse_rows_visited"] += searched.rows_visited
                stats["ellipse_coordinate_points"] += searched.coordinate_points_visited
                stats["ellipse_primitive_coprime_candidates"] += searched.primitive_coprime_candidates
                stats["ellipse_modular_strikes"] += searched.modular_strikes
                stats["ellipse_prime_divisions"] += searched.prime_divisions
                if searched.vector is not None:
                    found = candidate, searched.vector, searched.norm, "ellipse"
                    break
                if candidate["cap"] < candidate["allowed"]:
                    candidate["cap"] = min(candidate["allowed"], 4 * candidate["cap"])
                    next_active.append(candidate)
            active = next_active

        if found is None:
            failed.append({"generator_index": index, "prime": prime,
                           "maximum_exponent": maximum_exponent,
                           "eligible_target_powers": len(powers)})
        else:
            candidate, vector, norm, source = found
            exponent, target = candidate["exponent"], candidate["form"]
            assert norm < prime**exponent and norm <= height_bound
            word = norm_to_word(target, *vector, cutoff)
            if word is None or word.primitive_norm != norm:
                raise ArithmeticError("A searched norm failed direct word conversion.")
            encoded = exponent_vector(word, primes)
            row = tuple(exponent * (column == index) - entry
                        for column, entry in enumerate(encoded))
            if _vector_value(row, base, identity) != identity:
                raise ArithmeticError("An adaptive-power row failed its class identity.")
            stats["relations_verified"] += 1
            rows.append(row)
            pivots.append(index)
            witnesses.append({
                "generator_index": index, "generator_prime": prime,
                "target_exponent": exponent, "maximum_exponent": maximum_exponent,
                "target_form": target, "norm": norm, "primitive_vector": vector,
                "signed_prime_word": [(term.prime, term.exponent) for term in word.terms],
                "source": source,
            })
            if extraction_interval and len(rows) % extraction_interval == 0:
                factor = extract()
                if factor is not None:
                    return finish(factor, "adaptive powers and partial weighted relation halving")
        if progress is not None:
            progress({"generator_index": index, "generators": len(base),
                      "accepted_relations": len(rows), "failed_targets": len(failed),
                      "ellipse_coordinate_points": stats["ellipse_coordinate_points"]})
    if not extractions or extractions[-1]["rows"] != len(rows):
        factor = extract()
    else:
        factor = extractions[-1]["result"]["factor"]
    return finish(factor, "adaptive powers and weighted relation halving" if factor is not None
                  else "no useful top in the accepted adaptive relation torsion image")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("integer", type=int)
    parser.add_argument("--cutoff", type=int, default=127)
    parser.add_argument("--height-bits", type=int, required=True)
    parser.add_argument("--extraction-interval", type=int, default=0)
    parser.add_argument("--keep-certificates", action="store_true")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    if args.height_bits < 0:
        raise ValueError("Height bits must be nonnegative.")
    progress = (lambda value: print(json.dumps(value), file=sys.stderr, flush=True)) if args.progress else None
    print(json.dumps(factor_with_adaptive_height_norms(
        args.integer, args.cutoff, 1 << args.height_bits,
        extraction_interval=args.extraction_interval,
        keep_certificates=args.keep_certificates, progress=progress), indent=2))


if __name__ == "__main__":
    main()
