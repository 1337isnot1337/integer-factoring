#!/usr/bin/env python3
"""Bounded deterministic class-group factoring through actual norm search.

A successful factor is checked by integer division. Failure means only
that this finite search budget or generated subgroup did not split N.
The supplied practical prime cutoff is not asserted to be the uniform
analytic cutoff from the research theorem.
"""

import argparse
from dataclasses import asdict
from math import gcd
import json

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form
from early_ambiguous_factor import factor_from_annihilator
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from smooth_norm_sieve import SmoothNormSieve, primes_through
from smooth_norm_words import dominant_relation, exponent_vector, norm_to_word


def relation_determinant(matrix):
    rows = [list(row) for row in matrix]
    size, previous, sign = len(rows), 1, 1
    if size == 0:
        return 1
    for col in range(size - 1):
        pivot = next((row for row in range(col, size) if rows[row][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            rows[col], rows[pivot] = rows[pivot], rows[col]
            sign = -sign
        value = rows[col][col]
        for row in range(col + 1, size):
            for target in range(col + 1, size):
                numerator = rows[row][target] * value - rows[row][col] * rows[col][target]
                if numerator % previous:
                    raise ArithmeticError("inexact fraction-free elimination")
                rows[row][target] = numerator // previous
            rows[row][col] = 0
        previous = value
    return sign * rows[-1][-1]


def _vector_value(vector, base, identity):
    result = identity
    for exponent, (_, value) in zip(vector, base):
        if not exponent:
            continue
        if exponent < 0:
            value = inverse_form(value)
        result = compose_forms_hnf(
            result, power(value, abs(exponent), identity, compose_forms_hnf)
        )
    return result


def factor_with_short_norms(integer, smooth_bound=127, word_budget=None,
                           maximum_expansions=8, keep_relations=False):
    if integer < 3 or integer % 2 == 0:
        raise ValueError("this research prototype expects an odd integer >= 3")
    if smooth_bound < 3 or maximum_expansions < 0:
        raise ValueError("invalid finite search budget")
    if word_budget is None:
        word_budget = integer.bit_length()
    if word_budget < 0:
        raise ValueError("word budget must be nonnegative")
    record = {
        "integer": integer,
        "input_bits": integer.bit_length(),
        "discriminant": -4 * integer,
        "smooth_bound": smooth_bound,
        "word_budget": word_budget,
        "diagonal_parameter": word_budget + 1,
        "maximum_expansions_per_target": maximum_expansions,
        "complete_class_group_enumerated": False,
        "generator_orders_searched": False,
        "finite_cutoff_certified_to_generate_full_group": False,
    }
    for prime in primes_through(smooth_bound):
        divisor = gcd(prime, integer)
        if 1 < divisor < integer:
            record.update(status="factor", factor=divisor, method="factor-base gcd")
            return record

    discriminant = -4 * integer
    identity = (1, 0, integer)
    sieve = SmoothNormSieve(discriminant, smooth_bound)
    base = []
    for prime, root in sieve.prime_roots:
        root = root if root % 2 == 0 else root + prime
        root = min(root, 2 * prime - root)
        form = (prime, root, (root * root - discriminant) // (4 * prime))
        base.append((prime, reduce_with_transform(form)[0]))
    primes = [prime for prime, _ in base]
    record["factor_base_primes"] = primes
    if not base:
        record.update(status="unresolved", reason="empty finite factor base")
        return record

    rows, witnesses = [], []
    total_points = total_strikes = total_divisions = 0
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
        else:
            record.update(
                status="unresolved", reason="norm-search budget exhausted",
                failed_generator_index=index, failed_generator_prime=prime,
                failed_target_form=target, failed_attempts=attempts,
                relation_rows_completed=len(rows),
                coordinate_points_visited=total_points,
                modular_strikes=total_strikes, prime_divisions=total_divisions,
            )
            if keep_relations:
                record["witnesses"] = witnesses
                record["relation_matrix"] = rows
            return record
        word = norm_to_word(target, *found.vector, smooth_bound, word_budget)
        if word is None:
            raise ArithmeticError("sieve witness failed the independent direct converter")
        encoded = exponent_vector(word, primes)
        row = dominant_relation(index, [0] * len(base), encoded, word_budget, 0)
        if _vector_value(row, base, identity) != identity:
            raise ArithmeticError("constructed row is not a class-group relation")
        rows.append(row)
        witnesses.append({
            "generator_prime": prime,
            "target_form": target,
            "represented_norm": found.norm,
            "primitive_vector": found.vector,
            "multiplicity": found.multiplicity,
            "prime_word": [(term.prime, term.exponent) for term in word.terms],
            "search_attempts": attempts,
        })

    annihilator = abs(relation_determinant(rows))
    assert 0 < annihilator <= (2 * word_budget + 1) ** len(base)
    if not all(power(form, annihilator, identity, compose_forms_hnf) == identity for _, form in base):
        raise ArithmeticError("relation determinant failed to annihilate a generator")
    extraction = factor_from_annihilator(
        integer, [form for _, form in base], annihilator,
    )
    factor = extraction["factor"]
    if factor is not None:
        assert 1 < factor < integer and integer % factor == 0
    record.update(
        status="factor" if factor is not None else "unresolved",
        method="searched short norms, relation determinant, early ambiguous extraction",
        relation_rows_completed=len(rows),
        relation_determinant_bits=annihilator.bit_length(),
        ambiguous_extraction=extraction,
        coordinate_points_visited=total_points,
        modular_strikes=total_strikes,
        prime_divisions=total_divisions,
        maximum_witness_norm=max(witness["represented_norm"] for witness in witnesses),
    )
    if factor is not None:
        record["factor"] = factor
    else:
        record["reason"] = "no proper factor in the computed subgroup two-part"
    if keep_relations:
        record["witnesses"] = witnesses
        record["relation_matrix"] = rows
        record["relation_determinant"] = annihilator
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("integer", type=int)
    parser.add_argument("--smooth-bound", type=int, default=127)
    parser.add_argument("--word-budget", type=int)
    parser.add_argument("--maximum-expansions", type=int, default=8)
    parser.add_argument("--keep-relations", action="store_true")
    arguments = parser.parse_args()
    print(json.dumps(factor_with_short_norms(
        arguments.integer, arguments.smooth_bound, arguments.word_budget,
        arguments.maximum_expansions, arguments.keep_relations,
    ), indent=2))


if __name__ == "__main__":
    main()
