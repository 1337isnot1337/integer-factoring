#!/usr/bin/env python3
"""Exact controls for height-only, logarithm-weighted norm relations."""

import argparse
import json
from pathlib import Path
from random import Random

from arithmetic_short_norm_factor import relation_determinant
from arithmetic_two_group import power
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from relation_two_saturation import factor_from_partial_relations
from smooth_norm_sieve import SmoothNormSieve
from smooth_norm_words import exponent_vector, norm_to_word
from weighted_norm_relations import (
    check_weighted_dominance, factor_from_weighted_relations, least_power_above,
)


def searched_record(integer, cutoff):
    height = 2 * integer
    sieve = SmoothNormSieve(-4 * integer, cutoff)
    primes, generators = [], []
    for prime, root in sieve.prime_roots:
        root = root if root % 2 == 0 else root + prime
        root = min(root, 2 * prime - root)
        primes.append(prime)
        generators.append(reduce_with_transform(
            (prime, root, (root * root + 4 * integer) // (4 * prime)))[0])
    rows, witnesses, points, ordinary_failures = [], [], 0, 0
    for index, generator in enumerate(generators):
        exponent = least_power_above(primes[index], height)
        target = power(generator, exponent, (1, 0, integer), compose_forms_hnf)
        cap = min(height, 4 * max(target[0], target[2]))
        while True:
            result = sieve.search(target, cap, height.bit_length())
            points += result.coordinate_points_visited
            if result.vector is not None or cap == height:
                break
            cap = min(height, 4 * cap)
        assert result.vector is not None
        word = norm_to_word(target, *result.vector, cutoff)
        assert word is not None and word.primitive_norm <= height
        vector = exponent_vector(word, primes)
        row = tuple(exponent * (column == index) - value
                    for column, value in enumerate(vector))
        rows.append(row)
        ordinary = abs(row[index]) > sum(abs(value) for column, value in enumerate(row) if column != index)
        ordinary_failures += not ordinary
        witnesses.append({
            "prime": primes[index], "target_exponent": exponent,
            "target_form": target, "norm": word.primitive_norm,
            "primitive_vector": word.normalized_vector,
            "multiplicity": word.multiplicity,
            "signed_prime_exponents": vector,
            "ordinary_diagonal_dominance": ordinary,
        })
    extraction = factor_from_weighted_relations(
        integer, generators, rows, primes, height_bound=height,
        verify_relations=True,
    )
    assert extraction["factor"] is not None
    assert 1 < extraction["factor"] < integer and integer % extraction["factor"] == 0
    if ordinary_failures:
        try:
            factor_from_partial_relations(integer, generators, rows, verify_relations=False)
        except ValueError as error:
            assert "dominance" in str(error)
        else:
            raise AssertionError("The ordinary wrapper accepted a failing row.")
    # Nonconsecutive subsets retain their correctly selected dominant minor.
    partial_controls = []
    for offset in range(3):
        indices = list(range(offset, len(rows), 3))
        partial = factor_from_weighted_relations(
            integer, generators, [rows[index] for index in indices], primes,
            indices, height_bound=height, verify_relations=True,
        )
        if partial["factor"] is not None:
            assert integer % partial["factor"] == 0
        partial_controls.append({"pivots": indices, "factor": partial["factor"]})
    return {
        "integer": integer, "input_bits": integer.bit_length(),
        "prime_cutoff": cutoff, "height_bound": height,
        "factor_base_primes": primes, "relation_rows": rows,
        "witnesses": witnesses, "extraction": extraction,
        "coordinate_points_visited": points,
        "rows_failing_ordinary_dominance": ordinary_failures,
        "partial_controls": partial_controls,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    power_controls = 0
    for base in range(2, 130):
        for height in list(range(1, 258)) + [2**k - 1 for k in (32, 127, 512)]:
            exponent = least_power_above(base, height)
            assert base**(exponent - 1) <= height < base**exponent
            power_controls += 1

    random = Random(43117)
    matrix_controls = 0
    for size in range(1, 8):
        bases = [2, 3, 5, 7, 11, 13, 17][:size]
        for _ in range(20):
            rows = []
            for index in range(size):
                row = [random.randint(-3, 3) for _ in range(size)]
                row[index] = 0
                right = 1
                for base, entry in zip(bases, row):
                    right *= base**abs(entry)
                row[index] = least_power_above(bases[index], right)
                if random.getrandbits(1):
                    row[index] *= -1
                rows.append(row)
            check_weighted_dominance(rows, bases, tuple(range(size)), 2**256)
            assert relation_determinant(rows) != 0
            matrix_controls += 1

    rejection_controls = 0
    for rows, bases, pivots, height in (
        ([[1 << 10000]], [2], [0], 16),
        ([[1, -1], [-1, 1]], [2, 2], [0, 1], 16),
        ([[0, 1], [1, 0]], [2, 3], [0, 1], 16),
        ([[2, 0], [0, 2]], [2, 3], [0, 0], 16),
    ):
        try:
            check_weighted_dominance(rows, bases, pivots, height)
        except ValueError:
            rejection_controls += 1
        else:
            raise AssertionError("An invalid weighted minor was accepted.")

    records = [searched_record(integer, cutoff) for integer, cutoff in (
        (10807, 31), (13493, 13), (1022117, 127),
        (100160063, 127), (10002200057, 127),
        (1000000016000000063, 1021),
    )]
    assert sum(record["rows_failing_ordinary_dominance"] for record in records) > 0
    output = {
        "scope": "Bounded-height norm relations with exact logarithmic column weights; no fast norm-selector claim",
        "least_power_controls": power_controls,
        "independent_nonsingular_matrix_controls": matrix_controls,
        "guard_rejections": rejection_controls,
        "verified_class_group_relation_rows": sum(len(record["relation_rows"]) for record in records),
        "rows_failing_ordinary_dominance": sum(record["rows_failing_ordinary_dominance"] for record in records),
        "records": records,
    }
    encoded = json.dumps(output, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
