#!/usr/bin/env python3
"""Exact relation-halving controls, including the stored 60/92-bit rows."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd, lcm
from pathlib import Path
import random

from arithmetic_short_norm_factor import relation_determinant
from arithmetic_two_group import power
from binary_quadratic_forms import reduced_forms
from check_early_ambiguous_factor import has_multiple_prime_factors
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from relation_two_saturation import (
    factor_from_dominant_relations, factor_from_partial_relations,
    first_useful_relation_top, half_relation,
)
from smooth_norm_sieve import SmoothNormSieve
from smooth_norm_words import dominant_relation, exponent_vector, norm_to_word


def product_control(moduli, rows, bad_top):
    dimension = len(moduli)
    identity = (0,) * dimension
    generators = [tuple(int(i == j) % modulus for i, modulus in enumerate(moduli))
                  for j in range(dimension)]

    def multiply(left, right):
        return tuple((a + b) % modulus for a, b, modulus in zip(left, right, moduli))

    def inverse(value):
        return tuple(-value[index] % modulus for index, modulus in enumerate(moduli))

    top, certificate, stats = first_useful_relation_top(
        generators, rows, identity, multiply, inverse,
        lambda value: None if value == bad_top else "useful involution",
        keep_updates=True,
    )
    rank = sum(modulus % 2 == 0 for modulus in moduli)
    useful_exists = rank >= 2 or (rank == 1 and bad_top is None)
    assert (top is not None) == useful_exists
    if top is not None:
        assert top != identity and top != bad_top
        assert multiply(top, top) == identity and certificate is not None
    assert stats["tops_evaluated"] <= 2 * stats["kernel_passes"]
    assert stats["saturations"] <= stats["saturation_bound"]
    assert not stats["determinant_computed"] and not stats["annihilator_supplied"]

    # Independent integer determinants appear only in this checker.
    replayed = [list(row) for row in rows]
    determinant = relation_determinant(replayed)
    initial_determinant = determinant
    for update in stats["updates"]:
        vector = half_relation(replayed, update["dependency"])
        assert list(vector) == update["new_relation"]
        assert all(value % modulus == 0 for value, modulus in zip(vector, moduli))
        replayed[update["replaced_row"]] = list(vector)
        if dimension <= 3:
            next_determinant = relation_determinant(replayed)
            assert 2 * next_determinant == determinant
            determinant = next_determinant
    assert relation_determinant(replayed) * (1 << stats["saturations"]) == initial_determinant
    stats.pop("updates")
    return {
        "moduli": moduli, "two_rank": rank, "bad_top": bad_top,
        "returned_top": top, "stats": stats,
    }


def class_control(integer):
    forms = reduced_forms(-4 * integer)
    identity = (1, 0, integer)
    generators = [form for form in forms if form != identity]
    tau = (2, 2, (integer + 1) // 2) if integer % 4 == 1 else None
    if tau in generators:
        generators.remove(tau)
        generators.insert(0, tau)
    size = len(generators)
    # Exhaustive classes establish only this small input fixture. The
    # extraction routine receives generators and relations, not this order.
    rows = [[len(forms) * (1 << (index % 3)) * int(index == column)
             for column in range(size)] for index in range(size)]
    result = factor_from_dominant_relations(integer, generators, rows)
    factor = result["factor"]
    assert (factor is not None) == has_multiple_prime_factors(integer)
    if factor is not None:
        assert 1 < factor < integer and integer % factor == 0
    return {
        "integer": integer, "fixture_class_number": len(forms),
        "relation_dimension": size, **result,
    }


def partial_product_control(moduli, basis_rows, rank, shifts, bad_top):
    """The first rows of a unimodular matrix explicitly give L_sat."""
    dimension = len(moduli)
    identity = (0,) * dimension
    generators = [tuple(int(i == j) % modulus for i, modulus in enumerate(moduli))
                  for j in range(dimension)]
    annihilating_integer = lcm(*moduli) if moduli else 1
    rows = [[value * annihilating_integer * (1 << shifts[index]) for value in row]
            for index, row in enumerate(basis_rows[:rank])]

    def multiply(left, right):
        return tuple((a + b) % modulus for a, b, modulus in zip(left, right, moduli))

    def inverse(value):
        return tuple(-value[index] % modulus for index, modulus in enumerate(moduli))

    saturated_images = [tuple(value % modulus for value, modulus in zip(row, moduli))
                        for row in basis_rows[:rank]]
    subgroup, frontier = {identity}, [identity]
    while frontier:
        value = frontier.pop()
        for generator in saturated_images:
            candidate = multiply(value, generator)
            if candidate not in subgroup:
                subgroup.add(candidate)
                frontier.append(candidate)
    useful = {value for value in subgroup if value != identity and value != bad_top
              and multiply(value, value) == identity}
    top, certificate, stats = first_useful_relation_top(
        generators, rows, identity, multiply, inverse,
        lambda value: None if value == bad_top else "useful involution",
        keep_updates=True,
    )
    assert (top is not None) == bool(useful)
    if top is not None:
        assert top in useful and certificate is not None

    def maximal_minor_gcd(matrix):
        result = 0
        for columns in combinations(range(dimension), rank):
            minor = [[row[column] for column in columns] for row in matrix]
            result = gcd(result, abs(relation_determinant(minor)))
        return result

    replayed = [list(row) for row in rows]
    index = maximal_minor_gcd(replayed)
    assert index > 0
    for update in stats["updates"]:
        vector = half_relation(replayed, update["dependency"])
        assert list(vector) == update["new_relation"]
        assert all(value % modulus == 0 for value, modulus in zip(vector, moduli))
        replayed[update["replaced_row"]] = list(vector)
        new_index = maximal_minor_gcd(replayed)
        assert 2 * new_index == index
        index = new_index
    return {
        "generators": dimension, "relation_rows": rank,
        "torsion_image_size": len(subgroup), "useful_tops": len(useful),
        "saturations": stats["saturations"],
    }


def stored_control(record):
    integer, smooth_bound, budget = (record[key] for key in
                                     ("integer", "smooth_bound", "word_budget"))
    discriminant, identity = -4 * integer, (1, 0, integer)
    base = []
    for prime, root in SmoothNormSieve(discriminant, smooth_bound).prime_roots:
        root = root if root % 2 == 0 else root + prime
        root = min(root, 2 * prime - root)
        form = (prime, root, (root * root - discriminant) // (4 * prime))
        base.append((prime, reduce_with_transform(form)[0]))
    primes = [prime for prime, _ in base]
    assert primes == record["factor_base_primes"]
    assert len(record["witnesses"]) == len(base)
    rows = []
    for index, ((prime, generator), witness) in enumerate(zip(base, record["witnesses"])):
        assert witness["generator_prime"] == prime
        target = power(generator, budget + 1, identity, compose_forms_hnf)
        assert tuple(witness["target_form"]) == target
        vector = tuple(witness["primitive_vector"])
        word = norm_to_word(target, *vector, smooth_bound, budget)
        assert word is not None and word.normalized_vector == vector
        assert word.primitive_norm == witness["represented_norm"]
        assert word.multiplicity == witness["multiplicity"]
        assert [[term.prime, term.exponent] for term in word.terms] == witness["prime_word"]
        rows.append(dominant_relation(index, [0] * len(base),
                                      exponent_vector(word, primes), budget, 0))
    result = factor_from_dominant_relations(
        integer, [form for _, form in base], rows, verify_relations=True,
    )
    factor = result["factor"]
    assert factor is not None and 1 < factor < integer and integer % factor == 0
    assert result["stats"]["initial_relations_verified"] == len(base)
    prefix = {60: 89, 92: 280}[integer.bit_length()]
    partial = factor_from_partial_relations(
        integer, [form for _, form in base], rows[:prefix], verify_relations=False,
    )
    previous = factor_from_partial_relations(
        integer, [form for _, form in base], rows[:prefix - 1], verify_relations=False,
    )
    assert partial["factor"] is not None and previous["factor"] is None
    assert integer % partial["factor"] == 0
    return {
        "integer": integer, "input_bits": integer.bit_length(),
        "relation_dimension": len(base), "factor": factor,
        "cofactor": integer // factor, "ambiguous_form": result["ambiguous_form"],
        "stats": result["stats"],
        "first_successful_prefix": prefix,
        "prefix_result": partial,
        "preceding_prefix_result": previous,
        "source": "short-norm-factoring-certificates.json",
        "stored_annihilator_used": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-stored", action="store_true")
    args = parser.parse_args()
    rng = random.Random(537919)
    generic = []
    for dimension in range(13):
        for trial in range(12):
            moduli = [2 ** rng.randrange(0, 9) * rng.choice((1, 3, 5, 9, 15))
                      for _ in range(dimension)]
            rows = [[moduli[index] * 2 ** rng.randrange(0, 5) * int(index == column)
                     for column in range(dimension)] for index in range(dimension)]
            for _ in range(3 * dimension):
                if dimension >= 2:
                    left, right = rng.sample(range(dimension), 2)
                    sign = rng.choice((-1, 1))
                    rows[left] = [a + sign * b for a, b in zip(rows[left], rows[right])]
            even = [index for index, modulus in enumerate(moduli) if modulus % 2 == 0]
            bad_top = None
            if even and trial % 3:
                chosen = {rng.choice(even)} | {index for index in even if rng.randrange(2)}
                bad_top = tuple(modulus // 2 if index in chosen else 0
                                for index, modulus in enumerate(moduli))
            generic.append(product_control(moduli, rows, bad_top))
    paired = product_control([2, 2], [[2, 0], [2, 4]], (1, 0))
    assert paired["stats"]["paired_bad_top_updates"] == 1
    assert paired["returned_top"] == (1, 1)
    deep = product_control([3 * 2 ** 128], [[3 * 2 ** 256]], (3 * 2 ** 127,))
    assert deep["stats"]["saturations"] == 128 and deep["returned_top"] is None

    partial = []
    for dimension in range(5):
        for trial in range(12):
            moduli = [rng.choice((2, 3, 4, 6)) for _ in range(dimension)]
            unit_rows = [[int(index == column) for column in range(dimension)]
                         for index in range(dimension)]
            for _ in range(3 * dimension):
                if dimension >= 2:
                    left, right = rng.sample(range(dimension), 2)
                    sign = rng.choice((-1, 1))
                    unit_rows[left] = [a + sign * b for a, b in zip(unit_rows[left], unit_rows[right])]
            assert abs(relation_determinant(unit_rows)) == 1
            even = [index for index, modulus in enumerate(moduli) if modulus % 2 == 0]
            bad_top = None
            if even and trial % 2:
                selected = {rng.choice(even)}
                bad_top = tuple(modulus // 2 if index in selected else 0
                                for index, modulus in enumerate(moduli))
            for rank in range(dimension + 1):
                partial.append(partial_product_control(
                    moduli, unit_rows, rank, [rng.randrange(3) for _ in range(rank)], bad_top,
                ))

    integers = (3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 33, 35, 45, 49,
                55, 65, 75, 77, 81, 105, 121, 125, 165, 231, 385, 1155,
                1365, 10403, 10807, 16637, 15015, 19635)
    classes = [class_control(integer) for integer in integers]
    root = Path(__file__).resolve().parents[1]
    stored = []
    if not args.skip_stored:
        data = json.loads((root / "evidence" / "short-norm-factoring-certificates.json").read_text())
        stored = [stored_control(record) for record in data["records"]]
    report = {
        "status": "pass",
        "scope": "Exact determinant-free extraction from supplied relations; no new norm-search bound",
        "generic_product_controls": len(generic) + 2,
        "maximum_generic_two_rank": max(record["two_rank"] for record in generic),
        "generic_saturations": sum(record["stats"]["saturations"] for record in generic),
        "partial_presentation_controls": len(partial),
        "partial_saturation_updates": sum(record["saturations"] for record in partial),
        "partial_controls_with_useful_torsion": sum(record["useful_tops"] > 0 for record in partial),
        "paired_bad_top_control": paired, "deep_saturation_control": deep,
        "class_group_controls": len(classes),
        "successful_small_factor_targets": sum(record["factor"] is not None for record in classes),
        "class_group_records": classes,
        "stored_relation_rows_verified": sum(record["relation_dimension"] for record in stored),
        "stored_records": stored,
        "engine_sha256": sha256((root / "tools" / "relation_two_saturation.py").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == "__main__":
    main()
