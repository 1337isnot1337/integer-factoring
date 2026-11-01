#!/usr/bin/env python3
"""Replay the 60- and 92-bit norm-witness certificates without norm search.

This verifies every supplied representation, every relation row, the full
determinant, and the final proper factor. Exhausted search statistics are
reported by their original bounded runs and are not re-established here.
"""

import argparse
import json
from pathlib import Path

from arithmetic_short_norm_factor import relation_determinant
from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form
from early_ambiguous_factor import factor_from_annihilator
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from smooth_norm_sieve import SmoothNormSieve
from smooth_norm_words import dominant_relation, exponent_vector, norm_to_word


def replay(record):
    integer = record["integer"]
    discriminant = -4 * integer
    identity = (1, 0, integer)
    smooth_bound = record["smooth_bound"]
    budget = record["word_budget"]
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
        represented_class = identity
        for term in word.terms:
            value = term.form if term.exponent >= 0 else inverse_form(term.form)
            represented_class = compose_forms_hnf(
                represented_class,
                power(value, abs(term.exponent), identity, compose_forms_hnf),
            )
        assert represented_class == target
        encoded = exponent_vector(word, primes)
        rows.append(dominant_relation(index, [0] * len(base), encoded, budget, 0))
    annihilator = abs(relation_determinant(rows))
    assert annihilator == record["relation_determinant"]
    assert annihilator.bit_length() == record["relation_determinant_bits"]
    assert 0 < annihilator <= (2 * budget + 1) ** len(base)
    extraction = factor_from_annihilator(integer, [form for _, form in base], annihilator)
    assert extraction["factor"] == record["factor"]
    assert integer % extraction["factor"] == 0
    assert 1 < extraction["factor"] < integer
    assert json.loads(json.dumps(extraction)) == record["ambiguous_extraction"]
    assert sum(witness["coordinate_points_visited"] for witness in record["witnesses"]) == record["coordinate_points_visited"]
    return {
        "integer": integer,
        "input_bits": integer.bit_length(),
        "verified_rows": len(rows),
        "annihilator_bits": annihilator.bit_length(),
        "returned_factor": extraction["factor"],
        "cofactor": integer // extraction["factor"],
        "original_search_points": record["coordinate_points_visited"],
        "early_extraction_stats": extraction["stats"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", nargs="?", type=Path,
                        default=Path(__file__).resolve().parents[1] / "evidence" / "short-norm-factoring-certificates.json")
    arguments = parser.parse_args()
    data = json.loads(arguments.evidence.read_text())
    records = [replay(record) for record in data["records"]]
    print(json.dumps({
        "scope": "Exact certificate replay; no coordinate search or full class-group enumeration",
        "verified_rows": sum(record["verified_rows"] for record in records),
        "records": records,
    }, indent=2))


if __name__ == "__main__":
    main()
