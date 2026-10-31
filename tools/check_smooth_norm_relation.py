#!/usr/bin/env python3
"""Exact smooth-norm word and relation controls; witness search is exhaustive."""

from collections import deque
from functools import lru_cache
from math import gcd, isqrt, lcm
import json

from arithmetic_two_group import power, project_to_two_part, two_group_basis
from binary_quadratic_forms import (
    compose_forms, inverse_form, reduced_forms, reduce_form, smooth_prime_word,
)
from bounded_smooth_norm import smooth_factors
from check_positive_relation_kernel import determinant
from smooth_norm_relation import dominant_relation_row, norm_to_prime_word


def factor_base(discriminant, bound):
    output = []
    for prime in range(2, bound+1):
        if any(prime % divisor == 0 for divisor in range(2, isqrt(prime)+1)):
            continue
        if gcd(prime, discriminant) != 1:
            continue
        roots = [value for value in range(prime+1)
                 if (value*value-discriminant) % (4*prime) == 0]
        if roots:
            assert len(roots) == 1
            root = roots[0]
            output.append((prime, reduce_form((prime, root,
                                               (root*root-discriminant)//(4*prime)))))
    return output


def word_value(generators, vector, identity, multiply):
    value = identity
    for generator, exponent in zip(generators, vector):
        if exponent:
            base = generator if exponent > 0 else inverse_form(generator)
            value = multiply(value, power(base, abs(exponent), identity, multiply))
    return value


def exact_balance(integer, radius):
    return any(integer % divisor == 0 and integer//divisor <= radius
               for divisor in range(1, min(radius, isqrt(integer))+1))


def enumerate_witnesses(discriminant, radius, bound, forms, base, identity, multiply):
    primes = [prime for prime, _ in base]
    generators = [value for _, value in base]
    lookup = {}
    representations = negative_exponents = unbalanced = scaled = 0
    for form in forms:
        a, b, c = form
        y_bound = isqrt((4*a*radius*radius)//(-discriminant))
        for y in range(-y_bound, y_bound+1):
            root = isqrt(4*a*radius*radius+discriminant*y*y)
            for x in range(-((root+b*y)//(2*a)), (root-b*y)//(2*a)+1):
                if gcd(x, y) != 1:
                    continue
                norm = a*x*x+b*x*y+c*y*y
                if norm <= 0 or gcd(norm, discriminant) != 1 or smooth_factors(norm, bound) is None:
                    continue
                result = norm_to_prime_word(form, x, y, bound, radius*radius)
                assert result is not None and result["norm"] == norm
                vector = [0]*len(base)
                for entry in result["word"]:
                    position = primes.index(entry["prime"])
                    assert tuple(entry["form"]) == generators[position]
                    vector[position] += entry["exponent"]
                assert word_value(generators, vector, identity, multiply) == form
                assert sum(abs(value) for value in vector) <= (radius*radius).bit_length()-1
                assert norm_to_prime_word(form, 5*x, 5*y, bound, radius*radius) == result
                scaled += 1
                balanced = exact_balance(norm, radius)
                negative_exponents += any(value < 0 for value in vector)
                unbalanced += not balanced
                representations += 1
                candidate = {"form": form, "vector": (x, y), "word": result,
                             "balanced_at_radius": balanced}
                # Prefer an unbalanced witness to make the weaker contract visible.
                if (form not in lookup
                        or (balanced, norm) < (lookup[form]["balanced_at_radius"],
                                              lookup[form]["word"]["norm"])):
                    lookup[form] = candidate
    return lookup, {"representations_checked": representations,
                    "negative_orientation_controls": negative_exponents,
                    "unbalanced_norms_accepted": unbalanced,
                    "nonprimitive_scaled_controls": scaled}


def check(discriminant, radius=None, bound=7):
    if radius is None:
        radius = isqrt(-discriminant-1)//2
    assert 4*radius*radius < -discriminant
    forms = reduced_forms(discriminant)
    identity = next(form for form in forms if form[0] == 1)
    multiply = lru_cache(None)(compose_forms)
    base = factor_base(discriminant, bound)
    primes = [prime for prime, _ in base]
    generators = [value for _, value in base]
    words = {identity: [0]*len(base)}
    queue = deque([identity])
    while queue:
        value = queue.popleft()
        for i, generator in enumerate(generators):
            next_value = multiply(value, generator)
            if next_value not in words:
                vector = words[value][:]
                vector[i] += 1
                words[next_value] = vector
                queue.append(next_value)
    subgroup = set(words)
    for value, vector in words.items():
        assert word_value(generators, vector, identity, multiply) == value
    lookup, counts = enumerate_witnesses(discriminant, radius, bound, forms,
                                        base, identity, multiply)
    assert set(lookup) <= subgroup
    selected = {form for form in forms
                if form[0] <= radius and smooth_prime_word(form, bound) is not None}
    quotients = {multiply(left, inverse_form(right)) for left in selected for right in selected}
    assert quotients <= set(lookup)
    sample_bound = max(sum(abs(value) for value in word) for word in words.values())
    norm_word_bound = (radius*radius).bit_length()-1
    diagonal = sample_bound+norm_word_bound+1
    rows = []
    witnesses = []
    preferred = min(lookup.values(), key=lambda item: (
        item["balanced_at_radius"], item["word"]["norm"], item["form"]))
    for i, generator in enumerate(generators):
        target = power(generator, diagonal, identity, multiply)
        sample = multiply(inverse_form(target), preferred["form"])
        assert sample in words
        assert multiply(target, sample) == preferred["form"]
        row = dominant_relation_row(i, primes, words[sample], preferred["word"],
                                    sample_bound, norm_word_bound)
        assert word_value(generators, row, identity, multiply) == identity
        rows.append(row)
        witnesses.append({"generator_prime": primes[i], "sample": list(sample),
                          "target_after_sample": list(preferred["form"]),
                          "norm": preferred["word"]["norm"],
                          "balanced_at_radius": preferred["balanced_at_radius"]})
    relation_det = abs(determinant(rows))
    assert relation_det > 0 and relation_det % len(subgroup) == 0
    assert all(power(value, relation_det, identity, multiply) == identity for value in subgroup)
    index_bound = (len(forms)+len(selected)-1)//len(selected)
    correction = lcm(*range(1, index_bound+1))
    global_multiple = relation_det*correction
    assert global_multiple % len(forms) == 0
    projected = project_to_two_part(generators, relation_det, identity, multiply)
    basis = two_group_basis(projected, identity, multiply, inverse_form, len(forms).bit_length())
    tops = [power(value, 1 << (order-1), identity, multiply) for value, order in basis]
    span = {identity}
    for top in tops:
        span |= {multiply(value, top) for value in span}
    ambiguous = {value for value in subgroup if multiply(value, value) == identity}
    assert span == ambiguous
    return {"discriminant": discriminant, "radius": radius, "smooth_bound": bound,
            "class_number": len(forms), "factor_base_primes": primes,
            "factor_base_subgroup_size": len(subgroup), "smooth_set_size": len(selected),
            "quotient_set_size": len(quotients), "norm_word_target_classes": len(lookup),
            "norm_word_classes_outside_smooth_quotient": len(set(lookup)-quotients),
            "uniform_sample_size": len(words), "sample_word_bound": sample_bound,
            "norm_word_bound": norm_word_bound, "diagonal": diagonal,
            "relation_rows": rows, "relation_determinant": relation_det,
            "global_class_number_multiple": global_multiple,
            "recovered_subgroup_two_torsion_size": len(span),
            "relation_witnesses": witnesses, **counts}


def main():
    records = [check(value) for value in (
        -348, -420, -999, -1679, -2639, -9679,
        -10007, -41612, -43228, -66548,
    )]
    # A strict weakening: Q has a smooth norm 24 <= 5^2, but its minimum
    # primitive represented norm is 24 and 24 has no two factors <= 5.
    strict = check(-4799, radius=5, bound=3)
    form = (24, 1, 50)
    witness = norm_to_prime_word(form, 1, 0, 3, 25)
    assert witness is not None and witness["norm"] == 24
    assert not exact_balance(24, 5)
    assert strict["norm_word_classes_outside_smooth_quotient"] > 0
    records.append(strict)
    # Keep the established missing-global-direction example visible.
    missing = records[0]
    assert missing["factor_base_subgroup_size"] == 3 < missing["class_number"] == 6
    assert missing["recovered_subgroup_two_torsion_size"] == 1
    print(json.dumps({
        "scope": "Exact supplied smooth-norm words and diagonally dominant relation matrices; all witness and class enumeration is finite control only",
        "discriminants_checked": len(records),
        "representations_checked": sum(row["representations_checked"] for row in records),
        "unbalanced_norms_accepted": sum(row["unbalanced_norms_accepted"] for row in records),
        "relation_rows_checked": sum(len(row["relation_rows"]) for row in records),
        "relations_using_unbalanced_norms": sum(not witness["balanced_at_radius"]
            for row in records for witness in row["relation_witnesses"]),
        "strictly_weaker_oracle_control": {"form": list(form), "discriminant": -4799,
            "radius": 5, "smooth_bound": 3, "norm": 24,
            "word": witness["word"], "no_balanced_norm_at_most_radius_squared": True},
        "records": records,
        "fast_witness_search_supplied": False,
        "global_two_torsion_coverage_supplied": False,
    }, indent=2))


if __name__ == "__main__":
    main()
