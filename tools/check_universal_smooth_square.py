#!/usr/bin/env python3
"""Finite exhaustive certificates for fixed smooth-square coverage.

All group and shortest-word enumeration here is verification only.  Neither
enumeration is part of the conditional efficient representation oracle.
"""

from collections import deque
from functools import lru_cache
from math import gcd, isqrt, lcm
import json

from binary_quadratic_forms import reduced_forms
from bounded_smooth_norm import smooth_factors
from check_positive_relation_kernel import determinant
from quadratic_form_arithmetic import compose_primitive, inverse_primitive, power_primitive
from smooth_norm_relation import norm_to_prime_word
from universal_smooth_square import (
    dominant_square_relation, representation_from_square_word,
    universal_square_parameters, word_from_square_representation,
)


def check(discriminant, smooth_bound=7, radius=None):
    if radius is None:
        radius = isqrt(-discriminant-1)//2
    smooth_bound = min(smooth_bound, radius)
    assert radius >= 2 and 4*radius*radius < -discriminant
    forms = reduced_forms(discriminant)
    identity = next(form for form in forms if form[0] == 1)
    multiply = lru_cache(None)(compose_primitive)
    primes, roots, generators = [], {}, []
    for prime in range(2, smooth_bound+1):
        if (any(prime % divisor == 0 for divisor in range(2, isqrt(prime)+1))
                or gcd(prime, discriminant) != 1):
            continue
        root = next((value for value in range(prime+1)
                     if (value*value-discriminant) % (4*prime) == 0), None)
        if root is not None:
            primes.append(prime)
            roots[prime] = root
            generators.append((prime, root, (root*root-discriminant)//(4*prime)))
    size = len(primes)

    def word_value(word):
        result = identity
        for generator, exponent in zip(generators, word):
            result = multiply(result, power_primitive(generator, exponent))
        return result

    subgroup = {identity}
    queue = deque([identity])
    while queue:
        value = queue.popleft()
        for generator in generators:
            neighbor = multiply(value, generator)
            if neighbor not in subgroup:
                subgroup.add(neighbor)
                queue.append(neighbor)
    selected = {}
    for form in forms:
        leading = form[0]
        if (leading > radius or gcd(leading, discriminant) != 1
                or smooth_factors(leading, smooth_bound) is None):
            continue
        witness = norm_to_prime_word(form, 1, 0, smooth_bound)
        assert witness is not None
        word = [0]*size
        for entry in witness["word"]:
            word[primes.index(entry["prime"])] = entry["exponent"]
        assert word_value(word) == form
        selected[form] = word
    assert identity in selected and set(generators) <= set(selected) <= subgroup
    index_bound = (len(forms)+len(selected)-1)//len(selected)
    parameters = universal_square_parameters(primes, radius, index_bound)
    cap = parameters["exponent_cap"]
    quotient_words = {}
    for left, left_word in selected.items():
        for right, right_word in selected.items():
            quotient = multiply(left, inverse_primitive(right))
            candidate = [a-b for a, b in zip(left_word, right_word)]
            current = quotient_words.get(quotient)
            if current is None or (sum(map(abs, candidate)), candidate) < (sum(map(abs, current)), current):
                quotient_words[quotient] = candidate
    assert len(quotient_words) >= len(selected)
    shortest_words = {identity: [0]*size}
    depths = {identity: 0}
    queue = deque([identity])
    while queue:
        value = queue.popleft()
        for step, step_word in quotient_words.items():
            neighbor = multiply(value, step)
            if neighbor not in shortest_words:
                shortest_words[neighbor] = [a+b for a, b in zip(shortest_words[value], step_word)]
                depths[neighbor] = depths[value]+1
                queue.append(neighbor)
    assert set(shortest_words) == subgroup
    diameter = max(depths.values())
    assert diameter < 2*index_bound
    assert (diameter//2+1)*len(selected) <= len(subgroup)
    assert (diameter//3+1)*len(quotient_words) <= len(subgroup)
    maximum_word_length = max(sum(map(abs, word)) for word in shortest_words.values())
    assert maximum_word_length <= 2*(radius.bit_length()-1)*diameter <= cap
    for value, word in shortest_words.items():
        assert word_value(word) == value
    square_roots = {}
    for value in sorted(subgroup):
        square = multiply(value, value)
        word = shortest_words[value]
        old = square_roots.get(square)
        if old is None or (sum(map(abs, word)), word) < (sum(map(abs, old)), old):
            square_roots[square] = word
    nonprimitive = 0
    maximum_primitive_bits = 0
    for target, short_word in square_roots.items():
        witness = representation_from_square_word(target, primes, roots, short_word, cap)
        recovered = word_from_square_representation(target, *witness["vector"], primes, cap)
        assert word_value(recovered["word"]) == target
        assert recovered["half_word"] == short_word
        nonprimitive += witness["content"] > 1
        maximum_primitive_bits = max(maximum_primitive_bits, witness["primitive_norm"].bit_length())
    rows, witnesses = [], []
    for index, generator in enumerate(generators):
        target_exponent = parameters["target_exponent"]
        target = power_primitive(generator, target_exponent)
        assert target in square_roots
        witness = representation_from_square_word(target, primes, roots, square_roots[target], cap)
        recovered = word_from_square_representation(target, *witness["vector"], primes, cap)
        row = dominant_square_relation(index, primes, recovered["half_word"], cap)
        assert word_value(row) == identity
        assert row[index] > sum(abs(value) for position, value in enumerate(row) if position != index)
        assert target_exponent > 2*cap
        # Padding the supplied long word would need a negative rational content
        # exponent T-K at this prime.  It cannot represent the fixed norm.
        assert cap-target_exponent//2 < 0
        rows.append(row)
        witnesses.append({"generator_prime": primes[index], "target": target,
                          "primitive_norm": recovered["primitive_norm"],
                          "returned_word": recovered["word"],
                          "dominance_slack": row[index]-sum(abs(value) for position, value in enumerate(row)
                                                              if position != index)})
    relation_determinant = abs(determinant(rows))
    assert relation_determinant > 0 and relation_determinant % len(subgroup) == 0
    correction = lcm(*range(1, index_bound+1))
    global_multiple = relation_determinant*correction
    assert global_multiple % len(forms) == 0
    subgroup_two = {value for value in subgroup if multiply(value, value) == identity}
    global_two = {value for value in forms if multiply(value, value) == identity}
    return {"discriminant": discriminant, "radius": radius, "smooth_bound": smooth_bound,
            "class_number": len(forms), "factor_base_primes": primes,
            "subgroup_order": len(subgroup), "square_subgroup_order": len(square_roots),
            "smooth_set_size": len(selected), "quotient_set_size": len(quotient_words),
            "density_reciprocal_bound": index_bound, "quotient_cayley_diameter": diameter,
            "maximum_constructed_root_word_length": maximum_word_length,
            "exponent_cap_fixed_before_targets": cap,
            "fixed_norm_bits": parameters["fixed_norm"].bit_length(),
            "target_exponent": parameters["target_exponent"],
            "universal_representations_checked": len(square_roots),
            "nonprimitive_representations": nonprimitive,
            "largest_primitive_norm_bits": maximum_primitive_bits,
            "relation_rows": rows, "relation_determinant": relation_determinant,
            "global_class_number_multiple": global_multiple,
            "subgroup_two_torsion_size": len(subgroup_two),
            "global_two_torsion_size": len(global_two), "relation_witnesses": witnesses}


def main():
    records = [check(discriminant) for discriminant in range(-512, -19)
               if discriminant % 4 in (0, 1)]
    records += [check(value) for value in (-1679, -2639, -4799, -9679, -10007, -41612, -43228, -66548)]
    missing = next(record for record in records if record["discriminant"] == -348)
    assert missing["subgroup_order"] == 3 and missing["class_number"] == 6
    assert missing["subgroup_two_torsion_size"] == 1 < missing["global_two_torsion_size"] == 2
    print(json.dumps({
        "scope": "Fixed smooth-square coverage and dominant relations, using exhaustive finite word oracles only; no efficient representation search is supplied",
        "discriminants_checked": len(records),
        "universal_representations_checked": sum(record["universal_representations_checked"] for record in records),
        "relation_rows_checked": sum(len(record["relation_rows"]) for record in records),
        "maximum_fixed_norm_bits": max(record["fixed_norm_bits"] for record in records),
        "trivial_subgroup_controls": sum(record["subgroup_order"] == 1 for record in records),
        "missing_global_two_direction_controls": sum(record["subgroup_two_torsion_size"] < record["global_two_torsion_size"] for record in records),
        "records": records,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
