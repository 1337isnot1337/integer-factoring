#!/usr/bin/env python3
"""Exact supplied-norm controls; all witness searches here are finite scans."""

from collections import deque
from math import gcd, isqrt
import json

from arithmetic_two_group import power
from binary_quadratic_forms import (
    compose_forms, inverse_form, reduce_form, reduced_forms,
)
from bounded_smooth_norm import balanced_divisor, norm_to_pair, smooth_factors
from smooth_norm_words import (
    dominant_relation, exponent_vector, norm_to_word, word_value,
)


DISCRIMINANTS = (
    -23, -31, -47, -71, -80, -84, -120, -144, -180, -231,
    -348, -420, -999, -1007, -1679, -4620, -10007,
)


def identity_form(discriminant):
    parity = discriminant % 2
    return reduce_form((1, parity, (parity - discriminant) // 4))


def factor_base(discriminant, bound):
    result = []
    for prime in range(2, bound + 1):
        if any(prime % divisor == 0 for divisor in range(2, isqrt(prime) + 1)):
            continue
        if discriminant % prime == 0:
            continue
        roots = [root for root in range(prime + 1)
                 if (root * root - discriminant) % (4 * prime) == 0]
        if roots:
            root = min(roots)
            result.append((prime, reduce_form((
                prime, root, (root * root - discriminant) // (4 * prime)
            ))))
    return result


def evaluate_vector(vector, base, identity):
    result = identity
    for exponent, (_, value) in zip(vector, base):
        value = value if exponent >= 0 else inverse_form(value)
        result = compose_forms(result, power(value, abs(exponent), identity, compose_forms))
    return result


def finite_word_table(base, identity):
    """Exhaust the generated subgroup; no fast-search claim is attached."""
    words = {identity: (0,) * len(base)}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for index, (_, generator) in enumerate(base):
            for sign in (-1, 1):
                value = generator if sign == 1 else inverse_form(generator)
                target = compose_forms(current, value)
                if target not in words:
                    word = list(words[current])
                    word[index] += sign
                    words[target] = tuple(word)
                    queue.append(target)
    return words


def determinant(matrix):
    """Fraction-free elimination with exact pivoting."""
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
                assert numerator % previous == 0
                rows[row][target] = numerator // previous
            rows[row][col] = 0
        previous = value
    return sign * rows[-1][-1]


def relation_control(discriminant, bound, forms):
    identity = identity_form(discriminant)
    base = factor_base(discriminant, bound)
    if not base:
        return None
    samples = finite_word_table(base, identity)
    sample_budget = max(sum(map(abs, word)) for word in samples.values())
    radius = isqrt(-discriminant - 1) // 2
    word_budget = (radius * radius).bit_length() - 1
    dominant = sample_budget + word_budget + 1
    primes = [prime for prime, _ in base]
    rows, witnesses = [], []
    scan = [(x, y) for x in range(-3, 4) for y in range(-3, 4)
            if gcd(x, y) == 1]
    for index, (_, generator) in enumerate(base):
        target = power(generator, dominant, identity, compose_forms)
        found = None
        for sample, sample_word in sorted(samples.items(), key=lambda item: (sum(map(abs, item[1])), item[1])):
            form = compose_forms(target, sample)
            a, b, c = form
            for x, y in sorted(scan, key=lambda vector: (
                    a * vector[0] ** 2 + b * vector[0] * vector[1] + c * vector[1] ** 2,
                    vector)):
                word = norm_to_word(form, x, y, bound, word_budget)
                if word is not None:
                    found = sample_word, word
                    break
            if found is not None:
                break
        # The inverse of the target is in the exhaustively generated sample;
        # it reduces the represented class to the principal class, norm 1.
        assert found is not None
        sample_word, word = found
        encoded = exponent_vector(word, primes)
        row = dominant_relation(index, sample_word, encoded, word_budget, sample_budget)
        assert evaluate_vector(row, base, identity) == identity
        assert sum(map(abs, row)) <= 2 * dominant - 1
        rows.append(row)
        witnesses.append({
            "norm": word.primitive_norm,
            "multiplicity": word.multiplicity,
            "sample_word": sample_word,
            "target_word": encoded,
        })

        tautological = list(sample_word)
        tautological[index] += dominant
        try:
            dominant_relation(index, sample_word, tautological, word_budget, sample_budget)
        except ValueError:
            pass
        else:
            raise AssertionError("tautological input word passed the strict short-word budget")

    annihilator = abs(determinant(rows))
    assert 0 < annihilator <= (2 * dominant - 1) ** len(base)
    assert all(power(generator, annihilator, identity, compose_forms) == identity
               for _, generator in base)
    return {
        "factor_base_primes": primes,
        "generated_subgroup_size": len(samples),
        "full_class_group_generated": len(samples) == len(forms),
        "sample_word_budget": sample_budget,
        "target_word_budget": word_budget,
        "diagonal_parameter": dominant,
        "relation_matrix": rows,
        "annihilator": annihilator,
        "annihilator_bit_bound": 1 + len(base) * (2 * dominant - 2).bit_length(),
        "witnesses": witnesses,
    }


def main():
    smooth_bound = 7
    records = []
    totals = {
        "accepted_primitive_representations": 0,
        "word_replays": 0,
        "normalized_scaled_controls": 0,
        "strict_budget_rejections": 0,
        "coprimality_rejections": 0,
        "nonsmooth_rejections": 0,
        "accepted_beyond_old_norm_bound": 0,
        "accepted_without_old_balanced_split": 0,
    }
    two_signs = set()
    for discriminant in DISCRIMINANTS:
        forms = reduced_forms(discriminant)
        radius = isqrt(-discriminant - 1) // 2
        maximum_norm = radius * radius
        budget = maximum_norm.bit_length() - 1
        counts = dict.fromkeys(totals, 0)
        for form in forms:
            a, b, c = form
            for x in range(-4, 5):
                for y in range(-4, 5):
                    if gcd(x, y) != 1:
                        continue
                    norm = a * x * x + b * x * y + c * y * y
                    result = norm_to_word(form, x, y, smooth_bound, budget)
                    factors = smooth_factors(norm, smooth_bound)
                    if gcd(norm, discriminant) != 1:
                        assert result is None
                        counts["coprimality_rejections"] += 1
                        continue
                    if factors is None:
                        assert result is None
                        counts["nonsmooth_rejections"] += 1
                        continue
                    if len(factors) > budget:
                        assert result is None
                        continue
                    assert result is not None
                    assert result.primitive_norm == norm
                    assert result.multiplicity == len(factors)
                    assert word_value(result, discriminant) == form
                    counts["accepted_primitive_representations"] += 1
                    counts["word_replays"] += 1
                    for term in result.terms:
                        if term.prime == 2:
                            two_signs.add(1 if term.exponent > 0 else -1)
                    if result.multiplicity:
                        assert norm_to_word(form, x, y, smooth_bound, result.multiplicity - 1) is None
                        counts["strict_budget_rejections"] += 1
                    if norm > maximum_norm:
                        counts["accepted_beyond_old_norm_bound"] += 1
                    elif balanced_divisor(factors, radius) is None:
                        assert norm_to_pair(form, x, y, smooth_bound, radius) is None
                        counts["accepted_without_old_balanced_split"] += 1
                    for scale in (-3, 2):
                        scaled = norm_to_word(form, scale * x, scale * y, smooth_bound, budget)
                        assert scaled is not None
                        assert scaled.primitive_norm == result.primitive_norm
                        assert scaled.terms == result.terms
                        counts["normalized_scaled_controls"] += 1
            assert norm_to_word(form, 0, 0, smooth_bound, budget) is None
        for key, value in counts.items():
            totals[key] += value
        records.append({
            "discriminant": discriminant,
            "class_number": len(forms),
            "old_radius": radius,
            "old_norm_bound": maximum_norm,
            "multiplicity_budget": budget,
            "counts": counts,
            "finite_relation_control": relation_control(discriminant, smooth_bound, forms),
        })

    assert two_signs == {-1, 1}
    unbalanced = norm_to_word((2, -1, 125), 0, 1, 7, 7)
    assert unbalanced is not None and unbalanced.multiplicity == 3
    assert word_value(unbalanced, -999) == (2, -1, 125)
    assert norm_to_pair((2, -1, 125), 0, 1, 7, 12) is None
    assert totals["accepted_beyond_old_norm_bound"] > 0
    assert totals["accepted_without_old_balanced_split"] > 0
    assert norm_to_word(identity_form(-23), 1, 0, 7, 0).terms == ()
    print(json.dumps({
        "scope": "Direct supplied-smooth-norm to short-word conversion; finite witness and subgroup searches are exhaustive",
        "smooth_bound": smooth_bound,
        "coordinate_scan": [-4, 4],
        "totals": totals,
        "both_orientations_above_two_checked": True,
        "unbalanced_125_control": {
            "discriminant": -999,
            "form": [2, -1, 125],
            "primitive_norm": unbalanced.primitive_norm,
            "multiplicity": unbalanced.multiplicity,
            "word": [{"prime": term.prime, "exponent": term.exponent, "form": term.form}
                     for term in unbalanced.terms],
            "balanced_converter_accepts": False,
            "direct_converter_accepts": True,
        },
        "records": records,
        "fast_witness_search_supplied": False,
        "new_factoring_exponent_established": False,
    }, indent=2))


if __name__ == "__main__":
    main()
