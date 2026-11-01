#!/usr/bin/env python3
"""Complete small-class replay of powered prime absorption and decoding.

This checker deliberately enumerates the full class groups. It validates
the constructive identities, and is not a fast target-norm selector.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
import json
from math import gcd, isqrt, prod
from pathlib import Path

from binary_quadratic_forms import inverse_form, reduced_forms
from factored_norm_representation import represent_factored_norm
from prime_form_bias import prime_form_sample
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf, represented_norm_from_word


def check_discriminant(discriminant, cutoff, smooth_multiplicity):
    forms = tuple(reduced_forms(discriminant))
    order = len(forms)
    identity = forms[0]
    radius = (isqrt(-discriminant) - 1) // 2
    assert 4 * radius * radius < -discriminant
    base = [record for record in build_prime_form_factor_base(discriminant, cutoff)['records']
            if record.prime != 2]
    by_prime = {record.prime: record for record in base}

    @lru_cache(None)
    def multiply(left, right):
        return compose_forms_hnf(left, right)

    smooth = {}
    for chosen in combinations(base, smooth_multiplicity):
        norm = prod(record.prime for record in chosen)
        if norm > radius:
            continue
        for signs in product((-1, 1), repeat=smooth_multiplicity):
            word = tuple((record.prime, sign, record.root)
                         for record, sign in zip(chosen, signs))
            representation = represented_norm_from_word(discriminant, word)
            assert representation.norm == norm
            assert representation.form not in smooth
            smooth[representation.form] = word
    assert smooth

    differences = {}
    for left, left_word in smooth.items():
        for right, right_word in smooth.items():
            form = multiply(left, inverse_form(right))
            word = left_word + tuple((prime, -exponent, root)
                                     for prime, exponent, root in right_word)
            differences.setdefault(form, word)

    masses, source_records = prime_form_sample(discriminant, cutoff)
    assert masses
    prime_words = {}
    for prime, middle, form, _ in source_records:
        record = by_prime[prime]
        sign = 1 if middle == record.root else -1
        assert middle in (record.root, (-record.root) % (2 * prime))
        assert form == (record.positive if sign == 1 else record.negative)
        prime_words.setdefault(form, ((prime, sign, record.root),))
    assert set(prime_words) == set(masses)

    def covered_targets(words):
        result = {}
        for walk_form, walk_word in words.items():
            inverse_walk = inverse_form(walk_form)
            for difference_form, difference_word in differences.items():
                target = multiply(difference_form, inverse_walk)
                word = difference_word + tuple((prime, -exponent, root)
                                                for prime, exponent, root in walk_word)
                result.setdefault(target, word)
        return result

    first_coverage = len(covered_targets(prime_words))
    distribution = {identity: 1}
    words = {identity: ()}
    source_total = sum(masses.values())
    total = 1
    support_sizes = []
    for steps in range(1, 101):
        new_distribution = defaultdict(int)
        new_words = {}
        for current, current_mass in distribution.items():
            for prime_form, mass in masses.items():
                target = multiply(current, prime_form)
                new_distribution[target] += current_mass * mass
                new_words.setdefault(target, words[current] + prime_words[prime_form])
        distribution = dict(new_distribution)
        words = new_words
        total *= source_total
        assert sum(distribution.values()) == total
        support_sizes.append(len(distribution))
        error_numerator = sum(abs(order * distribution.get(form, 0) - total)
                              for form in forms)
        # l1 distance is error_numerator/(order*total), while delta=|S|/order.
        if 4 * error_numerator <= len(smooth) * total:
            break
    else:
        raise AssertionError("The chosen finite sample did not meet the exact bias budget.")

    witnesses = covered_targets(words)
    assert set(witnesses) == set(forms)
    multiplicity_bound = 2 * smooth_multiplicity + steps
    largest_norm_bits = 0
    largest_multiplicity = 0
    largest_distinct_support = 0
    decoder_multiplications = 0
    examples = []
    for target in forms:
        word = witnesses[target]
        assert len(word) <= multiplicity_bound
        exponents = defaultdict(int)
        for prime, exponent, root in word:
            assert root == by_prime[prime].root
            exponents[prime] += exponent
        compressed = [(prime, exponent, by_prime[prime].root)
                      for prime, exponent in sorted(exponents.items()) if exponent]
        factors = [(prime, abs(exponent), root) for prime, exponent, root in compressed]
        representation = represented_norm_from_word(discriminant, compressed)
        assert representation.form == target
        assert representation.norm == prod(prime ** exponent for prime, exponent, _ in factors)
        assert gcd(representation.norm, discriminant) == 1
        assert sum(exponent for _, exponent, _ in factors) <= multiplicity_bound
        assert representation.norm <= cutoff ** multiplicity_bound
        assert representation.norm <= cutoff ** steps * radius * radius
        decoded = represent_factored_norm(target, factors)
        actual = decoded['representation']
        assert actual is not None and actual.norm == representation.norm
        x, y = actual.vector
        a, middle, c = target
        assert a * x * x + middle * x * y + c * y * y == actual.norm
        assert gcd(x, y) == 1
        largest_norm_bits = max(largest_norm_bits, actual.norm.bit_length())
        largest_multiplicity = max(largest_multiplicity, sum(exponent for _, exponent, _ in factors))
        largest_distinct_support = max(largest_distinct_support, len(factors))
        decoder_multiplications += decoded['stats']['group_multiplications']
        if len(examples) < 3:
            examples.append({'target': target, 'norm': actual.norm,
                             'coordinates': actual.vector, 'factorization': factors})
    return {
        'discriminant': discriminant, 'class_number': order,
        'prime_cutoff': cutoff, 'radius': radius,
        'smooth_form_multiplicity': smooth_multiplicity,
        'smooth_classes': len(smooth), 'difference_classes': len(differences),
        'one_prime_covered_targets': first_coverage,
        'walk_steps': steps, 'walk_support_sizes': support_sizes,
        'exact_l1_below_density_over_four': True,
        'multiplicity_budget': multiplicity_bound,
        'primitive_norms_decoded': order,
        'largest_norm_bits': largest_norm_bits,
        'largest_actual_multiplicity': largest_multiplicity,
        'largest_distinct_prime_support': largest_distinct_support,
        'decoder_group_multiplications': decoder_multiplications,
        'examples': examples,
    }


def main():
    cases = [check_discriminant(-41612, 32, 1),
             check_discriminant(-300004, 32, 2)]
    assert any(case['one_prime_covered_targets'] < case['class_number'] for case in cases)
    evidence = {
        'claim_scope': 'Full small-class enumeration verifies absorption and decoding; no fast norm selector.',
        'cases': cases,
        'total_primitive_norms_decoded': sum(case['primitive_norms_decoded'] for case in cases),
    }
    destination = Path(__file__).resolve().parents[1] / 'evidence' / 'powered-prime-norm-witnesses.json'
    destination.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({'status': 'pass',
                      'primitive_norms_decoded': evidence['total_primitive_norms_decoded'],
                      'evidence': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
