#!/usr/bin/env python3
"""Independent finite checks for quadratic-code queries and annihilators."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import random

from hamming_restrictions import histogram_profile, word_histogram
from hamming_tensor_power import local_terms
from quadratic_crossing import (
    alternating_rank, binary_rank, binomial_profile_from_histogram,
    coefficient, expanded_profile, factored_local_terms, gauss_sum,
    joint_histogram, label_polynomial, polynomial_product,
    quadratic_annihilator_dimension, quadratic_monomials, truth_word,
)


def gauss_checks():
    cases = 0
    rng = random.Random(3102026)
    for variables in range(1, 9):
        monomials = quadratic_monomials(variables)
        labels = range(1 << len(monomials)) if variables <= 4 else [
            rng.getrandbits(len(monomials)) for _ in range(512)
        ]
        for label in labels:
            polynomial = label_polynomial(variables, len(monomials), label)
            expected = (1 << variables) - 2 * truth_word(variables, polynomial).bit_count()
            assert gauss_sum(variables, polynomial) == expected
            cases += 1
    return cases


def histogram_checks():
    cases = 0
    coefficient_checks = 0
    rng = random.Random(300102026)
    for variables in range(1, 8):
        width = len(quadratic_monomials(variables))
        if variables == 1:
            tuples = list(product(range(4), repeat=4))
        else:
            tuples = [tuple(rng.getrandbits(width) for _ in range(4)) for _ in range(200)]
            for _ in range(25):
                x, y, z = (rng.getrandbits(width) for _ in range(3))
                tuples += [(x, y, y, x), (x, y, z, y), (x, y, z, z), (x, x, y, y)]
        for labels in tuples:
            polys = tuple(label_polynomial(variables, width, label) for label in labels)
            words = tuple(truth_word(variables, p) for p in polys)
            histogram = joint_histogram(variables, polys)
            assert histogram == word_histogram(1 << variables, words)
            profile = binomial_profile_from_histogram(histogram)
            actual = expanded_profile(profile)
            assert actual == histogram_profile(histogram)
            assert actual.bit_count() % 2 == int(labels[0] == labels[3] and labels[1] == labels[2])
            for degree in range(actual.bit_length() + 2):
                assert coefficient(profile, degree) == (actual >> degree) & 1
                coefficient_checks += 1
            cases += 1
    return cases, coefficient_checks


def annihilator_checks():
    rng = random.Random(261030)
    cases = 0
    low_rank_cases = []
    for variables in range(2, 13):
        monomials = quadratic_monomials(variables)
        monomial_words = [truth_word(variables, frozenset([m])) for m in monomials]
        for pairs in range(1, variables // 2 + 1):
            top = frozenset((1 << (2*j)) | (1 << (2*j+1)) for j in range(pairs))
            affines = range(1 << (variables + 1)) if variables <= 7 else [
                0, 1, (1 << (variables + 1)) - 1,
                *(rng.getrandbits(variables + 1) for _ in range(24)),
            ]
            for affine in affines:
                polynomial = top.symmetric_difference(frozenset(
                    m for j, m in enumerate(monomials[:variables+1]) if affine >> j & 1
                ))
                rank = alternating_rank(variables, polynomial)
                assert rank == 2 * pairs
                word = truth_word(variables, polynomial)
                nullity = len(monomials) - binary_rank(word & w for w in monomial_words)
                if rank >= 6:
                    assert nullity == 1
                    complement = polynomial.symmetric_difference({0})
                    assert not polynomial_product(polynomial, complement)
                    cases += 1
                elif not affine and variables <= 6:
                    low_rank_cases.append({'variables': variables, 'rank': rank, 'nullity': nullity})
    # Rank four really cannot replace rank six in the uniform assertion.
    f = frozenset((3, 12))
    assert alternating_rank(4, f) == 4
    assert quadratic_annihilator_dimension(4, f) == 5
    return cases, low_rank_cases


def local_support_checks():
    forms = [sum((((a >> 4) & 1) ^ ((p & (a & 15)).bit_count() & 1)) << p
                 for p in range(16)) for a in range(32)]
    expressions = {
        (0, 5): ((5, 16),),
        (1, 6): ((1, 2), (4, 5)),
        (2, 7): ((1, 5), (2, 3)),
        (3, 8): ((1, 4), (2, 3)),
        (4, 20): ((6, 16), (1, 2), (1, 9)),
        (9, 21): ((6, 16), (2, 3), (8, 9)),
        (10, 19): ((6, 16), (1, 2), (1, 8)),
        (11, 12): ((1, 3), (4, 5)),
        (13, 14): ((3, 16),),
        (15, 18): ((6, 16), (1, 8), (1, 3)),
        (16, 17): ((6, 16), (1, 3), (1, 9)),
        (22,): ((6, 16), (9, 16)),
    }
    specified = {i: pairs for indices, pairs in expressions.items() for i in indices}
    result = []
    for index, (local, _) in enumerate(local_terms()):
        support = sum(bool(v) << p for p, v in enumerate(local))
        products = sorted(set(forms[a] & forms[b] for a in range(32) for b in range(a, 32)
                              if not (forms[a] & forms[b] & support)))
        zero_set = sum(all(not (f >> p & 1) for f in products) << p for p in range(16))
        assert zero_set == support
        table_zero_set = sum(all(not ((forms[a] & forms[b]) >> p & 1)
                                 for a, b in specified[index]) << p for p in range(16))
        assert table_zero_set == support
        result.append({'term': index, 'support_size': support.bit_count(),
                       'distinct_nonzero_product_relations': len(products) - 1})
    assert len(factored_local_terms()) == 23
    return result


def high_rank_support_checks():
    variables = 6
    rng = random.Random(610302026)
    q0 = frozenset((3, 12, 48))
    pairs = tuple((1 << i) | (1 << j) for i, j in combinations(range(variables), 2))
    for _ in range(10000):
        q1 = frozenset(m for m in pairs if rng.getrandbits(1))
        if min(alternating_rank(variables, q1),
               alternating_rank(variables, q0.symmetric_difference(q1))) >= 6:
            break
    else:
        raise AssertionError('deterministic finite high-rank example not found')
    polynomials = [frozenset()]
    for generator in (frozenset([0]), q0, q1):
        polynomials += [p.symmetric_difference(generator) for p in polynomials]
    assert all(not (p - {0}) or alternating_rank(variables, p) >= 6 for p in polynomials)
    words = [truth_word(variables, p) for p in polynomials]
    supports = [sum(bool(value) << i for i, value in enumerate(local)) for local, _ in local_terms()]
    cases, nonuniform = 0, 0
    for labels in product(range(8), repeat=4):
        histogram = word_histogram(1 << variables, tuple(words[i] for i in labels))
        image = [i for i, count in enumerate(histogram) if count]
        base = image[0]
        directions = {0}
        for point in image:
            directions |= {v ^ point ^ base for v in directions}
        hull = {base ^ v for v in directions}
        image_mask = sum(1 << i for i in image)
        hull_mask = sum(1 << i for i in hull)
        for support in supports:
            assert bool(image_mask & ~support) == bool(hull_mask & ~support)
        nonuniform += len({histogram[i] for i in image}) > 1
        cases += 1
    assert nonuniform
    return {'quartets': cases, 'nonuniform_histograms': nonuniform,
            'quadratic_generators': [sorted(q0), sorted(q1)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    gauss = gauss_checks()
    histograms, coefficients = histogram_checks()
    annihilators, controls = annihilator_checks()
    support = local_support_checks()
    high_rank_support = high_rank_support_checks()
    result = {
        'scope': 'Quadratic-code coefficient oracle and rank-six annihilator lemma; no improved factoring exponent',
        'gauss_sum_comparisons': gauss,
        'joint_histogram_and_source_profile_comparisons': histograms,
        'coefficient_parity_checks': coefficients,
        'rank_at_least_six_annihilator_cases': annihilators,
        'low_rank_controls': controls,
        'local_support_relations': support,
        'high_rank_support_checks': high_rank_support,
        'source_commit': '05d9c1b8aed4f8494ce230f0869f9048ef879474',
        'sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in [
            Path(__file__), Path(__file__).with_name('quadratic_crossing.py'),
        ]},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
