#!/usr/bin/env python3
"""Exact odd-residual certificates, quotient kernel and two static class groups."""
import argparse
from fractions import Fraction
from functools import cache
from hashlib import sha256
import json
from math import prod
from pathlib import Path

from arithmetic_two_group import power, order_log_two, project_to_two_part, two_group_basis
from binary_quadratic_forms import reduced_forms, smooth_prime_word, ambiguous_factor
from check_positive_relation_kernel import determinant
from odd_quotient_relations import odd_relation_annihilator, uniform_odd_multiplier, word_value
from quadratic_form_arithmetic import compose_primitive, inverse_primitive

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/odd-quotient-relations.json'


def product_control(moduli, generators, rows, witnesses, expected_logs):
    identity = (0,)*len(moduli)
    multiply = lambda a, b: tuple((x+y) % n for x, y, n in zip(a, b, moduli))
    inverse = lambda a: tuple(-x % n for x, n in zip(a, moduli))
    annihilator = odd_relation_annihilator(generators, rows, witnesses, identity, multiply, inverse)
    assert all(power(g, annihilator, identity, multiply) == identity for g in generators)
    projected = project_to_two_part(generators, annihilator, identity, multiply)
    basis = two_group_basis(projected, identity, multiply, inverse, prod(moduli).bit_length())
    assert sorted(k for _, k in basis) == expected_logs
    return dict(moduli=moduli, relation_matrix=rows, odd_witnesses=witnesses,
                determinant=determinant(rows), annihilator=annihilator,
                nonidentity_row_residues=sum(word_value(generators, row, identity, multiply, inverse)
                                            != identity for row in rows),
                extracted_two_order_logs=expected_logs)


def quotient_kernel_control():
    # Reuse the predecessor's exact quadratic small-bias set and five-set.
    sample = [x for x in range(16) if ((x & 1)*((x >> 1) & 1)
                                      ^ ((x >> 2) & 1)*((x >> 3) & 1)) == 1]
    selected = [x for x in range(16) if x not in sample][:5]
    identity = (0, 0)
    multiply = lambda a, b: (a[0] ^ b[0], (a[1]+b[1]) % 3)
    inverse = lambda a: (a[0], -a[1] % 3)
    generators = [(1 << i, 0) for i in range(4)]+[(0, 1)]
    word = lambda x: [(x >> i) & 1 for i in range(4)]+[0]
    diagonal = 1+max(x.bit_count() for x in sample)+2*max(x.bit_count() for x in selected)
    witnesses = {t: next((u, a, b) for u in sample for a in selected for b in selected
                         if t ^ u == a ^ b) for t in range(16)}
    rows = []
    for i, g in enumerate(generators):
        target = power(g, diagonal, identity, multiply)[0]
        u, a, b = witnesses[target]
        rows.append([diagonal*int(i == j)+word(u)[j]-word(a)[j]+word(b)[j]
                     for j in range(5)])
        assert rows[-1][i] > sum(abs(x) for j, x in enumerate(rows[-1]) if j != i)
    E = odd_relation_annihilator(generators, rows, [3]*5, identity, multiply, inverse)
    basis = two_group_basis(project_to_two_part(generators, E, identity, multiply),
                            identity, multiply, inverse, 6)
    assert len(basis) == 4 and all(k == 1 for _, k in basis)
    # For every target with nonzero odd coordinate there is NO exact pair
    # with these lifts, but there is a pair after the M=3 map.
    missed = 0
    for even in range(16):
        for odd in range(3):
            exact = any(even ^ u == a ^ b and odd == 0
                        for u in sample for a in selected for b in selected)
            missed += not exact
            u, a, b = witnesses[even]
            residual = (even ^ u ^ a ^ b, odd)
            assert power(residual, 3, identity, multiply) == identity
    assert missed == 32
    image_density = Fraction(len(selected), 16)
    original_density = Fraction(len(selected), 48)
    image_bias = max(abs(Fraction(sum(1 if (a & x).bit_count() % 2 == 0 else -1
                                      for x in sample), len(sample))) for a in range(1, 16))
    lower = image_density**2-image_bias*(image_density-image_density**2)
    correlation = {t: Fraction(sum((a ^ t) in selected for a in selected), 16)
                   for t in range(16)}
    averages = [sum(correlation[t ^ u] for u in sample)/len(sample) for t in range(16)]
    assert image_density >= original_density and image_bias == Fraction(1, 3)
    assert min(averages) >= lower > 0
    # The nontrivial character of the odd factor is identically one on
    # the lifted sample, proving full-group bias one without numerics.
    return dict(group='(Z/2)^4 x Z/3', multiplier=3, original_density='5/48',
                image_density=str(image_density), image_sample_bias=str(image_bias),
                correlation_lower_bound=str(lower),
                minimum_sampled_correlation=str(min(averages)),
                full_sample_bias='1', targets=48, exact_pair_targets_missed=missed,
                quotient_pair_targets=48, relation_matrix=rows, annihilator=E,
                extracted_two_rank=len(basis))


def class_control(discriminant, target):
    # Enumeration supplies ONLY these finite fixtures, never a uniform
    # generator or relation acquisition algorithm.
    forms = reduced_forms(discriminant)
    identity = next(f for f in forms if f[0] == 1)
    multiply = cache(compose_primitive)
    M = uniform_odd_multiplier((-discriminant).bit_length(), 11)
    image = lambda f: power(f, M, identity, multiply)
    images = {image(f) for f in forms}
    selected = [f for f in forms if smooth_prime_word(f, 11) is not None]
    selected_images = {image(f) for f in selected}
    assert len(selected_images)*len(forms) >= len(selected)*len(images)
    match = next(f for f in selected if image(f) == image(target))
    residual = multiply(target, inverse_primitive(match))
    assert residual != identity and power(residual, M, identity, multiply) == identity
    assert match[0] < target[2]
    rows = []
    for i, form in enumerate(forms):
        k = order_log_two(image(form), identity, multiply, (-discriminant).bit_length())
        rows.append([(1 << k)*int(i == j) for j in range(len(forms))])
    E = odd_relation_annihilator(forms, rows, [M]*len(forms), identity, multiply, inverse_primitive)
    projected = project_to_two_part(forms, E, identity, multiply)
    assert projected == [image(f) for f in forms]
    basis = two_group_basis(projected, identity, multiply, inverse_primitive,
                            (-discriminant).bit_length())
    assert 1 << sum(k for _, k in basis) == len(images)
    factors = [ambiguous_factor(power(f, 1 << (k-1), identity, multiply), -discriminant)
               for f, k in basis]
    assert any(factors)
    return dict(discriminant=discriminant, enumerated_classes=len(forms), smooth_bound=11,
                smooth_forms=len(selected), image_classes=len(images),
                image_smooth_forms=len(selected_images), multiplier=M,
                target=list(target), strict_smaller_quotient_preimage=list(match),
                nonidentity_odd_residual=list(residual), relation_determinant_bits=abs(determinant(rows)).bit_length(),
                annihilator_bits=E.bit_length(), extracted_two_order_logs=[k for _, k in basis],
                direct_filtered_two_group_test=True,
                recovered_factors=factors)


def run():
    products = [product_control([8, 2, 3], [(1, 0, 1), (1, 1, 1), (0, 0, 1)],
                               [[8, 0, 0], [-2, 2, 0], [0, 0, 1]], [3, 1, 3], [1, 3]),
                product_control([4, 3, 5], [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                                [[4, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 3, 5], [2])]
    multiplier_controls = 0
    for bits in range(9):
        for bound in (0, 2, 3, 5, 11, 17):
            M = uniform_odd_multiplier(bits, bound)
            assert M % 2 == 1
            for odd in range(1, (1 << bits)+1, 2):
                remainder = odd
                for p in range(3, bound+1, 2):
                    while remainder % p == 0:
                        remainder //= p
                if remainder == 1:
                    assert M % odd == 0
            multiplier_controls += 1
    rejected = []
    for rows, witnesses, label in (([[1]], [2], 'even_witness'),
                                    ([[1]], [3], 'false_odd_witness'),
                                    ([[0]], [1], 'singular_matrix')):
        try:
            odd_relation_annihilator([1], rows, witnesses, 0, int.__xor__, lambda x: x)
        except ValueError:
            rejected.append(label)
        else:
            raise AssertionError(label)
    # A large odd component is not mistaken for a two-group by the direct
    # smooth-odd filter. No class-number or odd-order oracle is used.
    M = uniform_odd_multiplier(5, 11)
    residual = power((1, 1), M, (0, 0),
                     lambda a, b: ((a[0]+b[0]) % 2, (a[1]+b[1]) % 13))
    try:
        order_log_two(residual, (0, 0),
                      lambda a, b: ((a[0]+b[0]) % 2, (a[1]+b[1]) % 13), 5)
    except ValueError:
        rejected.append('direct_filter_residual_odd_thirteen')
    else:
        raise AssertionError('Nonsmooth odd component passed the two-group test')
    names = ['tools/odd_quotient_relations.py', 'tools/check_odd_quotient_relations.py',
             'tools/arithmetic_two_group.py', 'tools/check_positive_relation_kernel.py',
             'tools/binary_quadratic_forms.py', 'tools/quadratic_form_arithmetic.py']
    assert (2, 1, 3) != (1, 1, 6)
    assert power((2, 1, 3), 3, (1, 1, 6), compose_primitive) == (1, 1, 6)
    assert 4*2 < 23
    return dict(scope='Exact odd-residual relation certificates and finite quotient controls; no fast quotient smooth-pair extractor, uniform generator claim or stronger factoring bound.',
                product_group_controls=products, quotient_kernel=quotient_kernel_control(),
                multiplier_controls=multiplier_controls, rejected_certificates=rejected,
                principal_gap_not_inherited=dict(discriminant=-23, prime_form=[2, 1, 3],
                    multiplier=3, nonzero_word_norm=2, old_safe_cap=[23, 4]),
                static_class_groups=[class_control(-7423, (29, 1, 64)),
                                     class_control(-26623, (13, 1, 512))],
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECEIPT.read_text()), 'receipt differs; inspect before regeneration'
    print('Odd quotient relations: certificate, image-density and static class controls PASS')
