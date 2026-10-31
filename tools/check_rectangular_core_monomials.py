#!/usr/bin/env python3
"""Endpoint-vanishing core monomials, including odd-degree perturbations.

Each specification is (x_or_None,y_or_None,u_or_None,v_or_None,terms,pole).
The constructor expands original Gaussian minors for mixed or full alphabets.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path

from check_rectangular_f4_repair import (
    ZERO, add, scale, multiply, times_binary, binary_mul, binary_power,
    c_order, fmul, one_particle, matrix_multiply, build,
)
from check_rectangular_refined_repairs import _verified_repair


HEIGHT = (0, 0, 1, 1)
OCCUPATION = (1, 0, 1, 2)


def core_height(index):
    return HEIGHT[index % 4]+HEIGHT[index//4]


def occupation_charge(modes):
    return int(modes[1] is not None)+int(modes[3] is not None) \
        -int(modes[0] is not None)-int(modes[2] is not None)


def positive_monomials(top=1, degree=None):
    return [modes for modes in product((None, 0, top), repeat=4)
            if (degree is None or sum(mode is not None for mode in modes) == degree)
            and sum(mode or 0 for mode in (modes[1], modes[3]))
            > sum(mode or 0 for mode in (modes[0], modes[2]))]


def validate_extra_monomials(top, extras, *, allow_marginal=False):
    if top <= 0:
        raise ValueError('The top mode must be positive.')
    result = []
    for extra in extras:
        if len(extra) != 6:
            raise ValueError('A monomial needs four optional modes, Laurent terms, and a c-pole order.')
        *modes, supplied, pole = extra
        if not any(mode is not None for mode in modes) or any(mode not in (None, 0, top) for mode in modes):
            raise ValueError('Use a nonempty extreme-mode monomial, with None for missing legs.')
        if not isinstance(pole, int) or pole < 0:
            raise ValueError('The c-pole order must be a nonnegative integer.')
        coefficients = {}
        for coefficient, exponent in supplied:
            if coefficient not in (0, 1, 2, 3) or not isinstance(exponent, int):
                raise ValueError('A Laurent term needs an F4 coefficient and an integer exponent.')
            coefficients[exponent] = coefficients.get(exponent, 0) ^ coefficient
        terms = [(coefficient, exponent) for exponent, coefficient in sorted(coefficients.items()) if coefficient]
        if not terms:
            continue
        offset = -min(0, min(exponent for _, exponent in terms))
        polynomial = tuple(sum(((coefficient >> bit)&1) << (exponent+offset)
                               for coefficient, exponent in terms) for bit in (0, 1))
        delta = sum(mode or 0 for mode in (modes[1], modes[3])) \
            -sum(mode or 0 for mode in (modes[0], modes[2]))
        excess = delta+(c_order(polynomial)-pole)*top
        if excess < 0 or (excess == 0 and not allow_marginal):
            raise ValueError('The extra monomial must vanish under the ramified energy gauge.'
                             if not allow_marginal else 'The residual monomial must be regular under the ramified energy gauge.')
        result.append((*modes, terms, pole))
    return result


def monomial_program_terms(M, place, specifications, *, allow_marginal=False):
    """Return exact CAR words and natural scalar shifts for validated extras."""
    if place not in ('infinity', 'zero'):
        raise ValueError('Choose one of the two valuation places.')
    infinity = place == 'infinity'
    result = []
    for xa, yb, uc, vd, terms, pole in validate_extra_monomials(M-1, specifications, allow_marginal=allow_marginal):
        if infinity:
            choices = ((vd, 'out_v', '_create'), (yb, 'out_y', '_create'),
                       (uc, 'in_u', '_annihilate'), (xa, 'in_x', '_annihilate'))
            gauge = (M if xa is not None else 0)-(M if yb is not None else 0) \
                +sum(2*mode-M+1 for mode in (uc, vd) if mode is not None)
        else:
            choices = ((xa, 'left_x', '_create'), (yb, 'left_y', '_create'),
                       (uc, 'right_u', '_annihilate'), (vd, 'right_v', '_annihilate'))
            gauge = (1-2*M+2*xa if xa is not None else 0)+(1+2*yb if yb is not None else 0)
        word = [prefix+str(mode)+suffix for mode, prefix, suffix in choices if mode is not None]
        result.append((word, gauge, terms, pole))
    return result


def core_monomial_seed(extras, indices=(0, 1, 2), ambient=4, alpha=1, beta=2,
                       full_alphabet_gamma=None, repair_kind='partial', check_types=True):
    indices = tuple(indices)
    extras = validate_extra_monomials(indices[-1], extras)
    full = full_alphabet_gamma is not None
    if full:
        from check_full_alphabet_quartic_repair import build as full_build
        seed = full_build(indices, ambient, full_alphabet_gamma, check_types=False,
                          independent_exterior=len(indices) <= 3, repair_kind=repair_kind)
        seed['ramified_endpoint_matrix'] = matrix_multiply(seed['repair_matrix'], seed['residual_matrix'])
    else:
        seed = build(indices, ambient, alpha, beta, check_types=False,
                     independent_exterior=len(indices) <= 3, repair=repair_kind)
    if not extras:
        return seed
    n, bulk = len(indices), seed['bulk_size']
    top = 1 << (n-1)
    A, B = seed['q_x'], seed['q_u']
    global_w = max(0, -min(exponent for extra in extras for _, exponent in extra[4]))
    denominator = max(seed['denominator_degree'], 2+max(extra[5] for extra in extras))
    factor = binary_power(3, denominator-seed['denominator_degree']) << global_w
    values = [times_binary(value, factor) for value in seed['values']]
    G = one_particle(indices, ambient)
    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i, answer = bit.bit_length()-1, 0
        remaining = columns
        while remaining:
            column = remaining & -remaining
            remaining ^= column
            j = column.bit_length()-1
            if G[i][j]:
                answer ^= binary_mul(G[i][j], minor(rows ^ bit, columns ^ column))
        return answer
    xmap = [((letter % bulk) << 1)|(1, 0, top, top | 1)[letter//bulk] for letter in range(A)]
    def code(coefficient):
        if full:
            return [[(mask, 1)] for mask in xmap]
        result = []
        for letter in range(B):
            middle = (letter % bulk) << 1
            entries = (((1, 1),), ((0, 1),), ((top, 1), (top | 1, coefficient)))[letter//bulk]
            result.append([(middle | mask, value) for mask, value in entries if value])
        return result
    umap, ymap = code(alpha), code(beta)
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    correction = seed['repair_matrix']
    for a, b, c, d, terms, pole in extras:
        xa, yb, uc, vd = (0 if mode is None else 1 << indices.index(mode) for mode in (a, b, c, d))
        polynomial = tuple(sum(((coefficient >> bit)&1) << (exponent+global_w)
                               for coefficient, exponent in terms) for bit in (0, 1))
        raw = [ZERO]*len(values)
        for x, y, u, v in product(range(A), range(B), range(B), range(A)):
            xx, vv = xmap[x], xmap[v]
            if xx & xa != xa or vv & vd != vd:
                continue
            answer = ZERO
            for (yy, cy), (uu, cu) in product(ymap[y], umap[u]):
                if yy & yb == yb and uu & uc == uc:
                    determinant = minor((vv ^ vd)|((yy ^ yb) << n),
                                        (uu ^ uc)|((xx ^ xa) << n))
                    answer = add(answer, scale(times_binary(polynomial, determinant), fmul(cy, cu)))
            raw[position(x, y, u, v)] = answer
        gates = [[(source, coefficient, binary_power(3, denominator-pole+core_height(source)-core_height(target)))
                  for source, coefficient in enumerate(row) if coefficient]
                 for target, row in enumerate(correction)]
        for x, y, u, v in product(range(A), range(B), range(B), range(A)):
            target = position(x, y, u, v)
            for source, coefficient, multiplier in gates[v//bulk+4*(y//bulk)]:
                old = raw[position(x, (source//4)*bulk+y % bulk, u, (source % 4)*bulk+v % bulk)]
                values[target] = add(values[target], scale(times_binary(old, multiplier), coefficient))
    provisional = dict(seed, denominator_degree=denominator)
    result = _verified_repair(provisional, values, seed['ramified_endpoint_matrix'],
                              [(core_height(i), 0) for i in range(len(correction))], check_types)
    result.update({'repair_kind': 'endpoint_vanishing_core_monomials', 'base_repair_kind': seed['repair_kind'],
                   'extra_monomials': extras, 'full_alphabet_gamma': full_alphabet_gamma,
                   'global_w_power': global_w, 'uncancelled_denominator_degree': denominator,
                   'cancelled_c_order': denominator-result['denominator_degree']})
    if 'non_gaussian_witness' in result:
        result['base_non_gaussian_witness'] = result.pop('non_gaussian_witness')
        result.pop('wick_minor_checks', None)
    return result


def typed_difference(base, changed, require_equal=False):
    """Compare complete two-copy paired-type tensors as rational polynomials."""
    A, B = base['q_x'], base['q_u']
    assert (A, B) == (changed['q_x'], changed['q_u'])
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    xtypes = [sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(A), 2)]
    ytypes = [sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(B), 2)]
    left_factor, right_factor = [binary_power(3, 2*seed['denominator_degree']) for seed in (base, changed)]
    global_w = 2*(changed.get('global_w_power', 0)-base.get('global_w_power', 0))
    assert global_w >= 0
    checks, witness = 0, None
    for X, Y in product(xtypes, ytypes):
        for x, v, y, u in product(X, X, Y, Y):
            positions = [position(x[i], y[i], u[i], v[i]) for i in (0, 1)]
            left = times_binary(multiply(*(changed['values'][i] for i in positions)), left_factor)
            right = times_binary(multiply(*(base['values'][i] for i in positions)), right_factor)
            right = tuple(plane << global_w for plane in right)
            if left != right and witness is None:
                witness = {'x': x, 'y': y, 'u': u, 'v': v,
                           'cross_multiplied_difference': [hex(plane) for plane in add(left, right)]}
            if require_equal:
                assert left == right
            checks += 1
    return {'paired_entries_checked': checks, 'nonzero_difference_witness': witness}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    patterns = positive_monomials()
    assert len(patterns) == 24
    cubics = positive_monomials(degree=3)
    assert len(cubics) == 10
    cases = []
    for pattern in patterns:
        for coefficient, ell in ((1, -2), (2, 0), (3, 3)):
            extra = (*pattern, [(coefficient, ell)], 0)
            seed = core_monomial_seed([extra], (0, 1), 2)
            cases.append({'modes': pattern, 'coefficient': coefficient, 'ell': ell,
                          'occupation_charge': occupation_charge(pattern),
                          'denominator_degree': seed['denominator_degree'],
                          'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                          'two_copy_type_checks': seed['two_copy_type_checks'],
                          'coefficient_sha256': seed['coefficient_sha256']})
    from check_full_alphabet_quartic_repair import build as full_build
    base = full_build((0, 1), 2, 2, check_types=False)
    for i, line in enumerate(base['repair_matrix']):
        for j, coefficient in enumerate(line):
            if coefficient:
                assert OCCUPATION[i % 4]+OCCUPATION[i//4] == OCCUPATION[j % 4]+OCCUPATION[j//4]
    negative = []
    for pattern in cubics:
        seed = core_monomial_seed([(*pattern, [(2, -1)], 0)], (0, 1), 2, full_alphabet_gamma=2)
        negative.append({'modes': pattern, 'occupation_charge': occupation_charge(pattern),
                         **typed_difference(base, seed, require_equal=True)})
    pairs = []
    for minus, plus in product([p for p in cubics if occupation_charge(p) == -1],
                               [p for p in cubics if occupation_charge(p) == 1]):
        extras = [(*minus, [(1, 0)], 0), (*plus, [(2, 0)], 0)]
        seed = core_monomial_seed(extras, (0, 1), 2, full_alphabet_gamma=2)
        pairs.append({'minus': minus, 'plus': plus, **typed_difference(base, seed),
                      'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                      'two_copy_type_checks': seed['two_copy_type_checks']})
    assert any(case['nonzero_difference_witness'] for case in pairs)
    higher = []
    for full in (None, 2):
        extras = [(0, 2, 0, None, [(2, -1)], 0),
                  (0, 0, None, 2, [(3, 1)], 0),
                  (None, 2, 0, 2, [(1, 0)], 1)]
        seed = core_monomial_seed(extras, (0, 1, 2), 4, full_alphabet_gamma=full)
        higher.append({key: value for key, value in seed.items()
                       if key not in ('values', 'repair_matrix', 'residual_matrix', 'c_order_histogram')})
    result = {'status': 'exact endpoint-vanishing monomials and full-alphabet single-sign control; no pressure claim',
              'positive_monomials': patterns, 'mixed_alphabet_cases': cases,
              'single_sign_full_alphabet_controls': negative,
              'opposite_charge_full_alphabet_pairs': pairs, 'higher_mode_cases': higher}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
