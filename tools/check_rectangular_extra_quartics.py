#!/usr/bin/env python3
"""Independent endpoint-vanishing extreme-mode quartic perturbations.

Each extra is (x_mode,y_mode,u_mode,v_mode,Laurent_terms,c_pole_order).
Strictly positive ramified energy preserves the partial endpoint. This
constructor uses original Gaussian minors, independently of active CARs.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from check_rectangular_f4_repair import (
    ZERO, add, scale, times_binary, binary_mul, binary_power, c_order,
    fmul, one_particle, build,
)
from check_rectangular_refined_repairs import height, _verified_repair


PATTERNS = ((0, 1, 0, 0), (0, 0, 0, 1), (1, 1, 0, 1), (0, 1, 1, 1))


def validate_extra_quartics(top, extras):
    """Normalize extras and certify strictly positive ramified energy."""
    if top <= 0:
        raise ValueError('The top mode must be positive.')
    result = []
    for extra in extras:
        if len(extra) != 6:
            raise ValueError('An extra quartic needs four modes, Laurent terms, and a c-pole order.')
        a, b, c, d, supplied, pole = extra
        if any(mode not in (0, top) for mode in (a, b, c, d)):
            raise ValueError('This bounded-core family uses only the two extreme modes.')
        if not isinstance(pole, int) or pole < 0:
            raise ValueError('A c-pole order must be a nonnegative integer.')
        coefficients = {}
        for coefficient, exponent in supplied:
            if coefficient not in (0, 1, 2, 3) or not isinstance(exponent, int):
                raise ValueError('A Laurent term needs an F4 coefficient and an integer exponent.')
            coefficients[exponent] = coefficients.get(exponent, 0) ^ coefficient
        terms = [(coefficient, exponent) for exponent, coefficient in sorted(coefficients.items()) if coefficient]
        if not terms:
            continue
        offset = -min(0, min(exponent for _, exponent in terms))
        polynomial = tuple(sum(((coefficient >> p)&1) << (exponent+offset)
                               for coefficient, exponent in terms) for p in (0, 1))
        excess = b+d-a-c+(c_order(polynomial)-pole)*top
        if excess <= 0:
            raise ValueError('The extra quartic must vanish after the ramified energy gauge.')
        result.append((a, b, c, d, terms, pole))
    return result


def extra_quartic_seed(extras, indices=(0, 1, 2), ambient=4,
                       alpha=1, beta=2, check_types=True):
    """Add validated quartics to the original Gaussian before the repair."""
    indices = tuple(indices)
    extras = validate_extra_quartics(indices[-1], extras)
    seed = build(indices, ambient, alpha, beta, check_types=False,
                 independent_exterior=len(indices) <= 3)
    if not extras:
        return seed
    n, top, bulk = len(indices), 1 << (len(indices)-1), seed['bulk_size']
    A, B = seed['q_x'], seed['q_u']
    global_w = max(0, -min(exponent for extra in extras for _, exponent in extra[4]))
    denominator = max(seed['denominator_degree'], 2+max(extra[5] for extra in extras))
    base_factor = binary_power(3, denominator-seed['denominator_degree']) << global_w
    values = [times_binary(value, base_factor) for value in seed['values']]
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
        xa, yb, uc, vd = (1 << indices.index(mode) for mode in (a, b, c, d))
        polynomial = tuple(sum(((coefficient >> p)&1) << (exponent+global_w)
                               for coefficient, exponent in terms) for p in (0, 1))
        raw = [ZERO]*len(values)
        for x, y, u, v in product(range(A), range(B), range(B), range(A)):
            xx, vv = xmap[x], xmap[v]
            if not xx & xa or not vv & vd:
                continue
            answer = ZERO
            for (yy, cy), (uu, cu) in product(ymap[y], umap[u]):
                if yy & yb and uu & uc:
                    determinant = minor((vv ^ vd)|((yy ^ yb) << n),
                                        (uu ^ uc)|((xx ^ xa) << n))
                    answer = add(answer, scale(times_binary(polynomial, determinant), fmul(cy, cu)))
            raw[position(x, y, u, v)] = answer
        gates = [[(source, coefficient, binary_power(3, denominator-pole+height(source)-height(target)))
                  for source, coefficient in enumerate(row) if coefficient]
                 for target, row in enumerate(correction)]
        for x, y, u, v in product(range(A), range(B), range(B), range(A)):
            target = position(x, y, u, v)
            for source, coefficient, factor in gates[v//bulk+4*(y//bulk)]:
                old = raw[position(x, (source//4)*bulk+y % bulk, u, (source % 4)*bulk+v % bulk)]
                values[target] = add(values[target], scale(times_binary(old, factor), coefficient))
    provisional = dict(seed, denominator_degree=denominator)
    result = _verified_repair(provisional, values, seed['ramified_endpoint_matrix'],
                              [(height(i), 0) for i in range(12)], check_types)
    result.update({'repair_kind': 'mixed_partial_extra_quartics', 'extra_quartics': extras,
                   'global_w_power': global_w, 'uncancelled_denominator_degree': denominator,
                   'cancelled_c_order': denominator-result['denominator_degree']})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    for label, pattern in enumerate(PATTERNS):
        for coefficient in (1, 2, 3):
            for ell in (-2, 0, 3):
                extra = (*pattern, [(coefficient, ell)], 0)
                seed = extra_quartic_seed([extra], (0, 1), 2)
                cases.append({'pattern': label, 'modes': list(pattern), 'coefficient': coefficient,
                              'ell': ell, 'denominator_degree': seed['denominator_degree'],
                              'global_w_power': seed['global_w_power'],
                              'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                              'two_copy_type_checks': seed['two_copy_type_checks'],
                              'coefficient_sha256': seed['coefficient_sha256']})
    mixed = []
    for indices, ambient in (((0, 1, 2), 4), ((0, 1, 2, 3), 4)):
        top = indices[-1]
        extras = [(*(mode*top for mode in pattern), [(1+i % 3, i-2)], 0)
                  for i, pattern in enumerate(PATTERNS)]
        seed = extra_quartic_seed(extras, indices, ambient)
        mixed.append({key: value for key, value in seed.items()
                      if key not in ('values', 'residual_matrix', 'repair_matrix', 'c_order_histogram')})
    try:
        validate_extra_quartics(3, [(0, 3, 0, 0, [(1, 0)], 1)])
    except ValueError:
        rejected = True
    else:
        raise AssertionError('A nonvanishing marginal endpoint was accepted.')
    result = {'status': 'exact endpoint-vanishing extreme quartics; no pressure claim',
              'single_extra_cases': cases, 'mixed_extra_cases': mixed,
              'marginal_endpoint_rejected': rejected}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
