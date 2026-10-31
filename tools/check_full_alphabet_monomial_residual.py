#!/usr/bin/env python3
"""Full-alphabet repairs with nonzero odd ramified residuals.

The explicit residual is computed in an eight-variable exterior algebra.
Every accepted correction is independently verified against that residual.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from check_full_alphabet_quartic_repair import CORE, HEIGHT, inverse
from check_rectangular_core_monomials import validate_extra_monomials, positive_monomials, occupation_charge
from check_rectangular_f4_repair import (
    ZERO, add, scale, times_binary, binary_mul, binary_power, divide_binary,
    c_order, evaluate_one, fmul, one_particle, matrix_multiply,
)
from check_rectangular_refined_repairs import _verified_repair


RANK_RESTORING_PAIRS = ((0, 2), (0, 7), (1, 3), (1, 6))


def height(index):
    return HEIGHT[index % 4]+HEIGHT[index//4]


def gaussian_residual():
    coefficients = {0: 1}
    for monomial in (65, 130, 20, 40, 96, 6):
        for mask, value in list(coefficients.items()):
            if not mask & monomial:
                coefficients[mask | monomial] = coefficients.get(mask | monomial, 0) ^ value
    return {mask: value for mask, value in coefficients.items() if value}


def critical_terms(top, specifications):
    normalized = validate_extra_monomials(top, specifications, allow_marginal=True)
    result = []
    for *modes, terms, pole in normalized:
        offset = -min(0, min(exponent for _, exponent in terms))
        polynomial = tuple(sum(((coefficient >> bit)&1) << (exponent+offset)
                               for coefficient, exponent in terms) for bit in (0, 1))
        order = c_order(polynomial)
        delta = sum(mode or 0 for mode in (modes[1], modes[3])) \
            -sum(mode or 0 for mode in (modes[0], modes[2]))
        if delta+(order-pole)*top == 0:
            coefficient = evaluate_one(divide_binary(polynomial, binary_power(3, order)))
            assert coefficient
            monomial = sum(1 << (2*leg+int(mode == top))
                           for leg, mode in enumerate(modes) if mode is not None)
            result.append((monomial, coefficient))
    return normalized, result


def residual(gamma, specifications, top=1):
    if gamma not in (0, 1, 2, 3):
        raise ValueError('The canonical quartic coefficient must lie in F4.')
    normalized, critical = critical_terms(top, specifications)
    gaussian = gaussian_residual()
    coefficients = dict(gaussian)
    for monomial, coefficient in [(153, gamma)]+critical:
        if not coefficient:
            continue
        for mask, value in gaussian.items():
            if not mask & monomial:
                coefficients[mask | monomial] = coefficients.get(mask | monomial, 0) ^ fmul(value, coefficient)
    matrix = [[coefficients.get(CORE[x] | (CORE[y] << 2) | (CORE[u] << 4) | (CORE[v] << 6), 0)
               for u in range(4) for x in range(4)] for y in range(4) for v in range(4)]
    return matrix, normalized, critical


def height_block_correction(matrix):
    """Canonical block-lower C with CF=I plus strictly height-lowering entries."""
    correction = [[0]*16 for _ in range(16)]
    blocks = []
    for level in range(3):
        selected = [i for i in range(16) if height(i) <= level]
        block = [[matrix[i][j] for j in selected] for i in selected]
        invert, rank, determinant = inverse(block)
        if invert is None:
            raise ValueError(f'The leading height block of order {len(selected)} is singular.')
        blocks.append({'order': len(selected), 'rank': rank, 'determinant': determinant})
        for position, row in enumerate(selected):
            if height(row) == level:
                for column, value in zip(selected, invert[position]):
                    correction[row][column] = value
    endpoint = matrix_multiply(correction, matrix)
    for i, line in enumerate(endpoint):
        assert line[i] == 1
        assert all(not value or i == j or height(i) < height(j) for j, value in enumerate(line))
    assert all(not value or height(j) <= height(i)
               for i, line in enumerate(correction) for j, value in enumerate(line))
    return correction, endpoint, blocks


def prepare_residual_repair(gamma, specifications, top=1, repair_kind='partial'):
    matrix, normalized, critical = residual(gamma, specifications, top)
    invert, rank, determinant = inverse(matrix)
    if invert is None:
        raise ValueError('The combined full-alphabet ramified residual is singular.')
    if repair_kind == 'partial':
        correction, endpoint, blocks = height_block_correction(matrix)
    elif repair_kind == 'inverse':
        correction, endpoint, blocks = invert, [[int(i == j) for j in range(16)] for i in range(16)], []
    else:
        raise ValueError('Use partial or inverse repair.')
    return {'residual_matrix': matrix, 'repair_matrix': correction,
            'ramified_endpoint_matrix': endpoint, 'residual_monomials': normalized,
            'critical_exterior_terms': critical, 'residual_rank': rank,
            'residual_determinant': determinant, 'leading_height_blocks': blocks}


def full_residual_monomial_seed(specifications, indices=(0, 1), ambient=2,
                                gamma=0, repair_kind='partial', check_types=True,
                                independent_exterior=False):
    indices = tuple(indices)
    if not indices or indices[0] != 0 or len(indices) < 2:
        raise ValueError('The mode set must contain zero and a positive top mode.')
    endpoint = prepare_residual_repair(gamma, specifications, indices[-1], repair_kind)
    extras = endpoint['residual_monomials']
    n, q, bulk = len(indices), 1 << len(indices), 1 << (len(indices)-2)
    top = 1 << (n-1)
    masks = [((letter % bulk) << 1) | (1, 0, top, top | 1)[letter//bulk] for letter in range(q)]
    gaussian = one_particle(indices, ambient)
    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i, result, remaining = bit.bit_length()-1, 0, columns
        while remaining:
            column = remaining & -remaining
            remaining ^= column
            j = column.bit_length()-1
            if gaussian[i][j]:
                result ^= binary_mul(gaussian[i][j], minor(rows ^ bit, columns ^ column))
        return result
    raw_terms = ([(0, indices[-1], 0, indices[-1], [(gamma, 0)], 2)] if gamma else [])+extras
    raw_denominator = max((term[-1] for term in raw_terms), default=0)
    global_w = max(0, -min((exponent for term in raw_terms for _, exponent in term[4]), default=0))
    gaussian_factor = binary_power(3, raw_denominator) << global_w
    prepared = []
    for *modes, terms, pole in raw_terms:
        removed = [0 if mode is None else 1 << indices.index(mode) for mode in modes]
        polynomial = tuple(sum(((coefficient >> bit)&1) << (exponent+global_w)
                               for coefficient, exponent in terms) for bit in (0, 1))
        polynomial = times_binary(polynomial, binary_power(3, raw_denominator-pole))
        prepared.append((removed, polynomial))
    def position(x, y, u, v):
        return x+q*(y+q*(u+q*v))
    raw = [ZERO]*(q**4)
    for x, y, u, v in product(range(q), repeat=4):
        xx, yy, uu, vv = (masks[letter] for letter in (x, y, u, v))
        answer = (binary_mul(gaussian_factor, minor(vv | (yy << n), uu | (xx << n))), 0)
        for (xa, yb, uc, vd), polynomial in prepared:
            if xx & xa == xa and yy & yb == yb and uu & uc == uc and vv & vd == vd:
                answer = add(answer, times_binary(polynomial, minor((vv ^ vd) | ((yy ^ yb) << n),
                                                                   (uu ^ uc) | ((xx ^ xa) << n))))
        raw[position(x, y, u, v)] = answer
    exterior_checks = 0
    if independent_exterior:
        support = {0: 1}
        for i, line in enumerate(gaussian):
            output = 3*n+i if i < n else i
            for j, coefficient in enumerate(line):
                if not coefficient:
                    continue
                input_index = 2*n+j if j < n else j-n
                monomial = (1 << output) | (1 << input_index)
                for mask, value in list(support.items()):
                    if not mask & monomial:
                        target = mask | monomial
                        support[target] = support.get(target, 0) ^ binary_mul(value, coefficient)
                        if not support[target]:
                            del support[target]
        exterior = {mask: (binary_mul(value, gaussian_factor), 0) for mask, value in support.items()}
        for removed, polynomial in prepared:
            monomial = sum(mask << (leg*n) for leg, mask in enumerate(removed))
            for mask, value in support.items():
                if not mask & monomial:
                    exterior[mask | monomial] = add(exterior.get(mask | monomial, ZERO), times_binary(polynomial, value))
        for x, y, u, v in product(range(q), repeat=4):
            mask = masks[x] | (masks[y] << n) | (masks[u] << (2*n)) | (masks[v] << (3*n))
            assert raw[position(x, y, u, v)] == exterior.get(mask, ZERO)
            exterior_checks += 1
    correction = endpoint['repair_matrix']
    gates = [[(source % 4, source//4, coefficient, binary_power(3, 2+height(source)-height(target)))
              for source, coefficient in enumerate(line) if coefficient]
             for target, line in enumerate(correction)]
    values = [ZERO]*len(raw)
    for x, y, u, v in product(range(q), repeat=4):
        answer = ZERO
        for source_v, source_y, coefficient, factor in gates[v//bulk+4*(y//bulk)]:
            old = raw[position(x, source_y*bulk+y % bulk, u, source_v*bulk+v % bulk)]
            answer = add(answer, scale(times_binary(old, factor), coefficient))
        values[position(x, y, u, v)] = answer
    heights = [sum(indices[j] for j in range(n) if mask >> j & 1) for mask in masks]
    denominator = raw_denominator+2
    provisional = {'field': 4, 'field_polynomial': 'omega^2+omega+1',
                   'indices': list(indices), 'ambient': ambient, 'gamma': gamma,
                   'dimensions': [q]*4, 'q_x': q, 'q_u': q, 'bulk_size': bulk,
                   'repair_kind': 'full_alphabet_monomial_residual_'+repair_kind,
                   'denominator_degree': denominator, 'strict_endpoint': False,
                   'typed_endpoint': 'independent matching full-letter types on X/V and Y/U',
                   'logical_heights_x_v': heights, 'logical_heights_u_y': heights,
                   'global_w_power': global_w, **endpoint}
    result = _verified_repair(provisional, values, endpoint['ramified_endpoint_matrix'],
                              [(height(i), 0) for i in range(16)], check_types)
    result['uncancelled_denominator_degree'] = denominator
    result['cancelled_c_order'] = denominator-result['denominator_degree']
    result['independent_exterior_checks'] = exterior_checks
    result['odd_degree_witness'] = None
    for index, value in enumerate(result['values']):
        if value == ZERO:
            continue
        rest, letters = index, []
        for _ in range(4):
            letters.append(rest % q)
            rest //= q
        if sum(masks[letter].bit_count() for letter in letters) % 2:
            result['odd_degree_witness'] = {'physical_letters_x_y_u_v': letters,
                                            'occupied_mode_masks': [masks[letter] for letter in letters],
                                            'nonzero_polynomial_numerator': [hex(plane) for plane in value]}
            break
    return result


def cubic_pair(top, minus, plus, coefficient=1, minus_power=0, plus_power=0):
    patterns = positive_monomials(top, degree=3)
    negative = [pattern for pattern in patterns if occupation_charge(pattern) == -1]
    positive = [pattern for pattern in patterns if occupation_charge(pattern) == 1]
    result = []
    for pattern, value, exponent in ((negative[minus], 1, minus_power), (positive[plus], coefficient, plus_power)):
        delta = sum(mode or 0 for mode in (pattern[1], pattern[3]))-sum(mode or 0 for mode in (pattern[0], pattern[2]))
        assert delta % top == 0
        result.append((*pattern, [(value, exponent)], delta//top))
    return result


def symbolic_rank_restoration():
    """Exact determinant identities in F2[a,b,gamma], via subset expansion."""
    gaussian = gaussian_residual()
    patterns = positive_monomials(degree=3)
    negative = [pattern for pattern in patterns if occupation_charge(pattern) == -1]
    positive = [pattern for pattern in patterns if occupation_charge(pattern) == 1]
    def determinant(matrix):
        rows = sorted(matrix, key=lambda line: sum(bool(value) for value in line))
        dp = {0: {0}}
        for line in rows:
            following = {}
            for mask, polynomial in dp.items():
                for j, coefficient in enumerate(line):
                    if mask >> j & 1 or not coefficient:
                        continue
                    target = following.setdefault(mask | (1 << j), set())
                    for left, right in product(polynomial, coefficient):
                        value = left+right
                        if value in target:
                            target.remove(value)
                        else:
                            target.add(value)
            dp = {mask: value for mask, value in following.items() if value}
        return dp.get((1 << len(matrix))-1, set())
    order = sorted(range(16), key=height)
    records = []
    for minus, plus in RANK_RESTORING_PAIRS:
        coefficients = {mask: {0} for mask in gaussian}
        for pattern, coefficient in [((0, 1, 0, 1), 1024), (negative[minus], 1), (positive[plus], 32)]:
            monomial = sum(1 << (2*leg+mode) for leg, mode in enumerate(pattern) if mode is not None)
            for mask in gaussian:
                if not mask & monomial:
                    coefficients.setdefault(mask | monomial, set()).symmetric_difference_update({coefficient})
        matrix = [[coefficients.get(CORE[x] | (CORE[y] << 2) | (CORE[u] << 4) | (CORE[v] << 6), set())
                   for u in range(4) for x in range(4)] for y in range(4) for v in range(4)]
        polynomials = []
        for size in (4, 12, 16):
            polynomial = determinant([[matrix[i][j] for j in order[:size]] for i in order[:size]])
            assert polynomial == ({0} if size < 16 else {0, 4+4*32, 4*1024})
            polynomials.append({'block_order': size,
                                'nonzero_monomial_exponents_a_b_gamma':
                                    [(term % 32, (term//32) % 32, term//1024) for term in sorted(polynomial)]})
        records.append({'minus': minus, 'plus': plus, 'determinants': polynomials})
    return {'ring': 'F2[a,b,gamma]', 'full_determinant': '(1+gamma+a*b)^4', 'pairs': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    symbolic = symbolic_rank_restoration()
    records = []
    for (minus, plus), coefficient, repair in product(RANK_RESTORING_PAIRS, (1, 2, 3), ('partial', 'inverse')):
        extras = cubic_pair(1, minus, plus, coefficient)
        seed = full_residual_monomial_seed(extras, gamma=1, repair_kind=repair, independent_exterior=True)
        assert seed['odd_degree_witness'] is not None
        records.append({key: value for key, value in seed.items()
                        if key not in ('values', 'residual_matrix', 'repair_matrix', 'ramified_endpoint_matrix', 'c_order_histogram')})
    for indices, ambient, minus, plus, coefficient in (((0, 1, 2), 4, 0, 2, 2),
                                                       ((0, 2, 5), 8, 1, 3, 3),
                                                       ((0, 1, 2, 3), 4, 0, 7, 1)):
        extras = cubic_pair(indices[-1], minus, plus, coefficient, minus_power=-1, plus_power=2)
        seed = full_residual_monomial_seed(extras, indices, ambient, gamma=1, independent_exterior=True)
        records.append({key: value for key, value in seed.items()
                        if key not in ('values', 'residual_matrix', 'repair_matrix', 'ramified_endpoint_matrix', 'c_order_histogram')})
    rejected = 0
    for coefficient in (1, 2, 3):
        try:
            prepare_residual_repair(1 ^ coefficient, cubic_pair(1, 0, 2, coefficient))
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('The singular parameter locus was accepted.')
    result = {'status': 'exact full-alphabet rank-restoring cubic residuals; no pressure claim',
              'symbolic_rank_restoration': symbolic, 'singular_cases_rejected': rejected, 'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
