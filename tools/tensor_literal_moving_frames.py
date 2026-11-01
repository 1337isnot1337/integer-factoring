#!/usr/bin/env python3
"""Larger literal codimension-one frames with parity-preserving Laurent flows.

The maps are L=[I,a], P=[I+a b^T;b^T], so LP=I identically.
Natural virtual weights remain explicit. Both weighted cuts use literal
Pfaffians and the infinity calculation reverses every map exponent.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

import tensor_hot_compounds_f16 as field
import tensor_literal_weighted_maps_f16 as weighted


def add_term(vector, key, coefficient):
    if coefficient:
        vector[key] = vector.get(key, 0) ^ coefficient
        if not vector[key]:
            del vector[key]


def frame_code(height_cap, anchor, support, a_mask, b_mask, flow=1, a_scalar=1, b_scalar=1):
    modes = (height_cap+1)//2
    words = sorted((w for w in range(1 << modes) if field.height(w) <= height_cap),
                   key=lambda w: (field.height(w), w))
    if anchor not in words or len(set(support)) != len(support) or any(w not in words or w == anchor for w in support):
        raise ValueError('Frame support must consist of distinct retained ambient words.')
    if any(w.bit_count() % 2 != anchor.bit_count() % 2 for w in support):
        raise ValueError('This constructor preserves physical occupation parity.')
    if not 0 <= a_mask < 1 << len(support) or not 0 <= b_mask < 1 << len(support):
        raise ValueError('Frame mask exceeds its declared support.')
    if not 0 < a_scalar < 16 or not 0 < b_scalar < 16:
        raise ValueError('Expected nonzero F16 frame scalars.')
    labels = [w for w in words if w != anchor]
    anchor_height = field.height(anchor)
    a, b = {}, {}
    for k, word in enumerate(support):
        difference = anchor_height-field.height(word)
        assert difference % 2 == 0
        shift = flow*difference//2
        if a_mask >> k & 1:
            a[word] = (a_scalar, shift)
        if b_mask >> k & 1:
            b[word] = (b_scalar, -shift)
    p, left = [], []
    for word in labels:
        row, column = {(word, 0): 1}, {(word, 0): 1}
        if word in a:
            coefficient, shift = a[word]
            add_term(row, (anchor, shift), coefficient)
        if word in b:
            coefficient, shift = b[word]
            add_term(column, (anchor, shift), coefficient)
            for other, (other_coefficient, other_shift) in a.items():
                add_term(column, (other, shift+other_shift), field.mul(coefficient, other_coefficient))
        left.append(row)
        p.append(column)
    # Independent Laurent-polynomial multiplication verifies the full split,
    # not just the scalar endpoint. Sparse dictionaries include negative powers.
    for i, row in enumerate(left):
        for j, column in enumerate(p):
            product = {}
            for (u, e), x in row.items():
                for (v, f), y in column.items():
                    if u == v:
                        add_term(product, e+f, field.mul(x, y))
            assert product == ({0: 1} if i == j else {})
    parities = [w.bit_count() % 2 for w in labels]
    assert all({w.bit_count() % 2 for w, _ in vector} == {parities[i]}
               for vectors in (p, left) for i, vector in enumerate(vectors))
    return {'height_cap': height_cap, 'anchor_word': anchor, 'support_words': support,
            'a_mask': a_mask, 'b_mask': b_mask, 'height_flow_in_w': flow,
            'a_scalar': a_scalar, 'b_scalar': b_scalar, 'q': len(labels), 'modes': modes,
            'labels': labels, 'heights': [field.height(w) for w in labels],
            'parities': parities, 'embedding': p, 'retraction': left}


def serialized_code(code):
    return {key: ([[[w, e, c] for (w, e), c in sorted(vector.items())] for vector in value]
                  if key in ('embedding', 'retraction') else value)
            for key, value in code.items()}


def cut_blocks(code, cut):
    p, left, heights, parities, q, modes = (code[k] for k in
        ('embedding', 'retraction', 'heights', 'parities', 'q', 'modes'))
    out = []
    for parity in (0, 1):
        pairs = [(i, j) for i in range(q) for j in range(q) if parities[i] ^ parities[j] == parity]
        if cut == 'zero':
            rows = [weighted.polynomial_tensor_terms(left[i], left[j], modes, heights[i], heights[j])
                    for i, j in pairs]
            columns = [weighted.polynomial_tensor_terms(p[j], p[i], modes) for i, j in pairs]
        elif cut == 'infinity':
            rows = [weighted.polynomial_tensor_terms(left[i], p[j], modes, None, heights[j], -1)
                    for i, j in pairs]
            columns = [weighted.polynomial_tensor_terms(left[j], p[i], modes, None, heights[i], -1)
                       for i, j in pairs]
        else:
            raise ValueError('Unknown place.')
        out.append({'parity': parity, 'pairs': pairs, 'rows': rows, 'columns': columns})
    clearing = max(weighted.clearing_power(block['rows'], block['columns']) for block in out)
    return out, clearing


def pressure_fit(zero, infinity, heights):
    zero, infinity, heights = map(lambda values: np.asarray(values, dtype=float), (zero, infinity, heights))
    queries = 0

    def objective(log_scale):
        nonlocal queries
        queries += 1
        scale = float(np.exp(log_scale))
        target = float(logsumexp(-scale*zero) + logsumexp(-scale*infinity))
        equation = lambda rho: target - 4*rho*float(logsumexp(-scale*heights/rho))
        if equation(1e-8) <= 0:
            raise ArithmeticError('No positive pressure root in the declared frame calculation.')
        high = 2
        while equation(high) > 0:
            high *= 2
            if high > 2**20:
                raise ArithmeticError('Pressure root exceeded the numerical search budget.')
        return float(brentq(equation, 1e-8, high, xtol=1e-13))

    grid = np.linspace(np.log(1e-4), np.log(2), 81)
    values = [objective(x) for x in grid]
    candidates = [{'rho': values[i], 'scale_t': float(np.exp(grid[i])), 'source': 'grid'}
                  for i in (0, len(grid)-1)]
    for i in range(1, len(grid)-1):
        if values[i] <= min(values[i-1], values[i+1]):
            fit = minimize_scalar(objective, bounds=(grid[i-1], grid[i+1]), method='bounded',
                                  options={'xatol': 1e-10, 'maxiter': 80})
            candidates.append({'rho': float(fit.fun), 'scale_t': float(np.exp(fit.x)),
                               'source': 'refined local minimum', 'completed': bool(fit.success)})
    best = min(candidates, key=lambda x: x['rho'])
    return {**best, 'queries': queries, 'grid_points': len(grid), 'scale_range': [1e-4, 2],
            'local_candidates': candidates, 'global_minimum_certified': False}


def profile(code, order_cap, cache_budget, callback):
    cuts, started = [], time.monotonic()
    for cut in ('zero', 'infinity'):
        blocks, clearing = cut_blocks(code, cut)
        orders, receipts = [], []
        for block in blocks:
            callback({'cut': cut, 'parity': block['parity'], 'block_size': len(block['pairs']),
                      'clearing_power': clearing})
            before = time.monotonic()
            result = weighted.block(code['modes'], block['rows'], block['columns'],
                                    order_cap=order_cap, clearing=clearing, cache_budget=cache_budget)
            orders.extend(result['orders_t'])
            receipts.append({key: value for key, value in result.items()
                             if key not in ('orders_t', 'cleared_orders_t')})
            receipts[-1].update(parity=block['parity'], seconds=time.monotonic()-before)
        assert len(orders) == code['q']**2
        cuts.append({'cut': cut, 'orders': sorted(orders), 'clearing_power': clearing, 'blocks': receipts})
    zero, infinity = (cut['orders'] for cut in cuts)
    fits = {'upper_counting_bound_fit': pressure_fit(zero, infinity, code['heights']),
            'lower_counting_bound_fit': pressure_fit([x for x in zero if x < order_cap],
                                                     [x for x in infinity if x < order_cap], code['heights'])}
    return {'q': code['q'], 'heights': code['heights'], 'order_cap': order_cap,
            'zero_counts': sorted(Counter(zero).items()), 'infinity_counts': sorted(Counter(infinity).items()),
            'computed_cuts': ['zero', 'infinity'], 'orders_include_all_laurent_shift_corrections': True,
            'blocks': [{k: v for k, v in cut.items() if k != 'orders'} for cut in cuts],
            'seconds': time.monotonic()-started}, fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--height', type=int, default=14)
    parser.add_argument('--anchor', type=int, default=9)
    parser.add_argument('--support', type=int, nargs='+', default=[0, 3, 5, 6])
    parser.add_argument('--a-masks', type=int, nargs='+', default=[0, 1, 3, 7, 15])
    parser.add_argument('--b-masks', type=int, nargs='+', default=[0, 1, 3, 7, 15])
    parser.add_argument('--height-flow', type=int, default=1)
    parser.add_argument('--a-scalar', type=int, default=1)
    parser.add_argument('--b-scalar', type=int, default=1)
    parser.add_argument('--order-cap', type=int, default=48)
    parser.add_argument('--cache-budget', type=int, default=8_000_000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    names = [Path(__file__).name, 'tensor_literal_weighted_maps_f16.py',
             'tensor_literal_weighted_maps_f16.cpp', 'tensor_literal_hot_pfaffians_f16.cpp',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'parameters': {k: v for k, v in vars(args).items() if k != 'output'},
              'completed': False, 'in_progress': None, 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':'))+'\n')
        temporary.replace(args.output)

    for a_mask in args.a_masks:
        for b_mask in args.b_masks:
            code = frame_code(args.height, args.anchor, args.support, a_mask, b_mask,
                              args.height_flow, args.a_scalar, args.b_scalar)

            def callback(block):
                report['in_progress'] = {'a_mask': a_mask, 'b_mask': b_mask, **block}
                save()

            print(json.dumps({'starting': [a_mask, b_mask], 'q': code['q'],
                              'height_flow': args.height_flow}), flush=True)
            data, fits = profile(code, args.order_cap, args.cache_budget, callback)
            report['records'].append({'code': serialized_code(code), 'profile': data, 'pressure_fits': fits})
            report['in_progress'] = None
            save()
            print(json.dumps({'completed': [a_mask, b_mask], 'fits': fits, 'seconds': data['seconds']}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
