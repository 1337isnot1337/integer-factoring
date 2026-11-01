#!/usr/bin/env python3
"""Finite-field rank-one physical retractions in a two-word height cell.

For every projective embedding line P this enumerates every retraction
L with LP=1, modulo scalar basis changes and field Frobenius. It retains
all other ambient height cells and contracts the literal Pfaffian tensor.
No deletion-closure shortcut is used. Both cuts are computed unless the
actual maps are orthogonal. Numerical pressure fits keep capped tails.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
import tensor_literal_hot_pfaffians_f16 as literal
from tensor_hot_boundary_retractions import validate_height_maps
from tensor_literal_boundary_profiles import pressure_fit


def square(value):
    return field.mul(value, value)


def orbit(case):
    infinity, a, b = case
    values, current = set(), case
    while current not in values:
        values.add(current)
        infinity, a, b = current
        current = (infinity, square(a), square(b))
    return sorted(values)


def all_cases(field_size):
    assert field_size in (2, 4, 16)
    cases = [(0, a, b) for a in range(field_size) for b in range(field_size)]
    cases.extend((1, 0, b) for b in range(field_size))
    representatives = [case for case in cases if case == min(orbit(case))]
    assert sum(len(orbit(case)) for case in representatives) == field_size * (field_size+1)
    return representatives


def cell_code(height_cap, boundary_height, case):
    modes = (height_cap + 1) // 2
    words = sorted((w for w in range(1 << modes) if field.height(w) <= height_cap),
                   key=lambda w: (field.height(w), w))
    cell = [w for w in words if field.height(w) == boundary_height]
    if len(cell) != 2:
        raise ValueError('The chosen height cell must contain exactly two words.')
    infinity, a, b = case
    if infinity:
        physical_p, physical_l = [0, 1], [b, 1]
    else:
        physical_p, physical_l = [1, a], [1 ^ field.mul(a, b), b]
    p_cell = {w: c for w, c in zip(cell, physical_p) if c}
    l_cell = {w: c for w, c in zip(cell, physical_l) if c}
    assert field.dot(p_cell, l_cell) == 1
    labels = [w for w in words if w != cell[1]]
    p = [p_cell if w == cell[0] else {w: 1} for w in labels]
    left = [l_cell if w == cell[0] else {w: 1} for w in labels]
    heights, diagnostics = validate_height_maps(p, left, modes)
    gram = field.dot(p_cell, p_cell)
    is_orthogonal = bool(gram) and all(field.mul(gram, physical_l[i]) == physical_p[i] for i in (0, 1))
    return {'height_cap': height_cap, 'boundary_height': boundary_height,
            'boundary_words': cell, 'case': list(case), 'frobenius_orbit': orbit(case),
            'physical_cell_embedding': physical_p, 'physical_cell_retraction': physical_l,
            'q': len(p), 'modes': modes, 'embedding': p, 'retraction': left,
            'heights': heights, 'orthogonal': is_orthogonal,
            'deletion_diagnostics': diagnostics}


def profile(code, order_cap, cache_budget, callback):
    p, left, heights, q, modes = (code[k] for k in ('embedding', 'retraction', 'heights', 'q', 'modes'))
    cuts = []
    started = time.monotonic()
    for cut in (('zero',) if code['orthogonal'] else ('zero', 'infinity')):
        orders, blocks = [], []
        for parity in (0, 1):
            pairs = [(i, j) for i in range(q) for j in range(q)
                     if (heights[i] + heights[j]) % 2 == parity]
            callback({'cut': cut, 'parity': parity, 'block_size': len(pairs)})
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else p[j], modes)
                    for i, j in pairs]
            columns = [field.tensor_vector(p[j] if cut == 'zero' else left[j], p[i], modes)
                       for i, j in pairs]
            before = time.monotonic()
            result = literal.compound_block(modes, rows, columns, order_cap=order_cap,
                                            cache_budget=cache_budget)
            orders.extend(result['orders_t'])
            blocks.append({key: value for key, value in result.items() if key != 'orders_t'})
            blocks[-1].update(parity=parity, seconds=time.monotonic()-before)
        assert len(orders) == q*q
        cuts.append({'cut': cut, 'orders': sorted(orders), 'blocks': blocks})
    zero = cuts[0]['orders']
    infinity = zero[:] if code['orthogonal'] else cuts[1]['orders']
    data = {'q': q, 'modes': modes, 'heights': heights, 'zero_orders_t': zero,
            'infinity_orders_t': infinity, 'order_cap': order_cap,
            'computed_cuts': [x['cut'] for x in cuts],
            'infinity_profile_source': ('literal orthogonal role-swap identity; not separately computed'
                                       if code['orthogonal'] else 'directly computed'),
            'blocks': [{'cut': x['cut'], 'blocks': x['blocks']} for x in cuts],
            'seconds': time.monotonic()-started}
    lower = {**data, 'zero_orders_t': [x for x in zero if x < order_cap],
             'infinity_orders_t': [x for x in infinity if x < order_cap]}
    fits = {'upper_counting_bound_fit': pressure_fit(data),
            'lower_counting_bound_fit': pressure_fit(lower)}
    data['zero_counts'] = sorted(Counter(data.pop('zero_orders_t')).items())
    data['infinity_counts'] = sorted(Counter(data.pop('infinity_orders_t')).items())
    return data, fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field-size', type=int, choices=(2, 4, 16), default=4)
    parser.add_argument('--height', type=int, default=14)
    parser.add_argument('--boundary-heights', type=int, nargs='+', default=[8, 9])
    parser.add_argument('--case', type=int, nargs=3)
    parser.add_argument('--order-cap', type=int, default=64)
    parser.add_argument('--cache-budget', type=int, default=8_000_000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = [tuple(args.case)] if args.case is not None else all_cases(args.field_size)
    if args.case is not None:
        assert args.case[0] in (0, 1) and 0 <= args.case[1] < args.field_size and 0 <= args.case[2] < args.field_size
        assert not args.case[0] or args.case[1] == 0
    names = [Path(__file__).name, 'tensor_literal_boundary_profiles.py',
             'tensor_literal_hot_pfaffians_f16.py', 'tensor_literal_hot_pfaffians_f16.cpp',
             'tensor_hot_boundary_retractions.py', 'tensor_hot_compounds_f16_fast.cpp',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'field_size': args.field_size, 'height_cap': args.height,
              'boundary_heights': args.boundary_heights, 'cases': cases, 'order_cap': args.order_cap,
              'complete_projective_catalogue': args.case is None, 'completed': False,
              'in_progress': None, 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':'))+'\n')
        temporary.replace(args.output)

    for boundary_height in args.boundary_heights:
        for case in cases:
            code = cell_code(args.height, boundary_height, case)

            def callback(block):
                report['in_progress'] = {'boundary_height': boundary_height, 'case': case, **block}
                save()

            print(json.dumps({'starting_boundary_height': boundary_height, 'case': case, 'q': code['q'],
                              'orthogonal': code['orthogonal']}), flush=True)
            data, fits = profile(code, args.order_cap, args.cache_budget, callback)
            report['records'].append({'code': code, 'profile': data, 'pressure_fits': fits})
            report['in_progress'] = None
            save()
            print(json.dumps({'completed_boundary_height': boundary_height, 'case': case,
                              'fits': fits, 'seconds': data['seconds']}), flush=True)
    assert len(report['records']) == len(cases)*len(args.boundary_heights)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
