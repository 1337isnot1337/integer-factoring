#!/usr/bin/env python3
"""Constructive F16 kernel profiles with the binary subset determinant.

Each output contains exact orders below a stated cap and certified lower
bounds on the remaining orders. The two floating pressure fits bracket
the counting-polynomial uncertainty; they are search diagnostics.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
import tensor_hot_compounds_f16_fast as fast
from search_coherent_hot_retractions import pressure


def profile(code, left, kind, order_cap, callback):
    embedding, modes = code['embedding'], code['modes']
    heights, charges = field.validate_maps(embedding, left, modes)
    q, records = len(embedding), []
    started = time.monotonic()
    for cut in (('zero',) if kind == 'orthogonal' else ('zero', 'infinity')):
        orders, blocks = [], []
        for charge in sorted({a + b for a in charges for b in charges}):
            pairs = [(i, j) for i in range(q) for j in range(q) if charges[i] + charges[j] == charge]
            callback({'cut': cut, 'charge': charge, 'block_size': len(pairs)})
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else embedding[j], modes)
                    for i, j in pairs]
            columns = [field.tensor_vector(embedding[j] if cut == 'zero' else left[j], embedding[i], modes)
                       for i, j in pairs]
            before = time.monotonic()
            block = fast.compound_block(modes, rows, columns, order_cap=order_cap)
            orders.extend(block['orders_t'])
            blocks.append({'charge': charge, 'size': len(pairs),
                           'exact_orders_below_cap': block['exact_orders_below_cap'],
                           'capped_tail_count': block['capped_tail_count'],
                           'precision_attempts': block['precision_attempts'],
                           'seconds': time.monotonic() - before})
        assert len(orders) == q * q
        records.append({'cut': cut, 'orders': sorted(orders), 'blocks': blocks})
    zero = records[0]['orders']
    infinity = zero[:] if kind == 'orthogonal' else records[1]['orders']
    data = {'q': q, 'modes': modes, 'heights': heights, 'charges': charges,
            'zero_orders_t': zero, 'infinity_orders_t': infinity, 'order_cap': order_cap,
            'computed_cuts': [x['cut'] for x in records],
            'infinity_profile_source': ('orthogonal identity; not separately computed'
                                       if kind == 'orthogonal' else 'directly computed'),
            'blocks': [{'cut': x['cut'], 'blocks': x['blocks']} for x in records],
            'seconds': time.monotonic() - started}
    lower = {**data, 'zero_orders_t': [x for x in zero if x < order_cap],
             'infinity_orders_t': [x for x in infinity if x < order_cap]}
    fits = {'upper_counting_bound_fit': pressure(data), 'lower_counting_bound_fit': pressure(lower)}
    data['zero_counts'] = sorted(Counter(data.pop('zero_orders_t')).items())
    data['infinity_counts'] = sorted(Counter(data.pop('infinity_orders_t')).items())
    return data, fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=('three-pair', 'four-pair', 'cubic'), required=True)
    parser.add_argument('--height', type=int, default=27)
    parser.add_argument('--coefficients', type=int, nargs='+')
    parser.add_argument('--retraction', choices=('coordinate', 'orthogonal'), default='orthogonal')
    parser.add_argument('--order-cap', type=int, default=64)
    parser.add_argument('--top-drop', type=int, default=0)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    count = 4 if args.family == 'four-pair' else 3
    coefficients = args.coefficients or [1 << i for i in range(count)]
    assert len(coefficients) == count and all(0 < x < 16 for x in coefficients)
    if args.family == 'cubic':
        operator = dict(zip((35, 21, 14), coefficients))
    else:
        operator = {(1 << i) | (1 << (2 * count - 1 - i)): coefficients[i] for i in range(count)}
    code = field.differential_kernel([operator], 16, args.height)
    if args.top_drop:
        assert 0 < args.top_drop < code['q']
        dropped = code['pivot_words'][-args.top_drop:]
        assert all(h == args.height for h in code['heights'][-args.top_drop:])
        for name in ('embedding', 'retraction', 'pivot_words', 'heights', 'charges'):
            code[name] = code[name][:-args.top_drop]
        code['q'] -= args.top_drop
        code['top_basis_words_removed'] = dropped
        field.validate_maps(code['embedding'], code['retraction'], code['modes'])
    left = code['retraction'] if args.retraction == 'coordinate' else field.orthogonal_retraction(code)
    assert left is not None, 'The selected kernel cell is orthogonally degenerate.'
    charge_counts = Counter(code['charges'])
    block_sizes = {s: sum(x * y for a, x in charge_counts.items() for b, y in charge_counts.items() if a + b == s)
                   for s in {a + b for a in charge_counts for b in charge_counts}}
    assert max(block_sizes.values()) <= 4608, block_sizes
    names = [Path(__file__).name, 'tensor_hot_compounds_f16_fast.py',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'search_coherent_hot_retractions.py']
    report = {'scope': __doc__, 'family': args.family, 'coefficients': coefficients,
              'height_cap': args.height, 'order_cap': args.order_cap,
              'retraction_kind': args.retraction, 'code': code, 'retraction': left,
              'completed': False, 'in_progress': None,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic() - started
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    def callback(block):
        report['in_progress'] = block
        save()

    print(json.dumps({'starting': args.family, 'coefficients': coefficients, 'q': code['q'],
                      'height_cap': args.height, 'block_sizes': block_sizes}), flush=True)
    data, fits = profile(code, left, args.retraction, args.order_cap, callback)
    report.update(profile=data, pressure_fits=fits, completed=True, in_progress=None)
    save()
    print(json.dumps({'completed': True, 'q': code['q'], 'fits': fits, 'seconds': data['seconds']}), flush=True)


if __name__ == '__main__':
    main()
