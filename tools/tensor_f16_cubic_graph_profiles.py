#!/usr/bin/env python3
"""Capped natural-pressure profiles of mixed-occupation F16 graph codes."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
import tensor_hot_compounds_f16_mixed as mixed
from tensor_f16_cubic_graph_codes import graph_code, orthogonal_retraction, validate_maps
from search_coherent_hot_retractions import pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coefficients', type=int, nargs=3, default=[2, 4, 8])
    parser.add_argument('--height', type=int, default=20)
    parser.add_argument('--order-cap', type=int, default=64)
    parser.add_argument('--top-drop', type=int, default=0)
    parser.add_argument('--retraction', choices=('orthogonal', 'coordinate'), default='orthogonal')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    code = graph_code(args.coefficients, args.height)
    if args.top_drop:
        assert 0 < args.top_drop < code['q']
        assert all(h == args.height for h in code['heights'][-args.top_drop:])
        code['top_basis_words_removed'] = code['pivot_words'][-args.top_drop:]
        for name in ('embedding', 'retraction', 'pivot_words', 'heights', 'parities'):
            code[name] = code[name][:-args.top_drop]
        code['q'] -= args.top_drop
    left = code['retraction'] if args.retraction == 'coordinate' else orthogonal_retraction(code)
    assert left is not None
    embedding, q, modes = code['embedding'], code['q'], code['modes']
    heights, parities = validate_maps(embedding, left, modes)
    counts = Counter(parities)
    sizes = [counts[0] ** 2 + counts[1] ** 2, 2 * counts[0] * counts[1]]
    assert max(sizes) <= 4608, sizes
    names = [Path(__file__).name, 'tensor_f16_cubic_graph_codes.py',
             'tensor_hot_compounds_f16_mixed.py', 'tensor_hot_compounds_f16_mixed.cpp',
             'tensor_hot_compounds_f16_fast.py', 'tensor_hot_compounds_f16_fast.cpp',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'search_coherent_hot_retractions.py']
    report = {'scope': __doc__, 'code': code, 'retraction': left,
              'retraction_kind': args.retraction, 'order_cap': args.order_cap,
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

    print(json.dumps({'starting': args.coefficients, 'height': args.height, 'q': q,
                      'parity_blocks': sizes, 'retraction_kind': args.retraction}), flush=True)
    cuts = []
    for cut in (('zero',) if args.retraction == 'orthogonal' else ('zero', 'infinity')):
        all_orders, blocks = [], []
        for parity in (0, 1):
            pairs = [(i, j) for i in range(q) for j in range(q)
                     if (parities[i] + parities[j]) % 2 == parity]
            report['in_progress'] = {'cut': cut, 'parity': parity, 'block_size': len(pairs)}
            save()
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else embedding[j], modes)
                    for i, j in pairs]
            columns = [field.tensor_vector(embedding[j] if cut == 'zero' else left[j], embedding[i], modes)
                       for i, j in pairs]
            before = time.monotonic()
            result = mixed.compound_block(modes, rows, columns, order_cap=args.order_cap)
            all_orders.extend(result['orders_t'])
            blocks.append({'parity': parity, 'size': len(pairs),
                           'exact_orders_below_cap': result['exact_orders_below_cap'],
                           'capped_tail_count': result['capped_tail_count'],
                           'precision_attempts': result['precision_attempts'],
                           'seconds': time.monotonic() - before})
        assert len(all_orders) == q * q
        cuts.append({'cut': cut, 'orders': sorted(all_orders), 'blocks': blocks})
    zero = cuts[0]['orders']
    infinity = zero[:] if args.retraction == 'orthogonal' else cuts[1]['orders']
    data = {'q': q, 'modes': modes, 'heights': heights, 'parities': parities,
            'zero_orders_t': zero, 'infinity_orders_t': infinity, 'order_cap': args.order_cap,
            'computed_cuts': [c['cut'] for c in cuts],
            'infinity_profile_source': ('orthogonal identity; not separately computed'
                                       if args.retraction == 'orthogonal' else 'directly computed'),
            'blocks': [{'cut': c['cut'], 'blocks': c['blocks']} for c in cuts],
            'seconds': time.monotonic() - started}
    lower = {**data, 'zero_orders_t': [x for x in zero if x < args.order_cap],
             'infinity_orders_t': [x for x in infinity if x < args.order_cap]}
    fits = {'upper_counting_bound_fit': pressure(data), 'lower_counting_bound_fit': pressure(lower)}
    data['zero_counts'] = sorted(Counter(data.pop('zero_orders_t')).items())
    data['infinity_counts'] = sorted(Counter(data.pop('infinity_orders_t')).items())
    report.update(completed=True, in_progress=None, profile=data, pressure_fits=fits)
    save()
    print(json.dumps({'completed': True, 'q': q, 'fits': fits, 'seconds': data['seconds']}), flush=True)


if __name__ == '__main__':
    main()
