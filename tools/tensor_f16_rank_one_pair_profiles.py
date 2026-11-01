#!/usr/bin/env python3
"""Constructive F16 rank-one three-/four-pair profile searches.

Binary-independent coefficients make every subset Gram nondegenerate.
Each completed profile stores exact capped Smith multiplicities. Replacing
larger orders by the cap gives an upper counting-polynomial bound; deleting
that tail gives a lower bound. Floating pressure fits guide further search.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
from search_coherent_hot_retractions import pressure


def binary_independent(values):
    basis = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit in basis:
                value ^= basis[bit]
            else:
                basis[bit] = value
                break
        if not value:
            return False
    return True


def canonical_frobenius(values):
    orbit, current = [], tuple(values)
    for _ in range(4):
        orbit.append(current)
        current = tuple(field.mul(x, x) for x in current)
    return min(orbit)


def pair_profile(code, left, kind, cap, callback):
    p, modes = code['embedding'], code['modes']
    heights, charges = field.validate_maps(p, left, modes)
    q, cuts = len(p), []
    started = time.monotonic()
    for cut in (('zero',) if kind == 'orthogonal' else ('zero', 'infinity')):
        orders, blocks = [], []
        for charge in sorted({a + b for a in charges for b in charges}):
            pairs = [(i, j) for i in range(q) for j in range(q) if charges[i] + charges[j] == charge]
            callback({'cut': cut, 'charge': charge, 'block_size': len(pairs)})
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else p[j], modes) for i, j in pairs]
            columns = [field.tensor_vector(p[j] if cut == 'zero' else left[j], p[i], modes) for i, j in pairs]
            before = time.monotonic()
            block = field.compound_block(modes, rows, columns, order_cap=cap)
            orders.extend(block['orders_t'])
            blocks.append({'charge': charge, 'size': len(pairs),
                           'exact_orders_below_cap': block['exact_orders_below_cap'],
                           'capped_tail_count': block['capped_tail_count'],
                           'precision_attempts': block['precision_attempts'],
                           'seconds': time.monotonic() - before})
        assert len(orders) == q * q
        cuts.append({'cut': cut, 'orders_t': sorted(orders), 'blocks': blocks})
    zero = cuts[0]['orders_t']
    infinity = zero[:] if kind == 'orthogonal' else cuts[1]['orders_t']
    result = {'q': q, 'modes': modes, 'heights': heights, 'charges': charges,
              'zero_orders_t': zero, 'infinity_orders_t': infinity,
              'order_cap': cap, 'computed_cuts': [cut['cut'] for cut in cuts],
              'infinity_profile_source': ('orthogonal identity; not separately computed'
                                           if kind == 'orthogonal' else 'directly computed'),
              'blocks': cuts, 'seconds': time.monotonic() - started}
    lower = dict(result)
    lower['zero_orders_t'] = [order for order in zero if order < cap]
    lower['infinity_orders_t'] = [order for order in infinity if order < cap]
    fits = {'upper_counting_bound_fit': pressure(result), 'lower_counting_bound_fit': pressure(lower)}
    return result, fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', type=int, choices=(3, 4), required=True)
    parser.add_argument('--caps', type=int, nargs='+', default=[20])
    parser.add_argument('--order-cap', type=int, default=64)
    parser.add_argument('--max-candidates', type=int, default=6)
    parser.add_argument('--coefficients', type=int, nargs='+')
    parser.add_argument('--coordinate', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    n = args.pairs
    if args.coefficients:
        assert len(args.coefficients) == n and binary_independent(args.coefficients)
        candidates = [tuple(args.coefficients)]
    else:
        pool = []
        bases = [(1, 2, 4), (1, 2, 8), (1, 4, 8)] if n == 3 else [(1, 2, 4, 8), (1, 2, 4, 9)]
        for base in bases:
            for tail in permutations(base[1:]):
                values = canonical_frobenius((1,) + tail)
                if values not in pool:
                    pool.append(values)
        candidates = pool[:args.max_candidates]
    assert candidates and all(binary_independent(values) for values in candidates)
    pairs = [(i, 2 * n - 1 - i) for i in range(n)]
    names = [Path(__file__).name, 'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'search_coherent_hot_retractions.py']
    report = {'scope': __doc__, 'field': 'F4[z]/(z^2+z+omega)', 'pairs': pairs, 'modes': 16,
              'caps': args.caps, 'order_cap': args.order_cap, 'coefficient_candidates': candidates,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'records': [], 'completed': False, 'in_progress': None}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic() - started
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    for values in candidates:
        coefficient_certificate = field.subset_gram_certificate([list(values)])
        assert coefficient_certificate['all_subsets_nondegenerate']
        operator = {(1 << a) | (1 << b): values[i] for i, (a, b) in enumerate(pairs)}
        for height_cap in args.caps:
            code = field.differential_kernel([operator], 16, height_cap)
            orthogonal = field.orthogonal_retraction(code)
            assert orthogonal is not None
            variants = [('orthogonal', orthogonal)] + ([('coordinate', code['retraction'])] if args.coordinate else [])
            for kind, left in variants:
                trial = {'coefficients': values, 'height_cap': height_cap, 'retraction_kind': kind, 'q': code['q']}

                def callback(block):
                    report['in_progress'] = {**trial, **block}
                    save()

                print(json.dumps({'starting': trial}), flush=True)
                profile, fits = pair_profile(code, left, kind, args.order_cap, callback)
                packed = {key: value for key, value in profile.items()
                          if key not in ('zero_orders_t', 'infinity_orders_t', 'blocks')}
                packed['zero_counts'] = sorted(Counter(profile['zero_orders_t']).items())
                packed['infinity_counts'] = sorted(Counter(profile['infinity_orders_t']).items())
                packed['blocks'] = [{'cut': cut['cut'], 'blocks': cut['blocks']} for cut in profile['blocks']]
                report['records'].append({**trial, 'coefficient_subset_certificate': coefficient_certificate,
                                          'code': code, 'retraction': left, 'profile': packed, 'pressure_fits': fits})
                report['in_progress'] = None
                save()
                print(json.dumps({'completed': trial, 'fits': fits, 'seconds': profile['seconds']}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
