#!/usr/bin/env python3
"""Exact source, independent contraction, and whole q12 profile controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np

from mixed_graph_retractions import WORDS as SMALL_WORDS, ambient_source
from physical_split_crossings import flatten
from q12_reflection_graphs import WORDS, graph_split, natural, source, tensor
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply, scalar_multiply
from windowed_f4_smith import pack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--full-independent', action='store_true',
                        help='Also recompute all 288 base orders with untruncated Python arithmetic.')
    args = parser.parse_args()
    begun = time.monotonic()
    ambient, small = source(), ambient_source()
    indices = [WORDS.index(w) for w in SMALL_WORDS]
    for index, value in enumerate(small['values']):
        target = sum(indices[index//6**leg % 6] * 16**leg for leg in range(4))
        assert ambient['values'][target] == (value[0] << 8, value[1] << 8)
    rng = np.random.default_rng(121620261101)
    coefficients = [0]*192
    for position in (0, 11, 17, 48, 57, 92, 97, 132, 145, 179):
        coefficients[position] = int(rng.integers(1, 4))
    examples = [(0,)*192, tuple(coefficients)]
    direct, independent_orders, records, derivatives = 0, 0, [], []
    for example, c in enumerate(examples):
        seed = tensor(c)
        lx, px = graph_split(c[:96])
        ly, py = graph_split(c[96:])
        maps = (lx, ly, py, px)
        for _ in range(96):
            coordinates = tuple(int(v) for v in rng.integers(0, 12, 4))
            terms = [[(j, value) for j, value in enumerate(matrix[i]) if value]
                     for matrix, i in zip(maps, coordinates)]
            value = 0
            for term in product(*terms):
                jj, cc = zip(*term)
                scalar = 1
                for a in cc:
                    scalar = polynomial_multiply(scalar, a)
                index = sum(j*16**k for k, j in enumerate(jj))
                value ^= scalar_multiply(pack(ambient['values'][index]), scalar)
            index = sum(j*12**k for k, j in enumerate(coordinates))
            assert value == pack(seed['values'][index])
            direct += 1
        conjugate = tensor(tuple((0, 1, 3, 2)[v] for v in c))
        assert conjugate['values'] == [(low ^ high, high) for low, high in seed['values']]
        model = ScaledAsymmetricPressure(seed)
        p = natural().astype(int)
        p[:-1] += rng.integers(-1, 2, len(p)-1)
        native = model.base.oracle.prefixes(p, 1)['orders']
        assert len(native) == 288
        point = np.r_[natural()[:-1]+rng.normal(0, .06, 88), math.log(.07), .3]
        _, gradient, _ = model.evaluate(point)
        errors = []
        for coordinate in range(90):
            delta = np.zeros(90)
            delta[coordinate] = 1e-6
            finite = (model.evaluate(point+delta, False)[0]-model.evaluate(point-delta, False)[0])/2e-6
            error = abs(float(finite-gradient[coordinate]))
            assert error < 2e-5, (example, coordinate, error)
            errors.append(error)
        derivatives.append({'example': example, 'point': point.tolist(), 'errors': errors})
        if args.full_independent and example == 0:
            print('Starting the 288-order independent untruncated base check.', flush=True)
            exact = []
            for side, cut in enumerate((1, 0)):
                matrix, rows, columns = flatten(seed, cut)
                h = [[0]+[int(v) for v in p[side*44+leg*11:side*44+(leg+1)*11]] for leg in range(4)]
                rw = [h[0][a]+h[2 if cut else 1][b] for a, b in rows]
                cw = [h[1 if cut else 2][a]+h[3][b] for a, b in columns]
                oracle = WeightedSmithOracle([[pack(v) for v in row] for row in matrix], infinity=bool(cut))
                profile, _ = oracle.profile(rw, cw)
                exact += list(profile)
                print(f'Completed independent cut {cut}.', flush=True)
            assert native == exact
            independent_orders += len(exact)
        # Independently compare a complete three-letter principal restriction of each graph.
        alphabet = (0, 2, 7)
        restricted = [seed['values'][alphabet[x]+12*(alphabet[y]+12*(alphabet[u]+12*alphabet[v]))]
                      for v, u, y, x in product(range(3), repeat=4)]
        tiny = {'values': restricted, 'q_x': 3, 'q_u': 3}
        tiny_model = ScaledAsymmetricPressure(tiny)
        tp = np.r_[np.asarray([p[side*44+leg*11+a-1] for side in range(2)
                              for leg in range(4) for a in alphabet[1:]]), 1]
        tiny_native = tiny_model.base.oracle.prefixes(tp, 1)['orders']
        exact = []
        for side, cut in enumerate((1, 0)):
            matrix, rows, columns = flatten(tiny, cut)
            h = [[0]+tp[side*8+2*leg:side*8+2*leg+2].astype(int).tolist() for leg in range(4)]
            rw = [h[0][a]+h[2 if cut else 1][b] for a, b in rows]
            cw = [h[1 if cut else 2][a]+h[3][b] for a, b in columns]
            profile, _ = WeightedSmithOracle([[pack(v) for v in row] for row in matrix], infinity=bool(cut)).profile(rw, cw)
            exact += list(profile)
        assert exact == tiny_native
        independent_orders += len(exact)
        records.append({'coefficients': c, 'degree': seed['degree'], 'parameters': p.tolist(),
                        'complete_native_orders': native, 'principal_alphabet': alphabet,
                        'independent_principal_orders': exact})
    names = ['q12_reflection_graphs.py', 'check_q12_reflection_graphs.py',
             'physical_split_crossings.py', 'reflection_dual_rail_gaussian.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'source_endpoints': 65536, 'source_degree': ambient['degree'],
              'small_source_entry_identities': 1296, 'compressed_endpoints': 2*12**4,
              'direct_packed_entries': direct, 'frobenius_entries': 2*12**4,
              'complete_native_orders': 576, 'independent_untruncated_orders': independent_orders,
              'finite_derivative_comparisons': 180, 'derivatives': derivatives,
              'full_288_independent_base_checked': args.full_independent,
              'records': records, 'seconds': time.monotonic()-begun,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Literal source and finite exact controls. Whole compressed nonbinary profiles use '
                       'the complete native window; independent full base replay is optional and explicitly flagged. '
                       'No favorable pressure, global optimum or stronger score.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'source_degree', 'direct_packed_entries',
                      'independent_untruncated_orders', 'seconds')}))


if __name__ == '__main__':
    main()
