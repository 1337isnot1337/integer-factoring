#!/usr/bin/env python3
"""Replay subalphabet cofactor controls and rational convex dual certificates."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_cell_minor_caps import (PRIME, Dyadic, all_caps, determinant_mod, literal_case,
    log_at_least_one, monomial_crossing, partition, point, poly_eval, polynomial_minor,
    random_crossing, smith_orders, tight_crossing)
from check_graded_adjugate_caps import characteristic_two_controls
from graded_adjugate_caps import inequalities as whole_alphabet_caps
from subalphabet_adjugate_caps import inequalities, linear_caps


def profile(matrix):
    bound = sum(max(len(p)-1 for p in row) for row in matrix)
    return smith_orders(matrix, bound+1)


def literal_controls():
    rng, records, checked = Random(0x01A195C2AB), [], 0
    for q, count in ((2, 6), (3, 8), (4, 4)):
        for trial in range(count):
            heights = tuple([[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2))
            family = trial % 3
            tensor = (monomial_crossing(rng, q, 2) if family == 0 else
                      tight_crossing(q, heights) if family == 1 else random_crossing(rng, q, 1, .5))
            rec = literal_case(rng, q, tensor, heights, 0)
            slacks = []
            for pair in rec['profiles'].values():
                r, s = pair['r'], pair['s']
                caps = inequalities(rec['deficit'], r, s,
                                    rec['paired_energies_X'], rec['paired_energies_Y'])
                assert all(cap['slack'] >= 0 for cap in caps)
                # Check conversion to nonnegative absolute-order coefficients.
                mu = [F(v)+F(17, 3) for v in r]
                nu = [F(v)-rec['deficit']-F(17, 3) for v in s]
                for cap, actual in zip(linear_caps(rec['deficit'], rec['paired_energies_X'],
                                                  rec['paired_energies_Y']), caps[::2]):
                    absolute = sum(a*v for a, v in zip(cap['absolute_first'], mu))
                    absolute += sum(a*v for a, v in zip(cap['absolute_second'], nu))
                    assert cap['physical_rhs']-absolute == actual['slack']
                checked += len(caps)
                slacks.append([str(cap['slack']) for cap in caps])
            records.append({'q': q, 'trial': trial, 'family': family, 'heights': heights,
                            'tensor': [[list(a), sorted(p.items())] for a, p in sorted(tensor.items())],
                            'normalization': rec, 'slacks': slacks})
    binary = characteristic_two_controls()
    for rec in binary:
        first, second = rec['complete_profiles']
        r, s = [v-first[0] for v in first], [v-second[0] for v in second]
        caps = inequalities(rec['deficit'], r, s, *rec['energies'])
        assert all(cap['slack'] >= 0 for cap in caps)
        rec['subalphabet_minimum_slack'] = str(min(cap['slack'] for cap in caps))
        checked += len(caps)
    return records, binary, checked


def interlacing_controls():
    from check_cell_minor_caps import series_multiply
    rng, records, comparisons = Random(0x5B1A), [], 0
    for trial in range(32):
        n = 3+trial % 4
        orders = sorted(rng.randrange(7) for _ in range(n))
        matrix = [[[0]*orders[i]+[1] if i == j else [0]
                   for j in range(n)] for i in range(n)]
        for _ in range(2*n):
            i, j = rng.sample(range(n), 2)
            scalar = [rng.randrange(PRIME), rng.randrange(PRIME)]
            if rng.randrange(2):
                for k in range(n):
                    term = series_multiply(scalar, matrix[j][k], len(matrix[j][k])+1)
                    width = max(len(term), len(matrix[i][k]))
                    matrix[i][k] = [((matrix[i][k][a] if a < len(matrix[i][k]) else 0)
                                     +(term[a] if a < len(term) else 0)) % PRIME for a in range(width)]
            else:
                for k in range(n):
                    term = series_multiply(scalar, matrix[k][j], len(matrix[k][j])+1)
                    width = max(len(term), len(matrix[k][i]))
                    matrix[k][i] = [((matrix[k][i][a] if a < len(matrix[k][i]) else 0)
                                     +(term[a] if a < len(term) else 0)) % PRIME for a in range(width)]
        assert profile(matrix) == orders
        size = n-1-trial % 2
        for _ in range(100):
            rows, columns = sorted(rng.sample(range(n), size)), sorted(rng.sample(range(n), size))
            sub = [[matrix[i][j] for j in columns] for i in rows]
            if determinant_mod([[poly_eval(p, 1) for p in row] for row in sub]):
                break
        else:
            raise AssertionError('No nonsingular submatrix selected for this control.')
        observed = profile(sub)
        assert all(a >= b for a, b in zip(observed, orders))
        comparisons += size
        records.append({'matrix': matrix, 'orders': orders, 'rows': rows,
                        'columns': columns, 'submatrix_orders': observed})
    return records, comparisons


def rectangular_channel_controls():
    rng, records, minor_count = Random(0x2EC7), [], 0
    for p, b, k in ((2, 3, 1), (2, 3, 2), (2, 4, 1), (1, 5, 2), (2, 5, 1)):
        n, h = p*b, min(p, b)**2
        gap = n-1-h
        a = [[rng.randrange(PRIME) for _ in range(b)] for _ in range(p)]
        bb = [[rng.randrange(PRIME) for _ in range(p)] for _ in range(b)]
        def entry(x, y, u, v):
            leading = a[x][u]*bb[y][v] % PRIME
            result = [leading, -leading % PRIME]+[0]*(k-1)
            result[k] = (result[k]+int(x == v and y == u)) % PRIME
            return result
        pairs = list(product(range(p), range(b)))
        row = [[entry(x, y, u, v) for v, u in pairs] for x, y in pairs]
        cross = [[entry(x, y, u, v) for y, v in product(range(b), range(p))]
                 for x, u in pairs]
        crossed = profile(cross)
        assert crossed[0] == 0 and crossed[1] == k
        mapping = {i: {j: row[i][j] for j in range(n)} for i in range(n)}
        selected = list(product(range(n), repeat=2)) if n <= 8 else [(i, j) for i, j in ((0,0),(1,3),(2,7),(5,2),(8,9),(9,4))]
        adj = [[[0] for _ in range(n)] for _ in range(n)]
        controls = []
        for i, j in selected:
            polynomial = polynomial_minor(mapping, [a for a in range(n) if a != j],
                                          [a for a in range(n) if a != i])
            if (i+j) % 2:
                polynomial = [-v % PRIME for v in polynomial]
            order = next((e for e, value in enumerate(polynomial) if value), None)
            assert order is None or order >= gap*k
            adj[i][j] = polynomial or [0]
            controls.append({'row': i, 'column': j, 'polynomial': polynomial, 'order': order})
        eta = None
        if n <= 8:
            eta = profile([[adj[v*b+u][x*b+y] for u, x in product(range(b), range(p))]
                           for v, y in pairs])
            assert all(v >= gap*k for v in eta)
        minor_count += len(selected)
        records.append({'alphabet': [p, b], 'jet_order': k, 'channel_rank_bound': h,
                        'cofactor_gap': gap, 'A': a, 'B': bb, 'crossed_profile': crossed,
                        'cofactor_controls': controls, 'adjugate_crossed_profile': eta})
    return records, minor_count


def dual_controls(root):
    packet = json.loads((root/'evidence/subalphabet-adjugate-duals.json').read_text())
    reports, arithmetic = [], Dyadic(90)
    for rec in packet:
        d, energy = F(rec['deficit']), list(map(F, rec['energies']))
        n, rows, rhs = len(energy)**2, {}, {}
        for k, (bound, _) in all_caps(d, energy, energy).items():
            rows[('prefix', k)] = [F(i < k) for i in range(n)]
            rhs[('prefix', k)] = bound
        for i in range(n-1):
            values = [F(0)]*n
            values[i], values[i+1] = 1, -1
            rows[('sort', i)], rhs[('sort', i)] = values, F(0)
        for i, cap in enumerate(linear_caps(d, energy, energy)):
            rows[('subalphabet', i)] = [a+b for a, b in zip(cap['first'], cap['second'])]
            rhs[('subalphabet', i)] = cap['rhs']
        probabilities = list(map(F, rec['probability']))
        assert len(probabilities) == n and min(probabilities) > 0 and sum(probabilities) == 1
        combined, bound = [F(0)]*n, F(0)
        for label, weight in rec['multipliers']:
            key, value = tuple(label), F(weight)
            assert value >= 0 and key in rows
            combined = [a+value*b for a, b in zip(combined, rows[key])]
            bound += value*rhs[key]
        assert all(combined[i] >= 2*probabilities[i] for i in range(1, n))
        assert bound == F(rec['dual_bound'])
        entropy = arithmetic.zero
        for p in probabilities:
            log = log_at_least_one(arithmetic, point(arithmetic, 1/p))
            entropy = arithmetic.add(entropy, arithmetic.multiply(point(arithmetic, 2*p), log))
        physical = arithmetic.multiply(arithmetic.root_two, arithmetic.multiply(
            point(arithmetic, 2), log_at_least_one(arithmetic, partition(arithmetic, energy, True))))
        lower = arithmetic.subtract(arithmetic.subtract(
            arithmetic.add(point(arithmetic, d), entropy), point(arithmetic, bound)), physical)
        assert lower[0] > arithmetic.rational(8, 100)[1]
        reports.append({'deficit': str(d), 'energies': list(map(str, energy)),
                        'lower_bound_interval': list(map(str, lower)), 'precision': 90,
                        'lower_bound_decimal': [float(F(v, arithmetic.scale)) for v in lower],
                        'nonnegative_dual_multipliers': len(rec['multipliers'])})
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    literal, binary, inequalities_checked = literal_controls()
    interlacing, comparisons = interlacing_controls()
    channels, minors = rectangular_channel_controls()
    target = (0,6,9,9,14,19,19,20,22,22,22,22,24,24,24,24,30,34,34,36,36,37,38,40,41)
    energy = (0,0,12,16,18)
    assert all(v['slack'] >= 0 for v in whole_alphabet_caps(6, target, target, energy, energy))
    violations = [v for v in inequalities(6, target, target, energy, energy) if v['slack'] < 0]
    assert any(v['alphabet'] == (2,2) and v['slack'] == -6 for v in violations)
    assert any(v['alphabet'] == (2,5) and v['size'] == 1 and v['slack'] == -20 for v in violations)
    names = ['tools/subalphabet_adjugate_caps.py', 'tools/check_subalphabet_adjugate_caps.py',
             'tools/check_cell_minor_caps.py', 'tools/check_graded_adjugate_caps.py',
             'tools/graded_adjugate_caps.py', 'tools/charged_gaussian_floor.py',
             'evidence/subalphabet-adjugate-duals.json']
    report = {'status': 'pass', 'literal_controls': literal, 'characteristic_two_controls': binary,
              'literal_subalphabet_inequalities': inequalities_checked,
              'individual_smith_interlacing_controls': interlacing,
              'individual_smith_interlacing_order_comparisons': comparisons,
              'rectangular_channel_controls': channels, 'exact_rectangular_cofactors': minors,
              'excluded_integer_target': {'orders': target, 'deficit': 6, 'energies': energy,
                  'violations': violations, 'whole_alphabet_graded_caps': 'pass'},
              'exact_fixed_parameter_pressure_lower_bounds': dual_controls(root),
              'source_sha256': {n: sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope': 'All-field necessary constraints have an ordinary proof in the note. '
              'Three fixed parameter points have rational dual and outward pressure certificates. '
              'No uniform energy-ray floor, realizing tensor, or improved factoring exponent.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True, default=str)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'literal_subalphabet_inequalities',
        'individual_smith_interlacing_order_comparisons', 'exact_rectangular_cofactors',
        'exact_fixed_parameter_pressure_lower_bounds')}))


if __name__ == '__main__':
    main()
