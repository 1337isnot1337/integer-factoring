#!/usr/bin/env python3
"""Check graded cofactor identities, literal profiles and old/new relaxation points."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_cell_minor_caps import (PRIME, Dyadic, literal_case, monomial_crossing,
    polynomial_minor, pressure, random_crossing, relaxed_orders, tight_crossing)
from graded_adjugate_caps import inequalities, prefix_form


def characteristic_two_controls():
    from check_joint_energy_conjugations import cut_matrix
    from search_f4_recentered_crossings import five_letter_source
    from weighted_smith_oracle_f4 import WeightedSmithOracle
    values, _, original_energy, _ = five_letter_source()
    rows = []
    rng = Random(0xAD205)
    for trial in range(6):
        if trial < 2:
            h = [v//2 for v in (original_energy if trial == 0 else [0,18,12,0,16])]
            zero = [0]*5
            heights = [[h, h, zero, zero], [zero, zero, h, h]]
        else:
            heights = [[[rng.randint(-3, 3) for _ in range(5)] for _ in range(4)] for _ in range(2)]
        h0, hi = heights
        first, _ = WeightedSmithOracle(cut_matrix(values, False)).profile(
            [h0[0][i % 5]+h0[1][i//5] for i in range(25)],
            [h0[3][i % 5]+h0[2][i//5] for i in range(25)])
        second, _ = WeightedSmithOracle(cut_matrix(values, True), infinity=True).profile(
            [hi[0][i % 5]+hi[2][i//5] for i in range(25)],
            [hi[1][i % 5]+hi[3][i//5] for i in range(25)])
        ex = [h0[0][j]+hi[0][j]+h0[3][j]+hi[3][j] for j in range(5)]
        ey = [h0[1][j]+hi[1][j]+h0[2][j]+hi[2][j] for j in range(5)]
        d = min(ex)+min(ey)-first[0]-second[0]
        ex, ey = [v-min(ex) for v in ex], [v-min(ey) for v in ey]
        r, s = [v-first[0] for v in first], [v-second[0] for v in second]
        caps = inequalities(d, r, s, ex, ey)
        assert all(cap['slack'] >= 0 for cap in caps)
        rows.append({'trial': trial, 'heights': heights, 'deficit': d,
                     'energies': [ex, ey], 'complete_profiles': [first, second],
                     'minimum_new_cap_slack': str(min(cap['slack'] for cap in caps))})
    return rows


def cofactor_controls(tensor, heights, original):
    q, n = 2, 4
    h0, hi = heights
    def shift(bits, place):
        return sum((h0 if place == 0 else hi)[l][bits[l]] for l in range(4))
    minimum0 = min(e+shift(a, 0) for a, p in tensor.items() for e in p)
    minimumi = min(-e+shift(a, 1) for a, p in tensor.items() for e in p)
    H = [[h0[l][s]+hi[l][s] for s in range(q)] for l in range(4)]
    EX, EY = [[H[a][s]+H[b][s] for s in range(q)] for a, b in ((0, 3), (1, 2))]
    bx, by = min(range(q), key=EX.__getitem__), min(range(q), key=EY.__getitem__)
    c = [[H[l][s]-H[l][bx if l in (0, 3) else by] for s in range(q)] for l in range(4)]
    d = original['deficit']
    constant = (n-1)*d+q*(sum(original['paired_energies_X'])+sum(original['paired_energies_Y']))
    def matrix(place):
        out = [[[0] for _ in range(n)] for _ in range(n)]
        for a, polynomial in tensor.items():
            x, y, u, v = a
            values = {(e if place == 0 else -e)+shift(a, place)
                      -(minimum0 if place == 0 else minimumi): value for e, value in polynomial.items()}
            row = [0]*(max(values)+1)
            for e, value in values.items():
                row[e] = value
            out[x+q*y][v+q*u] = row
        return out
    def adjugate(matrix):
        entries = {i: {j: matrix[i][j] for j in range(n)} for i in range(n)}
        result = []
        for i in range(n):
            row = []
            for j in range(n):
                minor = polynomial_minor(entries, [a for a in range(n) if a != j],
                                          [b for b in range(n) if b != i])
                row.append([(-value if (i+j) % 2 else value) % PRIME for value in minor] or [0])
            result.append(row)
        return result
    ar, ass = adjugate(matrix(0)), adjugate(matrix(1))
    for i, j in product(range(n), repeat=2):
        # Inverse row i is the original column (v,u); inverse column j
        # is the original row (x,y).
        degree = constant-c[3][i % q]-c[2][i//q]-c[0][j % q]-c[1][j//q]
        p = ar[i][j]
        if not any(p):
            assert not any(ass[i][j])
            continue
        assert degree >= len(p)-1
        reverse = (p+[0]*(degree+1-len(p)))[::-1]
        while len(reverse) > 1 and not reverse[-1]:
            reverse.pop()
        assert reverse == ass[i][j]
    new_tensor = {}
    for x, y, u, v in product(range(q), repeat=4):
        p = ass[x+q*y][v+q*u]
        if any(p):
            new_tensor[x, y, u, v] = {e: value for e, value in enumerate(p) if value}
    new_hi = [[-v for v in c[l]] for l in (3, 2, 1, 0)]
    new = literal_case(Random(1), q, new_tensor, ([[0]*q for _ in range(4)], new_hi), 0)
    new_min0 = min(e for p in new_tensor.values() for e in p)
    new_mini = min(-e+sum(new_hi[l][a[l]] for l in range(4))
                   for a, p in new_tensor.items() for e in p)
    opposite = new['profiles']['column cut at zero']
    eta = [v+new_min0 for v in opposite['r']]
    theta = [v+new_mini for v in opposite['s']]
    old = original['profiles']['row cut at zero']
    expected = [sum(old['r'])-v-constant for v in reversed(old['r'])]
    assert theta == expected
    assert all(v >= old['s'][1] for v in eta[1:])
    return {'graded_cofactor_identities': 16, 'eta': eta, 'theta': theta,
            'expected_theta': expected, 'reshuffled_second_order': old['s'][1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, literal, checks, identities = Random(0x01A195C2AD), [], 0, 0
    for q in (2, 3):
        for trial in range(20 if q == 2 else 8):
            heights = tuple([[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2))
            family = trial % 3
            tensor = (monomial_crossing(rng, q, 2) if family == 0 else
                      tight_crossing(q, heights) if family == 1 else random_crossing(rng, q, 2, .6))
            record = literal_case(rng, q, tensor, heights, 0)
            slacks = []
            for pair in record['profiles'].values():
                items = inequalities(record['deficit'], pair['r'], pair['s'],
                                     record['paired_energies_X'], record['paired_energies_Y'])
                assert all(item['slack'] >= 0 for item in items)
                absolute = prefix_form(record['deficit'], pair['r'], pair['s'],
                                       record['paired_energies_X'], record['paired_energies_Y'], F(17, 3))
                assert [a['slack'] for a in items] == [a['slack'] for a in absolute]
                slacks += [str(item['slack']) for item in items]
                checks += len(items)
            entry = {'q': q, 'trial': trial, 'family': family, 'heights': heights,
                     'normalization': record, 'slacks': slacks,
                     'tensor': [[list(a), sorted(p.items())] for a, p in sorted(tensor.items())]}
            if q == 2 and trial < 12:
                entry['cofactor_controls'] = cofactor_controls(tensor, heights, record)
                identities += 16
            literal.append(entry)
    excluded = []
    points = [(F(1), [0, 0, 2, F(8, 3), 3]),
              (F(7, 6), [0, 0, F(19, 9), F(23, 9), F(37, 9)]),
              (F(8, 7), [0, 0, F(15, 8), F(27, 10), 4, 4])]
    for deficit, energy in points:
        energy = tuple(map(F, energy))
        orders, _, _ = relaxed_orders(deficit, energy, energy)
        violations = [r for r in inequalities(deficit, orders, orders, energy, energy) if r['slack'] < 0]
        assert violations
        excluded.append({'deficit': str(deficit), 'energies': list(map(str, energy)),
                         'old_cell_cap_orders': list(map(str, orders)), 'violations': violations})
    energy = (0, 0, 12, 16, 18)
    old_integer = (0,6,9,9,14,19,19,20,22,22,22,22,24,24,24,24,30,34,34,36,36,38,40,40,42)
    assert any(r['slack'] < 0 for r in inequalities(6, old_integer, old_integer, energy, energy))
    new_integer = (0,6,9,9,14,19,19,20,22,22,22,22,24,24,24,24,30,34,34,36,36,37,38,40,41)
    new_caps = inequalities(6, new_integer, new_integer, energy, energy)
    assert all(r['slack'] >= 0 for r in new_caps)
    _, _, old_caps = relaxed_orders(F(6), tuple(map(F, energy)), tuple(map(F, energy)))
    assert all(sum(new_integer[:k]) <= bound for k, (bound, _) in old_caps.items())
    arithmetic = Dyadic(90)
    negative = pressure(arithmetic, F(1), [F(v, 6) for v in new_integer],
                        [F(v, 6) for v in energy], [F(v, 6) for v in energy])
    assert negative[1] < arithmetic.rational(-29, 10000)[0]
    binary_controls = characteristic_two_controls()
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_graded_adjugate_caps.py', 'tools/graded_adjugate_caps.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py']
    names += ['tools/check_joint_energy_conjugations.py', 'tools/search_f4_recentered_crossings.py',
              'tools/weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'literal_tensors': len(literal), 'prime': PRIME,
              'literal_graded_adjugate_inequalities': checks, 'exact_cofactor_identities': identities,
              'independent_characteristic_two_complete_orders': 300,
              'characteristic_two_controls': binary_controls,
              'literal_controls': literal, 'excluded_old_relaxation_points': excluded,
              'new_integer_relaxation': {'energies': energy, 'deficit': 6, 'orders': new_integer,
                  'scale': '1/6', 'pressure_interval': negative, 'precision': 90,
                  'all_new_caps': new_caps, 'scope': 'A necessary-inequality target, not a realizing tensor.'},
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names},
              'scope': 'All-field necessary inequalities proved in the note; finite controls include '
              'direct graded cofactor identities and independent prime-field complete profiles. '
              'Existing negative relaxed points are excluded; a nearby negative relaxation survives. '
              'No improved factoring exponent or general five-letter pressure floor.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True, default=str)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'literal_tensors', 'literal_graded_adjugate_inequalities', 'exact_cofactor_identities')}))


if __name__ == '__main__':
    main()
