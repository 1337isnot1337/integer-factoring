#!/usr/bin/env python3
"""Check exact leading obstructions for the full-alphabet F4 quartic repair."""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import check_rectangular_uniform_entry_bounds as symbolic
from check_rectangular_leading_face_lp import add_states, verify
from full_alphabet_quartic_profiles import FullAlphabetProgram
from rectangular_geometric_core import concentrated, MUL


def feature(place, mask, state, row, column):
    if place == 'infinity':
        rc, cc = row.bit_count()+mask.bit_count(), column.bit_count()+state.bit_count()
        v, y, u, x = mask & 3, mask >> 2, state & 3, state >> 2
        side = 0
    else:
        def charge(word):
            return (word & 63).bit_count()-(word >> 6).bit_count()
        rc = charge(row)-(mask & 3).bit_count()+(mask >> 2).bit_count()
        cc = charge(column)+(state & 3).bit_count()-(state >> 2).bit_count()
        y, x, v, u = mask & 3, mask >> 2, state & 3, state >> 2
        side = 1
    result = [0]*28
    result[2*side], result[2*side+1] = rc, cc
    for physical, position in ((v, 0), (y, 3), (u, 6), (x, 9)):
        if physical:
            result[4+12*side+position+physical-1] += 1
    return result


def physical():
    xs, ys = [[0]*28 for _ in range(4)], [[0]*28 for _ in range(4)]
    for side in (0, 1):
        for state in (1, 2, 3):
            for offset in (0, 9):
                xs[state][4+12*side+offset+state-1] = 1
            for offset in (3, 6):
                ys[state][4+12*side+offset+state-1] = 1
    return xs, [0, -2, 2, 0], ys, [0, -2, 2, 0]


def dual_checks(model):
    xs, xb, ys, yb = physical()
    rows, rhs = [], []
    for item in model['constraints']:
        row = [-value for value in item['feature']]+[0, 0]
        row[28+(item['place'] == 'zero')] = 1
        rows.append(row)
        rhs.append(item['slope'])
    equalities = [[0]*30 for _ in range(2)]
    equalities[0][0] = equalities[0][1] = 1
    equalities[1][2] = equalities[1][3] = 1
    table = [[None]*4 for _ in range(4)]
    for case in model['cases']:
        sx, sy = case['selected_x'], case['selected_y']
        a, b = add_states(rows, rhs, xs, xb, ys, yb, sx, sy)
        objective = [x+y for x, y in zip(xs[sx], ys[sy])]+[-1, -1]
        value = verify(case, objective, a, b, equalities, [0, 0], xb[sx]+yb[sy])
        assert value >= (0 if (sx, sy) == (1, 1) else 2)
        table[sx][sy] = str(value)
    a, b = add_states(rows, rhs, xs, xb, ys, yb, 1, 1)
    c_objective = [x+y for x, y in zip(xs[1], ys[1])]+[-1, -1]
    face = equalities+[c_objective]
    targets = [0, 0, 4]
    gaps = []
    for case in model['gaps']:
        group, offsets = (xs, xb) if case['letter'] == 'X' else (ys, yb)
        state = case['state']
        objective = [x-y for x, y in zip(group[state], group[1])]+[0, 0]
        value = verify(case, objective, a, b, face, targets, offsets[state]-offsets[1])
        multiplier = Fraction(dict(case['dual_equalities'])[2])
        assert value >= 2 and -2 <= multiplier <= -1
        gaps.append({'letter': case['letter'], 'state': state,
                     'constant_gap': str(value), 'off_face_C_multiplier': str(multiplier)})
    return {'leading_C_lower_bounds': table, 'gap_lower_bounds': gaps,
            'exact_primal_dual_certificates': len(model['cases'])+len(gaps)}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=root/'evidence/full-alphabet-leading-inputs.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    model = json.loads(args.input.read_text())
    assert (model['gamma'], model['conjugate_gamma'], model['repair_kind'], model['quartic_terms']) \
           == (2, 3, 'partial', [[1, 0]])
    witnesses = model['constraints']
    assert len(witnesses) == 34
    records = []
    original_build = symbolic.SymbolicProgram.build
    try:
        for gamma in (2, 3):
            def build(self, *unused, initial=None, **kwargs):
                self.gamma = gamma
                return FullAlphabetProgram.build(self, initial=initial)
            symbolic.SymbolicProgram.build = build
            for place in ('infinity', 'zero'):
                data = concentrated(place, gamma)
                program, outputs = symbolic.compile_program(data)
                for item in (w for w in witnesses if w['place'] == place):
                    mask, state = item['physical_masks']
                    row, column = item['vacuum_coordinates']
                    assert feature(place, mask, state, row, column) == item['feature']
                    evaluate = symbolic.evaluator(program, column)
                    value = evaluate(outputs[mask, state], row)
                    coefficient = (value.p[0] & 1)+2*(value.p[1] & 1)
                    expected = item['coefficient']
                    if gamma == 3:
                        expected = MUL[expected][expected]
                    assert coefficient == expected and coefficient
                    assert value.base == item['slope'] and value.known > value.base
                    records.append({'gamma': gamma, 'place': place,
                                    'physical_masks': item['physical_masks'],
                                    'vacuum_coordinates': item['vacuum_coordinates'],
                                    'certified_slope_upper_bound': value.base,
                                    'specialized_leading_coefficient': coefficient,
                                    'coefficients_known_below': value.known})
                    evaluate.cache_clear()
    finally:
        symbolic.SymbolicProgram.build = original_build
    duals = dual_checks(model)
    paths = [Path(__file__), root/'tools/full_alphabet_quartic_profiles.py',
             root/'tools/rectangular_geometric_core.py',
             root/'tools/check_rectangular_uniform_entry_bounds.py']
    result = {'status': 'All 34 uniform entry bounds pass for both conjugate quartics; 22 exact LP certificates pass.',
              'scope': 'The full-alphabet partial repair with gamma=omega or omega^2 and balanced aggregate charges.',
              'coefficient_checks': len(records), 'coefficients': records, 'duals': duals,
              'moving_gauge_constants': {'physical_core_cardinalities': [4, 4],
                  'coarse_constant': 'K0 = 20 log 2 + sqrt(2) log 16',
                  'coarse_pressure_bound': 'F >= T - (2B+8) ell - K0',
                  'nonpositive_pressure_gap_bound': 'gap >= 2 tau - 2K0 - (4B+20) ell'},
              'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in paths}}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
