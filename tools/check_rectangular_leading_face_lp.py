#!/usr/bin/env python3
"""Check rational LP certificates for a fitted rectangular entry-slope model.

This uses only the Python standard library. It certifies the displayed
finite model, not the unproved continuation of fitted entry slopes.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def dot(left, right):
    return sum((Fraction(a)*Fraction(b) for a, b in zip(left, right)), Fraction())


def physical(full):
    stride = 10 if full else 7
    d = 4+2*stride
    xs, ys = [[0]*d for _ in range(4)], [[0]*d for _ in range(3)]
    for side in (0, 1):
        for state in (1, 2, 3):
            xs[state][4+stride*side+state-1] = 1
            if full:
                xs[state][4+stride*side+7+state-1] = 1
        for state in (1, 2):
            ys[state][4+stride*side+3+state-1] = 1
            ys[state][4+stride*side+5+state-1] = 1
    return xs, [0, -2, 2, 0], ys, [0, -2, 2]


def add_states(base_rows, base_rhs, xs, xb, ys, yb, sx, sy):
    rows, rhs = list(base_rows), list(base_rhs)
    for group, constants, chosen in ((xs, xb, sx), (ys, yb, sy)):
        for k, comparison in enumerate(group):
            rows.append([a-b for a, b in zip(group[chosen], comparison)]+[0, 0])
            rhs.append(constants[k]-constants[chosen])
    return rows, rhs


def verify(certificate, objective, rows, rhs, equalities, targets, constant):
    assert certificate['success']
    x = [Fraction(v) for v in certificate['primal']]
    assert len(x) == len(objective)
    assert all(dot(row, x) <= bound for row, bound in zip(rows, rhs))
    assert all(dot(row, x) == bound for row, bound in zip(equalities, targets))
    y, z = [Fraction() for _ in rows], [Fraction() for _ in equalities]
    for i, value in certificate['dual_inequalities']:
        assert y[i] == 0
        y[i] = Fraction(value)
    for i, value in certificate['dual_equalities']:
        assert z[i] == 0
        z[i] = Fraction(value)
    assert all(value <= 0 for value in y)
    for j, coefficient in enumerate(objective):
        assert sum(row[j]*value for row, value in zip(rows, y)) \
             + sum(row[j]*value for row, value in zip(equalities, z)) == coefficient
    primal = dot(objective, x)+constant
    dual = dot(rhs, y)+dot(targets, z)+constant
    assert primal == dual == Fraction(certificate['optimum'])
    return primal


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=root/'evidence/rectangular-leading-face-lp.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    xs, xb, ys, yb = physical(data['full_core'])
    d = len(xs[0])
    rows, rhs = [], []
    for record in data['constraints']:
        assert len(record['feature']) == d
        row = [-value for value in record['feature']]+[0, 0]
        row[d+(record['place'] == 'zero')] = 1
        rows.append(row)
        rhs.append(record['fitted_slope'])
    eq = [[0]*(d+2) for _ in range(2)]
    eq[0][0] = eq[0][1] = eq[1][2] = eq[1][3] = 1
    matrix = [[None]*3 for _ in range(4)]
    for case in data['cases']:
        sx, sy = case['selected_x'], case['selected_y']
        inequalities, bounds = add_states(rows, rhs, xs, xb, ys, yb, sx, sy)
        objective = [a+b for a, b in zip(xs[sx], ys[sy])]+[-1, -1]
        value = verify(case, objective, inequalities, bounds, eq, [0, 0], xb[sx]+yb[sy])
        matrix[sx][sy] = str(value)
    assert all(value is not None for row in matrix for value in row)
    inequalities, bounds = add_states(rows, rhs, xs, xb, ys, yb, 1, 1)
    c_objective = [a+b for a, b in zip(xs[1], ys[1])]+[-1, -1]
    face_eq, targets = eq+[c_objective], [0, 0, -xb[1]-yb[1]]
    gaps = []
    for case in data['C0_physical_gaps']:
        group, offsets = (xs, xb) if case['letter'] == 'X' else (ys, yb)
        state = case['state']
        objective = [a-b for a, b in zip(group[state], group[1])]+[0, 0]
        value = verify(case, objective, inequalities, bounds, face_eq, targets,
                       offsets[state]-offsets[1])
        gaps.append({'letter': case['letter'], 'state': state, 'minimum_gap': str(value)})
    output = {'status': 'Exact rational primal-dual verification for the finite fitted model.',
              'full_core': data['full_core'], 'leading_C_optima': matrix,
              'C0_physical_gaps': gaps, 'certificates_checked': len(data['cases'])+len(gaps),
              'entry_constraints_used': len(data['used_entry_constraints']),
              'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    rendered = json.dumps(output, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
