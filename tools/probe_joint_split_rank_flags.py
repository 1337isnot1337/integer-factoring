#!/usr/bin/env python3
"""Finite extension-field queries for joint splits with both low ranks one.

Requires the optional z3-solver package. A solver result is a finite query,
not an all-field theorem or an independently checked unsatisfiability proof.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import z3

from check_fixed_field_descent import Field
from joint_energy_conjugations import model
from joint_five_letter_splits import HEIGHTS, coordinate_start, crossing, profiles
from search_joint_grid_compressions import numerical_cost

FIELD = Field(7)


def field_multiply(a, b):
    if isinstance(a, int) and isinstance(b, int):
        return FIELD.mul(a, b)
    if isinstance(a, int) and a in (0, 1):
        return 0 if a == 0 else b
    if isinstance(b, int) and b in (0, 1):
        return 0 if b == 0 else a
    degree = FIELD.degree
    a = z3.BitVecVal(a, degree) if isinstance(a, int) else a
    b = z3.BitVecVal(b, degree) if isinstance(b, int) else b
    result, shifted = z3.BitVecVal(0, degree), a
    for i in range(degree):
        result ^= shifted & (-(z3.LShR(b, i) & 1))
        carry = z3.LShR(shifted, degree-1)
        shifted = (shifted << 1) ^ ((FIELD.modulus ^ FIELD.order) & (-carry))
    return result


def add(values):
    result = 0
    for value in values:
        result ^= value
    return z3.simplify(result) if not isinstance(result, int) else result


def rank_flags(solver, first, second):
    solver.add(first[0][0] == 1)
    for i, j in product(range(1, 16), repeat=2):
        solver.add(first[i][j] == field_multiply(first[i][0], first[0][j]))
    residual = [[None]*16 for _ in range(16)]
    for i, j in product(range(1, 16), repeat=2):
        residual[i][j] = add((second[i][j], field_multiply(first[i][0], second[0][j]),
                              field_multiply(second[i][0], first[0][j]),
                              field_multiply(field_multiply(first[i][0], second[0][0]), first[0][j])))
    solver.add(residual[1][1] == 1)
    for i, j in product(range(1, 16), repeat=2):
        solver.add(residual[i][j] == field_multiply(residual[i][1], residual[1][j]))


def main():
    global FIELD
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--timeout-ms', type=int, default=600000)
    parser.add_argument('--field-degree', type=int, choices=(1, 2, 3, 4), default=2)
    parser.add_argument('--protect-below', type=int, choices=(3, 4), default=3)
    parser.add_argument('--tensor-core', action='store_true')
    parser.add_argument('--query', choices=('phase', 'zero-anchor', 'normalized'), default='phase')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.timeout_ms < 1:
        parser.error('Use a positive solver timeout.')
    FIELD = Field({1: 3, 2: 7, 3: 11, 4: 19}[args.field_degree])
    began = time.monotonic()
    names = (Path(__file__).name, 'joint_five_letter_splits.py', 'joint_energy_conjugations.py',
             'physical_split_crossings.py', 'search_f4_recentered_crossings.py',
             'search_joint_grid_compressions.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp', 'check_fixed_field_descent.py')
    source_hashes = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    a, b = z3.BitVecs('control_a control_b', FIELD.degree)
    test_product = field_multiply(a, b)
    for x, y in product(range(FIELD.order), repeat=2):
        evaluated = z3.simplify(z3.substitute(test_product, (a, z3.BitVecVal(x, FIELD.degree)),
                                              (b, z3.BitVecVal(y, FIELD.degree)))).as_long()
        assert evaluated == FIELD.mul(x, y)
    data, heights = model(), [a+b for a in HEIGHTS for b in HEIGHTS]
    initial_p, initial_l = coordinate_start()
    p, l, variables = [], [], []
    for kind, initial, target in (('p', initial_p, p), ('l', initial_l, l)):
        for j in range(16):
            row = []
            for r in range(25):
                value = 0
                if data['energies'][r] == heights[j]:
                    protected = heights[j] == 0 if args.tensor_core else heights[j] < args.protect_below
                    if protected:
                        value = initial[j][r]
                    else:
                        value = z3.BitVec(f'{kind}_{j}_{r}', FIELD.degree)
                        variables.append(value)
                row.append(value)
            target.append(row)
    solver = z3.SolverFor('QF_BV')
    solver.set(timeout=args.timeout_ms)
    for i, j in product(range(16), repeat=2):
        solver.add(add(field_multiply(l[i][r], p[j][r]) for r in range(25)) == int(i == j))
    raw = [[[0]*19 for _ in range(16)] for _ in range(16)]
    for i, j in product(range(16), repeat=2):
        terms = []
        for r in range(25):
            if isinstance(l[i][r], int) and l[i][r] == 0:
                continue
            for c in range(25):
                if isinstance(p[j][c], int) and p[j][c] == 0:
                    continue
                value, high = data['operator'][r][c]
                assert high == 0
                if value:
                    terms.append((value, field_multiply(l[i][r], p[j][c])))
        for k in range(19):
            raw[i][j][k] = add(value for mask, value in terms if mask >> k & 1)
    zero, infinity = [], []
    for jet in (0, 1):
        a = [[raw[i][j][8+jet-heights[i]] if 0 <= 8+jet-heights[i] < 19 else 0
              for j in range(16)] for i in range(16)]
        b = [[0]*16 for _ in range(16)]
        for x, y, u, v in product(range(4), repeat=4):
            i, j, degree = 4*x+y, 4*v+u, heights[4*v+u]+8-jet
            b[4*x+u][4*v+y] = raw[i][j][degree] if 0 <= degree < 19 else 0
        zero.append(a)
        infinity.append(b)
    rank_flags(solver, *zero)
    rank_flags(solver, *infinity)
    baseline = crossing(initial_p, initial_l)
    if args.tensor_core:
        for i, j in product((0, 1, 4, 5), repeat=2):
            for k in range(19):
                solver.add(raw[i][j][k] == ((baseline[i][j][0] >> k) & 1))
    # Four literal coefficients recover the independent physical phases.
    # Multiplying through by their characters avoids inverses and works
    # over each declared finite field, including when an anchor is zero.
    specifications = ((8, 1), (2, 4), (12, 5), (3, 5))
    anchors = [raw[i][j][8] for i, j in specifications]
    assert all((baseline[i][j][0] >> 8) & 1 for i, j in specifications)
    characters = []
    for x, y in product(range(4), repeat=2):
        first = anchors[0 if x == 2 else 2] if x in (2, 3) else 1
        second = anchors[1 if y == 2 else 3] if y in (2, 3) else 1
        characters.append(field_multiply(first, second))
    different = [value == 0 for value in anchors]
    for i, j in product(range(16), repeat=2):
        for k in range(19):
            expected = characters[i] if (baseline[i][j][0] >> k) & 1 else 0
            different.append(field_multiply(characters[j], raw[i][j][k]) != expected)
    phase_exclusion = z3.Or(*different)
    # Check feasibility and the exclusion predicate on independently supplied
    # physical diagonal phases. All map variables are fixed in these controls.
    controls = []
    phases = [(1, 1, 1, 1)]
    if FIELD.order > 2:
        phases.extend(((2, 1, 1, 1), (2, FIELD.order-1, FIELD.order-1, 2))
                      if args.tensor_core or args.protect_below == 3
                      else ((1, 1, 2, 1), (1, 1, FIELD.order-1, 2)))
    for phase in phases:
        phase_x, phase_y = (1, 1, phase[0], phase[2]), (1, 1, phase[1], phase[3])
        values = [FIELD.mul(phase_x[x], phase_y[y]) for x, y in product(range(4), repeat=2)]
        fixed = []
        for matrix, initial, inverse in ((p, initial_p, False), (l, initial_l, True)):
            for j in range(16):
                value = FIELD.power(values[j], FIELD.order-2) if inverse else values[j]
                for r in range(25):
                    expected = FIELD.mul(initial[j][r], value)
                    if isinstance(matrix[j][r], int):
                        assert matrix[j][r] == expected
                    else:
                        fixed.append(matrix[j][r] == expected)
        solver.push()
        solver.add(*fixed)
        feasible = solver.check()
        assert feasible == z3.sat, feasible
        solver.push()
        solver.add(phase_exclusion)
        excluded = solver.check()
        assert excluded == z3.unsat, excluded
        solver.pop(2)
        controls.append(dict(phase=list(phase), feasible=str(feasible), phase_exclusion=str(excluded)))
    if args.query == 'phase':
        solver.add(phase_exclusion)
    elif args.query == 'zero-anchor':
        solver.add(z3.Or(*(value == 0 for value in anchors)))
    else:
        solver.add(*(value == 1 for value in anchors))
        solver.add(z3.Or(*(raw[i][j][k] != ((baseline[i][j][0] >> k) & 1)
                           for i, j in product(range(16), repeat=2) for k in range(19))))
    smt_path = args.output.with_suffix('.smt2')
    smt_path.write_text(solver.to_smt2())
    build_seconds = time.monotonic()-began
    print(json.dumps(dict(field=FIELD.order, field_variables=len(variables), constraints=len(solver.assertions()),
                          tensor_core_only=args.tensor_core, build_seconds=build_seconds)), flush=True)
    answer = solver.check()
    report = dict(result=str(answer), field=FIELD.order, modulus=FIELD.modulus, field_variables=len(variables),
                  protected_joint_heights_below=None if args.tensor_core else args.protect_below,
                  tensor_core_only=args.tensor_core, query=args.query,
                  excluded_physical_phase_orbit=args.query == 'phase', fixed_phase_controls=controls,
                  multiplication_table_controls=FIELD.order**2,
                  z3_version=z3.get_version_string(), timeout_ms=args.timeout_ms,
                  build_seconds=build_seconds, total_seconds=time.monotonic()-began,
                  smt2_sha256=sha256(smt_path.read_bytes()).hexdigest(),
                  scope='Finite-field solver query for the declared protected splits; unsat requires trusting this solver run, and is not an all-field proof.')
    if answer == z3.sat:
        solution = solver.model()
        evaluate = lambda value: value if isinstance(value, int) else solution.eval(value, model_completion=True).as_long()
        pp, ll = ([[evaluate(value) for value in row] for row in matrix] for matrix in (p, l))
        assert FIELD.matmul(ll, list(zip(*pp))) == [[int(i == j) for j in range(16)] for i in range(16)]
        report.update(embedding=pp, retraction=ll, independent_retraction_replay=True)
        if FIELD.degree <= 2:
            operator = crossing(pp, ll)
            profile = profiles(operator)
            assert profile['zero_orders_t'].count(0) == profile['infinity_orders_t'].count(0) == 1
            assert profile['zero_orders_t'].count(1) == profile['infinity_orders_t'].count(1) == 1
            cost = numerical_cost(HEIGHTS, HEIGHTS, dict(zero_orders=profile['zero_orders_t'], infinity_orders=profile['infinity_orders_t']))
            report.update(profile=profile, numerical_cost=cost, native_model_replay=True)
            print(json.dumps(dict(result='sat', cost=cost)), flush=True)
        else:
            report['native_model_replay'] = False
            report['remaining_obligation'] = 'Recompute the complete profiles over this larger field before making any pressure claim.'
    elif answer == z3.unknown:
        report['reason_unknown'] = solver.reason_unknown()
    assert source_hashes == {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    report['source_sha256'] = source_hashes
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('embedding', 'retraction', 'profile', 'numerical_cost', 'source_sha256')}))


if __name__ == '__main__':
    main()
