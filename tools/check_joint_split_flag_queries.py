#!/usr/bin/env python3
"""Replay finite-query receipts and independent physical-phase jet controls.

This does not rerun Z3 or certify UNSAT independently. Use the query driver
with its recorded field and query arguments to rerun the solver itself.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_fixed_field_descent import Field
from check_joint_split_jet_guards import positive_certificate
from joint_energy_conjugations import model
from joint_five_letter_splits import HEIGHTS, coordinate_start, crossing, profiles

ROOT = Path(__file__).resolve().parents[1]


def jet_controls(field, phases):
    data = model()
    initial_p, initial_l = coordinate_start()
    base = crossing(initial_p, initial_l)
    heights = [x+y for x in HEIGHTS for y in HEIGHTS]
    orders, endpoint_entries, polynomial_entries = 0, 0, 0
    for phase in phases:
        dx, dy = (1, 1, phase[0], phase[2]), (1, 1, phase[1], phase[3])
        character = [field.mul(dx[x], dy[y]) for x, y in product(range(4), repeat=2)]
        p = [[field.mul(v, character[j]) for v in row] for j, row in enumerate(initial_p)]
        l = [[field.mul(v, field.power(character[j], field.order-2)) for v in row]
             for j, row in enumerate(initial_l)]
        assert field.matmul(l, list(map(list, zip(*p)))) == [[int(i == j) for j in range(16)] for i in range(16)]
        raw = [[[0]*19 for _ in range(16)] for _ in range(16)]
        for i, j in product(range(16), repeat=2):
            for r, a in enumerate(l[i]):
                if not a:
                    continue
                for c, b in enumerate(p[j]):
                    if not b:
                        continue
                    value, high = data['operator'][r][c]
                    assert high == 0
                    for k in range(19):
                        if value >> k & 1:
                            raw[i][j][k] ^= field.mul(a, b)
            scalar = field.mul(field.power(character[i], field.order-2), character[j])
            assert raw[i][j] == [scalar if base[i][j][0] >> k & 1 else 0 for k in range(19)]
            assert sum_xor(raw[i][j]) == int(i == j)
            endpoint_entries += 1
            polynomial_entries += 19
        for i, j in product((0, 1, 4, 5), repeat=2):
            assert raw[i][j] == [(base[i][j][0] >> k) & 1 for k in range(19)]
        anchors = [raw[i][j][8] for i, j in ((8, 1), (2, 4), (12, 5), (3, 5))]
        assert anchors == [field.power(v, field.order-2) for v in phase]
        for side in (0, 1):
            jets = []
            for jet in (0, 1):
                if side == 0:
                    matrix = [[raw[i][j][8+jet-heights[i]] if 0 <= 8+jet-heights[i] < 19 else 0
                               for j in range(16)] for i in range(16)]
                else:
                    matrix = [[0]*16 for _ in range(16)]
                    for x, y, u, v in product(range(4), repeat=4):
                        k = heights[4*v+u]+8-jet
                        matrix[4*x+u][4*v+y] = raw[4*x+y][4*v+u][k] if 0 <= k < 19 else 0
                jets.append(matrix)
            first, second = jets
            assert first[0][0] == 1 and field.rank(first) == 1
            toeplitz = [[0]*32 for _ in range(32)]
            for i, j in product(range(16), repeat=2):
                toeplitz[i][j] = toeplitz[i+16][j+16] = first[i][j]
                toeplitz[i+16][j] = second[i][j]
            assert field.rank(toeplitz) == 3
            # A separate Schur calculation checks the residual normalization.
            inverse = [[0]*16 for _ in range(16)]
            for i, j in product(range(1, 16), repeat=2):
                inverse[i][j] = sum_xor((second[i][j], field.mul(first[i][0], second[0][j]),
                                        field.mul(second[i][0], first[0][j]),
                                        field.mul(field.mul(first[i][0], second[0][0]), first[0][j])))
            assert inverse[1][1] == 1 and field.rank(inverse) == 1
            orders += 3
    return dict(field=field.order, phases=len(phases), rank_checks=orders,
                endpoint_entries=endpoint_entries, coefficient_entries=polynomial_entries)


def sum_xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def main():
    receipt_path = ROOT/'evidence/joint-split-flag-queries.json'
    receipt = json.loads(receipt_path.read_text())
    expected = {(q, query) for q in (2, 4, 8, 16) for query in ('zero-anchor', 'normalized')}
    seen, controls = set(), []
    for result in receipt['queries']:
        key = result['field'], result['query']
        assert key in expected and key not in seen
        seen.add(key)
        assert result['tensor_core_only'] and result['field_variables'] == 86
        assert result['multiplication_table_controls'] == result['field']**2
        assert result['result'] in ('unsat', 'unknown')
        if result['result'] == 'unknown':
            assert result['reason_unknown'] == 'timeout'
        for name, digest in result['source_sha256'].items():
            assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest, name
        assert len(result['smt2_sha256']) == 64
        fixed = result['fixed_phase_controls']
        assert all(c['feasible'] == 'sat' and c['phase_exclusion'] == 'unsat' for c in fixed)
        if result['query'] == 'zero-anchor':
            controls.append(jet_controls(Field(result['modulus']), [c['phase'] for c in fixed]))
    assert seen == expected
    classified = [q for q in (2, 4, 8, 16) if all(r['result'] == 'unsat' for r in receipt['queries'] if r['field'] == q)]
    assert classified == receipt['phase_classified_fields']
    baseline = profiles(crossing(*coordinate_start()))
    first, second = (baseline[key] for key in ('zero_orders_t', 'infinity_orders_t'))
    assert first == second == receipt['baseline_orders']
    floor = positive_certificate(Counter(first), {0: 1, 1: 1, 3: 1, 4: 1, 10: 12}, 5, 7)
    floor.pop('polynomial')
    assert floor == receipt['baseline_with_optimistic_other_cut_certificate']
    print(json.dumps(dict(scope='Receipt/source and positive jet replay only; UNSAT is still solver-trusted.',
                          queries=len(seen), phase_classified_fields=classified, independent_controls=controls,
                          baseline_bernstein=floor)))


if __name__ == '__main__':
    main()
