#!/usr/bin/env python3
"""Independent literal replay of the finite joint Laurent-retraction screen."""
from collections import Counter
from hashlib import sha256
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from check_joint_energy_conjugations import polynomial
from check_joint_split_jet_guards import positive_certificate
from joint_moving_retractions import moving_crossing, profiles as native_profiles
from physical_split_crossings import F4_PRODUCT
from search_f4_recentered_crossings import five_letter_source
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
HEIGHTS = (0, 1, 3, 4)


@lru_cache(None)
def literal_source():
    return five_letter_source()


def literal_crossing(specification):
    values, _, energies, _ = literal_source()
    physical = [energies.index(2*h) for h in HEIGHTS]
    indices = [x+5*y for x in physical for y in physical]
    i, out, k, scalar, side = (specification[key] for key in ('inside', 'outside', 'power', 'scalar', 'side'))
    assert 0 <= i < 16 and 0 <= out < 25 and out not in indices
    p = [[(j, 0, 1)] for j in indices]
    l = [[(j, 0, 1)] for j in indices]
    (p if side == 'column' else l)[i].append((out, k, scalar))
    for a, b in product(range(16), repeat=2):
        coefficients = {}
        for r, e, c in l[a]:
            for s, f, d in p[b]:
                if r == s:
                    coefficients[e+f] = coefficients.get(e+f, 0) ^ F4_PRODUCT[c][d]
        assert {e: c for e, c in coefficients.items() if c} == ({0: 1} if a == b else {})
    clearing = -min(0, k)
    matrix = [[None]*16 for _ in range(16)]
    for a, b in product(range(16), repeat=2):
        coefficients = {}
        for r, e, c in l[a]:
            for s, f, d in p[b]:
                pair = values[r+25*(s//5)+125*(s%5)]
                scale = F4_PRODUCT[c][d]
                for j in range(19):
                    coefficient = ((pair[0] >> j) & 1) | (((pair[1] >> j) & 1) << 1)
                    value = F4_PRODUCT[coefficient][scale]
                    if value:
                        degree = j+e+f+clearing
                        assert degree >= 0
                        coefficients[degree] = coefficients.get(degree, 0) ^ value
        matrix[a][b] = tuple(sum(((v >> plane) & 1) << e for e, v in coefficients.items()) for plane in range(2))
        assert (matrix[a][b][0].bit_count() % 2, matrix[a][b][1].bit_count() % 2) == (int(a == b), 0)
    return matrix, clearing


def python_profiles(matrix, clearing):
    first = [[polynomial(value) for value in row] for row in matrix]
    second = [[0]*16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        second[4*x+u][4*v+y] = first[4*x+y][4*v+u]
    zero, _ = WeightedSmithOracle(first).profile([x+y for x in HEIGHTS for y in HEIGHTS], [0]*16)
    infinity, _ = WeightedSmithOracle(second, infinity=True).profile([y for _ in HEIGHTS for y in HEIGHTS],
                                                                 [x for x in HEIGHTS for _ in HEIGHTS])
    return [v-8-clearing for v in zero], [v+8+clearing for v in infinity]


def main():
    path = ROOT/'evidence/joint-moving-monomials.json'
    report = json.loads(path.read_text())
    assert report['completed'] and report['finite_list_exhausted']
    assert (report['min_power'], report['max_power'], report['scalars']) == (-12, 12, [1, 2])
    for name, digest in report['source_sha256'].items():
        assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest, name
    values, _, energy, _ = five_letter_source()
    physical = [energy.index(2*h) for h in HEIGHTS]
    used = {x+5*y for x in physical for y in physical}
    expected = set(product(range(16), sorted(set(range(25))-used), ('row', 'column'), range(-12, 13), (1, 2)))
    seen, classes_used = set(), set()
    for inside, outside, side, power, scalar, cls in report['records']:
        key = inside, outside, side, power, scalar
        assert key in expected and key not in seen
        seen.add(key)
        assert cls in range(len(report['classes']))
        classes_used.add(cls)
        spec = dict(inside=inside, outside=outside, side=side, power=power, scalar=scalar)
        matrix, clearing = literal_crossing(spec)
        code = moving_crossing(**spec)
        assert matrix == code['matrix'] and clearing == code['clearing']
        rebuilt = native_profiles(dict(matrix=matrix, clearing=clearing))
        assert rebuilt['zero_orders'] == report['classes'][cls]['zero_orders']
        assert rebuilt['infinity_orders'] == report['classes'][cls]['infinity_orders']
    assert seen == expected and len(seen) == 14400
    assert classes_used == set(range(len(report['classes'])))
    assert len(report['classes']) == 1745
    input_provenance = {tuple(record[:5]): record[5] for record in report['records']}
    orders, entries, certificates, keys = 0, 0, [], set()
    for cls, result in enumerate(report['classes']):
        spec = result['representative']
        key = tuple(spec[name] for name in ('inside', 'outside', 'side', 'power', 'scalar'))
        assert input_provenance[key] == cls
        matrix, clearing = literal_crossing(spec)
        native_input = moving_crossing(**spec)
        assert native_input['matrix'] == matrix and native_input['clearing'] == clearing
        first, second = python_profiles(matrix, clearing)
        assert first == result['zero_orders'] and second == result['infinity_orders']
        assert (tuple(first), tuple(second)) not in keys
        keys.add((tuple(first), tuple(second)))
        floor = positive_certificate(Counter(first), Counter(second), 3, 4)
        floor.pop('polynomial')
        certificates.append(floor)
        orders += 32
        entries += 256
    controls = 0
    for side, power, scalar in product(('row', 'column'), (-12, -1, 0, 1, 12), (1, 2, 3)):
        spec = dict(inside=15, outside=1, power=power, scalar=scalar, side=side)
        matrix, clearing = literal_crossing(spec)
        assert matrix == moving_crossing(**spec)['matrix']
        # Nonbinary Galois conjugation preserves both full profiles.
        if scalar == 3:
            conjugate = moving_crossing(**dict(spec, scalar=2))['matrix']
            assert matrix == [[(a ^ b, b) for a, b in row] for row in conjugate]
        controls += 1
    print(json.dumps(dict(cases=len(seen), profile_classes=len(keys), independent_complete_orders=orders,
                          independently_constructed_native_orders=32*len(seen),
                          all_literal_endpoint_entries=256*len(seen),
                          independent_tensor_entries=entries, extra_split_controls=controls,
                          rho_floor=[3, 4], maximum_Bernstein_degree=max(c['Bernstein_degree'] for c in certificates),
                          minimum_integer_Bernstein_coefficient=min(c['minimum_integer_coefficient'] for c in certificates),
                          certificates_sha256=sha256(json.dumps(certificates, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                          scope='Complete declared F4 monomial family and exact natural balanced scale inequality; combinations, other coefficients and asymmetric charges open.')))


if __name__ == '__main__':
    main()
