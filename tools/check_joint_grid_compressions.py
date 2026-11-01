#!/usr/bin/env python3
"""Replay every joint grid through literal-source native weighted profiles."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import comb, factorial, prod
from pathlib import Path

from check_triangular_auxiliary_row_fiber import HEIGHT, LABELS, reflection_seed
from native_matrix_smith_f4 import MatrixOracle

ROOT = Path(__file__).resolve().parents[1]
BASELINE = [0, 1, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11, 12, 13, 15, 16]


def multiply(a, b):
    result = Counter()
    for i, x in a.items():
        for j, y in b.items():
            result[i+j] += x*y
    return {key: value for key, value in result.items() if value}


def power(value, exponent):
    result = {0: 1}
    for _ in range(exponent):
        result = multiply(result, value)
    return result


def bernstein_certificate(data):
    z = multiply(Counter(3*v for v in data['zero_orders']), Counter(3*v for v in data['infinity_orders']))
    f = multiply(Counter(4*v for v in data['heights_x']), Counter(4*v for v in data['heights_y']))
    difference = Counter(power(z, 2))
    for key, value in power(f, 3).items():
        difference[key] -= value
    difference = {key: value for key, value in difference.items() if value}
    order = min(difference)
    degree = max(difference)-order
    terms = {key-order: value for key, value in difference.items()}
    beta = [sum(value*comb(degree-key, i-key) for key, value in terms.items() if key <= i)
            for i in range(degree+1)]
    initial_degree = degree
    while min(beta) < 0 and degree < initial_degree+2048:
        beta = [beta[0]]+[beta[i]+beta[i-1] for i in range(1, len(beta))]+[beta[-1]]
        degree += 1
    assert min(beta) >= 0 and beta[0] > 0 and beta[-1] > 0
    return dict(monomial_order=order, polynomial_degree=initial_degree,
                Bernstein_degree=degree, minimum_integer_coefficient=min(beta),
                coefficient_sha256=sha256(json.dumps(beta, separators=(',', ':')).encode()).hexdigest())


def native_profiles(literal, hx, hy, mapping):
    qx, qy, size = len(hx), len(hy), len(mapping)
    zero = [[None]*size for _ in range(size)]
    infinity = [[None]*size for _ in range(size)]
    for x, y, u, v in product(range(qx), range(qy), range(qy), range(qx)):
        first, second = mapping[x*qy+y], mapping[v*qy+u]
        value = literal[(first//4, first%4, second%4, second//4)]
        assert value.bit_count() % 2 == int(x == v and y == u)
        zero[x*qy+y][v*qy+u] = (value, 0)
        infinity[x*qy+u][v*qy+y] = (value, 0)
    row_zero = [x+y for x in hx for y in hy]
    row_infinity = [y for _ in hx for y in hy]
    column_infinity = [x for x in hx for _ in hy]
    first = MatrixOracle(zero, size).profile(row_zero, [0]*size)['orders']
    second = MatrixOracle(infinity, size).profile(row_infinity, column_infinity, infinity=True)['orders']
    assert all(value == int(value) for value in first+second)
    return sorted(int(value)-4 for value in first), sorted(int(value)+4 for value in second)


def certificate():
    literal = dict(zip(LABELS, reflection_seed()))
    available = Counter(a+b for a in HEIGHT for b in HEIGHT)
    records, native_orders, endpoints, inputs = [], 0, 0, {}
    for qx, qy, expected_count in ((2, 3, 1232), (4, 4, 1152)):
        path = ROOT/f'evidence/joint-grid-{qx}x{qy}.json'
        report = json.loads(path.read_text())
        assert report['completed'] and report['finite_list_exhausted']
        assert len(report['records']) == expected_count
        for name, digest in report['source_sha256'].items():
            assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest, name
        expected = {}
        for a in combinations_with_replacement(range(1, 9), qx-1):
            for b in combinations_with_replacement(range(1, 9), qy-1):
                hx, hy = (0,)+a, (0,)+b
                needed = Counter(x+y for x in hx for y in hy)
                if all(count <= available[h] for h, count in needed.items()):
                    expected[(hx, hy)] = prod(factorial(available[h])//factorial(available[h]-count)
                                             for h, count in needed.items())
        assert sum(expected.values()) == expected_count
        seen, counts = set(), Counter()
        for record in report['records']:
            data = report['classes'][record['profile_class']]
            hx, hy = tuple(data['heights_x']), tuple(data['heights_y'])
            mapping = tuple(record['mapping'])
            assert (hx, hy) in expected and len(mapping) == qx*qy and len(set(mapping)) == len(mapping)
            assert all(0 <= image < 16 and HEIGHT[image//4]+HEIGHT[image%4] == hx[x]+hy[y]
                       for (x, y), image in zip(product(range(qx), range(qy)), mapping))
            key = hx, hy, mapping
            assert key not in seen
            seen.add(key)
            counts[(hx, hy)] += 1
            zero, infinity = native_profiles(literal, hx, hy, mapping)
            assert zero == data['zero_orders'] and infinity == data['infinity_orders']
            native_orders += 2*qx*qy
            endpoints += (qx*qy)**2
        assert dict(counts) == expected
        record = dict(dimensions=[qx, qy], height_grids=len(expected), coordinate_maps=len(seen),
                      complete_profile_classes=len(report['classes']))
        if qx == 2:
            record['balanced_rho_floor'] = [3, 4]
            record['certificates'] = [bernstein_certificate(data) for data in report['classes']]
        else:
            reference_partition = multiply(Counter(HEIGHT), Counter(HEIGHT))
            for data in report['classes']:
                assert multiply(Counter(data['heights_x']), Counter(data['heights_y'])) == reference_partition
                assert data['zero_orders'] == BASELINE
                a = b = 0
                for x, y in zip(data['infinity_orders'], BASELINE):
                    a += x
                    b += y
                    assert a <= b
            record['all_scale_no_gain_over_source_for_any_supporting_weights'] = True
        records.append(record)
        inputs[path.name] = sha256(path.read_bytes()).hexdigest()
    names = (Path(__file__).name, 'joint_grid_compressions.py', 'search_joint_grid_compressions.py',
             'check_triangular_auxiliary_row_fiber.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp', 'native_weighted_smith_f4.py')
    return dict(families=records, independent_native_weighted_orders=native_orders,
                literal_strict_endpoint_entries=endpoints, retained_evidence_sha256=inputs,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Complete declared coordinate-grid families over every characteristic-two extension; 2x3 bound is balanced only, 4x4 dominance is relative to the original source. Arbitrary joint field-valued projections remain open.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-grid-compression-controls.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({key: value for key, value in result.items() if key not in ('families', 'source_sha256')}))


if __name__ == '__main__':
    main()
