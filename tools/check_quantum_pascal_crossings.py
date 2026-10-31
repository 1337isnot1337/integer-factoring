#!/usr/bin/env python3
"""Exact endpoints and bare-profile obstructions for quantum-Pascal crossings."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

from check_fixed_field_descent import Field
from quantum_pascal_crossing import construct, cut_blocks, evaluate, quantum_binomials
from weighted_smith_oracle import WeightedSmithOracle, polynomial_multiply


def determinant(matrix):
    answer = 0
    for permutation in permutations(range(len(matrix))):
        term = 1
        for i, j in enumerate(permutation):
            term = polynomial_multiply(term, matrix[i][j])
        answer ^= term
    return answer


def primitive_root(field, order):
    for a in range(1, field.order):
        if field.power(a, order) == 1 and all(field.power(a, j) != 1 for j in range(1, order)):
            return a
    raise ValueError('the field has no root of the requested order')


def check():
    factorials = [1]
    for i in range(1, 16):
        factorials.append(polynomial_multiply(factorials[-1], (1 << i)-1))
    qbinomials = quantum_binomials(15)
    identities = 0
    for n in range(16):
        for k in range(n+1):
            assert polynomial_multiply(qbinomials[n][k], polynomial_multiply(factorials[k], factorials[n-k])) == factorials[n]
            identities += 1
    records, endpoint_entries, trace_entries, determinants, scalar_tests = [], 0, 0, 0, 0
    for d, modulus in ((3, 7), (5, 19), (7, 11), (9, 67), (15, 19)):
        field = Field(modulus)
        assert all(field.power(a, field.order-1) == 1 for a in range(1, field.order))
        root = primitive_root(field, d)
        tensor, degree, ladder = construct(d)
        assert degree == d*(d-1)
        cache = {value: evaluate(value, root, field) for value in set(tensor.values())}
        for x, y, u, v in product(range(d), repeat=4):
            assert cache.get(tensor.get((x, y, u, v), 0), 0) == int(x == v and y == u)
            endpoint_entries += 1
        # The high total-grade blocks each have a literal constant-one entry.
        for y in range(d):
            assert tensor[d-1, y, y, d-1] == 1
        raw_infinity, zero = [], []
        if d <= 5:
            for cut in (0, 1):
                for _, _, _, matrix in cut_blocks(tensor, d, cut):
                    n = len(matrix)
                    orders, minors = WeightedSmithOracle(matrix, infinity=cut == 0).profile([0]*n, [0]*n)
                    (raw_infinity if cut == 0 else zero).extend(orders)
                    if cut == 0:
                        assert determinant(matrix) == 1 == minors[-1].polynomial
                        determinants += 1
            assert sum(raw_infinity) == 0
            assert zero.count(0) >= d
            for scale in map(Fraction, ('1', '5/4', '3/2', '2', '3')):
                cost = sum(scale**(-v) for v in raw_infinity)*sum(scale**(-v) for v in zero)
                assert cost >= d**3
                scalar_tests += 1
            # Literal operator matrices provide a second trace construction.
            operators = {}
            for x, y in product(range(d), repeat=2):
                operator = [[0]*d for _ in range(d)]
                for k in range(y, d):
                    target = k-y+x
                    if target < d:
                        operator[target][k] = ladder[k, y]
                operators[x, y] = operator
            for x, y, u, v in product(range(d), repeat=4):
                left, right = operators[x, y], operators[u, v]
                value = 0
                for i in range(d):
                    for j in range(d):
                        value ^= polynomial_multiply(left[i][j], right[j][i])
                assert value == tensor.get((x, y, u, v), 0)
                trace_entries += 1
        records.append(dict(alphabet=d, degree=degree, field_modulus=modulus,
                            root=root, nonzero_polynomials=sum(bool(v) for v in tensor.values()),
                            raw_row_infinity_orders=sorted(raw_infinity), column_zero_orders=sorted(zero)))
    return dict(q_factorial_identities=identities, complete_endpoint_entries=endpoint_entries,
                independent_trace_entries=trace_entries, exact_unit_block_determinants=determinants,
                rational_bare_comparisons=scalar_tests, records=records,
                uniform_bare_rho_floor='3/4', uniform_multipoint_rho_floor='3/4',
                scope='Quantum-Pascal trace basis with arbitrary polynomial diagonal vacuum, P0=1; physical word-type gauges are outside the obstruction.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    result['source_sha256'] = {p.name: sha256(p.read_bytes()).hexdigest()
                              for p in (Path(__file__), Path(__file__).with_name('quantum_pascal_crossing.py'))}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('complete_endpoint_entries', 'independent_trace_entries', 'exact_unit_block_determinants', 'factoring_exponent')}))


if __name__ == '__main__':
    main()
