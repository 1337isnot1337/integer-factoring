#!/usr/bin/env python3
"""Independently check singular Clifford profiles against literal Fock matrices.

The reference uses coefficient-list F4 polynomials and reduced rational
functions. It neither calls the Clifford construction nor its Bareiss engine.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from random import Random
from time import process_time

from tensor_singular_affine_clifford import (
    singular_profile, polynomial_smith_orders, exact_divide, multiply,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE_NAMES = ('tensor_singular_affine_clifford.py',
                'check_tensor_singular_affine_clifford.py',
                'tensor_affine_clifford_profiles.py')


def field_mul(a, b):
    value = 0
    for i in range(2):
        for j in range(2):
            if (a >> i) & (b >> j) & 1:
                value ^= 1 << (i+j)
    if value & 4:
        value ^= 7
    return value


TABLE = [[field_mul(a, b) for b in range(4)] for a in range(4)]


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return tuple(a)


def padd(a, b):
    return trim([(a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def pmul(a, b):
    if not a or not b:
        return ()
    result = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] ^= TABLE[x][y]
    return trim(result)


def pdiv(a, b):
    if not b:
        raise ZeroDivisionError
    result = [0]*max(0, len(a)-len(b)+1)
    left = list(a)
    inverse = TABLE[b[-1]][b[-1]]
    while len(left) >= len(b):
        power = len(left)-len(b)
        c = TABLE[left[-1]][inverse]
        result[power] = c
        for j, value in enumerate(b):
            left[power+j] ^= TABLE[c][value]
        while left and not left[-1]:
            left.pop()
    return trim(result), tuple(left)


def pgcd(a, b):
    while b:
        a, b = b, pdiv(a, b)[1]
    if not a:
        return ()
    inverse = TABLE[a[-1]][a[-1]]
    return tuple(TABLE[x][inverse] for x in a)


def rat(a, b=(1,)):
    if not b:
        raise ZeroDivisionError
    if not a:
        return (), (1,)
    common = pgcd(a, b)
    a, ra = pdiv(a, common)
    b, rb = pdiv(b, common)
    assert not ra and not rb
    inverse = TABLE[b[-1]][b[-1]]
    return tuple(TABLE[x][inverse] for x in a), tuple(TABLE[x][inverse] for x in b)


def radd(a, b):
    return rat(padd(pmul(a[0], b[1]), pmul(b[0], a[1])), pmul(a[1], b[1]))


def rmul(a, b):
    return rat(pmul(a[0], b[0]), pmul(a[1], b[1]))


def rdiv(a, b):
    return rat(pmul(a[0], b[1]), pmul(a[1], b[0]))


def porder(a):
    return next(i for i, value in enumerate(a) if value)


def rorder(a):
    return porder(a[0])-porder(a[1])


def reference_orders(matrix):
    current = [[rat(value) for value in row] for row in matrix]
    orders = []
    while current and current[0]:
        found = min(((rorder(x), i, j) for i, row in enumerate(current)
                     for j, x in enumerate(row) if x[0]), default=None)
        if found is None:
            break
        order, i, j = found
        orders.append(order)
        current[0], current[i] = current[i], current[0]
        for row in current:
            row[0], row[j] = row[j], row[0]
        pivot = current[0][0]
        current = [[radd(row[k], rmul(rdiv(row[0], pivot), current[0][k]))
                    for k in range(1, len(current[0]))] for row in current[1:]]
    assert orders == sorted(orders)
    return orders


def coefficients(planes):
    return trim([((planes[0] >> i) & 1) | (((planes[1] >> i) & 1) << 1)
                 for i in range(max(x.bit_length() for x in planes))])


def planes(coeffs):
    a = b = 0
    for i, c in enumerate(coeffs):
        a |= (c & 1) << i
        b |= ((c >> 1) & 1) << i
    return a, b


def ident(size):
    return [[(1,) if i == j else () for j in range(size)] for i in range(size)]


def product(A, B):
    result = [[() for _ in B[0]] for _ in A]
    for i, row in enumerate(A):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(B[k]):
                    if y:
                        result[i][j] = padd(result[i][j], pmul(x, y))
    return result


def literal_fock(n, factors):
    size = 1 << n
    result, clearing = ident(size), 0
    for factor in factors:
        kind = factor['kind']
        if kind == 'local':
            mode = factor['mode']
            gate = [[() for _ in range(size)] for _ in range(size)]
            for source in range(size):
                b = (source >> mode) & 1
                for a in (0, 1):
                    target = (source & ~(1 << mode)) | (a << mode)
                    gate[target][source] = coefficients(factor['matrix'][a][b])
        elif kind in ('scalar', 'diagonal'):
            if kind == 'scalar':
                powers = [factor['power']]*size
            else:
                powers = [sum(value for i, value in enumerate(factor['weights'])
                              if (state >> i) & 1) for state in range(size)]
            shift = max(0, -min(powers))
            clearing += shift
            gate = [[() for _ in range(size)] for _ in range(size)]
            for i, power in enumerate(powers):
                gate[i][i] = (0,)*(power+shift)+(1,)
        else:
            gate = ident(size)
            i, j = factor['modes']
            c = coefficients(factor['coefficient'])
            for source in range(size):
                bi, bj = (source >> i) & 1, (source >> j) & 1
                if kind == 'create_pair' and not bi and not bj:
                    target = source | (1 << i) | (1 << j)
                elif kind == 'annihilate_pair' and bi and bj:
                    target = source & ~(1 << i) & ~(1 << j)
                elif kind == 'mode_shear' and not bi and bj:
                    target = (source | (1 << i)) & ~(1 << j)
                else:
                    continue
                gate[target][source] = padd(gate[target][source], c)
        result = product(gate, result)
    return result, clearing


def matrix_rank(matrix):
    rows = [list(row) for row in matrix]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = TABLE[rows[rank][column]][rows[rank][column]]
        rows[rank] = [TABLE[x][inverse] for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            c = rows[i][column]
            if c:
                rows[i] = [x ^ TABLE[c][y] for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def toeplitz_rank(matrix, precision):
    size = len(matrix)
    block = [[0]*(size*precision) for _ in range(size*precision)]
    for a in range(precision):
        for b in range(a+1):
            for i, row in enumerate(matrix):
                for j, polynomial in enumerate(row):
                    if a-b < len(polynomial):
                        block[a*size+i][b*size+j] = polynomial[a-b]
    return matrix_rank(block)


def local(mode, matrix):
    return {'kind': 'local', 'mode': mode, 'matrix': matrix}


def case_inputs():
    zero, one = (0, 0), (1, 0)
    cases = []
    for mask in range(256):
        values = [planes([(mask >> (2*i)) & 3]) for i in range(4)]
        cases.append((f'F4-local-{mask}', 1, [local(0, [values[:2], values[2:]])]))
    for a in (1, 3, 9):
        cases.append((f'nilpotent-order-{a}', 1,
                      [local(0, [[zero, (1 << a, 0)], [zero, zero]])]))
    cases += [
        ('opposing-projectors', 1, [local(0, [[one, zero], [zero, zero]]),
                                   local(0, [[zero, zero], [zero, one]])]),
        ('negative-Laurent-scalar', 1, [local(0, [[one, zero], [zero, zero]]),
                                       {'kind': 'scalar', 'power': -7}]),
        ('free-three-mode-identity', 3, []),
    ]
    rng = Random(113071)
    for n in (2, 3, 4):
        for index in range(12 if n == 4 else 24):
            factors = []
            for step in range(2):
                mode = rng.randrange(n)
                if index % 3 == 0:
                    u = [trim([rng.randrange(4), rng.randrange(2)]) for _ in range(2)]
                    v = [trim([rng.randrange(4)]) for _ in range(2)]
                    matrix = [[planes(pmul(a, b)) for b in v] for a in u]
                elif index % 3 == 1:
                    power = rng.randrange(2)
                    matrix = [[zero, (1 << power, 0)], [zero, zero]]
                else:
                    matrix = [[planes(trim([rng.randrange(4), rng.randrange(2)]))
                               for _ in range(2)] for _ in range(2)]
                factors.append(local(mode, matrix))
                i, j = rng.sample(range(n), 2)
                factors.append({'kind': ('create_pair', 'annihilate_pair', 'mode_shear')[(index+step) % 3],
                                'modes': [i, j], 'coefficient': planes([rng.randrange(1, 4), 1])})
            factors.insert(2, {'kind': 'diagonal', 'weights': [rng.randrange(-1, 2) for _ in range(n)]})
            factors.append({'kind': 'scalar', 'power': rng.randrange(-2, 3)})
            cases.append((f'coupled-{n}-{index}', n, factors))
    return cases


def compact_large_controls():
    """Known tensor profiles under integral invertible factors at both ends."""
    zero, one = (0, 0), (1, 0)
    output = []
    for n in (8, 16, 32):
        weights = [0, 0] + [i % 3 for i in range(n-2)]
        factors = [
            {'kind': 'create_pair', 'modes': [0, n-1], 'coefficient': (3, 0)},
            local(0, [[zero, (2, 0)], [zero, zero]]),
            local(1, [[one, zero], [zero, zero]]),
            {'kind': 'diagonal', 'weights': weights},
            {'kind': 'annihilate_pair', 'modes': [1, n-2], 'coefficient': (1, 1)},
        ]
        profile = singular_profile(n, factors)
        assert 'nonzero_orders' not in profile
        assert profile['rank'] == 1 << (n-2)
        assert profile['minimum_nonzero_order'] == 1
        assert profile['nonzero_binomial_gaps'] == sorted(weights[2:])
        output.append({'physical_modes': n, 'factors': factors, 'profile': profile,
                       'literal_Fock_matrix_materialized': False,
                       'reference': 'Disjoint local tensor factors with integral unimodular Gaussian factors at both ends.'})
    bad = [(0, []), (1, [local(1, [[one, zero], [zero, one]])]),
           (1, [local(0, [[one]])]),
           (1, [local(0, [[(-1, 0), zero], [zero, one]])]),
           (2, [{'kind': 'create_pair', 'modes': [0, 0], 'coefficient': one}]),
           (2, [{'kind': 'diagonal', 'weights': [1]}]),
           (1, [{'kind': 'scalar', 'power': 0.5}]),
           (1, [{'kind': 'unknown'}])]
    for n, factors in bad:
        try:
            singular_profile(n, factors)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid circuit input was accepted.')
    return output, len(bad)


def arithmetic_controls():
    rng = Random(43319)
    divisions = comparisons = 0
    for _ in range(128):
        a = trim([rng.randrange(4) for _ in range(rng.randrange(1, 7))])
        b = trim([rng.randrange(4) for _ in range(rng.randrange(1, 6))]) or (1,)
        assert exact_divide(multiply(planes(a), planes(b)), planes(b)) == planes(a)
        divisions += 1
    for _ in range(96):
        h, w = rng.randrange(1, 5), rng.randrange(1, 5)
        matrix = [[trim([rng.randrange(4) for _ in range(rng.randrange(4))])
                   for _ in range(w)] for _ in range(h)]
        assert polynomial_smith_orders([[planes(x) for x in row] for row in matrix]) == reference_orders(matrix)
        comparisons += 1
    rejected = 0
    for numerator, denominator in [((1, 0), (2, 0)), ((3, 0), (4, 0))]:
        try:
            exact_divide(numerator, denominator)
        except ArithmeticError:
            rejected += 1
        else:
            raise AssertionError('A nonexact division was accepted.')
    assert rejected == 2
    return {'exact_polynomial_divisions': divisions,
            'independent_rectangular_Smith_comparisons': comparisons,
            'nonexact_divisions_rejected': rejected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes = {f'tools/{name}': sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
              for name in SOURCE_NAMES}
    start = process_time()
    arithmetic = arithmetic_controls()
    results = []
    toeplitz_checks = literal_entries = nonzero_orders = 0
    for label, n, factors in case_inputs():
        profile = singular_profile(n, factors, expand=True)
        literal, clearing = literal_fock(n, factors)
        independent = [value-clearing for value in reference_orders(literal)]
        assert profile['nonzero_orders'] == independent, label
        assert profile['rank'] == len(independent), label
        literal_entries += len(literal)**2
        nonzero_orders += len(independent)
        if label.startswith(('coupled-', 'nilpotent-', 'opposing-')):
            for precision in (1, 2, 4):
                expected = sum(max(0, precision-value-clearing) for value in independent)
                assert toeplitz_rank(literal, precision) == expected, (label, precision)
                toeplitz_checks += 1
        if label.startswith('nilpotent-order-'):
            a = int(label.rsplit('-', 1)[1])
            assert profile['regularized_binomial_gaps'] == [[2, -2*a]]
        if label == 'opposing-projectors':
            assert profile['rank'] == 0 and profile['regularized_minimum_pair'][0] == 1
        results.append({'case': label, 'physical_modes': n, 'factors': factors,
                        'profile': profile, 'literal_clearing': clearing,
                        'independent_nonzero_orders': independent})
    compact, invalid = compact_large_controls()
    assert hashes == {f'tools/{name}': sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                      for name in SOURCE_NAMES}
    report = {'claim': 'Complete nonzero Smith profiles of singular one-mode/Gaussian circuits',
              'scope': 'Uniform characteristic-two theorem with finite F4 implementation controls; no pressure or factoring-exponent improvement claimed.',
              'source_sha256': hashes, 'arithmetic_controls': arithmetic,
              'case_count': len(results), 'literal_matrix_entries': literal_entries,
              'independent_nonzero_orders': nonzero_orders,
              'independent_Toeplitz_prefix_ranks': toeplitz_checks,
              'compact_large_controls': compact, 'invalid_inputs_rejected': invalid,
              'cpu_seconds': process_time()-start, 'cases': results}
    if args.output:
        args.output.write_text(json.dumps(report, sort_keys=True, indent=2)+'\n')
    summary = {k: v for k, v in report.items() if k not in ('cases', 'compact_large_controls')}
    summary['compact_large_controls'] = [
        {'physical_modes': item['physical_modes'],
         'clifford_matrix_order': item['profile']['clifford_matrix_order'],
         'rank': item['profile']['rank'], 'literal_Fock_matrix_materialized': False}
        for item in compact]
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
