#!/usr/bin/env python3
"""Fixed exact high-root projection controls; no root acquisition search."""
import argparse
from hashlib import sha256
import json
from math import prod
from pathlib import Path

from arithmetic_two_group import two_group_basis
from check_positive_relation_kernel import determinant
from high_root_projection import adjugate, project_from_high_roots, two_power
from odd_quotient_relations import word_value

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/high-root-projection.json'


def product_control(moduli, generators, rows, roots, expected, expected_logs):
    identity = (0,)*len(moduli)
    multiply = lambda a, b: tuple((x+y) % n for x, y, n in zip(a, b, moduli))
    inverse = lambda a: tuple(-x % n for x, n in zip(a, moduli))
    b = (prod(moduli)-1).bit_length()
    det = determinant(rows)
    s = (abs(det) & -abs(det)).bit_length()-1
    depth = b+s
    projected = project_from_high_roots(generators, rows, roots, depth, b,
                                         identity, multiply, inverse)
    assert projected == expected
    basis = two_group_basis(projected, identity, multiply, inverse, b)
    assert sorted(k for _, k in basis) == expected_logs
    matrix = adjugate(rows)
    assert [[sum(matrix[i][k]*rows[k][j] for k in range(len(rows)))
             for j in range(len(rows))] for i in range(len(rows))] == [
                 [det*int(i == j) for j in range(len(rows))] for i in range(len(rows))]
    return dict(moduli=moduli, generators=[list(g) for g in generators], row_matrix=rows,
                supplied_roots=[list(z) for z in roots], ambient_order_bit_bound=b, root_depth=depth,
                determinant=det, determinant_two_valuation=s, adjugate=matrix,
                projected=[list(w) for w in projected], extracted_two_order_logs=expected_logs)


def guards():
    multiply = lambda a, b: (a+b) % 4
    inverse = lambda a: -a % 4
    # All algebraic root identities hold, but low precision destroys C4.
    assert two_power(1, 1, multiply) == word_value([1], [2], 0, multiply, inverse)
    unchecked = (1-two_power(1, 0, multiply)) % 4
    assert unchecked == 0
    rejected = []
    cases = [dict(rows=[[2]], roots=[1], depth=1, b=2, label='insufficient_precision'),
             dict(rows=[[4]], roots=[1], depth=4, b=2, label='accepted_identity_control'),
             dict(rows=[[2]], roots=[0], depth=3, b=2, label='false_root'),
             dict(rows=[[0]], roots=[0], depth=2, b=2, label='singular_rows'),
             dict(rows=[[1 << 5]], roots=[1], depth=5, b=2, label='tautological_rows')]
    for case in cases:
        try:
            answer = project_from_high_roots([1], case['rows'], case['roots'],
                        case['depth'], case['b'], 0, multiply, inverse)
        except ValueError:
            assert case['label'] != 'accepted_identity_control'
            rejected.append(case['label'])
        else:
            assert case['label'] == 'accepted_identity_control' and answer == [1]
    # A subgroup-only bound is a FALSE ambient promise and cannot be verified
    # from these black-box identities: in C16, H=<4> has size4, but z=1 is outside H.
    mul16 = lambda a, b: (a+b) % 16
    wrong = project_from_high_roots([4], [[2]], [1], 3, 2, 0, mul16, lambda a: -a % 16)
    assert wrong == [0]
    return dict(rejected_certificates=rejected, low_precision_unchecked_projection=unchecked,
                false_ambient_promise=dict(modulus=16, generator=4, root=1,
                    row=[2], depth=3, falsely_declared_b=2, actual_ambient_b=4,
                    lost_projection=wrong[0]))


def charged_depth_control():
    counts = []
    for depth in [64, 2048]:
        operations = 0

        def multiply(a, b):
            nonlocal operations
            operations += 1
            return (a+b) % 12

        # 2**depth == 4 (mod12); root1 is supplied for row4 of generator1.
        result = project_from_high_roots([1], [[4]], [1], depth, 4,
                                          0, multiply, lambda a: -a % 12)
        assert result == [9]
        counts.append(dict(root_depth=depth, group_multiplications=operations))
    assert counts[1]['group_multiplications']-counts[0]['group_multiplications'] == 2*(2048-64)
    return counts


def run():
    controls = [
        product_control([12, 40], [(1, 0), (0, 1)], [[4, 8], [12, -8]],
                        [(1, 3), (3, 2)], [(3, 0), (0, 15)], [2, 3]),
        product_control([16, 9], [(4, 3)], [[4]], [(1, 3)], [(4, 0)], [2]),
        product_control([15], [(1,)], [[3]], [(3,)], [(0,)], []),
    ]
    # The second root is outside the supplied subgroup: all its first
    # coordinates are divisible by4, whereas that root's first coordinate is1.
    assert controls[1]['supplied_roots'][0][0] % 4 != 0
    names = ['tools/high_root_projection.py', 'tools/check_high_root_projection.py',
             'tools/odd_quotient_relations.py', 'tools/check_positive_relation_kernel.py',
             'tools/arithmetic_two_group.py']
    return dict(scope='Ordinary high-root certificate projection and fixed abstract controls; no root acquirer or complete factoring improvement.',
                product_controls=controls, guards=guards(), charged_depth=charged_depth_control(),
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECEIPT.read_text()), 'receipt differs; inspect before regeneration'
    print('High-root projection: full two-part, precision and ambient guards PASS')
