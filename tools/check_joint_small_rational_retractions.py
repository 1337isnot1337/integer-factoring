#!/usr/bin/env python3
"""Exact row-module certificate and a fixed third-order bordered minor."""
import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

from check_joint_hot_retraction_module import product as matrix_product, at_two, rank_f4
from joint_five_letter_splits import coordinate_start, crossing
from joint_hot_retraction_module import degree, dependency, polynomial_right_inverse, polynomial_multiply as mul
from joint_small_rational_retractions import H, jet, project, regular_infinity, regular_zero, row_space, saturated_basis, source_rows, two_pivot_flag

ROOT = Path(__file__).resolve().parents[1]


def patterns(matrix):
    """Cofactor-of-order-one formula for the coefficient-three minors."""
    output = []
    for i, j in product(range(2, 16), repeat=2):
        rr, cc = (0, 1, i), (0, 1, j)
        terms = []
        for a, b in product(range(3), repeat=2):
            r, c = [rr[k] for k in range(3) if k != a], [cc[k] for k in range(3) if k != b]
            cofactor = mul(matrix[r[0]][c[0]], matrix[r[1]][c[1]]) ^ mul(matrix[r[0]][c[1]], matrix[r[1]][c[0]])
            if cofactor >> 1 & 1:
                terms.append((rr[a], cc[b]))
        output.append(terms)
    return output


def minor_image(matrix, formula):
    return sum((sum(matrix[i][j] for i, j in terms) % 2) << k for k, terms in enumerate(formula))


def determinant(matrix, rows=(0, 1, 2), columns=(0, 1, 2)):
    result = 0
    for p in permutations(range(3)):
        value = 1
        for i in range(3):
            value = mul(value, matrix[rows[i]][columns[p[i]]])
        result ^= value
    return result


def symbolic_minor(base, directions):
    """Actual polynomial-ring expansion; repeated indices mean powers."""
    ring = [[[set() for _ in range(5)] for _ in range(3)] for _ in range(3)]
    for i, j, k in product(range(3), range(3), range(5)):
        if base[i][j] >> k & 1:
            ring[i][j][k].add(())
        for q, direction in enumerate(directions):
            if direction[i][j] >> k & 1:
                ring[i][j][k].add((q,))
    def multiply(a, b):
        result = set()
        for p in a:
            for q in b:
                monomial = tuple(sorted(p+q))
                if monomial in result:
                    result.remove(monomial)
                else:
                    result.add(monomial)
        return result
    result = set()
    for p in permutations(range(3)):
        for a in range(5):
            for b in range(5-a):
                c = 4-a-b
                result ^= multiply(multiply(ring[0][p[0]][a], ring[1][p[1]][b]), ring[2][p[2]][c])
    assert result == {()}, result
    return [list(m) for m in sorted(result)]


def certificate():
    basis, steps = saturated_basis()
    assert len(basis) == 9 and [degree(r) for r in basis] == [3, 3, 0, 0, 0, 9, 2, 12, 9]
    inverse, bound = polynomial_right_inverse(basis)
    identity = [[int(i == j) for j in range(9)] for i in range(9)]
    assert bound == 7 and matrix_product(basis, inverse) == identity
    conversion = matrix_product(source_rows(), inverse)
    assert matrix_product(conversion, basis) == list(map(list, source_rows()))
    assert rank_f4(at_two(conversion)) == 9
    assert not dependency([sum(((p >> degree(row)) & 1) << j for j, p in enumerate(row)) for row in basis])
    dims = [[len(row_space(h, k)) for h in range(9)] for k in range(7)]
    assert dims[-1] == [0]*9
    for k, h in product(range(7), range(9)):
        rr = row_space(h, k)
        assert all(p.bit_count() % 2 == 0 for row in rr for p in row)
        assert not dependency([sum(p << (17*j) for j, p in enumerate(row)) for row in rr])
    literal = crossing(*coordinate_start())
    base = [[p for p, high in row] for row in literal]
    assert all(high == 0 for row in literal for _, high in row)
    normalized = [regular_zero(base), regular_infinity(base)]
    first_minor_checks = 0
    for row, h in enumerate(H):
        for raw in row_space(h, 1):
            changed = [r[:] for r in base]
            changed[row] = [a ^ b for a, b in zip(changed[row], raw)]
            for m in (regular_zero(changed), regular_infinity(changed)):
                minor = mul(m[0][0], m[1][1]) ^ mul(m[0][1], m[1][0])
                assert minor & 3 == 2
                first_minor_checks += 1
    assert first_minor_checks == 958
    gauges, formula = [two_pivot_flag(m) for m in normalized], [patterns(m) for m in normalized]
    variables, images, rawimages = [], [], []
    for row, h in enumerate(H):
        for raw in row_space(h, 2):
            jets = [jet(row, raw, 2, side) for side in (0, 1)]
            pair = [minor_image(jets[s], formula[s]) for s in (0, 1)]
            assert pair == [project(jets[s], gauges[s]) for s in (0, 1)]
            images.append(pair[0] | (pair[1] << 196))
            rawimages.append(sum(jets[s][i][j] << (256*s+16*i+j)
                                 for s in (0, 1) for i, j in product(range(16), repeat=2)))
            variables.append((row, raw))
    pivots, nullspace = {}, []
    for i, image in enumerate(images):
        combo = 1 << i
        while image:
            k = image.bit_length()-1
            if k in pivots:
                image ^= pivots[k][0]
                combo ^= pivots[k][1]
            else:
                pivots[k] = image, combo
                break
        if not image:
            nullspace.append(combo)
    quotient, directions = {}, []
    for combo in nullspace:
        image = 0
        for i, v in enumerate(rawimages):
            if combo >> i & 1:
                image ^= v
        while image:
            k = image.bit_length()-1
            if k in quotient:
                image ^= quotient[k][0]
                combo ^= quotient[k][1]
            else:
                quotient[k] = image, combo
                direction = [[0]*16 for _ in range(16)]
                for i, (row, raw) in enumerate(variables):
                    if combo >> i & 1:
                        direction[row] = [a ^ b for a, b in zip(direction[row], raw)]
                directions.append(direction)
                break
    higher = sum(len(row_space(h, 3)) for h in H)
    assert (len(variables), len(pivots), len(nullspace), len(directions), higher) == (284, 140, 144, 6, 138)
    assert len(nullspace)-len(directions) == higher
    assert determinant(normalized[0]) & 15 == 0 and determinant(normalized[0]) >> 4 & 1
    higher_checks = 0
    for row, h in enumerate(H):
        for raw in row_space(h, 3):
            changed = [r[:] for r in base]
            changed[row] = [a ^ b for a, b in zip(changed[row], raw)]
            # A direction beginning at order three enters coefficient four
            # only linearly: every pair of such factors begins at order six.
            assert determinant(regular_zero(changed)) >> 4 & 1
            higher_checks += 1
    symbolic = symbolic_minor(normalized[0], [regular_zero(d) for d in directions])
    names = (Path(__file__).name, 'joint_small_rational_retractions.py', 'joint_hot_retraction_module.py',
             'check_joint_hot_retraction_module.py', 'joint_five_letter_splits.py', 'joint_energy_conjugations.py',
             'search_f4_recentered_crossings.py', 'physical_split_crossings.py')
    return dict(scope='All characteristic-two fields: complete fixed-embedding rational row module; under the stated degree bound, frozen first two raw jets and no order two at either cut, the third zero Smith order is exactly three.',
                basis=basis, polynomial_right_inverse=inverse, source_to_basis=conversion,
                right_inverse_degree=bound, source_rank_at_omega=9, saturation=steps,
                frozen_leading_space_dimension=479, first_minor_coefficient=1,
                first_minor_direction_checks=first_minor_checks,
                bounded_row_dimensions=dims, k2_variables=len(variables), no_order_two_rank=len(pivots),
                no_order_two_dimension=len(nullspace), second_jet_quotient_dimension=len(directions),
                higher_jet_dimension=higher, higher_jet_minor_checks=higher_checks,
                symbolic_minor_coefficient=symbolic, principal_zero_minor_rows=[0, 1, 2],
                principal_zero_minor_columns=[0, 1, 2], coefficient_degree=4,
                second_jet_directions=directions,
                linear_constraint_sha256=sha256(json.dumps(images, separators=(',', ':')).encode()).hexdigest(),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/joint-small-rational-retractions.json').read_text())
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('basis', 'polynomial_right_inverse', 'source_to_basis', 'saturation',
                                   'second_jet_directions', 'source_sha256', 'bounded_row_dimensions')}))


if __name__ == '__main__':
    main()
