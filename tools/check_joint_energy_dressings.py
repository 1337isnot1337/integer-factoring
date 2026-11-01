#!/usr/bin/env python3
"""All-field polynomial dressing reduction by two explicit leading minors."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_joint_leading_tangent import realign
from joint_energy_conjugations import IDENTITY, Space, binary_product, kernel, model, unit_characters

ROOT = Path(__file__).resolve().parents[1]


def polynomial_product(a, b):
    # Encode xi^i eta^j by bit 5*i+j. Inputs have total degree <=1,
    # so the products used in the two minors have no exponent aliasing.
    result = 0
    while a:
        x = a & -a
        second = b
        while second:
            y = second & -second
            result ^= 1 << (x.bit_length()+y.bit_length()-2)
            second ^= y
        a ^= x
    return result


def certificate():
    data = model()
    leading_square = binary_product(data['weighted'][-8], data['weighted'][-8])
    assert binary_product(leading_square, leading_square) == IDENTITY
    characters, full_characters = unit_characters(data)
    radical = [data['leading'].matrix((v, 0))[0] for v in kernel(characters, 10)]
    trivial = [data['full'].matrix((v, 0))[0] for v in kernel(full_characters, 4)]
    assert len(radical) == 4 and len(trivial) == 2
    assert all(binary_product(a, b) == 0 for a, b in product(radical, repeat=2))
    for value in trivial:
        Space(radical).binary_coordinates(value)
    quadratic_identities = 0
    for coefficient in data['weighted'].values():
        for i, a in enumerate(radical):
            assert binary_product(binary_product(a, coefficient), a) == 0
            quadratic_identities += 1
            for b in radical[i+1:]:
                assert binary_product(binary_product(a, coefficient), b) == binary_product(
                    binary_product(b, coefficient), a)
                quadratic_identities += 1
    quotient = Space(trivial+radical)
    assert len(quotient.basis) == 4
    first, second = quotient.basis[2:]
    next_coefficient = data['weighted'][-7]
    d1 = binary_product(first, next_coefficient) ^ binary_product(next_coefficient, first)
    d2 = binary_product(second, next_coefficient) ^ binary_product(next_coefficient, second)
    assert binary_rank([d1, d2]) == 2
    constant, p, q = map(realign, (data['weighted'][-8], d1, d2))
    entries = [((constant >> i) & 1) | (((p >> i) & 1) << 5) | (((q >> i) & 1) << 1)
               for i in range(625)]
    minors = []
    for rows, columns, expected in (((0, 3), (0, 9), 1 << 5), ((0, 1), (0, 23), 1 << 1)):
        r, s = rows
        c, d = columns
        value = polynomial_product(entries[25*r+c], entries[25*s+d]) ^ polynomial_product(
            entries[25*r+d], entries[25*s+c])
        assert value == expected
        minors.append(dict(rows=rows, columns=columns, polynomial=value))
    return dict(nilpotent_sector_dimension=4, trivial_sector_dimension=2,
                leading_operator_fourth_power_identity=True,
                effective_polynomial_parameters=2, sector_square_zero_products=16,
                quadratic_matrix_identities=quadratic_identities,
                quotient_generators=[first, second], leading_minor_witnesses=minors,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_energy_dressings.py', 'joint_energy_conjugations.py',
                                            'check_joint_leading_tangent.py', 'search_f4_recentered_crossings.py')},
                scope='Every characteristic-two field and every polynomial degree, for radical-valued energy-preserving unimodular dressings retaining the original worst infinity pole and rank-one leading coefficient.',
                conclusion='The transformed seed equals a constant conjugate; if the dressing is identity at zero, the seed is unchanged.',
                rank_one_extension='Every nonconstant scalar-polynomial rank-one energy dressing retaining that pole and leading rank is trivial, even without assuming its generator lies in the leading centralizer.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-energy-dressings.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == json.loads(json.dumps(result))
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'quotient_generators')}))


if __name__ == '__main__':
    main()
