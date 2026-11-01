#!/usr/bin/env python3
"""Exact controls for the paired-spectrum average face-rank floor."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from cube_symmetric_exterior_blocks import multiplication_matrix
from high_axis_plucker import elementary, power, MUL
from recursive_symmetric_exterior_blocks import recursive_blocks
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / 'evidence/multiaxis-face-rank-controls.json'


def binary_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            i = value.bit_length()-1
            if i not in pivots:
                pivots[i] = value
                break
            value ^= pivots[i]
    return len(pivots)


def field_rank(matrix):
    rows = []
    for row in matrix:
        low = high = 0
        for j, value in enumerate(row):
            for bit, scalar in enumerate((value, MUL[value][2])):
                low |= (scalar & 1) << (2*j+bit)
                high |= ((scalar >> 1) & 1) << (2*j+bit)
        rows.extend((low, high))
    rank = binary_rank(rows)
    assert rank % 2 == 0
    return rank // 2


def selection_controls():
    rng = Random(20261101)
    records = []
    for d in (2, 3, 4, 5, 9, 16, 47):
        for trial in range(8):
            m, D = 16, Q(37, 3)
            low = []
            ranks = []
            for a in range(d):
                r = rng.randrange(m//2+1)
                ranks.append(r)
                half = [Q(0)]*r + [D*Q(rng.randrange(1,101),100)
                                         for _ in range(m//2-r)]
                profile = sorted(half+[2*D-s for s in half])
                assert all(x+y == 2*D for x,y in zip(profile,reversed(profile)))
                assert sum(s == 0 for s in profile) == r
                low.extend(profile[:m//2])
            c = Q(sum(ranks),d*m)
            rho = min(Q(1,2),c+Q(1,d))
            count = d*m*rho
            assert sum(ranks) <= count <= Q(d*m,2)
            ordered = sorted(low)
            whole = count.numerator//count.denominator
            cost = sum(ordered[:whole],Q(0))
            if count != whole:
                cost += (count-whole)*ordered[whole]
            assert cost <= (count-sum(ranks))*D <= m*D
            records.append(dict(d=d,m=m,D=str(D),face_ranks=ranks,c=str(c),
                                rho=str(rho),retained_count=str(count),
                                retained_order=str(cost),budget=str(m*D)))
    comparisons = []
    for d in range(4,257):
        rho = min(Q(1,2),Q(1,4)+Q(1,d))
        assert d*rho*rho >= 1
        comparisons.append(dict(d=d,floor=str(rho),squared_margin=str(d*rho*rho-1)))
    thresholds = []
    for k, start in ((4,47),(5,223),(6,959)):
        c = Q(1,1 << (k-1))
        assert (start-1)*(c+Q(1,start-1))**2 < 1
        assert start*(c+Q(1,start))**2 >= 1
        thresholds.append(dict(subset_size=k,first_dimension=start,
                               floor=str(c+Q(1,start))))
    return dict(unequal_profile_controls=records,quarter_rank_comparisons=comparisons,
                higher_subset_thresholds=thresholds)


def fixture(values):
    coefficients = []
    for mask in range(256):
        selected = [power(v,8//((mask.bit_count()+1)&-(mask.bit_count()+1)))
                    for i,v in enumerate(values) if mask >> i & 1]
        k = mask.bit_count()+1
        r = k & -k
        coefficients.append({0:elementary(selected,r),
                             8//r:elementary(selected,r-1)})
    return coefficients


def pack_binary(polynomial):
    return sum(1 << (2*i) for i in range(polynomial.bit_length())
               if polynomial >> i & 1)


def fixture_controls():
    ones = fixture([1]*8)
    exponents = [8//((j+1)&-(j+1)) for j in range(9)]
    for mask, poly in enumerate(ones):
        assert poly[0] == 0 and poly[exponents[mask.bit_count()]] == 1
    histogram, blocks = Counter(), []
    for multiplicity, local, kinds in recursive_blocks([1 << e for e in exponents]):
        matrix = [[pack_binary(entry) for entry in row]
                  for row in multiplication_matrix(local)]
        n = len(local)
        orders, _ = WeightedSmithOracle(matrix).profile([0]*n,[0]*n)
        orders = list(orders)
        assert all(type(s) is int for s in orders)
        for s in orders:
            histogram[s] += multiplicity
        blocks.append(dict(type=list(kinds),multiplicity=multiplicity,
                           dimension=n,orders=orders))
    flat = sorted(histogram.elements())
    assert len(flat) == 256 and sum(flat) == 256*8
    assert all(a+b == 16 for a,b in zip(flat,reversed(flat)))
    assert sum(flat[:86]) == 258 and 3*258 > 8*86
    assert all(3*sum(flat[:m]) > 8*m for m in range(86,257))
    coefficients = [poly[0] for poly in fixture([2]+[1]*7)]
    matrix = multiplication_matrix(coefficients)
    support = 3
    cols = [i for i in range(256) if not i & support]
    rows = [support | i for i in cols]
    minor = [[matrix[i][j] for j in cols] for i in rows]
    coefficient = coefficients[support]
    assert coefficient == power(2,8)^power(1,8) == 2
    assert all(minor[i][i] == coefficient for i in range(64))
    assert all(minor[i][j] == 0 for i in range(64) for j in range(i+1,64))
    determinant = 1
    for i in range(64):
        determinant = MUL[determinant][minor[i][i]]
    rank = field_rank(matrix)
    assert determinant and rank >= 64 and rank <= 128
    return dict(degree=8,axes=9,all_one_coefficient_exponents=exponents,
                all_one_complete_profile=flat,symmetric_blocks=blocks,
                all_one_face_rank=0,first86_sum=258,
                nonbinary_other_values=[2]+[1]*7,nonbinary_face_rank=rank,
                supported_mask=support,minor_rows=rows,minor_columns=cols,
                literal_minor=minor,minor_determinant=determinant,
                generic_floor='13/36')


def report():
    sources = ['tools/check_multiaxis_face_rank_floor.py',
               'tools/cube_symmetric_exterior_blocks.py',
               'tools/recursive_symmetric_exterior_blocks.py',
               'tools/high_axis_plucker.py','tools/high_axis_cayley.py',
               'tools/weighted_smith_oracle_f4.py']
    return dict(scope='Paired-spectrum face-rank floor and exact specialization controls; no stronger factoring bound.',
                selection=selection_controls(),fixture=fixture_controls(),
                source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    result = report()
    if args.write:
        RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        assert result == json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],unequal_profiles=len(result['selection']['unequal_profile_controls']),
                          complete_small_block_orders=sum(b['dimension'] for b in result['fixture']['symmetric_blocks']),
                          nonbinary_face_rank=result['fixture']['nonbinary_face_rank'],
                          literal_minor_size=64,all_one_first86=258)))


if __name__ == '__main__':
    main()
