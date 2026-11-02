#!/usr/bin/env python3
"""Exact leading-jet controls for one homogeneous additive-charge cube."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import argparse
import json
import cube_homogeneous_gaussian as cube

MODULUS = 0x11B
LENGTH = 5
ELL2 = frozenset(((2,0,0),(0,2,0),(0,0,2)))
QUARTIC = frozenset(((8,4,4),(4,8,4),(4,4,8),
                     (6,5,5),(5,6,5),(5,5,6)))


def multiply(a, b):
    out = 0
    while b:
        if b & 1: out ^= a
        a <<= 1
        if a & 256: a ^= MODULUS
        b >>= 1
    return out


def power(a, exponent):
    out = 1
    for _ in range(exponent): out = multiply(out, a)
    return out


def evaluate(poly):
    out = [0]*LENGTH
    for a,b,c in poly:
        if a < LENGTH: out[a] ^= multiply(power(2,b), power(3,c))
    return out


def polynomial_product(a,b):
    out = [0]*LENGTH
    for i,x in enumerate(a):
        for j,y in enumerate(b[:LENGTH-i]): out[i+j] ^= multiply(x,y)
    return out


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            i = row.bit_length()-1
            if i in pivots: row ^= pivots[i]
            else:
                pivots[i] = row
                break
    return len(pivots)


def restriction_blocks(blocks):
    # Multiplication by each F256 coefficient is independently represented over F2.
    values = {a for block in blocks for row in block for a in row}
    planes = {a:[sum(((multiply(a,1<<k)>>bit)&1)<<k for k in range(8))
                 for bit in range(8)] for a in values}
    return [[sum(planes[a][bit]<<(8*j) for j,a in enumerate(row))
             for row in block for bit in range(8)] for block in blocks]


def toeplitz_ranks(blocks):
    binary = restriction_blocks(blocks)
    width = 8*len(blocks[0]);ranks = [0]
    for h in range(1,LENGTH+1):
        rows = [sum(binary[d][i] << (width*(r-d)) for d in range(r+1))
                for r in range(h) for i in range(width)]
        rank = binary_rank(rows)
        assert rank % 8 == 0
        ranks.append(rank//8)
    return ranks


def formal_controls(kernel):
    gaussian = cube.gaussian_tensor(cube.kernel())
    for mask,entry in enumerate(gaussian):
        paired = all(bool(mask>>v&1)==bool(mask>>(v^7)&1) for v in range(8))
        assert (len(entry)%2)==paired
    assert QUARTIC and all(sum(e)==16 for e in QUARTIC) and len(QUARTIC)%2==0
    assert all(sum(e)==8 for row in kernel for entry in row for e in entry)
    records = []
    for axis in range(3):
        left,right = cube.axis_parts(axis)
        def leading(poly,degree):
            return frozenset(tuple(v if j!=axis else 0 for j,v in enumerate(e))
                             for e in poly if e[axis]==degree)
        a = leading(kernel[left[0]][right[0]],2)
        assert all(cube.valuation(kernel[i][j],axis)==2 and leading(kernel[i][j],2)==a
                   for i in left for j in right)
        assert cube.valuation(QUARTIC,axis)==4 and leading(QUARTIC,4)==cube.multiply(a,a)
        spatial = [[kernel[i][j] for j in right] for i in left]
        two = [cube.minor(spatial,r,c) for r in combinations(range(4),2)
               for c in combinations(range(4),2)]
        three = [cube.minor(spatial,r,c) for r in combinations(range(4),3)
                 for c in combinations(range(4),3)]
        assert min(cube.valuation(p,axis) for p in two if p)==6
        assert all(not p or cube.valuation(p,axis)>=11 for p in three)
        records.append(dict(axis=axis,cross_entries=16,two_minors=36,three_minors=16,
                            first_spatial_orders=[2,4],later_spatial_orders_at_least=5))
    return records


def low_jet_controls(kernel):
    left,right = cube.axis_parts(0)
    words = []
    for word in product(range(8),repeat=4):
        colors = [a for a in word if a];charge = 0
        for a in colors: charge ^= a
        if charge==1 and 1<=len(colors)<=3 and len(set(colors))==len(colors): words.append(word)
    assert Counter(sum(a!=0 for a in w) for w in words)=={1:4,2:36,3:96}
    n = len(words)
    pair = [[evaluate(kernel[a][b]) for b in right] for a in left]
    quartic = evaluate(QUARTIC)
    assert quartic[4]==multiply(pair[0][0][2],pair[0][0][2]) != 0
    records = {}
    for fusion in (False,True):
        blocks = [[[0]*n for _ in range(n)] for _ in range(LENGTH)]
        for i,a in enumerate(words):
            aa = [(j,c) for j,c in enumerate(a) if c]
            for j,b in enumerate(words):
                bb = [(k,c) for k,c in enumerate(b) if c]
                entry = [0]*LENGTH
                if len(aa)==len(bb)<=2 and sorted(c for _,c in aa)==sorted(c for _,c in bb):
                    entry[0] = 1
                    for place,color in aa:
                        other = next(k for k,c in bb if c==color)
                        entry = polynomial_product(entry,pair[place][other])
                colors = [c for _,c in aa+bb]
                if fusion and len(colors)==len(set(colors))==4:
                    assert colors[0]^colors[1]^colors[2]^colors[3]==0
                    entry = [a^b for a,b in zip(entry,quartic)]
                for degree,value in enumerate(entry): blocks[degree][i][j] = value
        ranks = toeplitz_ranks(blocks)
        cumulative = [ranks[k+1]-ranks[k] for k in range(LENGTH)]
        profile = {k:cumulative[k]-(cumulative[k-1] if k else 0) for k in range(LENGTH)}
        profile = {k:v for k,v in profile.items() if v}
        assert profile==({2:1,4:2} if fusion else {2:1,4:4})
        records['fusion' if fusion else 'gaussian'] = dict(ranks=ranks,profile=profile)
    pairs = list(combinations(range(1,8),2))
    fused = [sum(((a^b)==(c^d))<<j for j,(c,d) in enumerate(pairs)) for a,b in pairs]
    assert binary_rank(fused)==7
    return dict(block_size=n,coefficient_field_modulus=MODULUS,y=2,z=3,
                first_layers=records,color_pair_ranks=[len(pairs),binary_rank(fused)],
                full_prefix_from_sector_proof=dict(gaussian=[[0,1],[2,7],[4,28]],
                                                  fusion=[[0,1],[2,7],[4,14]]))


def main():
    parser = argparse.ArgumentParser();parser.add_argument('--output',type=Path);args = parser.parse_args()
    assert all(power(a,255)==1 for a in range(1,256))
    kernel = [[cube.multiply(ELL2,a) for a in row] for row in cube.kernel()]
    result = dict(scope='Formal kernel/endpoint/minor identities and one q8 low-jet charge block, independently restricted from F256 to binary Toeplitz ranks. The generic prefix uses the ordinary sector argument; no full profile, favorable normalized pressure or stronger factoring exponent is established.',
                  physical_degree_per_occupied_vertex=4,quartic_terms=sorted(QUARTIC),
                  formal=formal_controls(kernel),low_jets=low_jet_controls(kernel),
                  source_sha256={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in ('check_quartic_charge_fusion.py','cube_homogeneous_gaussian.py')})
    if args.output: args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',axes=3,charge_block_size=136,
                          gaussian_prefix=[[0,1],[2,7],[4,28]],fusion_prefix=[[0,1],[2,7],[4,14]],
                          full_profile_known=False,factoring_exponent='1/2; unchanged')))


if __name__=='__main__': main()
