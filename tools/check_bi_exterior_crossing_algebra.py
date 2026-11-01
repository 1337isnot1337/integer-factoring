#!/usr/bin/env python3
"""Exact controls for the bi-exterior crossing representation."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
from random import Random


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def transpose(rows, size):
    return [sum(((rows[i] >> j) & 1) << i for i in range(size)) for j in range(size)]


def multiply(a, b):
    columns = transpose(b, len(b))
    return [sum(((row & col).bit_count() % 2) << j for j, col in enumerate(columns)) for row in a]


def reshuffle(rows, d):
    answer = [0] * (d*d)
    for x, y, u, v in product(range(d), repeat=4):
        answer[x*d+u] |= ((rows[x*d+y] >> (u*d+v)) & 1) << (y*d+v)
    return answer


LOCAL = ([1,2,4,8], [6,8,8,0], [0,1,1,6], [0,6,6,0])


def local_controls():
    identity, a, b, n = LOCAL
    zero = [0]*4
    assert multiply(a,a) == multiply(b,b) == zero
    assert multiply(a,b) == multiply(b,a) == n
    assert rank([2,1,8,6]) == 4
    r, j = reshuffle(identity,2), reshuffle(n,2)
    aa, bb = multiply(j,reshuffle(a,2)), multiply(j,reshuffle(b,2))
    assert (r,j,aa,bb) == ([9,0,0,9],[8,4,2,1],[2,0,9,2],[4,9,0,4])
    assert multiply(j,j) == identity and multiply(j,r) == r
    assert multiply(aa,aa) == multiply(bb,bb) == zero
    assert multiply(aa,bb) == multiply(bb,aa) == r
    assert rank([1,4,2,9]) == 4
    assert [x^y for x,y in zip(identity,n)] == [1,4,2,8]


def support(polynomial):
    result = []
    while polynomial:
        bit = polynomial & -polynomial
        result.append(bit.bit_length()-1)
        polynomial ^= bit
    return result


def exterior_product(a, b):
    result = 0
    for x in support(a):
        for y in support(b):
            if not x & y:
                result ^= 1 << (x | y)
    return result


def exterior_rank(width, polynomial):
    terms = support(polynomial)
    return rank([sum(1 << (a|b) for a in terms if not a & b) for b in range(1 << (2*width))])


def dual(width, polynomial):
    result = 0
    for monomial in support(polynomial):
        target = 0
        for j in range(width):
            digit = (monomial >> (2*j)) & 3
            target |= (3,1,2,0)[digit] << (2*j)
        result ^= 1 << target
    return result


def crossing_polynomial(width):
    result = 1
    for j in range(width):
        result = exterior_product(result, 1 | (1 << (3 << (2*j))))
    return result


def tensor(width, polynomial):
    d = 1 << width
    rows = [0] * (d*d)
    terms = support(polynomial)
    for x,y,u,v in product(range(d),repeat=4):
        entry = 0
        for monomial in terms:
            value = 1
            for j in range(width):
                matrix = LOCAL[(monomial >> (2*j)) & 3]
                row = 2*((x>>j)&1)+((y>>j)&1)
                column = 2*((u>>j)&1)+((v>>j)&1)
                value &= (matrix[row] >> column) & 1
            entry ^= value
        rows[x*d+y] |= entry << (u*d+v)
    return rows


def substitute(polynomial, linear_forms):
    result = 0
    for monomial in support(polynomial):
        value = 1
        for j in range(len(linear_forms)):
            if monomial >> j & 1:
                value = exterior_product(value, sum(1 << (1 << k) for k in range(len(linear_forms)) if linear_forms[j] >> k & 1))
        result ^= value
    return result


def boolean_complement_controls():
    """Independent Boolean-selector lift of the same six-term polynomial.

    Write N_i=1+g_i in F_2[(F_2)^4].  The coefficient of g_w in
    N_B is one exactly when w is a subset of B.  Thus evaluating
    sum_{|A|=2} S_A N_{A^c} at Boolean S gives the kernel below.
    The two crossing cuts split into convolution blocks along the
    rows and columns of this kernel, respectively.
    """
    d = 16
    kernel = [sum(((sum(1 for a in range(d)
                             if a.bit_count() == 2 and a & s == a
                             and not a & w) % 2) << w)
                  for w in range(d)) for s in range(d)]
    columns = transpose(kernel, d)

    def convolution_rank(word):
        return rank([sum(((word >> (a ^ b)) & 1) << b
                         for b in range(d)) for a in range(d)])

    row_blocks = [convolution_rank(word) for word in kernel]
    column_blocks = [convolution_rank(word) for word in columns]
    assert sum(row_blocks) == sum(column_blocks) == 46
    assert sorted(row_blocks) == sorted(column_blocks) == [0]*5 + [4]*10 + [6]
    return {'row_rank':sum(row_blocks), 'column_rank':sum(column_blocks),
            'row_block_ranks':row_blocks, 'column_block_ranks':column_blocks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    local_controls()
    rng = Random(20261030075)
    dense = changes = 0
    for width in range(1,4):
        d = 1 << width
        endpoint = tensor(width,crossing_polynomial(width))
        assert endpoint == [1 << (y*d+x) for x in range(d) for y in range(d)]
        for _ in range(10):
            f = rng.getrandbits(1 << (2*width))
            rows = tensor(width,f)
            assert rank(rows) == exterior_rank(width,f)
            assert rank(reshuffle(rows,d)) == exterior_rank(width,dual(width,f))
            dense += 1
            while True:
                forms = [rng.randrange(1 << (2*width)) for _ in range(2*width)]
                if rank(forms) == 2*width:
                    break
            g = substitute(f,forms)
            assert exterior_rank(width,g) == exterior_rank(width,f)
            assert exterior_rank(width,dual(width,g)) == exterior_rank(width,dual(width,f))
            changes += 1
    f = 0
    for subset in range(16):
        if subset.bit_count() == 2:
            monomial = sum((1 if subset >> j & 1 else 2) << (2*j) for j in range(4))
            f ^= 1 << monomial
    assert dual(4,f) == f and exterior_rank(4,f) == 36
    rows = tensor(4,f)
    assert rank(rows) == rank(reshuffle(rows,16)) == 36
    assert 36*36 < 16*16*6
    boolean_control = boolean_complement_controls()
    result = {'status':'PASS','scope':'Exact characteristic-two bi-exterior algebra; no improved factoring exponent',
              'dense_random_tensor_controls':dense,'invertible_substitution_controls':changes,
              'cheap_direction':{'labels':16,'row_rank':36,'column_rank':36,'coefficient_rank':6},
              'boolean_selector_lift_of_same_polynomial':boolean_control,
              'endpoint_dimensions':[2,4,8],
              'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
