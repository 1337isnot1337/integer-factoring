#!/usr/bin/env python3
"""Standard-library proof checks for supplied tiny polynomial factorizations."""
import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import product
import json
import lzma
from math import gcd, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def multiply(a, b):
    c = [0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):
            c[i+j] += u*v
    return trim(c)


def divisors(n):
    assert n > 0
    out = []
    for a in range(1,isqrt(n)+1):
        if n % a == 0:
            out.append(a)
            if a*a != n:
                out.append(n//a)
    return sorted(out)


def evaluate(a, radix):
    return sum(c*radix**i for i,c in enumerate(a))


def in_box(a, box):
    return len(a) <= len(box) and all(lo <= (a[i] if i<len(a) else 0) <= hi
                                    for i,(lo,hi) in enumerate(box))


def irreducible_tiny(a):
    """Over Q, degree two/three is irreducible iff it has no rational root."""
    assert a == trim(a) and 1 <= len(a)-1 <= 3 and a[-1] > 0 and gcd(*a) == 1
    if len(a) == 2:
        return
    assert a[0] != 0
    degree = len(a)-1
    for u,v,sign in product(divisors(abs(a[0])),divisors(a[-1]),(1,-1)):
        assert sum(c*(sign*u)**i*v**(degree-i) for i,c in enumerate(a)) != 0


def compute():
    path = ROOT/'evidence/polynomial-digit-fibres.factorizations.json.xz'
    compressed = path.read_bytes()
    raw = lzma.decompress(compressed)
    acquisition = json.loads(raw)
    certificates = {}
    factors_checked = 0
    for W,c,fs in acquisition['factorizations']:
        W = tuple(W)
        assert W == trim(W) and any(W) and W not in certificates
        assert abs(c) == gcd(*W) and c != 0
        factors = [(tuple(a),e) for a,e in fs]
        assert len(set(a for a,e in factors)) == len(factors)
        reconstruction = (c,)
        for a,e in factors:
            irreducible_tiny(a)
            assert isinstance(e,int) and e >= 1
            factors_checked += 1
            for _ in range(e):
                reconstruction = multiply(reconstruction,a)
        assert reconstruction == W
        certificates[W] = c,factors
    seen = set()

    @cache
    def candidates(W):
        seen.add(W)
        c,fs = certificates[W]
        result = []
        for selected in product(*(range(e+1) for a,e in fs)):
            P,Q = (1,),(1,)
            for (a,e),chosen in zip(fs,selected):
                for _ in range(chosen): P = multiply(P,a)
                for _ in range(e-chosen): Q = multiply(Q,a)
            for a,sign in product(divisors(abs(c)),(1,-1)):
                a *= sign
                A,B = tuple(a*v for v in P),tuple((c//a)*v for v in Q)
                assert multiply(A,B) == W
                result.append((A,B))
        assert len(result) == len(set(result))
        assert len(result) <= 2*len(divisors(abs(c)))*2**(len(W)-1)
        return result

    def count(W,A,B):
        if not any(W):
            a0,b0 = in_box((0,),A),in_box((0,),B)
            return prod(hi-lo+1 for lo,hi in A)*b0+prod(hi-lo+1 for lo,hi in B)*a0-a0*b0
        return sum(in_box(a,A) and in_box(b,B) for a,b in candidates(W))

    records = []
    for radix,length in ((2,4),(4,3),(8,2)):
        words = list(product(range(radix),repeat=length))
        hist = [Counter(),Counter(),Counter()]
        content_hist = Counter()
        rough,carry_controls = 0,0
        for A,B in product(words,repeat=2):
            W = multiply(A,B)
            N = evaluate(A,radix)*evaluate(B,radix)
            assert evaluate(W,radix) == N
            hist[0][W] += 1
            if A[-1] == 0: hist[1][W] += 1
            if A[-1] == radix-1: hist[2][W] += 1
            carry,number,cs,ns = 0,N,[],[]
            for j in range(2*length):
                value = (W[j] if j<len(W) else 0)+carry
                digit,number = number % radix,number//radix
                assert value % radix == digit
                carry = value//radix
                assert 0 <= carry <= length*(radix-1)
                cs.append(carry); ns.append(digit)
            assert number == carry == 0
            assert trim([ns[j]+radix*cs[j]-(cs[j-1] if j else 0) for j in range(2*length)]) == W
            carry_controls += 1
            if N > 0:
                content_hist[gcd(*W)] += 1
                if all(N % a for a in range(2,radix)):
                    assert gcd(*A) == gcd(*B) == gcd(*W) == 1
                    rough += 1
        box = ((0,radix-1),)*length
        for i,H in enumerate(hist):
            A = box if i==0 else box[:-1]+(((0,0) if i==1 else (radix-1,radix-1)),)
            for W in hist[0]:
                assert count(W,A,box) == H[W]
        records.append(dict(radix=radix,operand_digits=length,literal_pairs=len(words)**2,
            distinct_product_polynomials=len(hist[0]),three_box_fibre_checks=3*len(hist[0]),
            carry_identity_controls=carry_controls,rough_pairs_primitive_controls=rough,
            nonzero_fibre_content_histogram=dict(sorted(content_hist.items())),zero_fibre=hist[0][(0,)]))
    signed = list(product(range(-1,2),repeat=2))
    signed_hist = Counter(multiply(A,B) for A,B in product(signed,repeat=2))
    for W,amount in signed_hist.items():
        assert count(W,((-1,1),)*2,((-1,1),)*2) == amount
    assert seen == set(certificates)
    names = ['tools/check_polynomial_digit_fibres.py','research/polynomial-digit-fibres.md']
    return dict(scope='Classical known-polynomial fibre acquisition, exact carry identity and roughness/primitive consequence. No fast carry selection, complete-bound improvement or LLL implementation.',
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in names},
        acquisition=acquisition['acquisition'],factorizations=len(certificates),
        primitive_irreducible_factor_controls=factors_checked,
        factor_certificate_raw_sha256=sha256(raw).hexdigest(),factor_certificate_archive_sha256=sha256(compressed).hexdigest(),
        sizes=records,signed_literal_pairs=len(signed)**2,signed_fibres_checked=len(signed_hist))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(compute()))
    expected = ROOT/'evidence/polynomial-digit-fibres.json'
    if args.write:
        expected.write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert result == json.loads(expected.read_text())
    print('Polynomial digit fibres: exact factor certificates, box counts, carries and rough primitive controls PASS')
