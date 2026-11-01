#!/usr/bin/env python3
"""Exact checks for the shift/cross-count lemma behind Theorem T1.

Setting: x bits 0..n-1, z bits 0..2n-1, each coloured A/B (side of a
balanced separation).  Cross(j) = #{i<n : colour(x_i) != colour(z_{i+j})}.
Sum_j Cross(j) = Sum_i O_i, O_i = #z in window [i,i+n) of colour != colour(x_i),
so min over x-colourings of Sum_j Cross(j) = Sum_i min(#A, #B in window i).

(L1) Lemma: if 6|n and both z colour classes have >= 19n/30 elements then
     Sum_i min(...) >= n^2/120.  Checked exhaustively over all z colourings.
(L2) Exact min over balanced (x,z) colourings of max_j Cross(j), small n.
(L3) Odd variant: y = 2^j+1, x < 2^j: Cross'(j) = #{1<=i<j: x_i, z_i, z_{i+j}
     not monochromatic}; exact min of max_j Cross'(j), small n.
(L4) Arithmetic facts used by the restriction: x*2^j and x*(2^j+1) (x<2^j)
     have the claimed bit patterns (no carries), all x for small n.
"""
import math, itertools, sys, time

def L1(n):
    need = math.ceil(19*n/30); full = (1 << n) - 1; worst = None
    for d in range(1 << (2*n)):
        a = d.bit_count()
        if a < need or 2*n - a < need: continue
        s = 0
        for i in range(n):
            w = ((d >> i) & full).bit_count(); s += min(w, n - w)
        if worst is None or s < worst: worst = s
    return worst

def L2(n, frac=19/30):
    need = math.ceil(frac*n); best = None
    for d in range(1 << (2*n)):
        a = d.bit_count()
        if a < need or 2*n - a < need: continue
        for c in range(1 << n):
            m = max(sum(((c >> i) & 1) != ((d >> (i+j)) & 1) for i in range(n)) for j in range(n))
            if best is None or m < best: best = m
            if best == 0: return 0
    return best

def L3(n, frac=19/30):
    need = math.ceil(frac*n); best = None
    for d in range(1 << (2*n)):
        a = d.bit_count()
        if a < need or 2*n - a < need: continue
        for c in range(1 << n):
            m = 0
            for j in range(2, n):
                cnt = 0
                for i in range(1, j):
                    cols = {(c >> i) & 1, (d >> i) & 1, (d >> (i+j)) & 1}
                    cnt += len(cols) > 1
                m = max(m, cnt)
            if best is None or m < best: best = m
    return best

def L4(n):
    for j in range(n):
        for x in range(1 << n):
            assert x * (1 << j) == x << j
        for x in range(1 << j):
            p = x * ((1 << j) + 1)
            for i in range(j):
                b = (x >> i) & 1
                assert (p >> i) & 1 == b and (p >> (i+j)) & 1 == b
            assert p < (1 << (2*j)) or j == 0
    return True

if __name__ == '__main__':
    ok = True
    for n in (6, 12):
        t = time.time(); w = L1(n); bound = n*n/120
        print(f"L1 n={n:2d}: min Sum_i min(#A,#B) over balanced z = {w}  "
              f"(lemma bound n^2/120 = {bound:.2f}; ratio {w/(n*n):.4f})  {time.time()-t:.1f}s", flush=True)
        ok &= w >= bound
    for n in (3, 4, 5, 6):
        t = time.time()
        print(f"L2 n={n}: min over balanced colourings of max_j Cross(j) = {L2(n)}  ({time.time()-t:.1f}s)", flush=True)
    for n in (4, 5, 6):
        t = time.time()
        print(f"L3 n={n}: odd variant min max_j Cross'(j) = {L3(n)}  ({time.time()-t:.1f}s)", flush=True)
    for n in (2, 4, 6, 8):
        ok &= L4(n)
    print("L4 bit-pattern facts: OK")
    print("PASS" if ok else "FAIL")
