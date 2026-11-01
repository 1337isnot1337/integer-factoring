"""Path B, checker 2: the leading-power window as a p-adic lattice chain.

Lemma B2 (see note.md).  Let D<0 be a discriminant, p an odd prime with
(D/p)=1, and b_k the lift of sqrt(D) in Z/2p^kZ with b_k = D mod 2 and
b_k = b_1 mod p (Hensel).  Put
    L_k = {(u, y) in Z^2 : u = b_k*y (mod 2p^k)}   (index 2p^k in Z^2)
Then L_{k+1} is an index-p sublattice of L_k, and the reduced leading
coefficient of g_p^k is
    A_{p,k} = min{ (u^2 + |D| y^2) / (4 p^k) : (u, y) in L_k, (u, y) != 0 }.
So the window {A_{p,k} : K <= k < K+W} is the sequence of normalized first
minima of a descending chain of index-p lattices fixed by the p-adic digits of
sqrt(D) -- a p-adic analogue of the CFRAC denominators Q_k.

The script checks the identity exactly (composition vs. brute-force lattice
minimum vs. Gauss reduction of (p^k, b_k, c_k)), checks the chain inclusion,
and records finite smoothness statistics of the window (VERIFIED data only).
Run:  python3 check_padic_window.py
"""
import random
import sys
from math import gcd, isqrt
import bqf

random.seed(1101)


def hensel_sqrt(D, p, k, r1):
    """x with x^2 = D mod p^k and x = r1 mod p (p odd, p not dividing D)."""
    x, pk = r1 % p, p
    for _ in range(1, k):
        pk2 = pk * p
        # Newton step: x <- x - (x^2 - D)/(2x) mod pk2
        x = (x - (x * x - D) * pow(2 * x, -1, pk2)) % pk2
        pk = pk2
    assert (x * x - D) % pk == 0
    return x


def lift_b(D, p, k, r1):
    pk = p**k
    x = hensel_sqrt(D, p, k, r1)
    if (x - D) % 2:
        x += pk
    # x mod 2p^k, x = D mod 2, x^2 = D mod 4p^k
    assert (x * x - D) % (4 * pk) == 0
    return x % (2 * pk)


def lattice_min(D, p, k, b):
    pk = p**k
    m2 = 2 * pk
    best = None
    ymax = 1
    y = 1
    while True:
        u = (b * y) % m2
        if u > pk:
            u -= m2
        val = u * u + (-D) * y * y
        if best is None or val < best:
            best = val
        # stop once |D| y^2 alone exceeds the best value
        if (-D) * (y + 1) * (y + 1) > best:
            break
        y += 1
    best = min(best, 4 * pk * pk)  # y = 0, u = +-2p^k
    assert best % (4 * pk) == 0
    return best // (4 * pk), y


def smooth(n, Y):
    for p in bqf.primes_upto(Y):
        while n % p == 0:
            n //= p
    return n == 1


def main():
    ok = True
    checked = 0
    brute = 0
    for trial in range(40):
        while True:
            N = random.randrange(10**9, 10**11) | 1
            if N % 4 == 3 and not bqf.is_prime(N):
                break
        D = -N
        H = isqrt(-D // 3)
        for p in [q for q in bqf.primes_upto(60) if q > 2]:
            g = bqf.prime_form(D, p)
            if g is None:
                continue
            r1 = g[1] % p
            K = 1
            while p**K <= H:
                K += 1
            f = bqf.power(g, K)
            prev_b = None
            for k in range(K, K + 25):
                b = lift_b(D, p, k, r1)
                pk = p**k
                A_gauss = bqf.reduce_form((pk, b, (b * b - D) // (4 * pk)))[0]
                ok &= (f[0] == A_gauss)
                if k < K + 3:
                    # brute force enumerates y up to ~sqrt(p^k / sqrt|D|): only feasible for small k
                    A_lat, ystop = lattice_min(D, p, k, b)
                    ok &= (A_lat == A_gauss)
                    brute += 1
                if prev_b is not None:
                    # L_{k} subset L_{k-1}: b_k = b_{k-1} mod 2p^{k-1}
                    ok &= (b - prev_b) % (2 * p**(k - 1)) == 0
                prev_b = b
                f = bqf.compose(f, g)
                checked += 1
    print(f"[1] A_(p,k) from composition == Gauss reduction of (p^k,b_k,c_k) and L_(k+1) subset L_k: "
          f"{checked} (N,p,k) cases (k=K..K+24); brute-force lattice minimum also equal in {brute} cases "
          f"(k=K..K+2): {'OK' if ok else 'FAIL'}", flush=True)

    # [2] finite statistics: Y-smooth rate in windows vs. leading coefficients of random classes
    Y = 200
    rows = []
    for trial in range(6):
        while True:
            N = random.randrange(10**13, 10**14) | 1
            if N % 4 == 3 and not bqf.is_prime(N):
                break
        D = -N
        H = isqrt(-D // 3)
        gens = [bqf.prime_form(D, p) for p in bqf.primes_upto(Y) if p > 2]
        gens = [g for g in gens if g is not None]
        win_hits = win_tot = 0
        first_hits = []
        for g in gens[:8]:
            p = g[0]
            K = 1
            while p**K <= H:
                K += 1
            f = bqf.power(g, K)
            first = None
            for k in range(K, K + 400):
                win_tot += 1
                if smooth(f[0], Y):
                    win_hits += 1
                    if first is None:
                        first = k - K
                f = bqf.compose(f, g)
            first_hits.append(first)
        rnd_hits = rnd_tot = 0
        for _ in range(3200):
            f = bqf.identity(D)
            for g in random.sample(gens, 6):
                f = bqf.compose(f, bqf.power(g, random.randrange(1, 10**6)))
            rnd_tot += 1
            rnd_hits += smooth(f[0], Y)
        rows.append((N, win_hits, win_tot, rnd_hits, rnd_tot, first_hits))
        print(f"[2] N={N}: window smooth rate {win_hits}/{win_tot}={win_hits/win_tot:.4f}; "
              f"random-class leading coeff rate {rnd_hits}/{rnd_tot}={rnd_hits/rnd_tot:.4f}; "
              f"first hit offsets {first_hits}", flush=True)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
