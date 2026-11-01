#!/usr/bin/env python3
"""Path C: explicit N on which one target's power window misses for ~2W steps.

Construction (proof in note.md, claim P3):  pick split prime p, exponents K>W>=1,
a prime A > Y with (3/4)p^(K-W-2) <= A < (3/4)p^(K-W), and odd x, 0<x<A, p∤x.
Put Delta = x^2 - 4 p^(K+W) A  (so N = -Delta = 3 mod 4 and the repo uses Delta=-N).
Then K_p = K, and for 0 <= j < W the reduced leading coefficient of g_p^(K+j)
is A*p^(W-j); at j = W it is A; for j = W+i it is A*p^i while that is reduced.
None of these is Y-smooth because A > Y is prime.

This script builds such N (optionally requiring N to be a semiprime, factored
by PARI only to *select* x), then runs the repo's exact orbit for target p and
checks each predicted leading coefficient.

Usage: python3 adversary.py --p 3 --bits 64 96 128 160 --cutoff 4093 [--semiprime]
"""
import argparse, json, math, os, sys
sys.path.insert(0, os.environ.get('IF_TOOLS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))
from math import isqrt, prod
from arithmetic_two_group import power
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above
from leading_power_smooth_relations import compose_with_prime_form
from window_scan import is_probable_prime


def build(p, bits, cutoff, semiprime, max_tries=4000):
    # b ~ (K+W)log2 p + log2 A + 2 ~ 2K log2 p ; choose K, then max W with A > cutoff.
    K = max(3, round((bits - 2) / (2 * math.log2(p))))
    W = K - 2
    while W >= 1 and (3 * p ** (K - W)) // 4 <= cutoff:  # need room for prime A > Y
        W -= 1
    lo = max(cutoff + 1, -(-3 * p ** (K - W - 2) // 4))
    hi = (3 * p ** (K - W)) // 4 - 1              # A < (3/4) p^(K-W)
    A = lo | 1
    while not (is_probable_prime(A) and A % p):
        A += 2
    assert A <= hi, (p, bits, K, W, lo, hi)
    M = p ** (K + W)
    fac = None
    if semiprime:
        import cypari2
        pari = cypari2.Pari()
    for x in range(1, min(A, 2 * max_tries + 1), 2):
        if x % p == 0:
            continue
        n = 4 * M * A - x * x
        if semiprime:
            f = pari.factor(n)
            ps = [int(v) for v in f[0]]
            es = [int(v) for v in f[1]]
            if len(ps) != 2 or es != [1, 1] or min(ps).bit_length() < n.bit_length() // 4:
                continue
            fac = ps
        return dict(p=p, K=K, W=W, A=A, x=x, N=n, factors=fac)
    raise RuntimeError('no x found')


def verify(c, cutoff):
    p, K, W, A, n = c['p'], c['K'], c['W'], c['A'], c['N']
    assert n % 4 == 3
    delta = -n
    identity = (1, 1, (1 - delta) // 4)
    height = isqrt(n // 3)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    assert p in records, 'p must be a base (split) prime'
    primes = sorted(records)
    P = prod(primes)
    sexp = height.bit_length()
    Kp = least_power_above(p, height)
    form = power(records[p][1], Kp, identity, compose_forms_hnf)
    predicted, observed, j = [], [], 0
    first_hit = None
    limit = 2 * W + 8
    while True:
        a = form[0]
        observed.append(a)
        smooth = pow(P, sexp, a) == 0
        if smooth and first_hit is None:
            first_hit = j
            if j >= limit:
                break
        if j < W:
            predicted.append(A * p ** (W - j))
        elif A * p ** (j - W) * 2 < isqrt(n):         # norm < sqrt|Delta|/2 => reduced
            predicted.append(A * p ** (j - W))
        else:
            predicted.append(None)
        j += 1
        if first_hit is not None and j >= limit:
            break
        if j > limit + 20000:   # tail only for information; censored as null
            break
        form = compose_with_prime_form(form, p, *records[p])
    checks = [(i, pr, ob) for i, (pr, ob) in enumerate(zip(predicted, observed)) if pr is not None]
    ok = all(pr == ob for _, pr, ob in checks)
    guaranteed = len(checks)
    return dict(bits=n.bit_length(), p=p, K=K, K_p=Kp, W=W, A=A, x=c['x'], N=str(n),
                factors=[str(f) for f in c['factors']] if c['factors'] else None,
                base_size=len(primes), cutoff=cutoff,
                predicted_misses=guaranteed, predictions_match=ok, K_p_equals_K=(Kp == K),
                observed_first_hit=first_hit,
                lower_bound_formula=round(n.bit_length() / math.log2(p) - 2 * math.log(cutoff, p), 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p', type=int, nargs='+', default=[3])
    ap.add_argument('--bits', type=int, nargs='+', default=[64, 96, 128, 160])
    ap.add_argument('--cutoff', type=int, default=4093)
    ap.add_argument('--semiprime', action='store_true')
    ap.add_argument('--out')
    a = ap.parse_args()
    rows = []
    for p in a.p:
        for b in a.bits:
            c = build(p, b, a.cutoff, a.semiprime)
            r = verify(c, a.cutoff)
            assert r['predictions_match'] and r['K_p_equals_K'], r
            rows.append(r)
            print(json.dumps({k: r[k] for k in ('bits', 'p', 'K', 'K_p', 'W', 'A', 'x',
                  'predicted_misses', 'predictions_match', 'observed_first_hit',
                  'lower_bound_formula', 'factors')}), flush=True)
    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(rows, fh, indent=1)


if __name__ == '__main__':
    main()
