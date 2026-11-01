#!/usr/bin/env python3
"""Path C empirics: first-hit windows of the prime-form power search.

For each N (deterministic seeded semiprime) and each sampled target prime p,
this walks the exact same orbit the repo selector walks,
    Q_{p,k} = red(g_p^k),  k = K_p, K_p+1, ...,  K_p = least k with p^k > H,
using the repo's own composition and base code, and records the offsets j
(0-based, so w_p = j+1 in the repo's notation) at which the reduced leading
coefficient A is base-smooth (P^bits(H) == 0 mod A, P = product of the base
primes, exactly the repo filter).  It keeps scanning until `hits` smooth
offsets are seen or `cap` steps are spent (censored).

Usage:  python3 window_scan.py --bits 60 80 ... --per-size 8 --schedule fixed4093 \
            --targets 16 --hits 3 --cap 3000000 --jobs 8 --out data/scan.jsonl
"""
import argparse, hashlib, json, math, os, sys, time
from multiprocessing import Pool

REPO_TOOLS = os.environ.get('IF_TOOLS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, REPO_TOOLS)
from math import gcd, isqrt, prod                                  # noqa: E402
from arithmetic_two_group import power                              # noqa: E402
from binary_quadratic_forms import reduce_form                      # noqa: E402
from prime_form_factor_base import build_prime_form_factor_base     # noqa: E402
from quadratic_ideal_hnf import compose_forms_hnf                   # noqa: E402
from weighted_norm_relations import least_power_above               # noqa: E402
from leading_power_smooth_relations import compose_with_prime_form  # noqa: E402


def is_probable_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def seeded_prime(bits, label):
    h = int.from_bytes(hashlib.sha256(label.encode()).digest() * 4, 'big')
    x = (h % (1 << (bits - 2))) | (1 << (bits - 1)) | (1 << (bits - 2)) | 1
    while not is_probable_prime(x):
        x += 2
    assert x.bit_length() == bits
    return x


def seeded_semiprime(bits, index):
    a = bits // 2
    p = seeded_prime(a, f'pathC|{bits}|{index}|P')
    q = seeded_prime(bits - a, f'pathC|{bits}|{index}|Q')
    n = p * q
    return n, p, q


def cutoff_for(schedule, bits):
    if schedule.startswith('fixed'):
        return int(schedule[5:])
    if schedule == 'bsq':          # Y = b^2  (polylog schedule, c = 2)
        return bits * bits
    raise ValueError(schedule)


def scan_one(job):
    bits, index, schedule, ntargets, hits_wanted, cap = job
    t0 = time.time()
    n, p_fac, q_fac = seeded_semiprime(bits, index)
    delta = -n if n % 4 == 3 else -4 * n
    parity = delta % 2
    identity = (1, parity, (parity - delta) // 4)
    height = isqrt((-delta) // 3)
    cutoff = cutoff_for(schedule, bits)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    if n % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (n + 1) // 2)))
    primes = sorted(records)
    prime_product = prod(primes)
    sexp = height.bit_length()
    m = len(primes)
    # deterministic spread of targets across the base (always includes smallest)
    if ntargets >= m:
        targets = primes
    else:
        targets = sorted({primes[(i * (m - 1)) // (ntargets - 1)] for i in range(ntargets)})
    out_targets = []
    steps_total = 0
    for p in targets:
        first = least_power_above(p, height)
        form = power(records[p][1], first, identity, compose_forms_hnf)
        hit_offsets = []
        j = 0
        gcd_hit = None
        while j < cap and len(hit_offsets) < hits_wanted:
            a = form[0]
            g = gcd(a, n)
            if 1 < g < n:
                gcd_hit = j
                hit_offsets.append(j)
            elif pow(prime_product, sexp, a) == 0:
                hit_offsets.append(j)
            j += 1
            if len(hit_offsets) < hits_wanted and j < cap:
                form = compose_with_prime_form(form, p, *records[p])
        steps_total += j
        out_targets.append({'p': p, 'K_p': first, 'hits': hit_offsets,
                            'scanned': j, 'censored': len(hit_offsets) < hits_wanted,
                            'gcd_hit': gcd_hit})
    return {'bits': bits, 'index': index, 'schedule': schedule, 'N': str(n),
            'P': str(p_fac), 'Q': str(q_fac), 'delta': str(delta),
            'H_bits': height.bit_length(), 'cutoff': cutoff, 'base_size': m,
            'targets': out_targets, 'steps': steps_total,
            'seconds': round(time.time() - t0, 2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bits', type=int, nargs='+', required=True)
    ap.add_argument('--per-size', type=int, default=8)
    ap.add_argument('--first-index', type=int, default=0)
    ap.add_argument('--schedule', default='fixed4093')
    ap.add_argument('--targets', type=int, default=16)
    ap.add_argument('--hits', type=int, default=3)
    ap.add_argument('--cap', type=int, default=3_000_000)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    jobs = [(b, i, a.schedule, a.targets, a.hits, a.cap)
            for b in a.bits for i in range(a.first_index, a.first_index + a.per_size)]
    with Pool(a.jobs) as pool, open(a.out, 'a') as fh:
        for rec in pool.imap_unordered(scan_one, jobs):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            fs = [t['hits'][0] if t['hits'] else None for t in rec['targets']]
            print(rec['bits'], rec['index'], rec['schedule'], 'm=%d' % rec['base_size'],
                  'steps=%d' % rec['steps'], '%.1fs' % rec['seconds'],
                  'max_first=%s' % max((f for f in fs if f is not None), default=None),
                  'censored=%d' % sum(t['censored'] for t in rec['targets']), flush=True)


if __name__ == '__main__':
    main()
