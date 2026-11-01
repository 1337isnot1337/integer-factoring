"""Exact class-group arithmetic for the duality experiments (Delta < 0).

Group structure and generators come from PARI quadclassunit (GRH-conditional for
large |Delta|; certified by an exhaustive reduced-form count when |Delta| <= CERT_MAX).
Every class is enumerated as a mixed-radix product of the PARI generators, so each
reduced form gets an exact discrete-log vector x in prod Z/c_i, and every character
is chi_k(x) = e(sum_i k_i x_i / c_i).  The factor base, prime forms, K_p and the
smoothness test are the repo's own (tools/ of 1337isnot1337/integer-factoring),
found relative to this file (override with IF_TOOLS).
"""
import math, os, sys
from math import isqrt, gcd, prod
HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.environ.get('IF_TOOLS', os.path.join(HERE, '..'))
sys.path.insert(0, TOOLS)
from binary_quadratic_forms import reduce_form                   # noqa: E402
from prime_form_factor_base import build_prime_form_factor_base  # noqa: E402
from weighted_norm_relations import least_power_above            # noqa: E402
import numpy as np                                               # noqa: E402
import cypari2                                                   # noqa: E402
pari = cypari2.Pari()
pari.allocatemem(1 << 30, silent=True)
CERT_MAX = 10 ** 7


def disc_of(N):
    """Repo convention: Delta = -N if N = 3 mod 4 else -4N."""
    return -N if N % 4 == 3 else -4 * N


def count_reduced_forms(D):
    """Exhaustive count of primitive reduced forms of discriminant D (certifies h)."""
    n = 0
    amax = isqrt(-D // 3)
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0):
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            n += 1
    return n


class ClassGroup:
    def __init__(self, N, certify=True):
        self.N = N
        D = disc_of(N)
        self.D = D
        pari.setrand(1)  # deterministic generators -> reproducible outputs
        r = pari.quadclassunit(D)
        self.h = int(r[0])
        self.cyc = [int(c) for c in r[1]]
        gens = list(r[2])
        if not self.cyc:  # trivial group
            self.cyc, gens = [1], [pari.qfbred(pari.Qfb(1, D % 2, (D % 2 - D) // 4))]
        self.shape = tuple(self.cyc)
        h = self.h
        ident = pari.qfbred(pari.Qfb(1, D % 2, (D % 2 - D) // 4))
        index = {}
        A = np.zeros(h, dtype=np.int64)
        # mixed radix enumeration, last axis fastest (numpy C order)
        strides = [int(np.prod(self.cyc[i + 1:])) for i in range(len(self.cyc))]
        self.strides = strides
        def rec(level, form, flat):
            if level == len(self.cyc):
                v = pari.Vec(form)
                key = (int(v[0]), int(v[1]), int(v[2]))
                index[key] = flat
                A[flat] = key[0]
                return
            f = form
            for x in range(self.cyc[level]):
                rec(level + 1, f, flat + x * strides[level])
                f = pari.qfbcomp(f, gens[level])
        rec(0, ident, 0)
        assert len(index) == h, ('generators do not give h distinct classes', len(index), h)
        self.index, self.A = index, A
        self.certified = None
        if certify and -D <= CERT_MAX:
            self.certified = (count_reduced_forms(D) == h)
            assert self.certified, 'class number certification failed'

    def dlog(self, form):
        flat = self.index[tuple(int(t) for t in reduce_form(form))]
        return np.array(np.unravel_index(flat, self.shape), dtype=np.int64)

    def flat(self, xvec):
        xvec = np.asarray(xvec) % np.array(self.cyc)
        return int(np.ravel_multi_index(tuple(xvec), self.shape))

    def char_grid(self):
        """All characters k (h x r int array, C order matching fftn output)."""
        g = np.indices(self.shape).reshape(len(self.shape), -1).T
        return g

    def phases(self, tvec):
        """theta(k) = sum_i k_i t_i / c_i mod 1 for every character k (length h)."""
        K = self.char_grid()
        num = np.zeros(self.h)
        for i, c in enumerate(self.cyc):
            num = num + (K[:, i] * int(tvec[i]) % c) / c
        return np.mod(num, 1.0)


def repo_base(N, D, Y):
    base = build_prime_form_factor_base(D, Y)
    recs = {r.prime: r.positive for r in base['records']}
    if N % 4 == 1:   # repo adds the ramified 2 when Delta = -4N, N = 1 mod 4
        recs[2] = reduce_form((2, 2, (N + 1) // 2))
    return dict(sorted(recs.items()))


def smooth_mask(G, primes):
    H = isqrt(-G.D // 3)
    P = prod(primes)
    e = H.bit_length()
    return np.array([pow(P, e, int(a)) == 0 for a in G.A]), H


def dist1(theta):
    """||theta|| distance to nearest integer."""
    return np.abs(theta - np.round(theta))
