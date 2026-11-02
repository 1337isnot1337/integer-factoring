#!/usr/bin/env python3
"""Exact finite controls for adjacent CM points with opposite smooth flags."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path

from bounded_smooth_norm import smooth_factors
from prime_form_factor_base import build_prime_form_factor_base, primes_through
from quadratic_form_arithmetic import form_discriminant, reduce_primitive
from smooth_norm_relation import norm_to_prime_word

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / 'evidence/cm-smooth-flag-geometry.json'


def compute():
    controls = []
    for exponent, cutoff in ((3, 5), (8, 5), (10, 11), (12, 19), (16, 31)):
        start = 2 ** exponent
        q = next(p for p in primes_through(2*cutoff) if p > cutoff)
        rough_multiple = start + q - start % q
        a = start
        while smooth_factors(a + 1, cutoff) is not None:
            a += 1
        assert a < rough_multiple <= start + 2*cutoff
        d = 20*a*(a + 1) - 1
        radius = isqrt(d)//4
        forms = ((a, 1, 5*(a + 1)), (a + 1, 1, 5*a))
        assert a + 1 <= radius
        for form in forms:
            assert form_discriminant(form) == -d and reduce_primitive(form) == form
            assert gcd(form[0], d) == 1
        left = norm_to_prime_word(forms[0], 1, 0, cutoff)
        right = norm_to_prime_word(forms[1], 1, 0, cutoff)
        assert left is not None and right is None
        records = build_prime_form_factor_base(-d, cutoff)['records']
        candidates = [r.prime for r in records if a % r.prime]
        pfree_prime = min(candidates) if candidates else None
        if pfree_prime is not None:
            assert all(t['prime'] != pfree_prime for t in left['word'])
        # cosh(distance)-1 = |z_a-z_(a+1)|^2/(2 Im(z_a) Im(z_(a+1))).
        cosh_minus_one = Fraction(d + 1, 2*d*a*(a + 1))
        assert cosh_minus_one == Fraction(10, d)
        controls.append(dict(exponent=exponent, cutoff=cutoff, discriminant=-d,
                             forms=[list(f) for f in forms], radius=radius,
                             searched_integers=a - start + 1, raw_flags=[1, 0],
                             p_free_prime=pfree_prime,
                             cosh_lifted_distance_minus_one=str(cosh_minus_one),
                             lifted_distance_squared_upper_bound=str(Fraction(20, d)),
                             exact_interpolation_lipschitz_squared_lower_bound=str(Fraction(d, 20))))
    # Small-cutoff guard: there need not be a split target prime avoiding a.
    assert controls[0]['p_free_prime'] is None
    assert all(c['p_free_prime'] is not None for c in controls[1:])
    return dict(scope='Two constructed reduced forms per fixed control; no class-cycle enumeration, analytic theorem verification or complete bound change.',
                controls=controls,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.write:
        EXPECTED.write_text(json.dumps(result, indent=2) + '\n')
    else:
        assert result == json.loads(EXPECTED.read_text())
    print('CM smooth-flag geometry: five exact adjacent-form controls and the small-cutoff guard passed')
