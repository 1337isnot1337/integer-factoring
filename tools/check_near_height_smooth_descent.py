#!/usr/bin/env python3
"""Fixed target controls for the near-height ellipse theorem; no orbit sweep."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt, prod
from pathlib import Path

from bounded_smooth_norm import smooth_factors
from near_height_smooth_descent import prime_power_target, smaller_norm_vectors
from prime_form_factor_base import build_prime_form_factor_base, primes_through
from quadratic_form_arithmetic import compose_primitive, inverse_primitive, power_primitive, reduce_primitive
from smooth_norm_relation import norm_to_prime_word
from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "evidence/near-height-smooth-descent.json"


def audit(discriminant, prime, exponent, cutoff, root):
    cap = prime ** exponent
    form, known = prime_power_target(discriminant, prime, exponent, root)
    seed = reduce_primitive((prime, root, (root*root - discriminant)//(4*prime)))
    assert power_primitive(seed, exponent) == form
    a, b, c = form
    d = -discriminant
    assert cap > isqrt(d // 3)
    leading_success = gcd(a, d) == 1 and smooth_factors(a, cutoff) is not None
    points = list(smaller_norm_vectors(form, cap))
    # Independent signed-box exhaustion using direct Q, not the row roots.
    ymax = isqrt(4*a*(cap - 1) // d)
    xmax = isqrt((cap - 1) // a) + (ymax + 1)//2 + 1
    direct = {(x, y, a*x*x + b*x*y + c*y*y)
              for y in range(ymax + 1) for x in range(-xmax, xmax + 1)
              if (y > 0 or x > 0) and 0 < a*x*x + b*x*y + c*y*y < cap}
    assert set(points) == direct and len(points) == len(direct)
    if not leading_success:
        assert known[1] != 0 and 4*a*cap >= d
        integer_ratio = (cap + isqrt(d) - 1) // isqrt(d)
        assert len(points) <= 17 * integer_ratio
    successes = []
    parity = discriminant % 2
    identity = (1, parity, (parity - discriminant)//4)
    for x, y, norm in points:
        if gcd(x, y) != 1 or gcd(norm, d) != 1:
            continue
        decoded = norm_to_prime_word(form, x, y, cutoff, cap - 1)
        if decoded is None:
            continue
        value = identity
        for term in decoded['word']:
            base = term['form'] if term['exponent'] > 0 else inverse_primitive(term['form'])
            value = compose_primitive(value, power_primitive(base, abs(term['exponent'])))
        assert value == form
        assert prod(term['prime'] ** abs(term['exponent']) for term in decoded['word']) == norm < cap
        # This exact product inequality proves weighted dominance of k*e_p-v.
        p_coordinate = next((t['exponent'] for t in decoded['word'] if t['prime'] == prime), 0)
        assert exponent - p_coordinate > 0
        assert prime ** (exponent - p_coordinate) > norm // prime ** abs(p_coordinate)
        successes.append(dict(vector=[x, y], norm=norm,
                              word=[dict(prime=t['prime'], exponent=t['exponent'])
                                    for t in decoded['word']]))
    return dict(discriminant=discriminant, prime=prime, exponent=exponent,
                cutoff=cutoff, cap=cap, reduced_form=list(form),
                known_primitive_vector=list(known), leading_success=leading_success,
                smaller_vectors=len(points), smooth_successes=successes)


def rough_composite_control(exponent, ratio, cutoff, sieve_cutoff):
    """One bounded CRT/sieve construction, not a factoring or orbit search."""
    cap = 3 ** exponent
    q = next(p for p in primes_through(2*sieve_cutoff) if p > sieve_cutoff)
    residue = (q + 1) * pow(4*cap, -1, q*q) % (q*q)
    lower = cap // (8*ratio*ratio)
    first = lower + (residue - lower) % (q*q)
    primorial = prod(primes_through(sieve_cutoff))
    for rejected in range(4096):
        rough = first + rejected*q*q
        assert rough < 2*lower
        d = 4*cap*rough - 1
        if gcd(rough*d, primorial) == 1:
            break
    else:
        raise AssertionError('The fixed CRT/sieve construction exhausted its budget.')
    assert d > q and d % (q*q) == q  # A known proper prime factor of valuation one.
    trial_cutoff = 1 << (isqrt(d.bit_length() - 1) + 1)
    assert trial_cutoff <= sieve_cutoff
    record = audit(-d, 3, exponent, cutoff, 1)
    assert record['reduced_form'] == [rough, -1, cap]
    assert not record['smooth_successes']
    assert sum(gcd(x, y) == 1 for x, y, _ in smaller_norm_vectors((rough, -1, cap), cap)) == 1
    return dict(**record, sieve_cutoff=sieve_cutoff, proper_prime_factor=q,
                prime_factor_valuation=1, rejected_crt_candidates=rejected,
                trial_division_cutoff=trial_cutoff)


def bonferroni_control():
    """Exact product/intersection controls for the elementary lower-sieve proof."""
    primes = primes_through(19)
    m, q, length = 81, 23, 257
    residue = (q + 1) * pow(4*m, -1, q*q) % (q*q)
    roots = {p: [t for t in range(p)
                 if ((residue + q*q*t) * (4*m*(residue + q*q*t) - 1)) % p == 0]
             for p in primes}
    assert all(len(roots[p]) == (1 if p in (2, 3) else 2) for p in primes)
    # Group intersection counts by number of selected prime conditions.
    counts, masses, errors = [0]*9, [Fraction(0)]*9, [0]*9
    states = [(1, [0], 0)]
    for p in primes:
        extended = []
        for modulus, residues, order in states:
            next_roots = [(r + modulus*((s - r)*pow(modulus, -1, p) % p))
                          for r in residues for s in roots[p]]
            extended.append((modulus*p, next_roots, order + 1))
        states += extended
    for modulus, residues, order in states:
        count = sum((length - 1 - r)//modulus + 1 for r in residues if r < length)
        mass = Fraction(len(residues), modulus)
        assert abs(count - length*mass) <= len(residues)
        counts[order] += count
        masses[order] += mass
        errors[order] += len(residues)
    product_mass = prod(1 - Fraction(len(roots[p]), p) for p in primes)
    survivors = sum(all(t % p not in roots[p] for p in primes) for t in range(length))
    assert sum((-1)**j * counts[j] for j in range(9)) == survivors
    for order in (1, 3, 5, 7):
        lower = sum((-1)**j * masses[j] for j in range(order + 1))
        lower_count = sum((-1)**j * counts[j] for j in range(order + 1))
        assert lower <= product_mass <= lower + masses[order + 1]
        assert lower_count <= survivors <= lower_count + counts[order + 1]
        assert abs(lower_count - length*lower) <= sum(errors[:order + 1])
    return dict(primes=primes, length=length, survivors=survivors,
                independent_product=str(product_mass),
                bonferroni_odd_orders=[1, 3, 5, 7], exact_crt_intersections=len(states))


def compute():
    controls = []
    # Existing finite inputs, only two predeclared powers per fixed prime.
    for integer, cutoff, prime in (
            (1000000016000000063, 1000, 7),
            (906694351050364540708897, 1000, 7),
            (4951760154835678088235319297, 1021, 13)):
        discriminant = -integer if integer % 4 == 3 else -4*integer
        records = {r.prime: r for r in build_prime_form_factor_base(discriminant, cutoff)['records']}
        first = least_power_above(prime, isqrt(-discriminant // 3))
        for offset in (0, 1):
            controls.append(audit(discriminant, prime, first + offset, cutoff,
                                  min(records[prime].root, 2*prime - records[prime].root)))
    counterexamples = []
    # Roughness is checked by division; no primality claim about a is needed.
    for exponent in (6, 14, 32, 62):
        cap, rough = 2 ** exponent, 2 ** (exponent - 1) - 1
        cutoff = min(1000, 2 ** (exponent - 3))
        discriminant = 1 - 4*cap*rough
        assert smooth_factors(rough, cutoff) is None
        assert least_power_above(2, isqrt(-discriminant // 3)) == exponent
        record = audit(discriminant, 2, exponent, cutoff, 1)
        assert record['reduced_form'] == [rough, -1, cap]
        assert record['smaller_vectors'] == 1 and not record['smooth_successes']
        counterexamples.append(record)
    odd = audit(-7451, 3, 4, 19, 1)  # m=81, a=23; the same obstruction at an odd prime.
    assert least_power_above(3, isqrt(7451 // 3)) == 4
    assert odd['reduced_form'] == [23, -1, 81]
    assert odd['smaller_vectors'] == 1 and not odd['smooth_successes']
    counterexamples.append(odd)
    assert any(c['smooth_successes'] and not c['leading_success'] for c in controls)
    sources = ['tools/near_height_smooth_descent.py', 'tools/check_near_height_smooth_descent.py']
    return dict(scope='Complete smaller-norm ellipses at fixed near-height prime powers; no uniform smooth-existence or score claim.',
                controls=controls, first_power_counterexamples=counterexamples,
                rough_composite_counterexamples=[rough_composite_control(24, 16, 19, 1009),
                                                 rough_composite_control(40, 64, 43, 10007)],
                bonferroni_control=bonferroni_control(),
                source_sha256={name: sha256((ROOT / name).read_bytes()).hexdigest() for name in sources})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.write:
        EXPECTED.write_text(json.dumps(result, indent=2) + '\n')
    else:
        assert result == json.loads(EXPECTED.read_text())
    print(f"near-height descent: {len(result['controls'])} fixed targets, "
          f"{len(result['first_power_counterexamples'])} first-power failures, "
          "two rough composite controls and 256 CRT intersections passed")
