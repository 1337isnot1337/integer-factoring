#!/usr/bin/env python3
"""Literal cofactors, normalization controls and a strict relaxation separation."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path
from random import Random

from check_cell_minor_caps import (PRIME, Dyadic, all_caps, literal_case,
                                   polynomial_minor, pressure, random_crossing,
                                   smith_orders, valuation)
from hereditary_adjugate_caps import inequalities
from multichannel_cofactor_caps import capacity_weights, cofactor_order_lower, linear_caps
from rectangular_adjugate_jet_caps import linear_caps as jet_caps

ROOT = Path(__file__).resolve().parents[1]


def constant_basis(rng, n, prime=PRIME):
    result = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(5*n):
        i, j = rng.sample(range(n), 2)
        scalar = rng.randrange(1, prime)
        result[i] = [(a+scalar*b) % prime for a, b in zip(result[i], result[j])]
    return result


def cofactor_controls():
    rng, records, count = Random(0xCA9AC17), [], 0
    for p, b, orders in ((1, 5, [0, 1, 2, 3, 4]),
                         (2, 5, [0, 1]+[3]*8),
                         (2, 6, [0, 1]+[2]*10)):
        n = p*b
        u, v = constant_basis(rng, n), constant_basis(rng, n)
        source = [[[sum(u[i][k]*v[k][j] for k in range(n) if orders[k] == e) % PRIME
                    for e in range(max(orders)+1)] for j in range(n)] for i in range(n)]
        observed = smith_orders(source, n*max(orders)+1)
        assert observed == orders
        row = [[source[x*b+z][y*p+w] for w, z in product(range(p), range(b))]
               for x, y in product(range(p), range(b))]
        entries = {i: dict(enumerate(line)) for i, line in enumerate(row)}
        cofactors = [[polynomial_minor(entries, [a for a in range(n) if a != j],
                                      [a for a in range(n) if a != i])
                      for j in range(n)] for i in range(n)]
        actual = min(valuation(poly) for line in cofactors for poly in line if any(poly))
        lower = cofactor_order_lower(p, b, observed)
        old = (n-1-min(p, b)**2)*orders[1]
        assert actual >= lower > old
        records.append(dict(alphabet=[p, b], crossed_source=source,
                            crossed_source_orders=observed, unsigned_cofactors=cofactors,
                            minimum_cofactor_order=actual, multichannel_lower=lower,
                            old_rank_one_lower=old))
        count += n*n
    return records, count


def binary_controls():
    from check_adjugate_seven_letter_floor import adjugate
    from check_joint_energy_conjugations import polynomial
    from weighted_smith_oracle_f4 import WeightedSmithOracle
    rng, records, count = Random(0xB1CA9), [], 0
    for p, b, orders in ((1, 5, [0, 1, 2, 3, 4]), (2, 3, [0, 1]+[3]*4)):
        n = p*b
        u, v = constant_basis(rng, n, 2), constant_basis(rng, n, 2)
        source = [[sum((sum(u[i][k]*v[k][j] for k in range(n) if orders[k] == e) % 2) << e
                       for e in range(max(orders)+1)) for j in range(n)] for i in range(n)]
        observed, _ = WeightedSmithOracle([[polynomial((a, 0)) for a in line] for line in source]).profile([0]*n, [0]*n)
        assert list(observed) == orders
        row = [[source[x*b+z][y*p+w] for w, z in product(range(p), range(b))]
               for x, y in product(range(p), range(b))]
        cofactors = adjugate(row)
        actual = min((a & -a).bit_length()-1 for line in cofactors for a in line if a)
        lower = cofactor_order_lower(p, b, orders)
        assert actual >= lower
        records.append(dict(alphabet=[p, b], crossed_source=source, source_orders=orders,
                            adjugate=cofactors, minimum_cofactor_order=actual, lower=lower))
        count += n*n
    return dict(permutation_expanded_cofactors=count, cases=records)


def literal_controls():
    rng, count, profiles, summaries = Random(0xC0FAC), 0, 0, []
    for q in (2, 3, 4):
        for trial in range(3):
            tensor = random_crossing(rng, q, 2, .65)
            heights = [[[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2)]
            data = literal_case(rng, q, tensor, heights, cell_samples=0)
            caps = linear_caps(data['deficit'], data['paired_energies_X'], data['paired_energies_Y'])
            slacks = []
            for pair in data['profiles'].values():
                profiles += 2
                for r, s in ((pair['r'], pair['s']), (pair['s'], pair['r'])):
                    mu = [F(v)+F(17, 3) for v in r]
                    nu = [F(v)-data['deficit']-F(17, 3) for v in s]
                    for cap in caps:
                        value = lambda a, b: (sum(c*v for c, v in zip(cap['first'], a))+
                                              sum(c*v for c, v in zip(cap['second'], b)))
                        slack = cap['rhs']-value(r, s)
                        assert slack >= 0
                        assert cap['physical_rhs']-value(mu, nu) == slack
                        slacks.append(slack)
            count += len(slacks)
            summaries.append(dict(alphabet=q, inequalities=len(slacks), minimum_slack=str(min(slacks))))
    return dict(tensors=9, complete_profiles=profiles, inequalities=count, summaries=summaries)


def separation():
    path = ROOT/'evidence/multichannel-q8-separation.json'
    data = json.loads(path.read_text())
    scale = F(data['denominator'])
    d, r, e = F(data['deficit'], scale), [F(v, scale) for v in data['orders']], [F(v, scale) for v in data['energies']]
    assert len(e) == 8 and len(r) == 64 and r == sorted(r) and r[0] == 0
    assert all(sum(r[:k]) <= bound for k, (bound, _) in all_caps(d, e, e).items())
    square = inequalities(d, r, r, e, e)
    assert all(cap['slack'] >= 0 for cap in square)
    previous = jet_caps(d, e, e)
    assert all(sum(c*v for c, v in zip(cap['first'], r))+
               sum(c*v for c, v in zip(cap['second'], r)) <= cap['rhs'] for cap in previous)
    arithmetic = Dyadic(100)
    interval = pressure(arithmetic, d, r, e, e)
    assert interval[1]*500 < -arithmetic.scale
    rejected = []
    for cap in linear_caps(d, e, e):
        violation = sum(c*v for c, v in zip(cap['first'], r))+sum(c*v for c, v in zip(cap['second'], r))-cap['rhs']
        if violation > 0:
            rejected.append(dict(alphabet=cap['alphabet'], family=cap['family'], size=cap['size'], violation=str(violation)))
    assert rejected
    return dict(scope='Exact negative point of the previous necessary relaxation, excluded by multichannel caps; no tensor.',
                denominator=data['denominator'], square_caps=len(square), rectangular_jet_caps=len(previous),
                pressure_interval=list(interval), pressure_precision=100, pressure_upper_bound='-1/500',
                multichannel_violations=len(rejected), largest_violations=sorted(rejected, key=lambda c: F(c['violation']), reverse=True)[:4],
                input_sha256=sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records, count = cofactor_controls()
    names = ['tools/multichannel_cofactor_caps.py', 'tools/check_multichannel_cofactor_caps.py',
             'tools/rectangular_adjugate_jet_caps.py', 'tools/hereditary_adjugate_caps.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py',
             'tools/check_adjugate_seven_letter_floor.py', 'tools/check_joint_energy_conjugations.py',
             'tools/weighted_smith_oracle_f4.py']
    result = dict(scope='All-field necessary multichannel rectangular cofactor constraints; no stronger factoring algorithm.',
                  prime=PRIME, cofactor_controls=records, interpolated_cofactors=count,
                  binary_controls=binary_controls(),
                  literal_controls=literal_controls(), strict_relaxation_separation=separation(),
                  source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    if args.output:
        args.output.write_bytes(lzma.compress((json.dumps(result, sort_keys=True, default=str)+'\n').encode()))
    else:
        expected = json.loads(lzma.decompress((ROOT/'evidence/multichannel-cofactor-controls.json.xz').read_bytes()))
        assert json.loads(json.dumps(result, default=str)) == expected
    print(json.dumps({key: result[key] for key in ('scope', 'interpolated_cofactors', 'literal_controls', 'strict_relaxation_separation')}))


if __name__ == '__main__':
    main()
