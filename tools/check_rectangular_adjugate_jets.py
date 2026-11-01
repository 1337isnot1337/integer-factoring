#!/usr/bin/env python3
"""Check rectangular first possible cofactor jets and a nine-letter witness."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path
from random import Random

from check_cell_minor_caps import (PRIME, Dyadic, log_at_least_one, point,
                                   polynomial_minor, smith_orders)
from check_subalphabet_bare_entropy import independent_constraint
from rectangular_adjugate_jet_caps import linear_caps


def profile(matrix):
    precision = 1+sum(max(len(p)-1 for p in row) for row in matrix)
    return smith_orders(matrix, precision)


def rectangular_controls():
    rng, records, minors, counterexample = Random(0x2EC7AA), [], 0, None
    cases = [(2,3,1,0),(2,3,2,1),(2,3,1,2),(2,4,1,0),
             (2,4,2,1),(2,4,1,4),(2,5,1,6),(1,4,2,1)]
    for trial, (p,b,k,kernel_rank) in enumerate(cases):
        n, h = p*b, p*p
        nullity, gap = n-h, n-h-1
        pairs = list(product(range(p), range(b)))
        u = [[rng.randrange(PRIME) for _ in range(n)] for _ in range(n)]
        kernel = [i for i,(_,y) in enumerate(pairs) if y >= p]
        if kernel_rank == 0:
            for i,j in product(kernel, repeat=2):
                u[i][j] = 0
        elif kernel_rank == 1:
            left, right = [[rng.randrange(PRIME) for _ in kernel] for _ in range(2)]
            for a,i in enumerate(kernel):
                for c,j in enumerate(kernel):
                    u[i][j] = left[a]*right[c] % PRIME
        row, crossing = [], []
        for i,(x,y) in enumerate(pairs):
            line = []
            for j,(v,z) in enumerate(pairs):
                channel = int(x == z and y == v and z < p and y < p)
                polynomial = [channel, -channel % PRIME]+[0]*k
                polynomial[k] = (polynomial[k]+u[i][j]) % PRIME
                polynomial[k+1] = (int(i == j)-u[i][j]) % PRIME
                line.append(polynomial)
            row.append(line)
        cross = [[row[x*b+y][v*b+z] for y,v in product(range(b),range(p))]
                 for x,z in pairs]
        source_orders = profile(cross)
        assert source_orders[0] == 0 and source_orders[1] == k
        entries = {i:{j:row[i][j] for j in range(n)} for i in range(n)}
        adj = []
        for i in range(n):
            line = []
            for j in range(n):
                polynomial = polynomial_minor(entries,[a for a in range(n) if a != j],
                                               [a for a in range(n) if a != i])
                line.append([(-v if (i+j) % 2 else v) % PRIME for v in polynomial] or [0])
            adj.append(line)
        adj_cross = [[adj[v*b+z][x*b+y] for z,x in product(range(b),range(p))]
                     for v,y in pairs]
        observed = profile(adj_cross)
        lower = [gap*k]*nullity+[(gap+1)*k]*(n-nullity)
        assert all(a >= b for a,b in zip(observed,lower))
        # This tempting second refinement fails when the kernel block has a
        # singular leading term; retain the full counterexample, not a claim.
        proposed = [gap*k]*nullity+[(gap+1)*k]+[(gap+2)*k]*(h-1)
        stronger_pass = all(a >= b for a,b in zip(observed,proposed))
        rec = {'alphabet':[p,b], 'jet_order':k, 'kernel_rank_policy':kernel_rank,
               'U':u, 'row_tensor':row, 'adjugate':adj,
               'crossed_source_orders':source_orders,'crossed_adjugate_orders':observed,
               'proved_lower_profile':lower,'unproved_stronger_profile':proposed,
               'stronger_profile_pass':stronger_pass}
        if not stronger_pass and counterexample is None:
            counterexample = rec
        records.append(rec)
        minors += n*n
    assert counterexample is not None
    return records, minors, counterexample


def literal_cap_controls(root):
    previous = json.loads(lzma.decompress((root/'evidence/subalphabet-adjugate-controls.json.xz').read_bytes()))
    checked, summaries = 0, []
    cases = []
    for rec in previous['literal_controls']:
        data = rec['normalization']
        for pair in data['profiles'].values():
            cases.append((data['deficit'],pair['r'],pair['s'],data['paired_energies_X'],data['paired_energies_Y']))
    for rec in previous['characteristic_two_controls']:
        first,second = rec['complete_profiles']
        cases.append((rec['deficit'],[v-first[0] for v in first],
                      [v-second[0] for v in second],*rec['energies']))
    for d,r,s,x,y in cases:
        caps = linear_caps(d,x,y)
        slacks = []
        for cap in caps:
            for a,b in ((r,s),(s,r)):
                lhs = sum(c*v for c,v in zip(cap['first'],a))+sum(c*v for c,v in zip(cap['second'],b))
                slack = cap['rhs']-lhs
                assert slack >= 0
                mu,nu = [F(v)+F(17,3) for v in a], [F(v)-d-F(17,3) for v in b]
                absolute = sum(c*v for c,v in zip(cap['absolute_first'],mu))
                absolute += sum(c*v for c,v in zip(cap['absolute_second'],nu))
                assert cap['physical_rhs']-absolute == slack
                slacks.append(slack)
        checked += len(slacks)
        summaries.append({'alphabet':len(x),'inequalities':len(slacks),'minimum_slack':str(min(slacks))})
    return checked,summaries


def entropy_control(root):
    data = json.loads((root/'evidence/rectangular-jet-entropy-witness.json').read_text())
    q, rho = data['alphabet'], F(data['rho_floor'])
    coefficients, rhs = [F(0)]*(q*q), F(0)
    reference = {(c['family'],*c['alphabet'],c['size']):c for c in linear_caps(1,[0]*q,[0]*q)}
    for label,weight in data['integer_constraint_weights']:
        assert type(weight) is int and weight > 0
        if label[0] == 'rectangular jet cell':
            _,p,b,k = label
            n, nullity = p*b, p*b-min(p,b)**2
            assert k == n
            row = [F(n) if i < n-1 else F(0) for i in range(q*q)]
            row[1] += (n-1)*nullity
            bound = F(n*(n-1))
            cap = reference[tuple(label)]
            assert row == [a+b for a,b in zip(cap['first'],cap['second'])] and bound == cap['rhs']
        else:
            row,bound = independent_constraint(q,label)
        coefficients = [a+weight*b for a,b in zip(coefficients,row)]
        rhs += weight*bound
    probabilities = [v/(2*rhs) for v in coefficients]
    probabilities[0] = 1-sum(probabilities[1:])
    assert min(probabilities)>0 and sum(probabilities)==1
    arithmetic,entropy = Dyadic(100),(0,0)
    for p,multiplicity in Counter(probabilities).items():
        entropy = arithmetic.add(entropy,arithmetic.multiply(point(arithmetic,multiplicity*p),
            log_at_least_one(arithmetic,point(arithmetic,1/p))))
    required = arithmetic.multiply(point(arithmetic,2*rho),log_at_least_one(arithmetic,point(arithmetic,q)))
    margin = arithmetic.subtract(entropy,required)
    assert margin[0]>0 and 2*rho.numerator**2>rho.denominator**2
    return {'alphabet':q,'rho_floor':str(rho),'complete_probability':list(map(str,probabilities)),
            'unnormalized_rhs':str(rhs),'normalized_weighted_rhs':'1/2','precision':100,
            'entropy_margin_interval':list(map(str,margin)),
            'entropy_margin_decimal':[float(F(v,arithmetic.scale)) for v in margin]}


def binary_tail_control(root):
    from check_adjugate_seven_letter_floor import adjugate
    from weighted_smith_oracle_f4 import WeightedSmithOracle
    from windowed_f4_smith import pack
    data = json.loads((root/'evidence/rectangular-binary-tail-counterexample.json').read_text())
    row = data['row_tensor_binary_polynomials']
    assert len(row)==6 and all(len(line)==6 for line in row)
    assert all(v.bit_count()%2 == int(i==j)
               for i,line in enumerate(row) for j,v in enumerate(line))
    computed = adjugate(row)
    assert computed == data['adjugate_binary_polynomials']
    pairs = list(product(range(2),range(3)))
    cross = [[row[x*3+y][v*3+z] for y,v in product(range(3),range(2))] for x,z in pairs]
    adj_cross = [[computed[v*3+z][x*3+y] for z,x in product(range(3),range(2))] for v,y in pairs]
    observed = []
    for matrix in (cross,adj_cross):
        orders,_ = WeightedSmithOracle([[pack((v,0)) for v in line] for line in matrix]).profile([0]*6,[0]*6)
        observed.append(list(orders))
    assert observed == [data['crossed_source_orders'],data['crossed_adjugate_orders']]
    assert observed == [[0,2,2,2,2,2],[2,3,5,5,6,7]]
    assert all(a>=b for a,b in zip(observed[1],[2,2,4,4,4,4]))
    assert any(a<b for a,b in zip(observed[1],data['failed_stronger_lower']))
    return {'field':2,'exact_endpoint_entries':36,'permutation_expanded_cofactors':36,
            'independent_fraction_free_orders':12,'source_profiles':observed,
            'failed_stronger_lower':data['failed_stronger_lower']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    rectangles,minors,counterexample = rectangular_controls()
    checked,summaries = literal_cap_controls(root)
    names = ['tools/rectangular_adjugate_jet_caps.py','tools/check_rectangular_adjugate_jets.py',
             'tools/check_subalphabet_bare_entropy.py','tools/check_cell_minor_caps.py',
             'tools/charged_gaussian_floor.py','evidence/rectangular-jet-entropy-witness.json',
             'evidence/subalphabet-adjugate-controls.json.xz',
             'evidence/rectangular-binary-tail-counterexample.json',
             'tools/check_adjugate_seven_letter_floor.py','tools/weighted_smith_oracle_f4.py',
             'tools/windowed_f4_smith.py']
    result = {'status':'pass','literal_rectangular_controls':rectangles,
              'independently_interpolated_cofactors':minors,'literal_lifted_inequalities':checked,
              'literal_cap_summaries':summaries,'counterexample_to_stronger_tail':counterexample,
              'binary_counterexample_replay':binary_tail_control(root),
              'nine_letter_entropy_certificate':entropy_control(root),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'First possible rectangular cofactor jet has bounded reshuffled rank. '
              'Ordinary all-field necessary inequalities and all-scale bare nine-letter floor. '
              'No stronger factoring exponent or global physical-energy exclusion.'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({k:result[k] for k in ('status','independently_interpolated_cofactors',
        'literal_lifted_inequalities','binary_counterexample_replay')}
        | {'nine_letter_floor':result['nine_letter_entropy_certificate']['rho_floor'],
           'positive_entropy_margin':result['nine_letter_entropy_certificate']['entropy_margin_decimal']}))


if __name__ == '__main__':
    main()
