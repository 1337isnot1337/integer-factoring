#!/usr/bin/env python3
"""Independently verify a complete crossing and a cofactor nonexistence proof."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import flatten,rank
from degree_five_matching import KEYS,candidate_problem,singular_script

ROOT = Path(__file__).resolve().parents[1]


def import_certificate(log,variables):
    lines = log.read_text().splitlines()
    assert 'BASIS_SIZE=1' in lines and 'CAS_IDENTITY=1' in lines
    terms = []
    start = lines.index('CERTIFICATE_BEGIN')+1
    end = lines.index('CERTIFICATE_END')
    for line in lines[start:end]:
        index,coefficient,exponents = line.split(':')
        powers = list(map(int,exponents.split(',')))
        assert len(powers) == variables and all(x >= 0 for x in powers)
        monomial = [j for j,n in enumerate(powers) for _ in range(n)]
        terms.append([int(index),int(coefficient) % 5,monomial])
    return terms


def intertwiner_constraint_rank(seed,pair):
    rows = []
    for a,b,k in product(range(2),range(2),range(6)):
        f = [seed[(a,y,u,b) if pair == 'yu' else (y,a,b,u)][k]
             for y,u in product(range(2),repeat=2)]
        a0,b0,c0,d0 = f
        rows.extend(([0,-c0,c0,0],[-b0,a0-d0,0,c0],
                     [b0,0,d0-a0,-c0],[0,b0,-b0,0]))
    return rank(rows,5)


def check(candidate):
    rebuilt = json.loads(json.dumps(candidate_problem()))
    for key,value in rebuilt.items():
        assert candidate[key] == value,key
    assert rebuilt['linear_rank'] == 17 and len(rebuilt['free_indices']) == 15
    assert len(rebuilt['quadratic_equations']) == 72
    assert rebuilt['entry_span_dimension'] == 6 and rebuilt['determinant_gcd_degree'] == 0
    seed = {key:value for key,value in zip(KEYS,rebuilt['seed'])}
    assert all(sum(seed[key]) % 5 == int(key[0] == key[3] and key[1] == key[2]) for key in KEYS)
    ranks = []
    for point,value in enumerate(rebuilt['values']):
        table = {key:[a] for key,a in zip(KEYS,value)}
        pair = [rank(flatten(table,2,cut),5) for cut in (0,1)]
        assert pair == ([4,4] if point == 1 else [1,4] if point in (0,2,3) else [4,1])
        ranks.append(pair)
    intertwiner_ranks = {pair:intertwiner_constraint_rank(seed,pair) for pair in ('yu','xv')}
    assert intertwiner_ranks == {'yu':4,'xv':4}
    identity = {}
    terms = candidate['nullstellensatz_certificate']
    maximum = 0
    for index,coefficient,monomial in terms:
        assert 0 <= index < 72 and 0 < coefficient < 5
        assert monomial == sorted(monomial) and all(0 <= i < 15 for i in monomial)
        maximum = max(maximum,len(monomial))
        for c,other in rebuilt['quadratic_equations'][index]:
            product_monomial = tuple(sorted(monomial+other))
            identity[product_monomial] = (identity.get(product_monomial,0)+coefficient*c) % 5
    identity = {m:c for m,c in identity.items() if c}
    assert identity == {():1}
    assert maximum == candidate['certificate_multiplier_degree']
    assert candidate['certificate_identity'] == {'constant':1,'other_terms':0}
    return {'status':'pass','prime':5,'degree':5,'strict_endpoint_entries':16,
            'projective_point_cut_ranks':ranks,'entry_span_dimension':6,
            'cut_determinant_common_degree':0,'intertwiner_constraint_ranks':intertwiner_ranks,
            'cofactor_linear_rank':17,'remaining_variables':15,'rank_one_minor_equations':72,
            'nullstellensatz_multiplier_monomials':len(terms),
            'maximum_multiplier_degree':maximum,'verified_polynomial_identity':'sum multipliers times equations = 1 mod 5',
            'scope':'Absence of a shared projected cofactor over every F5 extension; no factoring improvement or uniform low-rank upper bound.',
            'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('matching_curve.py','degree_five_matching.py',
                                           'check_degree_five_matching_candidate.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,default=ROOT/'evidence/degree-five-matching-candidate.json')
    parser.add_argument('--write-singular',type=Path)
    parser.add_argument('--import-certificate',type=Path)
    parser.add_argument('--candidate-output',type=Path)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    problem = candidate_problem()
    if args.write_singular:
        args.write_singular.write_text(singular_script(problem))
    if args.import_certificate:
        assert args.candidate_output
        problem['nullstellensatz_certificate'] = import_certificate(args.import_certificate,15)
        problem['certificate_multiplier_degree'] = max(len(x[2]) for x in problem['nullstellensatz_certificate'])
        problem['certificate_identity'] = {'constant':1,'other_terms':0}
        args.candidate_output.write_text(json.dumps(problem,separators=(',',':'))+'\n')
        args.candidate = args.candidate_output
    candidate = json.loads(args.candidate.read_text())
    result = check(candidate)
    result['candidate_sha256'] = sha256(args.candidate.read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
