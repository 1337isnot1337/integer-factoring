#!/usr/bin/env python3
"""Exact endpoint majorization, shared-prefix and saturated variance controls."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path
from random import Random

from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import pack


def var(values):
    return Fraction(sum(v*v for v in values),len(values))-Fraction(sum(values),len(values))**2


def brauer(a,b,power,reverse=False):
    degree = power+2
    values = []
    for v,u,y,x in product(range(a),range(b),range(b),range(a)):
        e = int(x==v and y==u)
        # Swapping the middle legs exchanges the cuts while preserving E.
        r = int(x==u and y==v) if reverse else int(x==y and u==v)
        change = e^r
        pair = (r ^ ((1<<power) if change else 0),
                ((1<<(power+1))^(1<<(power+2))) if change else 0)
        if reverse:
            pair = tuple(sum(((p>>i)&1)<<(degree-i) for i in range(degree+1)) for p in pair)
        values.append(pair)
    return values


def profiles(values,a,b,h):
    n,result = a*b,[]
    for place in range(2):
        matrix = [[0]*n for _ in range(n)]
        for v,u,y,x in product(range(a),range(b),range(b),range(a)):
            r,c = (x+a*u,y+b*v) if place==0 else (x+a*y,v+a*u)
            matrix[r][c] = pack(values[x+a*(y+b*(u+b*v))])
        row = [h[place][0][i%a]+h[place][2 if place==0 else 1][i//a] for i in range(n)]
        col = ([h[place][1][i%b]+h[place][3][i//b] for i in range(n)] if place==0 else
               [h[place][3][i%a]+h[place][2][i//a] for i in range(n)])
        orders,_ = WeightedSmithOracle(matrix,infinity=place==0).profile(row,col)
        result.append(list(orders))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--recorded-principal',action='store_true')
    args = parser.parse_args()
    rng,records = Random(0xA815C0),[]
    entry_count,order_count,prefix_count,shared_count,variance_count = 0,0,0,0,0
    cases = []
    for a,b in ((2,2),(2,3),(3,3),(4,4)):
        for power in (1,2,3):
            for strong in (0,1):
                h = [[[0]*q for q in (a,b,b,a)] for _ in range(2)]
                for leg in range(4):
                    h[strong][leg][1:] = [rng.randrange(-3,4) for _ in h[strong][leg][1:]]
                cases.append((a,b,brauer(a,b,power,reverse=strong==1),h,
                              {'kind':'nonbinary_brauer','power':power,'target_strong_place':strong}))
    for a,b in ((2,2),(2,3),(3,3)):
        for _ in range(2):
            values = []
            for v,u,y,x in product(range(a),range(b),range(b),range(a)):
                low,high = rng.getrandbits(3),rng.getrandbits(3)
                values.append((int(x==v and y==u)^low^(low<<1),high^(high<<1)))
            h = [[[0]+[rng.randrange(-2,3) for _ in range(q-1)] for q in (a,b,b,a)] for _ in range(2)]
            cases.append((a,b,values,h,{'kind':'random_strict_polynomial'}))
    for a,b,values,h,label in cases:
        n = a*b
        for v,u,y,x in product(range(a),range(b),range(b),range(a)):
            low,high = values[x+a*(y+b*(u+b*v))]
            assert (low.bit_count()%2 | (high.bit_count()%2)<<1)==int(x==v and y==u)
            entry_count += 1
        cuts = profiles(values,a,b,h)
        order_count += 2*n
        gx = [sum(h[p][leg][i] for p in range(2) for leg in (0,3)) for i in range(a)]
        gy = [sum(h[p][leg][i] for p in range(2) for leg in (1,2)) for i in range(b)]
        energies = sorted(x+y for x in gx for y in gy)
        assert var(energies)==var(gx)+var(gy)
        checks = []
        for strong in (0,1):
            r,s = cuts[strong],cuts[1-strong]
            m,delta = s[0],s[1]-s[0]
            lam = [v+m for v in r]
            for k in range(1,n+1):
                assert sum(lam[:k])<=sum(energies[:k])
                prefix_count += 1
            assert sum(lam[:2])<=sum(energies[:2])-delta
            shared_count += 1
            gap = Fraction(sum(energies)-sum(lam),n)
            assert gap>=0
            variance = var(lam)-var(energies)
            if gap==0:
                assert variance>=0
                if delta>0:
                    assert variance>=Fraction(delta*delta,2*(n-2))
                    variance_count += 1
            if label.get('target_strong_place')==strong:
                assert (gap==0 and delta>0) if a==b else gap>0
            checks.append({'strong_place':strong,'first_order_gap':str(gap),
                           'weak_gap':delta,'variance_difference':str(variance)})
        records.append(dict(label,alphabet_sizes=[a,b],heights=h,profiles=cuts,
                            physical_energies=[gx,gy],checks=checks))
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_polynomial_crossing_axis.py','tools/weighted_smith_oracle_f4.py',
             'tools/windowed_f4_smith.py']
    report = {'status':'pass','literal_endpoints':entry_count,'independent_complete_orders':order_count,
              'endpoint_prefix_inequalities':prefix_count,'shared_second_prefix_inequalities':shared_count,
              'saturated_simple_variance_controls':variance_count,'cases':records,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Exact finite controls for the ordinary all-field endpoint argument. '
                      'No numerical scale mesh is promoted to the uniform theorem, and no '
                      'interior-angle, moving-height or factoring improvement is claimed.'}
    if args.recorded_principal:
        archive = root/'evidence/pencil-principal-screen-complete.json.xz'
        data = json.loads(lzma.decompress(archive.read_bytes()))
        assert data['summary']['status']=='complete_registered_screen'
        pc,sc,vc,states = 0,0,0,0
        for line in data['profiles_jsonl'].splitlines():
            record = json.loads(line)
            gx,gy = record['physical_energies']
            assert all(int(v)==v for v in gx+gy)
            energies = sorted(int(x+y) for x in gx for y in gy)
            n = len(energies)
            for strong in (0,1):
                r,s = record['orders'][strong],record['orders'][1-strong]
                lam = [v+s[0] for v in r]
                lp,ep = 0,0
                for x,y in zip(lam,energies):
                    lp,ep = lp+x,ep+y
                    assert lp<=ep
                    pc += 1
                delta = s[1]-s[0]
                assert sum(lam[:2])<=sum(energies[:2])-delta
                sc += 1
                if sum(lam)==sum(energies) and delta>0:
                    assert var(lam)-var(energies)>=Fraction(delta*delta,2*(n-2))
                    vc += 1
                states += 1
        assert states==90584 and pc==753750 and sc==90584 and vc==939
        report['recorded_principal_replay'] = {'endpoint_states':states,'prefix_inequalities':pc,
            'shared_second_prefix_inequalities':sc,'saturated_variance_floors':vc,
            'archive_sha256':sha256(archive.read_bytes()).hexdigest(),
            'scope':'Replays new exact inequalities on every retained profile; this does not independently recompute all native profiles.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','source_sha256')}))


if __name__=='__main__':
    main()
