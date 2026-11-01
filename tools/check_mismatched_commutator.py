#!/usr/bin/env python3
"""Check the two-mode commutator projector, with an optional full-shell replay."""
import argparse
import hashlib
import json
from pathlib import Path
from time import monotonic

from check_mismatched_weyl_spectrum import product, rank, read_dump


def local_check(bound):
    entries = 0
    spaces = 0
    for ha in range(bound+1):
        first = [(ha-3*j,j) for j in range(ha//3+1)]
        for hb in range(bound+1):
            second = [(hb-3*j,j) for j in range(hb//3+1)]
            states = [(a,y) for a in first for y in second]
            index = {state:i for i,state in enumerate(states)}
            H = [0]*len(states)
            for j,(beta,u) in enumerate(states):
                input_d = (beta[0]&1)*(u[0]&1)
                moves = []
                if beta[1]&1 and u[1]>=1 and u[0]%4==0:
                    moves.append(((beta[0]+3,beta[1]-1),(u[0]+3,u[1]-1)))
                if beta[0]%4==3 and u[0]>=3 and u[1]%2==0:
                    moves.append(((beta[0]-3,beta[1]+1),(u[0]-3,u[1]+1)))
                for alpha,y in moves:
                    output_d = (alpha[0]&1)*(y[0]&1)
                    if input_d ^ output_d:
                        H[index[(alpha,y)]] ^= 1 << j
            H2 = product(H,H)
            H3 = product(H2,H)
            assert product(H2,H2) == H2
            for i,(alpha,y) in enumerate(states):
                zero = alpha[0]%4==y[0]%4==0 and alpha[1]%2==y[1]%2==1
                one = alpha[0]%4==y[0]%4==3 and alpha[1]%2==y[1]%2==0
                assert H2[i] == int(zero or one) << i
                if zero:
                    target = ((alpha[0]+3,alpha[1]-1),(y[0]+3,y[1]-1))
                    assert H3[i] == 1 << index[target]
                elif one:
                    target = ((alpha[0]-3,alpha[1]+1),(y[0]-3,y[1]+1))
                    assert H3[i] == 1 << index[target]
                else:
                    assert not H3[i]
            entries += len(states)**2
            spaces += 1
    return {"local_energy_bound":bound,"independent_energy_spaces":spaces,"local_entries":entries}


def source_check(prefix):
    h,d,weights,labels,blocks = read_dump(str(prefix)+".row.bin")
    assert weights[:2] == (1,3) and all(w>3 for w in weights[2:])
    f0 = {i for i,a in enumerate(labels) if a[0]%4==0 and a[1]%2==1}
    f1 = {i for i,a in enumerate(labels) if a[0]%4==3 and a[1]%2==0}
    assert len(f0) == len(f1)
    ids = {a:i for i,a in enumerate(labels)}
    tau = {}
    for i in f0:
        a = list(labels[i]); a[0]+=3; a[1]-=1
        j = ids[tuple(a)]
        assert j in f1
        tau[i]=j; tau[j]=i
    mask1 = sum(1<<j for j in range(h+1) if j&1)
    mask3 = sum(1<<j for j in range(h+1) if j&3==3)
    rh = rp = rr = entries = 0
    for pairs,polys in blocks:
        n = len(pairs)
        local = {code:i for i,code in enumerate(pairs)}
        diagonal = [(labels[code//d][0]&1)*(labels[code%d][0]&1) for code in pairs]
        B3 = []
        for i in range(n):
            row = 0
            for j in range(n):
                p = polys[i*n+j]
                assert (p&mask1).bit_count()%2 == diagonal[i]*(i==j)
                row |= ((p&mask3).bit_count()%2) << j
                entries += 1
            B3.append(row)
        H = [sum(((B3[i]>>j)&1)*(diagonal[i]^diagonal[j]) << j for j in range(n)) for i in range(n)]
        P = product(H,H)
        H3 = product(P,H)
        R = [a^b for a,b in zip(P,H3)]
        for i,code in enumerate(pairs):
            alpha,y = divmod(code,d)
            active = (alpha in f0 and y in f0) or (alpha in f1 and y in f1)
            assert P[i] == int(active) << i
            expected = (1<<local[tau[alpha]*d+tau[y]]) if active else 0
            assert H3[i] == expected
        rh += rank(H); rp += rank(P); rr += rank(R)
    m = len(f0)
    assert rp == 2*m*m and rr == m*m
    result = {"energy":h,"labels":d,"sector_labels":m,
            "source_entries_checked":entries,"rank_H":rh,"rank_projector":rp,
            "pairing_row_rank":rr,"pairing_column_rank_by_verified_permutation":4*m*m}
    if h == 63 and weights == (1,3,9,27):
        alpha = (57,2,0,0)
        beta = (3,20,0,0)
        partner = (0,21,0,0)
        assert ids[alpha] not in f0|f1
        assert ids[beta] in f1 and tau[ids[beta]] == ids[partner]
        def entry(a,b):
            return int(all(x & ~(x+y) == 0 for x,y in zip(a,b)))
        assert entry(alpha,beta) == 0 and entry(alpha,partner) == 1
        result["noncommutation_witness"] = {"row_pair":[alpha,alpha],"column_pair":[beta,beta],"AR":1,"RA":0}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound",type=int,default=24)
    parser.add_argument("--prefix")
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    result = local_check(args.bound)
    if args.prefix:
        result["source"] = source_check(args.prefix)
    result["source_sha256"] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(__file__).with_name("check_mismatched_weyl_spectrum.py"))}
    result["seconds"] = round(monotonic()-started,3)
    text = json.dumps(result,indent=2)+"\n"
    if args.output:
        Path(args.output).write_text(text)
    print(text,end="")


if __name__ == "__main__":
    main()
