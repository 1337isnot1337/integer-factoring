#!/usr/bin/env python3
"""Replay all-field quadratic-output elimination and literal Smith controls."""
from pathlib import Path
from hashlib import sha256
import json,random
from joint_hot_jet_rankone import analyze,lift,regular_matrices,leading_residual,project,pack
from native_matrix_smith_f4 import MatrixOracle


def binary_rank(rows):
    P={}
    for r in rows:
        while r:
            k=r.bit_length()-1
            if k in P:r^=P[k]
            else:P[k]=r;break
    return len(P)


def main():
    records=[];controls=[];rng=random.Random(202611010602)
    for k,expected in ((8,(14,46,60)),(10,(12,36,48)),(12,(12,36,48)),(14,(10,25,35))):
        result=analyze(k)
        assert (result['anchor_dimension'],result['quadratic_output_rank'],result['constraint_output_dimension'])==expected
        assert result['linear_rank_mod_quadratics']==result['anchor_dimension']
        records.append(result)
        d=lift(k);letters,hs,A,M,N=regular_matrices();U0,V0,C0,r0=leading_residual(M,k);Ui,Vi,Ci,ri=leading_residual(N,k)
        for trial in range(8):
            delta=[[0]*64 for _ in range(64)]
            for index in rng.sample(range(len(d['variables'])),min(12,len(d['variables']))):
                i,row=d['variables'][index];delta[i]=[a^b for a,b in zip(delta[i],row)]
            AA=[[a^b for a,b in zip(row,dr)] for row,dr in zip(A,delta)]
            assert all(v.bit_count()%2==int(i==j) for i,row in enumerate(AA) for j,v in enumerate(row))
            Z=[[(v<<hs[i])>>32 for v in row] for i,row in enumerate(AA)]
            def rev(v,degree):
                assert v.bit_length()-1<=degree
                return sum(1<<(degree-j) for j in range(v.bit_length()) if v>>j&1)
            F=[[rev(AA[x+8*y][v+8*u],hs[v+8*u]+32) for v in range(8) for y in range(8)] for u in range(8) for x in range(8)]
            for cut,base,new,U,V,C,r in [('zero',M,Z,U0,V0,C0,r0),('infinity',N,F,Ui,Vi,Ci,ri)]:
                diff=[[x^y for x,y in zip(a,b)] for a,b in zip(base,new)]
                assert all(not(v&((1<<k)-1)) for row in diff for v in row)
                jet=[sum(((v>>k)&1)<<j for j,v in enumerate(row)) for row in diff]
                projected=project(jet,U,V);predicted=binary_rank([x^y for x,y in zip(C,projected)])
                orders=MatrixOracle([[(v,0) for v in row] for row in new],64,precision=128).profile()['orders']
                assert sum(o<k for o in orders)==r
                assert sum(o==k for o in orders)==predicted
                controls.append({'jet':k,'trial':trial,'cut':cut,'earlier_order_count':r,'predicted_first_rank':predicted})
    names=[Path(__file__).name,'joint_hot_jet_rankone.py','joint_hot_retraction_module.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp']
    report={'status':'passed','records':records,'literal_controls':controls,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'All-field first-residual rigidity under simultaneous raw-jet protection for the fixed coordinate embedding; not full-profile or factoring rigidity.'}
    Path('evidence/joint-hot-jet-rankone.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','jets':[8,10,12,14],'all_field_linear_quotient_ranks':[r['linear_rank_mod_quadratics'] for r in records],'literal_full_smith_controls':len(controls)}))

if __name__=='__main__':main()
