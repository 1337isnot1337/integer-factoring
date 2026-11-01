#!/usr/bin/env python3
"""All-algebraic-place bare-pressure exclusion for the retained Hecke seeds."""
import argparse,json,time
from hashlib import sha256
from pathlib import Path
import numpy as np
from three_row_hecke_crossing import coefficient_tensor,example_paths
from three_row_hecke_paths import M,F
from check_three_row_hecke_crossing import common_factor
from binary_extension_smith import ExtensionMatrixOracle
from finite_field_place_budget import all_place_bound,toeplitz_order_sum

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--alphabet',type=int,choices=[5,9],default=5);args=ap.parse_args();start=time.monotonic();q=args.alphabet;n=q*q
    T,construction=coefficient_tensor(example_paths(q));D=len(T)-1;assert common_factor(T)==[1]
    old=json.loads(Path('evidence/three-row-hecke-q'+str(q)+'-replay.json').read_text());assert construction['coefficient_sha256']==old['construction']['coefficient_sha256']
    mt=np.array(M,dtype=np.uint8);places=[]
    for point in list(range(8))+[None]:
        if point is None:A=T[::-1].copy()
        else:
            A=np.zeros_like(T)
            for k in range(D+1):
                for j in range(k+1):
                    if j&k==j:A[j]^=mt[F.power(point,k-j),T[k]]
        profiles=[];controls=[]
        for cut in (0,1):
            C=A.reshape(D+1,n,n)if cut==0 else A.transpose(0,1,3,2,4).reshape(D+1,n,n)
            oracle=ExtensionMatrixOracle([[C[:,i,j].tolist()for j in range(n)]for i in range(n)],11,n)
            orders=oracle.profile()['orders'];profiles.append(orders);controls.append(toeplitz_order_sum(C,orders))
            if point==0:assert orders==old['profiles'][cut]['zero']['orders']
            if point is None:assert orders==old['profiles'][cut]['regular_infinity_orders']
        if point==1:assert all(not any(p)for p in profiles)
        places.append({'point':point,'profiles':profiles,'independent_toeplitz_controls':controls})
        print(json.dumps({'alphabet':q,'point':point,'order_sums':list(map(sum,profiles))}),flush=True)
    bounds=[all_place_bound([p['profiles'][cut]for p in places],D,8)for cut in(0,1)]
    margin=D*n-sum(b['all_algebraic_point_order_bound']for b in bounds);assert margin>=0
    # A scalar polynomial with a nonrational root saturates the orbit bound.
    for k in range(1,6):assert all_place_bound([[0],[0],[0]],2*k,2)['all_algebraic_point_order_bound']==k
    bad=0
    for p,d,f in [([[0],[0]],2,2),([[3],[3],[3]],2,2),([[0],[-1],[0]],2,2),([[0],[0],[0]],2,3)]:
        try:all_place_bound(p,d,f)
        except ValueError:bad+=1
        else:raise AssertionError('Invalid degree table accepted')
    names=['check_hecke_all_algebraic_places.py','finite_field_place_budget.py','three_row_hecke_crossing.py','three_row_hecke_paths.py','check_three_row_hecke_crossing.py','binary_extension_smith.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp']
    out={'completed':True,'alphabet':q,'matrix_order':n,'degree':D,'primitive_scalar_polynomial':[1],'coefficient_sha256':construction['coefficient_sha256'],'places':places,'degree_budget_bounds':bounds,'moment_margin':margin,'exact_bare_two_place_infimum':1,'unweighted_supporting_cost_lower_bound':1,'all_positive_supporting_angles_excluded':True,'invalid_controls':bad,'source_sha256':{s:sha256(Path('tools',s).read_bytes()).hexdigest()for s in names},'scope':'Every pair of distinct algebraic projective places over an algebraic closure ofF8, either cut orientation, for the primitive retained tensor. Unweighted uniform-type two-place pressure at every positive supporting angle. Arbitrary physical gauges and mixed-place allocations are outside scope.'}
    Path('evidence/hecke-q'+str(q)+'-all-algebraic-places.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'alphabet':q,'all_place_order_bounds':[b['all_algebraic_point_order_bound']for b in bounds],'moment_margin':margin,'seconds':time.monotonic()-start}))
if __name__=='__main__':main()
