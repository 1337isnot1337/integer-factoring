#!/usr/bin/env python3
"""Exact deletion-module checks and independent weighted-oracle controls."""
import json
from pathlib import Path
from hashlib import sha256
import numpy as np
from coherent_hot_compounds_f4 import mul,dot,axpy
from check_coherent_hot_f4 import maps_q8
from weighted_coherent_hot_f4 import CoherentWeightedOracle
from weighted_hot_code_oracle import HotWeightedOracle
from search_coherent_hot_retractions import pressure


def matmul(a,b):
    out=[[0]*len(b[0]) for _ in a]
    for i,row in enumerate(a):
        for k,x in enumerate(row):
            for j,y in enumerate(b[k]):out[i][j]^=mul(x,y)
    return out


def main():
    p,l,_=maps_q8();counts={'intertwining_vectors':0,'square_zero_operators':0,'commuting_pairs':0}
    for vectors,duals in ((p,l),(l,p)):
        operators=[]
        for bit in range(4):
            deleted=[{w^(1<<bit):c for w,c in v.items() if w>>bit&1} for v in vectors]
            matrix=[[dot(dual,d) for d in deleted] for dual in duals]
            for j,d in enumerate(deleted):
                rebuilt={}
                for i,v in enumerate(vectors):rebuilt=axpy(rebuilt,matrix[i][j],v)
                assert rebuilt==d;counts['intertwining_vectors']+=1
            assert not any(any(row) for row in matmul(matrix,matrix));counts['square_zero_operators']+=1
            operators.append(matrix)
        for i in range(4):
            for j in range(i):assert matmul(operators[i],operators[j])==matmul(operators[j],operators[i]);counts['commuting_pairs']+=1
    o=CoherentWeightedOracle(p,l,4);base=pressure(o.raw);value=o.evaluate(np.zeros(32),base['scale_t'])
    assert abs(value[0]-base['rho'])<1e-12
    A,b=o.cone();g=np.tile(np.array(o.heights)*.113,4)+np.random.default_rng(732).normal(0,.003,32);scale=.337
    assert min(A@g-b*scale)>0
    rho,dh,ds,_=o.evaluate(g,scale);eps=1e-6;errors=[]
    for j in range(32):
        step=np.zeros(32);step[j]=eps
        errors.append(abs((o.evaluate(g+step,scale)[0]-o.evaluate(g-step,scale)[0])/(2*eps)-dh[j]))
    errors.append(abs((o.evaluate(g,scale+eps)[0]-o.evaluate(g,scale-eps)[0])/(2*eps)-ds));assert max(errors)<2e-6
    # Compare to the previous independent coordinate-code wrapper, accounting
    # for its original-w scale versus this module's deflated-t scale.
    letters=sorted((w for w in range(16) if sum(2*i+1 for i in range(4) if w>>i&1)<=8),key=lambda w:(sum(2*i+1 for i in range(4) if w>>i&1),w))[:8]
    coordinate=[{w:1} for w in letters];ours=CoherentWeightedOracle(coordinate,coordinate,4);old=HotWeightedOracle(letters)
    checks=0
    for sw in (.17,.3):
        newv,newg,news,_=ours.log_profiles(g,2*sw);oldv,oldg,olds,_=old.log_profile_with_gradient(g,sw)
        assert abs(newv-oldv)<1e-11 and max(abs(newg-oldg))<1e-11 and abs(2*news-olds)<1e-10;checks+=1
    bad=-100*np.tile(np.array(o.heights),4);assert min(A@bad-b)<0
    try:o.evaluate(bad,1.)
    except ValueError:pass
    else:raise AssertionError('Outside-cone cost accepted.')
    names=[Path(__file__).name,'weighted_coherent_hot_f4.py','coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','weighted_hot_code_oracle.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp']
    report={'status':'passed','exact_deletion_module_controls':counts,'cone_inequalities':A.shape[0],'binary_weighted_wrapper_comparisons':checks,'smooth_chamber_gradient_coordinates':33,'maximum_gradient_error':max(errors),'outside_cone_control_rejected':True,'source_sha256':{n:sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},'scope':'Exact module identities; numerical weighted-wrapper and derivative controls. No exact pressure improvement claimed.'}
    Path('evidence/weighted-coherent-hot-f4-controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
