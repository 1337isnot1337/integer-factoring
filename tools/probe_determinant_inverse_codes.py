#!/usr/bin/env python3
"""Registered fixed determinant-code panel; exact profiles, numerical fits.

Run from a clean frozen checkout with --output outside the checkout.
No approximate profile or numerical fit is a factoring bound certificate.
"""
import argparse
import ctypes
import datetime
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import platform
import subprocess
import time

from binary_extension_smith import ExtensionMatrixOracle
from coherent_f8_maps import mul,tensor_vector
from coherent_hot_compounds import height,one_particle
from coherent_hot_compounds_f8 import library,partial_profiles
from determinant_inverse_codes import code,serializable

ROOT=Path(__file__).resolve().parents[1]


def pmul(first,second):
    value=0
    while second:
        bit=second & -second;value^=first << (bit.bit_length()-1);second^=bit
    return value


@lru_cache(None)
def determinant(modes,rowmask,colmask):
    rows=[i for i in range(2*modes) if rowmask >> i & 1]
    cols=[i for i in range(2*modes) if colmask >> i & 1]
    assert len(rows)==len(cols);matrix=one_particle(modes);values={0:1}
    for row in rows:
        updated={}
        for mask,value in values.items():
            for j,col in enumerate(cols):
                if not(mask >> j & 1) and matrix[row][col]:
                    target=mask | (1 << j)
                    updated[target]=updated.get(target,0)^pmul(value,matrix[row][col])
        values=updated
    return values.get((1 << len(cols))-1,0)


def entry_controls(data):
    k=data['k'];modes=k*k;p=data['embedding'];left=data['retraction'];q=len(p)
    charges=[x['charge'] for x in data['labels']];records=[];entries=nonbinary=0
    for charge in sorted({a+b for a in charges for b in charges}):
        available=[(i,j) for i in range(q) for j in range(q) if charges[i]+charges[j]==charge]
        indices=list(range(len(available))) if k==2 else sorted({0,len(available)//2,len(available)-1})
        chosen=[available[i] for i in indices];n=len(chosen)
        for cut in ('zero','infinity'):
            rows=[tensor_vector(left[x],left[y] if cut=='zero' else p[y],modes) for x,y in chosen]
            cols=[tensor_vector(p[y] if cut=='zero' else left[y],p[x],modes) for x,y in chosen]
            low=(1 << modes)-1
            rowh=[height(next(iter(v)) & low)+height(next(iter(v)) >> modes) for v in rows]
            colh=[height(next(iter(v)) & low)+height(next(iter(v)) >> modes) for v in cols]
            bound=max(rowh)+max(colh);cutoff=max(127,bound)
            words=1
            while 64*words<=cutoff:words*=2
            def pack(vectors):
                offsets=[0];masks=[];scalars=[]
                for v in vectors:
                    for w,c in sorted(v.items()):masks.append(w);scalars.append(c)
                    offsets.append(len(masks))
                return tuple((ctypes.c_uint32*len(v))(*v) for v in (offsets,masks,scalars))
            orders=(ctypes.c_int64*n)();out=(ctypes.c_uint64*(3*n*n*words))()
            status=library(words,cutoff).coherent_hot_compound_f8(modes,n,*pack(rows),*pack(cols),orders,out)
            assert status==0,status
            matrix=[];polynomials=[]
            for i,row in enumerate(rows):
                line=[];stored=[]
                for j,col in enumerate(cols):
                    planes=[0,0,0]
                    for a,ca in row.items():
                        for b,cb in col.items():
                            term=determinant(modes,a,b);scalar=mul(ca,cb)
                            for bit in range(3):
                                if scalar >> bit & 1:planes[bit]^=term
                    native=[sum(int(out[(3*(i*n+j)+bit)*words+w]) << (64*w)
                                for w in range(words)) for bit in range(3)]
                    assert native==planes,(k,data['kind'],cut,charge,i,j)
                    endpoint=sum((v.bit_count()%2) << bit for bit,v in enumerate(planes))
                    assert endpoint==int(i==j)
                    assert max(v.bit_length()-1 for v in planes)<=bound
                    line.append([sum(((v >> degree)&1) << bit for bit,v in enumerate(planes))
                                 for degree in range(max(v.bit_length() for v in planes))])
                    stored.append(planes);entries+=1;nonbinary+=bool(planes[1] or planes[2])
                matrix.append(line);polynomials.append(stored)
            independent=ExtensionMatrixOracle(matrix,11,n).profile()['orders']
            assert sorted(int(v) for v in orders if v<=cutoff)==[v for v in independent if v<=cutoff]
            records.append(dict(cut=cut,charge=charge,selected_pairs=chosen,
                all_charge_pairs_checked=k==2,degree_bound=bound,complete_polynomial_planes=polynomials,
                independent_binary_restriction_orders=independent,
                native_retained_orders=sorted(int(v) for v in orders if v<=cutoff)))
    return dict(entries=entries,nonbinary_entries=nonbinary,records=records)


def lse(values):
    top=max(values)
    return top+math.log(sum(math.exp(v-top) for v in values))


def fit(profiles,heights):
    orders=[c['retained_orders_t'] for c in profiles['cuts']];queries=0
    def evaluate(logscale):
        nonlocal queries
        queries+=1;scale=math.exp(logscale)
        raw=sum(lse([-scale*v for v in cut]) for cut in orders)
        def equation(rho):return raw-4*rho*lse([-scale*h/rho for h in heights])
        low,high=1e-10,2.
        assert equation(low)>0
        while equation(high)>0:high*=2
        for _ in range(64):
            middle=(low+high)/2
            if equation(middle)>0:low=middle
            else:high=middle
        return (low+high)/2
    start,end=-9.,math.log(4.)
    grid=[start+(end-start)*i/80 for i in range(81)]
    values=[evaluate(v) for v in grid]
    candidates=[(values[0],grid[0]),(values[-1],grid[-1])]
    ratio=(math.sqrt(5)-1)/2
    for i in range(1,80):
        if not(values[i]<=values[i-1] and values[i]<=values[i+1]):continue
        a,b=grid[i-1],grid[i+1];c=b-ratio*(b-a);d=a+ratio*(b-a)
        fc,fd=evaluate(c),evaluate(d)
        for _ in range(64):
            if fc<fd:b,d,fd=d,c,fc;c=b-ratio*(b-a);fc=evaluate(c)
            else:a,c,fc=c,d,fd;d=a+ratio*(b-a);fd=evaluate(d)
        candidates.append(min((fc,c),(fd,d)))
    rho,logscale=min(candidates);scale=math.exp(logscale)
    z=lse([-scale*h/rho for h in heights]);pi=[math.exp(-scale*h/rho-z) for h in heights]
    return dict(rho=rho,scale=scale,physical_probabilities=pi,
        entropy=-sum(p*math.log(p) for p in pi),queries=queries,
        logscale_grid=[start,end,81],global_minimum_certified=False,
        strictly_below_critical_numerically=rho<1/math.sqrt(2)-1e-7)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    def git(*parts):return subprocess.check_output(['git',*parts],cwd=ROOT,text=True).strip()
    paths=['tools/probe_determinant_inverse_codes.py','tools/determinant_inverse_codes.py',
           'tools/coherent_f8_maps.py','tools/coherent_hot_compounds.py',
           'tools/coherent_hot_compounds_f8.py','tools/coherent_hot_compounds_f8.cpp',
           'tools/selected_gaussian_coefficients.cpp','tools/native_rectangular_smith_f4.cpp',
           'tools/binary_extension_smith.py','tools/native_matrix_smith_f4.py',
           'tools/native_matrix_smith_f4.cpp']
    data=dict(status='running',scope='Six prescribed determinant-code cases; complete exact native profiles and selected independent entries, numerical scale fits only. No score claim.',
        source_commit=git('rev-parse','HEAD'),source_tree=git('rev-parse','HEAD^{tree}'),
        pre_status=git('status','--porcelain'),started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        python=platform.python_version(),source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        cases=[],factoring_exponent_improved=False)
    assert not data['pre_status'];begun=time.monotonic()
    def save():
        temporary=args.output.with_suffix('.checkpoint.tmp')
        temporary.write_text(json.dumps(serializable(data),indent=2,sort_keys=True)+'\n');temporary.replace(args.output)
    save()
    try:
        for k in (2,3):
            for kind in ('diagonal','antidiagonal','cauchy'):
                acquired=code(k,kind);record=dict(k=k,kind=kind,maps=acquired,attempts=[])
                data['cases'].append(record);save()
                hs=[x['height'] for x in acquired['labels']];ch=[x['charge'] for x in acquired['labels']]
                detbound=max(2*sum(hs[x]+hs[y] for x in range(len(hs)) for y in range(len(hs))
                    if ch[x]+ch[y]==charge) for charge in set(a+b for a in ch for b in ch))
                record['determinant_degree_bound']=detbound
                cutoff=min(127,detbound)
                while True:
                    def progress(event):print(json.dumps(dict(k=k,kind=kind,**event)),flush=True)
                    result=partial_profiles(acquired['embedding'],acquired['retraction'],k*k,cutoff,progress)
                    record['attempts'].append(result);save()
                    if result['all_orders_resolved']:break
                    assert cutoff<detbound,'Endpoint/determinant-degree bound failed to resolve all orders'
                    cutoff=min(2*cutoff+1,detbound)
                record['complete_profile']=result
                record['independent_controls']=entry_controls(acquired)
                if k==2:
                    for checked in record['independent_controls']['records']:
                        cut=next(c for c in result['cuts'] if c['cut']==checked['cut'])
                        charge=checked['charge']
                        block=cut['blocks'][charge]
                        assert block['retained_orders_t']==checked['independent_binary_restriction_orders']
                record['numerical_fit']=fit(result,hs);save()
                print(json.dumps(dict(k=k,kind=kind,complete=True,
                    rho=record['numerical_fit']['rho'],controls=record['independent_controls']['entries'])),flush=True)
        data['status']='complete'
    except Exception as error:
        data['status']='failed';data['failure']=repr(error);save();raise
    finally:
        data.update(post_status=git('status','--porcelain'),seconds=time.monotonic()-begun,
                    finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    distinct_original_matrix_minors=determinant.cache_info().currsize)
        save()
    assert not data['post_status']
    print(json.dumps(dict(status=data['status'],seconds=data['seconds'],cases=len(data['cases']))),flush=True)


if __name__=='__main__':main()
