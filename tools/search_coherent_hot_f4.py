#!/usr/bin/env python3
"""Finite F4 rank-codimension-one boundary search, with complete paired profiles.

Numerical pressure is a discovery diagnostic, never an exact rate certificate.
"""
import argparse,json,time,random
from itertools import product
from pathlib import Path
from hashlib import sha256
from coherent_hot_compounds_f4 import profiles,mul,height
from search_coherent_hot_retractions import pressure


def dot(a,b):
    z=0
    for x,y in zip(a,b):z^=mul(x,y)
    return z


def catalogue(d):
    for k in product(range(4),repeat=d):
        if not any(k) or next(v for v in k if v)!=1:continue
        for c in product(range(4),repeat=d):
            if dot(k,c)==1:yield k,c


def maps(k,c,boundary,base):
    pivot=next(i for i,v in enumerate(c) if v);inv=mul(c[pivot],c[pivot]);keep=[i for i in range(len(k)) if i!=pivot]
    ps=[];ls=[]
    for i in keep:
        p={boundary[i]:1};v=mul(c[i],inv)
        if v:p[boundary[pivot]]=v
        l={boundary[j]:mul(k[i],c[j]) for j in range(len(k)) if mul(k[i],c[j])}
        l[boundary[i]]=l.get(boundary[i],0)^1
        if not l[boundary[i]]:del l[boundary[i]]
        ps.append(p);ls.append(l)
    return [{w:1} for w in base]+ps,[{w:1} for w in base]+ls


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--family',choices=['q8','q32'],default='q32');ap.add_argument('--trials',type=int,default=100);ap.add_argument('--seconds',type=float,default=1800);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--seed',type=int,default=20261101);args=ap.parse_args()
    modes,level=(4,8) if args.family=='q8' else (8,16)
    lower=sorted((w for w in range(1<<modes) if height(w)<level),key=lambda w:(height(w),w))
    boundary=[w for w in range(1<<modes) if height(w)==level and w.bit_count()==2]
    base=lower+([] if args.family=='q8' else [15])
    queue=list(catalogue(len(boundary)));random.Random(args.seed).shuffle(queue)
    # First retain an explicit coordinate projector as the in-process baseline.
    k=(1,)+(0,)*(len(boundary)-1);c=k;queue.remove((k,c));queue.insert(0,(k,c))
    names=['coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp',Path(__file__).name,'coherent_hot_compounds_fast.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp']
    report={'family':args.family,'field':'F4: omega^2=omega+1','catalogue_size':len(queue),'source_sha256':{n:sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},'records':[],'completed':False,'scope':__doc__}
    start=time.monotonic();best=10
    for k,c in queue[:args.trials]:
        if time.monotonic()-start>=args.seconds:break
        p,l=maps(k,c,boundary,base);result=profiles(p,l,modes);fit=pressure(result)
        row={'kernel':k,'covector':c,'zero_orders_t':result['zero_orders_t'],'infinity_orders_t':result['infinity_orders_t'],'heights':result['heights'],'pressure':fit,'seconds':result['seconds']}
        report['records'].append(row)
        if fit['rho']<best:
            best=fit['rho'];report['best_index']=len(report['records'])-1;print(json.dumps({'trial':len(report['records']),'best':best,'kernel':k,'covector':c,'seconds':time.monotonic()-start}),flush=True)
        report['elapsed_seconds']=time.monotonic()-start
        temp=args.output.with_suffix('.tmp');temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(args.output)
    report['completed']=True;report['catalogue_exhausted']=len(report['records'])==len(queue)
    report['elapsed_seconds']=time.monotonic()-start
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'trials':len(report['records']),'best':best,'catalogue_exhausted':report['catalogue_exhausted'],'seconds':report['elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
