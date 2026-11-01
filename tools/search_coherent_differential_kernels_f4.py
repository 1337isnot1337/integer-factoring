#!/usr/bin/env python3
"""Complete-profile F4 interior differential-code experiments.

Two quadratic equations use four distinct F4 projective columns; this
rank-two four-column MDS pattern has no binary counterpart. Every completed
trial retains both full profiles. Numerical pressure is not a certificate.
"""
import argparse,json,time
from pathlib import Path
from hashlib import sha256
from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds_f4 import profiles
from search_coherent_hot_retractions import pressure


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--caps',type=int,nargs='+',default=[20,22,24]);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    words=[24,36,66,129]
    families={'mds_two':[dict.fromkeys(words,1),{words[i]:i for i in range(1,4)}], 'weighted_one':[{words[i]:c for i,c in enumerate([1,1,2,3])}]}
    names=[Path(__file__).name,'coherent_differential_kernels_f4.py','coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','search_coherent_hot_retractions.py']
    report={'scope':__doc__,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'records':[],'completed':False,'in_progress':None};start=time.monotonic()
    def save():
        temp=args.output.with_suffix('.tmp');temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(args.output)
    for cap in args.caps:
        for name,ops in families.items():
            code=differential_kernel(ops,12,cap);report['in_progress']={'family':name,'height_cap':cap,'code':code};save()
            print(json.dumps({'starting':name,'cap':cap,'q':code['q']}),flush=True)
            result=profiles(code['embedding'],code['retraction'],12);fit=pressure(result)
            report['records'].append({'family':name,'height_cap':cap,'code':code,'profile':result,'pressure':fit});report['in_progress']=None;report['elapsed_seconds']=time.monotonic()-start;save()
            print(json.dumps({'completed':name,'cap':cap,'q':code['q'],'rho':fit['rho'],'seconds':report['elapsed_seconds']}),flush=True)
    report['completed']=True;save()

if __name__=='__main__':main()
