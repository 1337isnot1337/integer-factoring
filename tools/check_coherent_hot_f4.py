#!/usr/bin/env python3
"""Independent literal F4 coefficient controls for coherent compound profiles."""
import json
from hashlib import sha256
from pathlib import Path
from coherent_hot_compounds_f4 import profiles,mul,height,validate_maps
from coherent_hot_compounds_fast import profiles as binary_profiles


def maps_q8(a=2,c=1):
    lower=sorted((w for w in range(16) if height(w)<8),key=lambda w:(height(w),w))
    boundary=[w for w in range(16) if height(w)==8 and w.bit_count()==2]
    assert len(lower)==7 and len(boundary)==2
    p={boundary[0]:1};l={}
    if a:p[boundary[1]]=a
    b=1^mul(a,c)
    if b:l[boundary[0]]=b
    if c:l[boundary[1]]=c
    return [{w:1} for w in lower]+[p],[{w:1} for w in lower]+[l],lower+boundary


def main():
    p,l,ambient=maps_q8();full=binary_profiles([1<<w for w in ambient],[1<<w for w in ambient],4,keep_blocks=True)
    out=profiles(p,l,4,keep_blocks=True);count=0;native_nonbinary=0
    for cut_index,cut in enumerate(('zero','infinity')):
        source_by_charge={sum(w.bit_count() for w in (ambient[b['pairs'][0][0]],ambient[b['pairs'][0][1]])):b for b in full['cuts'][cut_index]['blocks']}
        for block in out['cuts'][cut_index]['blocks']:
            # Match each charge block by its output size and order in the profile.
            index=out['cuts'][cut_index]['blocks'].index(block)
            charge=sorted({a+b for a in out['charges'] for b in out['charges']})[index]
            old=source_by_charge[charge];lookup={(ambient[x],ambient[y]):i for i,(x,y) in enumerate(old['pairs'])}
            pairs=[(x,y) for x in range(8) for y in range(8) if out['charges'][x]+out['charges'][y]==charge]
            for i,(x,y) in enumerate(pairs):
                for j,(u,v) in enumerate(pairs):
                    expected=[0,0]
                    r1,r2=l[x],l[y] if cut=='zero' else p[y]
                    c1,c2=p[u],p[v] if cut=='zero' else l[v]
                    for a,ca in r1.items():
                        for b,cb in r2.items():
                            for d,cd in c1.items():
                                for e,ce in c2.items():
                                    f=mul(mul(ca,cb),mul(cd,ce));value=old['matrix'][lookup[a,b]][lookup[d,e]]
                                    if f&1:expected[0]^=value
                                    if f&2:expected[1]^=value
                    assert tuple(expected)==block['matrix'][i][j],(cut,charge,i,j)
                    native_nonbinary+=bool(expected[1]);count+=1
    conjugate=lambda maps:[{w:mul(c,c) for w,c in v.items()} for v in maps]
    frob=profiles(conjugate(p),conjugate(l),4)
    scaledp=[{w:mul(c,2) for w,c in v.items()} for v in p]
    scaledl=[{w:mul(c,3) for w,c in v.items()} for v in l]
    scaled=profiles(scaledp,scaledl,4)
    for name in ('zero_orders_t','infinity_orders_t'):
        assert out[name]==frob[name]==scaled[name]
    binaryp=[{w:1} for w in ambient];fp=profiles(binaryp,binaryp,4)
    for name in ('zero_orders_t','infinity_orders_t'):assert fp[name]==full[name]
    rejected=0
    for pp,ll in [(p,p),([{0:1,1:2}], [{0:1}]),([{3:1}],[{3:1}])]:
        try:validate_maps(pp,ll,4)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid map accepted')
    report={'literal_F4_polynomial_entries':count,'nonbinary_entries':native_nonbinary,'binary_cut_profile_controls':2,'frobenius_and_rescaling_profile_controls':4,'invalid_maps_rejected':rejected,'q8_zero':out['zero_orders_t'],'q8_infinity':out['infinity_orders_t'],'status':'passed','scope':'Complete two-cut F4 coefficient and profile controls; no numerical improvement claimed.'}
    names=['coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp',Path(__file__).name,'coherent_hot_compounds_fast.py','coherent_hot_compounds_fast.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp']
    report['source_sha256']={n:sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}
    print(json.dumps(report));Path('evidence/coherent-hot-f4-controls.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
