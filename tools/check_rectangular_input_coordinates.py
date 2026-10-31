#!/usr/bin/env python3
"""Replay input repairs and full inverses against independent original tensors."""
import argparse
import hashlib
import json
from pathlib import Path
from check_rectangular_f4_repair import build
from check_rectangular_refined_repairs import single_shear,precorrect_seed,postcorrect_seed
from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import point_control,weighted_controls
from rectangular_defect_coordinates import coordinate_matrix

def check():
    base=build(tuple(range(4)),4,1,2,check_types=False)
    cases=[]
    for row,column in ((3,8),(8,3),(10,11),(1,0)):
        spec=single_shear(row,column,2)+(column!=0,)
        seed=precorrect_seed(base,*spec[:3],check_types=False,keep_background=spec[3])
        cases.append((f'input_{row}_{column}',seed,{'input_spec':spec}))
    spec=single_shear(3,8,2)
    joint=precorrect_seed(postcorrect_seed(base,*spec,check_types=False),*spec,check_types=False)
    cases.append(('common_flag_both_sides',joint,{'repair_spec':spec,'input_spec':spec}))
    for alpha,beta in ((1,1),(1,2)):
        seed=build(tuple(range(4)),4,alpha,beta,repair='inverse',check_types=False)
        cases.append((f'full_inverse_{alpha}_{beta}',seed,{'alpha':alpha,'beta':beta,'repair_kind':'inverse'}))
    records=[]
    for name,seed,parameters in cases:
        for place in ('infinity','zero'):
            data=data_for(place,4,full_core=True)
            compiled=coordinate_matrix(data,**parameters)
            entry=point_control(data,compiled,2,seed)
            profiles=weighted_controls(data,compiled,2,seed)
            records.append({'case':name,'entrywise':entry,'weighted_profiles':profiles})
            print(json.dumps({'case':name,'place':place,'entries':entry['entries_checked'],'status':'pass'}),flush=True)
    names=['rectangular_defect_profiles.py','rectangular_defect_coordinates.py','check_rectangular_refined_repairs.py','check_rectangular_f4_repair.py','check_rectangular_defect_fock.py','check_rectangular_defect_coordinates.py',Path(__file__).name]
    return {'status':'pass','scope':'Independent full M4 tensor entries at an F16 point and complete arbitrary-core weighted profiles for four input shears, a two-sided common flag, and two full inverses. No new rate claimed.','cases':records,'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=check();rendered=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered)

if __name__=='__main__':main()
