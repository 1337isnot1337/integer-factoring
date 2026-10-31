#!/usr/bin/env python3
"""Cache locally constructed coordinate matrices under their source digest.

Cache files contain Python pickle data and must be generated locally by this
module. Never load a cache supplied by an untrusted source.
"""
from __future__ import annotations
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import tempfile
from check_rectangular_defect_core_basis import data_for
from rectangular_defect_coordinates import coordinate_matrix

SOURCES=['rectangular_defect_profiles.py','rectangular_defect_coordinates.py','native_fock_f4.cpp','check_rectangular_defect_core_basis.py','check_rectangular_defect_support.py','check_rectangular_full_core_support.py','check_rectangular_f4_repair.py','check_rectangular_refined_repairs.py','check_rectangular_identity_transports.py']

def cached_coordinates(directory,M=8,alpha=1,beta=2,quartic_terms=None,full_core=False,repair_kind='partial',profile_precision=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    terms=[tuple(term) for term in (quartic_terms or [(1,0)])]
    source_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in SOURCES}
    inputs={'M':M,'alpha':alpha,'beta':beta,'quartic_terms':terms,'full_core':full_core,'repair_kind':repair_kind,'profile_precision':profile_precision,'source_sha256':source_hashes}
    digest=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    filename=directory/(digest+'.pickle.gz')
    if filename.exists():
        with gzip.open(filename,'rb') as stream:result=pickle.load(stream)
        if result['digest']!=digest or result['inputs']!=inputs:raise ValueError('Coordinate cache provenance mismatch.')
        return result['data'],result['compiled'],{'digest':digest,'cache_hit':True,'path':str(filename),'inputs':inputs}
    data=[data_for(place,M,full_core=full_core) for place in ('infinity','zero')]
    compiled=[coordinate_matrix(d,alpha,beta,terms,repair_kind=repair_kind,profile_precision=profile_precision) for d in data]
    result={'digest':digest,'inputs':inputs,'data':data,'compiled':compiled}
    fd,temporary=tempfile.mkstemp(prefix=digest+'-',suffix='.tmp',dir=directory);os.close(fd)
    try:
        with gzip.open(temporary,'wb',compresslevel=1) as stream:pickle.dump(result,stream,protocol=5)
        os.replace(temporary,filename)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)
    return data,compiled,{'digest':digest,'cache_hit':False,'path':str(filename),'inputs':inputs}

def divided_difference(terms):
    """Return Q with (1+w)Q=P+1 for a Laurent P satisfying P(1)=1."""
    values={}
    for coefficient,exponent in terms:
        if coefficient not in (0,1,2,3) or not isinstance(exponent,int):raise ValueError('Invalid Laurent numerator.')
        values[exponent]=values.get(exponent,0)^coefficient
    endpoint=0
    for coefficient in values.values():endpoint^=coefficient
    if endpoint!=1:raise ValueError('The quartic numerator must have P(1)=1.')
    values[0]=values.get(0,0)^1;values={k:v for k,v in values.items() if v}
    if not values:return []
    low,high=min(values),max(values);previous=0;answer=[]
    for exponent in range(low,high):
        previous^=values.get(exponent,0)
        if previous:answer.append((previous,exponent))
    if previous!=values[high]:raise ArithmeticError('Laurent divided difference has a nonzero remainder.')
    return answer

def quartic_from_two_bases(base_one,base_w,terms):
    from rectangular_coordinate_transforms import linear_combination_coordinate_matrices
    quotient=divided_difference(terms)
    if not quotient:return base_one
    results=[]
    for one,w in zip(base_one,base_w):
        sign=-1 if one[2]['place']=='infinity' else 1
        q=[(coefficient,sign*power) for coefficient,power in quotient]
        results.append(linear_combination_coordinate_matrices([(one,[(1,0)]+q),(w,q)]))
    return results
