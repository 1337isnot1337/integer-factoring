#!/usr/bin/env python3
"""Export complete coordinate-entry valuations for exploratory scale fitting.

Every exported matrix is a complete polynomial numerator with a valuation-unit
denominator. No absent entry is inferred from a truncated coefficient window.
The export itself is not a bivariate or all-length profile certificate.
"""
import argparse
import hashlib
import json
import gzip
import pickle
from pathlib import Path
import time
from check_rectangular_defect_core_basis import data_for
from rectangular_defect_coordinates import coordinate_matrix

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes',nargs='+',type=int,default=[16,32,64])
    parser.add_argument('--full-core',action='store_true')
    parser.add_argument('--save-compiled',action='store_true')
    parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args();args.output_directory.mkdir(parents=True,exist_ok=True)
    names=['rectangular_defect_profiles.py','rectangular_defect_coordinates.py','check_rectangular_defect_core_basis.py','check_rectangular_defect_support.py','native_fock_f4.cpp',Path(__file__).name]
    sources={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    for M in args.modes:
        for place in ('infinity','zero'):
            started=time.monotonic();data=data_for(place,M,full_core=args.full_core)
            precision=max(512,1<<(64*M-1).bit_length())
            if precision>8192:raise ValueError('Requested exact export exceeds the supported window.')
            while True:
                matrix,shifts,metadata=coordinate_matrix(data,profile_precision=precision)
                if metadata['exact_polynomial_matrix']:break
                if precision==8192:raise ArithmeticError('The maximum window did not produce a complete numerator.')
                precision*=2
            values=[[shifts[i][j] if value[0]|value[1] else None for j,value in enumerate(row)] for i,row in enumerate(matrix)]
            result={'status':'complete finite polynomial entry valuations; scale fitting remains exploratory','M':M,'place':place,'full_core':args.full_core,'valuations':values,'metadata':metadata,'active_orders':data['active_orders'],'active_indices':data['active_indices'],'spectator_orders':data['spectator_orders'],'source_sha256':sources}
            filename=args.output_directory/f'entries-{place}-M{M}.json'
            filename.write_text(json.dumps(result,separators=(',',':'))+'\n')
            if args.save_compiled:
                with gzip.open(args.output_directory/f'compiled-{place}-M{M}.pickle.gz','wb',compresslevel=1) as stream:
                    pickle.dump({'compiled':(matrix,shifts,metadata),'data':data,'source_sha256':sources},stream,protocol=5)
            print(json.dumps({'M':M,'place':place,'matrix_order':len(matrix),'nonzero_entries':sum(v is not None for row in values for v in row),'precision':precision,'seconds':time.monotonic()-started,'path':str(filename)}),flush=True)

if __name__=='__main__':main()
