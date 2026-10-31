#!/usr/bin/env python3
"""Original Fock and rational controls for growing physical Gaussian frames."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_rectangular_defect_support import Z
from check_reflection_unrestricted_gaussian import (
    conjugation_matrix, from_polynomial, rational_smith,
)
from reflection_gaussian_physical_frames import (
    frame_matrix, literal_tensor, literal_profile, predicted,
)


def check():
    rng=np.random.default_rng(55631)
    records=[]
    independent=0
    conjugations=0
    cases=[[{'leg':leg,'kind':kind,'powers':(-3,2),'gamma':gamma,'matching':'reverse'}]
           for leg in range(4) for kind in ('create','annihilate') for gamma in (1,2)]
    cases += [[{'leg':0,'kind':'create','powers':(-5,1),'gamma':1},
               {'leg':0,'kind':'annihilate','powers':(-1,4),'gamma':1}],
              [{'leg':3,'kind':'annihilate','powers':(-4,2),'gamma':3},
               {'leg':3,'kind':'create','powers':(0,5),'gamma':2},
               {'leg':1,'kind':'create','powers':(-2,3),'gamma':1}]]
    for M in (4,8):
        selected=cases if M==4 else [cases[-2],cases[-1],
            [{'leg':2,'kind':'create','powers':(-3,5),'gamma':2,'matching':'all'}]]
        for index,gates in enumerate(selected):
            values,clearing=literal_tensor(M,gates)
            degree=max(max(x).bit_length()-1 for x in values)
            negative=M//4*sum(max(0,-min(gate['powers'])) for gate in gates)
            positive=M//4*sum(max(0,max(gate['powers'])) for gate in gates)
            assert clearing<=negative and degree<=3*M*M//2+negative+positive
            profiles=[]
            for side in (0,1):
                heights=rng.integers(-2*M,2*M+1,2*M).tolist()
                actual=literal_profile(M,side,heights,values)
                expected=[x+(clearing if side==0 else -clearing)
                          for x in predicted(M,side,heights,gates)]
                assert actual==expected,(M,side,gates,actual,expected)
                record={'side':side,'heights':heights,'fock_orders':actual}
                if M==4 and index in (0,12,16):
                    assert all(high==0 for low,high in values)
                    b,rows,columns=frame_matrix(M,side,gates)
                    # This control uses the Clifford mode order directly;
                    # the production full-tensor flattening puts the first
                    # physical letter in the more significant block.
                    row_masks=[sum(((state>>k)&1)<<physical for k,physical in enumerate(rows))
                               for state in range(1<<M)]
                    col_masks=[sum(((state>>k)&1)<<physical for k,physical in enumerate(columns))
                               for state in range(1<<M)]
                    raw=[[values[r|c] for c in col_masks] for r in row_masks]
                    fock=[[from_polynomial(low,infinity=bool(side)) for low,high in row]
                          for row in raw]
                    row_weights=[sum(heights[physical] for k,physical in enumerate(rows)
                                     if state>>k&1) for state in range(1<<M)]
                    col_weights=[sum(heights[physical] for k,physical in enumerate(columns)
                                     if state>>k&1) for state in range(1<<M)]
                    weighted=[[value*Z**(row_weights[i]+col_weights[j])
                               for j,value in enumerate(row)] for i,row in enumerate(fock)]
                    assert rational_smith(weighted)==actual
                    independent+=1
                    # No scalar clearing affects this conjugation identity.
                    exact_b=[[from_polynomial(low,shift,bool(side)) for low,high,shift in row]
                             for row in b]
                    assert conjugation_matrix(fock,M)==exact_b
                    conjugations+=1
                    record['independent_rational_fock']=True
                    record['full_conjugation_identity']=True
                profiles.append(record)
            records.append({'M':M,'gates':gates,'clearing':clearing,'degree':degree,
                            'profiles':profiles})
    return {'status':'exact original-tensor controls passed',
            'frame_words':len(records),'complete_profiles':2*len(records),
            'independent_rational_fock_profiles':independent,
            'full_conjugation_identities':conjugations,
            'records':records}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=check()
    root=Path(__file__).resolve().parents[1]
    sources=['research/growing-gaussian-physical-frames.md',
             'tools/reflection_gaussian_physical_frames.py',
             'tools/check_reflection_gaussian_physical_frames.py',
             'tools/reflection_unrestricted_gaussian.py',
             'tools/check_reflection_unrestricted_gaussian.py']
    result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in sources}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items()
                      if key not in ('records','source_sha256')},sort_keys=True),flush=True)


if __name__=='__main__':
    main()
