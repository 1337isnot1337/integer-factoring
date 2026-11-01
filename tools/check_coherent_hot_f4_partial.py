#!/usr/bin/env python3
"""Compare truncated coherent profiles with complete exact calculations."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
import json
from coherent_differential_kernels_f4 import differential_kernel
from orthogonal_differential_pairs_f4 import orthogonal_retraction
from coherent_hot_compounds_f4 import profiles
from coherent_hot_f4_partial import partial_profiles


def main():
    cases=[];rejected=0
    for modes,cap,ops in [(4,16,[{9:1,6:2}]),(4,16,[{9:1,6:3}]),(8,20,[{129:1,66:1,36:2,24:2},{66:1,36:3,24:1}])]:
        code=differential_kernel(ops,modes,cap);P=code['embedding'];L=orthogonal_retraction(code);assert L is not None
        full=profiles(P,L,modes)
        for cutoff in (0,7,20,40):
            partial=partial_profiles(P,L,modes,cutoff);checks=[]
            for cut in partial['cuts']:
                values=full[cut['cut']+'_orders_t'];retained=[v for v in values if v<=cutoff]
                assert cut['retained_orders_t']==retained
                assert cut['tail_count']==sum(v>cutoff for v in values)
                for t in (Fraction(1,2),Fraction(9,10)):
                    actual=sum(t**v for v in values);lower=sum(t**v for v in retained);upper=lower+cut['tail_count']*t**(cutoff+1)
                    assert lower<=actual<=upper
                checks.append({'cut':cut['cut'],'retained_count':len(retained),'tail_count':cut['tail_count']})
            cases.append({'modes':modes,'q':code['q'],'cutoff':cutoff,'cuts':checks})
        for invalid in (-1,1.5):
            try:partial_profiles(P,L,modes,invalid)
            except ValueError:rejected+=1
            else:raise AssertionError('Invalid cutoff accepted.')
    names=[Path(__file__).name,'coherent_hot_f4_partial.py','coherent_hot_f4_partial.cpp','coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','orthogonal_differential_pairs_f4.py']
    report={'status':'passed','cases':cases,'complete_reference_constructions':3,'cut_comparisons':2*len(cases),'rational_partition_comparisons':4*len(cases),'invalid_cutoffs_rejected':rejected,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Exact orders through the chosen cutoff and rigorous channel-sum envelopes; omitted orders are not claimed known.'}
    Path('evidence/coherent-hot-f4-partial-controls.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','source_sha256')}))

if __name__=='__main__':main()
