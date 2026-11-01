#!/usr/bin/env python3
"""Compute the complete coordinate core shared by height-25 boundary retracts.

Only one cut is computed: the coordinate embedding equals its retraction,
so the two coherent compound formulas are literally identical. No score is
claimed for the larger boundary class by this construction alone.
"""
import json,time
from pathlib import Path
from hashlib import sha256
from coherent_hot_compounds import height
from coherent_hot_compounds_fast import compound_block


def main():
    modes=12;words=sorted((w for w in range(1<<modes) if height(w)<25),key=lambda w:(height(w),w));assert len(words)==92
    pairs=[(x,y) for x in words for y in words];blocks=[];orders=[];start=time.monotonic();path=Path('evidence/hot-height25-common-core.json')
    names=[Path(__file__).name,'coherent_hot_compounds_fast.py','coherent_hot_compounds_fast.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp']
    report={'completed':False,'words':words,'heights':[height(w) for w in words],'q':92,'modes':modes,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':__doc__,'blocks':blocks}
    for charge in sorted({x.bit_count()+y.bit_count() for x,y in pairs}):
        chosen=[(x,y) for x,y in pairs if x.bit_count()+y.bit_count()==charge]
        rows=[(x|(y<<modes),) for x,y in chosen];cols=[(y|(x<<modes),) for x,y in chosen]
        b=compound_block(modes,rows,cols);b['charge']=charge;blocks.append(b);orders+=b['orders_t']
        report['elapsed_seconds']=time.monotonic()-start;path.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'charge':charge,'size':len(chosen),'seconds':report['elapsed_seconds']}),flush=True)
    assert len(orders)==92**2;report.update(completed=True,orders_t=sorted(orders),paired_cuts_identical=True);path.write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
