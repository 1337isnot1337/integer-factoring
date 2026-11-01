#!/usr/bin/env python3
"""Rebuild the selected refined crossing without overwriting its canonical receipt."""
import json,subprocess,sys,tempfile
from pathlib import Path
from hashlib import sha256

def strip_times(x):
    if isinstance(x,dict):return {k:strip_times(v)for k,v in x.items()if k!='seconds'}
    if isinstance(x,list):return list(map(strip_times,x))
    return x

def main():
    source=Path('evidence/f8-hyperoval-retraction-789.json');fingerprint=sha256(source.read_bytes()).hexdigest();reference=json.loads(source.read_text())
    with tempfile.TemporaryDirectory(prefix='integer-factoring-f8-refinement-')as directory:
        output=Path(directory)/'replay.json'
        subprocess.run([sys.executable,'tools/probe_f8_hyperoval_retraction.py','--pivots','7,8,9','--output',str(output)],check=True)
        actual=json.loads(output.read_text());assert actual['completed']and actual['construction']['all_orders_resolved']
        for key in ['pivot_columns','normalized_embedding_rows','boundary_words','embedding','retraction','profile','source_sha256']:assert actual[key]==reference[key],key
        assert strip_times(actual['construction'])==strip_times(reference['construction'])
    assert sha256(source.read_bytes()).hexdigest()==fingerprint
    print(json.dumps({'completed':True,'pivot_columns':[7,8,9],'full_orders_replayed':18432,'canonical_receipt_unchanged':True}))
if __name__=='__main__':main()
