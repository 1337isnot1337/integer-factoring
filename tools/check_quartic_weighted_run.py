#!/usr/bin/env python3
"""Check a floating-search archive's integrity, not its mathematical optimality."""
import collections
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def lse(values):
    m=max(values)
    return m+math.log(sum(math.exp(x-m) for x in values))


def main():
    data=json.loads((ROOT/'evidence/quartic-weighted-run.json').read_text())
    assert data['source_commit']=='3ddcf36f9ab335893f0662010e5b5783340dc37f'
    assert data['source_tree']=='e119e09624810d814a9f91588f715a53c3805fee'
    assert data['status']=='complete_no_witness'
    assert data['pre_status']==data['post_status']==''
    assert data['seed']==0x51554152544943
    assert data['query_budget']==50000 and data['max_iterations']==250
    assert data['witness_threshold']==-1e-7
    for path,digest in data['source_sha256'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    compressed=(ROOT/'evidence/quartic-weighted-run.queries.jsonl.gz').read_bytes()
    assert sha256(compressed).hexdigest()==data['compressed_query_sha256']
    plain=gzip.decompress(compressed)
    assert sha256(plain).hexdigest()==data['uncompressed_query_sha256']
    queries=[json.loads(line) for line in plain.splitlines()]
    assert data['attempted_queries']==data['completed_queries']==len(queries)==5995
    assert data['rejected_queries']==0
    assert len(data['runs'])==64 and len(data['gradient_controls'])==4
    assert max(c['maximum_error'] for c in data['gradient_controls'])<2e-9
    rate=2*math.sqrt(2);counts=collections.Counter()
    for index,q in enumerate(queries,1):
        assert q['query']==index
        x=q['params'];assert len(x)==25
        assert all(lo<=v<=hi for v,(lo,hi) in zip(x,data['bounds']))
        w=[0.]+[sum(x[12*e+3*l+a-1] for e in range(2) for l in range(4))
                  for a in range(1,4)]
        assert q['weights']==w and len(q['spectra'])==32 and len(q['gradient'])==25
        assert all(math.isfinite(v) for v in x+q['spectra']+q['gradient'])
        raw=lse([-v for v in q['spectra'][:16]])+lse([-v for v in q['spectra'][16:]])
        assert abs(raw-q['raw_pressure'])<1e-10
        z=lse([-v/rate for v in w]);margin=raw-rate*z
        assert abs(margin-q['margin'])<1e-10
        assert margin>=data['witness_threshold']
        assert all(abs(p-math.exp(-v/rate-z))<1e-14
                   for p,v in zip(q['probabilities'],w))
        counts[q['tensor'],q['cut'],q['start']]+=1
    for run in data['runs']:
        assert counts[run['tensor'],run['cut'],run['start']]==run['completed_queries']
        assert run['status'] in (0,1,2)
    for tensor in ('Gaussian baseline','S149:x4+y4'):
        for cut in (0,1):
            assert counts[tensor,cut,'gradient control']==51
            assert {r['start'] for r in data['runs']
                    if r['tensor']==tensor and r['cut']==cut}==set(range(16))
    best=min(queries,key=lambda q:q['margin'])
    assert abs(best['margin']-data['best']['margin'])<1e-14
    print('PASS: complete 64-start/5995-query floating-search archive and hashes; no strict witness.')
    print('This checks archive integrity and recorded arithmetic, not an all-gauge floor or new score.')


if __name__=='__main__':main()
