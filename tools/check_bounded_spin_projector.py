#!/usr/bin/env python3
"""Exact full-word projector matrices versus independent integral path vectors.

NumPy is used only by this focused checker, not by the projector constructor.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths, recursive_vector, path_norm, binary_path
from bounded_spin_projector import STATES, projector_entry, binary_projector_entry
from check_bounded_spin_rank_algebra import rank

ROOT = Path(__file__).resolve().parents[1]


def direct_projector(n, ps):
    vectors = [recursive_vector(p) for p in ps]
    frame = np.array([[v.get(w,0) % 5 for w in range(1 << n)] for v in vectors],
                     dtype=np.int64)
    inverse = np.array([pow(path_norm(p),-1,5) for p in ps],dtype=np.int64)
    return (frame.T*inverse[None,:]) @ frame % 5


def check():
    assert len(STATES) == 30
    entries, cases = 0, []
    for n in (0,2,4,6,8):
        ps = paths(n)
        expected = direct_projector(n,ps)
        actual = np.array([[projector_entry(n,x,y) for y in range(1 << n)]
                           for x in range(1 << n)],dtype=np.int64)
        assert np.array_equal(actual,expected)
        assert np.array_equal(actual,actual.T)
        assert np.array_equal(actual @ actual % 5,actual)
        assert rank(actual) == len(ps)
        entries += actual.size
        cases.append({'length':n,'ambient_dimension':1 << n,'projector_rank':len(ps)})
    selected = []
    for width in (1,2,3):
        last = binary_path(width,(1 << width)-1)
        n = len(last)
        ps = paths(n)[:1 << width]
        assert last == ps[-1]
        expected = direct_projector(n,ps)
        actual = np.array([[binary_projector_entry(width,x,y) for y in range(1 << n)]
                           for x in range(1 << n)],dtype=np.int64)
        assert np.array_equal(actual,expected)
        assert np.array_equal(actual @ actual % 5,actual)
        assert rank(actual) == 1 << width
        entries += actual.size
        selected.append({'binary_width':width,'length':n,'projector_rank':1 << width})
    random = Random(1337)
    last = binary_path(4,15)
    n = len(last)
    vectors = [recursive_vector(p) for p in paths(n)[:16]]
    norms = [path_norm(p) for p in paths(n)[:16]]
    for _ in range(128):
        x,y = random.randrange(1 << n),random.randrange(1 << n)
        expected = sum(v.get(x,0)*v.get(y,0)*pow(norm,-1,5)
                       for v,norm in zip(vectors,norms)) % 5
        assert projector_entry(n,x,y,last) == expected
    vector = recursive_vector(paths(2)[0])
    raw = np.array([[vector.get(x,0)*vector.get(y,0) % 5 for y in range(4)]
                    for x in range(4)],dtype=np.int64)
    assert not np.array_equal(raw @ raw % 5,raw)
    return {'status':'pass','prime':5,'full_auxiliary_states':30,'selected_auxiliary_states':60,
            'full_projector_cases':cases,'binary_projector_cases':selected,
            'complete_word_pair_controls':int(entries),'additional_width_four_controls':128,
            'missing_norm_negative_control':'unnormalized two-site matrix is not idempotent',
            'scope':'Uniform finite-state projectors and exact finite identities; no favorable crossing cut-rank bound or factoring improvement.',
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ('tools/bounded_spin_crossing.py','tools/bounded_spin_projector.py',
                              'tools/check_bounded_spin_projector.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
