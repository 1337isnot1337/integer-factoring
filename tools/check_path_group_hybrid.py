#!/usr/bin/env python3
"""Exact finite budget/separation controls for one fixed path/group dictionary."""
import argparse
from hashlib import sha256
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from check_fixed_field_descent import Field
from check_quantum_braid_path_profiles import cut_ranks
from quantum_braid_nilpotent_seam import FIELD, coefficient_tensor

ROOT = Path(__file__).resolve().parents[1]


def dense_rank(matrix):
    """Independent dense elimination, rather than the packed rank backend."""
    a = matrix.copy()
    rank = 0
    for j in range(a.shape[1]):
        found = np.flatnonzero(a[rank:, j])
        if not len(found):
            continue
        p = rank + int(found[0])
        a[[rank, p]] = a[[p, rank]]
        a[rank] = FIELD[(0, 1, 3, 2)[a[rank, j]], a[rank]]
        factors = a[rank+1:, j].copy()
        a[rank+1:] ^= FIELD[factors[:, None], a[rank]]
        rank += 1
        if rank == len(a):
            break
    return rank


def pairing(tensor, witness):
    flat = tensor.reshape(-1)
    value = 0
    for index, coefficient in witness:
        assert 0 <= index < len(flat) and 1 <= coefficient <= 3
        value ^= int(FIELD[coefficient, flat[index]])
    return value


def check(receipt):
    field = Field(0b111)
    assert all(FIELD[a, b] == field.mul(a, b) for a in range(4) for b in range(4))
    # sqrt(2)<10/7 and sqrt(2)<99/70; both are strict integer tests.
    assert 2*7**2 < 10**2 and 5**10 < 10**7
    assert 2*70**2 < 99**2 and 13**99 < 38**70
    counts = []
    for saved in receipt['sizes']:
        n, d = saved['length'], saved['labels']
        assert (n, d) in ((6, 5), (8, 13))
        tensor = coefficient_tensor(n)
        assert tensor.shape[2:] == (d,)*4
        assert sha256(tensor.tobytes()).hexdigest() == saved['source_tensor_sha256']
        p, q, r, s = np.indices((d,)*4)
        identity = ((p == r) & (q == s)).astype(np.uint8)
        crossing = ((p == s) & (q == r)).astype(np.uint8)
        ones = np.ones((d,)*4, dtype=np.uint8)
        assert np.array_equal(np.bitwise_xor.reduce(tensor.reshape(-1,d,d,d,d), axis=0), crossing)
        assert cut_ranks(identity) == [d*d, 1]
        assert cut_ranks(ones) == [1, 1]
        candidates, profiles = [], []
        for a in range(tensor.shape[0]):
            row = []
            for b in range(tensor.shape[1]):
                t = tensor[a, b]
                ranks = cut_ranks(t)
                independent = [dense_rank(t.reshape(d*d, d*d)),
                               dense_rank(t.transpose(0,2,1,3).reshape(d*d, d*d))]
                assert ranks == independent
                row.append(ranks)
                if ranks[0]:
                    candidates.extend(((ranks[0]*ranks[1], t),
                                       (ranks[0]*ranks[1], t.transpose(0,2,1,3))))
            profiles.append(row)
        assert profiles == saved['complete_cut_ranks']
        lines = [((p+k*q-s-k*r) % d == 0).astype(np.uint8) for k in range(d)]
        lines.append((q == r).astype(np.uint8))
        assert all(cut_ranks(t) == [d,d] for t in lines)
        assert np.array_equal(np.bitwise_xor.reduce(lines), crossing ^ ones)
        expected = [(25, None)] if d == 5 else [(625, None)]+[(144,k) for k in range(14)]
        assert [(c['path_product_strict_cutoff'], c['line']) for c in saved['cases']] == expected
        unreserved = ([(100,())]+[(25,(i,)) for i in range(6)] if d == 5 else
                      [(1444,())]+[(625,(i,)) for i in range(14)]+
                      [(144,pair) for pair in combinations(range(14),2)])
        assert [(c['path_product_strict_cutoff'], tuple(c['lines']))
                for c in saved['unreserved_cases']] == unreserved
        plans = [(c, [] if c['line'] is None else [c['line']]) for c in saved['cases']]
        plans += [(c,c['lines']) for c in saved['unreserved_cases']]
        for case, selected_lines in plans:
            cutoff = case['path_product_strict_cutoff']
            entries = [identity, ones]+[t for cost,t in candidates if cost < cutoff]
            for line in selected_lines:
                entries.append(lines[line])
            witness = case['separator']
            assert len(entries) == case['candidate_count']
            assert all(pairing(t, witness) == 0 for t in entries)
            assert pairing(crossing, witness) == case['separator_target_value'] != 0
            assert dense_rank(np.array([t.reshape(-1) for t in entries])) == case['span_dimension']
        counts.append({'length':n, 'labels':d, 'full_coefficient_rank_pairs':len(profiles)*len(profiles[0]),
                       'separators':len(expected)+len(unreserved), 'line_directions':len(lines)})
    assert [x['length'] for x in counts] == [6,8]
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    receipt = json.loads((ROOT/'evidence/path-group-hybrid.json').read_text())
    counts = check(receipt)
    sources = ['tools/check_path_group_hybrid.py', 'tools/quantum_braid_nilpotent_seam.py',
               'tools/check_quantum_braid_path_profiles.py', 'tools/quantum_braid_rank_one.py',
               'tools/quantum_braid_path.py', 'tools/bounded_spin_crossing.py',
               'tools/search_auxiliary_extension_span.py', 'tools/check_fixed_field_descent.py',
               'research/path-group-hybrid.md']
    hashes = {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}
    if args.output:
        receipt['source_sha256'] = hashes
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert receipt['source_sha256'] == hashes
    print(json.dumps({'status':'pass','sizes':counts,'scope':'Fixed individual-piece dictionary only; no score gain.'}))


if __name__ == '__main__':
    main()
