#!/usr/bin/env python3
"""Complete F3 quadratic seeds with two-cut rank product at most two."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

from search_multipoint_rank_one_seeds import (
    KEYS, profiles, representatives, seed_from_values, uniform_prefix_threshold,
)


def matrix_rank(value):
    if not any(value):
        return 0
    return 1 if (value[0]*value[3]-value[1]*value[2]) % 3 == 0 else 2


def low_rank_tensors():
    vectors = list(product(range(3), repeat=4))
    projective = [v for v in vectors if any(v) and next(x for x in v if x) == 1]
    values = {(0,)*16}
    for a,b in product(projective,vectors[1:]):
        if matrix_rank(a)*matrix_rank(b) > 2:
            continue
        tensor = tuple(a[2*x+y]*b[2*u+v] % 3 for x,y,u,v in KEYS)
        values.add(tensor)
        values.add(tuple(tensor[KEYS.index((x,u,y,v))] for x,y,u,v in KEYS))
    assert len(values) == 3585
    return sorted(values)


def low_rank_representatives():
    answer = representatives()
    for beta in ((1,0,0,1),(2,0,0,1),(0,1,1,0)):
        answer.append(tuple(int(x == y == 0)*beta[2*u+v] for x,y,u,v in KEYS))
    assert len(answer) == 9
    return answer


def run():
    signatures,witnesses = Counter(),{}
    for first in low_rank_representatives():
        for second in low_rank_tensors():
            signature = profiles(seed_from_values(first,second))
            signatures[signature] += 1
            witnesses.setdefault(signature,(first,second))
    thresholds = {pair: uniform_prefix_threshold(pair)
                  for signature in signatures for pair in signature}
    failures = sum(sum(thresholds[pair][0] for pair in signature) > 2
                   for signature in signatures)
    return {'scope': 'F3 binary-alphabet quadratic crossings with two-cut rank product at most two at both 0 and 2; representatives allow a global cut swap.',
            'representatives': 9, 'second_values': 3585,
            'profiled_seeds': sum(signatures.values()),
            'signatures': len(signatures), 'local_profile_pairs': len(thresholds),
            'prefix_certificate_failures': failures,
            'uniform_prefix_rho_floor': '3/4' if failures == 0 else None,
            'catalogue': [{'profiles': signature,'count': count,
                           'values': witnesses[signature],
                           'prefix_certificates': [(str(thresholds[pair][0]),
                                                   *thresholds[pair][1:]) for pair in signature]}
                          for signature,count in sorted(signatures.items())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'catalogue'},indent=2))


if __name__ == '__main__':
    main()
