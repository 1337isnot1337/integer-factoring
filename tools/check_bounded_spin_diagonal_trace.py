#!/usr/bin/env python3
"""Trace-form control for the bounded-spin primary decomposition.

For a commutative coefficient algebra acting on the path-pair module, a
projection has rank equal to its trace.  In the normalized cut basis the
diagonal entry at the pair (a,b) is the original-pencil evaluation on the
four paths (a,b,b,a).  This checker records those trace forms at lengths
ten and twelve, and verifies the published ten-site primary dimensions
without rebuilding a cut operator.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import product
from math import comb, sqrt
from pathlib import Path

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry

ROOT = Path(__file__).resolve().parents[1]
PRIME = 5


def original_form(four_paths):
    """Normalized original-pencil coefficients at one path quadruple."""
    n = len(four_paths[0])
    result = [0]*(n+1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2*b+1):
            result[2*c+j] = (result[2*c+j]
                             + value*comb(2*b, j)*(-1)**j) % PRIME
    return tuple(result)


def trace_forms(length):
    ps = paths(length)
    counts = Counter(original_form((a, b, b, a)) for a in ps for b in ps)
    return counts


def check():
    model = json.loads((ROOT/'evidence/bounded-spin-ten-algebra-model.json').read_text())
    basis = model['basis_indices']
    ten = trace_forms(10)
    traces = []
    for e in model['idempotents']:
        value = sum(sum(count*form[basis[j]] for form,count in ten.items())*int(e[j])
                    for j in range(len(basis))) % PRIME
        traces.append(value)
    assert traces == [d % PRIME for d in model['dimensions']]
    assert sum(ten.values()) == 34*34

    twelve = trace_forms(12)
    identity_trace = sum(count*sum(form) for form,count in twelve.items()) % PRIME
    assert identity_trace == 89*89 % PRIME
    assert sum(twelve.values()) == 89*89
    exact_identity = all(sum(form) % PRIME == 1 for form in twelve)
    assert exact_identity

    def summary(counts):
        multiplicities = sorted(counts.values(), reverse=True)
        return {'distinct_forms': len(multiplicities),
                'multiplicities': multiplicities,
                'square_root_sum': sum(sqrt(m) for m in multiplicities),
                'class_size_sum': sum(multiplicities)}

    return {
        'status': 'pass',
        'prime': PRIME,
        'length_ten': {'primary_dimensions': model['dimensions'],
                       'trace_values': traces,
                       'forms': summary(ten)},
        'length_twelve': {'labels': 89,
                          'identity_trace': identity_trace,
                          'identity_diagonal_one': exact_identity,
                          'forms': summary(twelve)},
        'scope': ('Exact diagonal trace control for the ten-site primary '
                  'decomposition; the twelve-site forms reduce the open '
                  'primary-dimension computation to the coefficient '
                  'idempotents. No score change.'),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
