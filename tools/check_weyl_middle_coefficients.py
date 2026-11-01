#!/usr/bin/env python3
"""Complete F4 middle-coefficient tests for two ordered Weyl boundaries."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from weyl_ordered_crossing import (MULTIPLICATION, ROOT_POWERS, add, canonical_weights,
                                    coefficient_cuts, field_rank, scale)

ROOT = Path(__file__).resolve().parents[1]
POSITIONS = (0, 1, 2, 4, 5, 6)
SOURCES = ('check_weyl_middle_coefficients.py', 'weyl_ordered_crossing.py',
           'check_steiner_cubic_crossings.py')


def directions(boundary):
    cuts = coefficient_cuts(2, canonical_weights(boundary))
    current = [side[3] for side in cuts]
    changes = [[scale(side[k], 1 << b) for side in cuts]
               for k in POSITIONS for b in range(2)]
    previous = 0
    for index in range(4096):
        code = index ^ (index >> 1)
        if index:
            changed = (code ^ previous).bit_length()-1
            current = [add(a, b) for a, b in zip(current, changes[changed])]
        previous = code
        yield code, current


def coefficient_code(code):
    return (1 << 6) | sum(((code >> (2*i)) & 3) << (2*k)
                         for i, k in enumerate(POSITIONS))


def generate():
    records = {}
    for boundary in (0, 1694):
        pairs = [None]*4096
        for code, matrices in directions(boundary):
            pairs[code] = [field_rank(a) for a in matrices]
        records[str(boundary)] = pairs
    return dict(records=records, middle_position=3,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in SOURCES},
                scope='Complete F4 middle-detecting functionals at length two, for boundary codes 0 and 1694 only.',
                factoring_exponent='1/2; unchanged')


def verify(data):
    assert data['middle_position'] == 3
    assert set(data['records']) == {'0', '1694'}
    # Independent pure-boundary formula on all nine Weyl basis matrices.
    psi, phi = (1, 2, 2), (1, 3, 3)
    for a in range(3):
        for b in range(3):
            value = 0
            for j in range(3):
                value ^= MULTIPLICATION[phi[(j+a) % 3]][
                    MULTIPLICATION[ROOT_POWERS[b*j % 3]][psi[j]]]
            assert value == canonical_weights(0)[3*a+b]
    minima = {}
    witnesses = {}
    for boundary, expected in ((0, 576), (1694, 1950)):
        saved = data['records'][str(boundary)]
        assert len(saved) == 4096
        minimum = 81**2
        witness = None
        for code, matrices in directions(boundary):
            pair = [field_rank(a, low=True) for a in matrices]
            assert pair == saved[code]
            cost = pair[0]*pair[1]
            if cost < minimum:
                minimum, witness = cost, coefficient_code(code)
        assert minimum == expected
        minima[str(boundary)] = minimum
        witnesses[str(boundary)] = witness
    # rho > 5/7 > 1/sqrt(2), without floating-point logarithms.
    assert 576**7 > 9**20 and 2*5**2 > 7**2
    for name, digest in data['source_sha256'].items():
        assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest
    return dict(complete_normalized_middle_functionals=8192,
                exact_cut_ranks=16384, minimum_rank_products=minima,
                attaining_coefficient_codes=witnesses,
                strict_finite_rho_floor='5/7 > 1/sqrt(2)',
                scope=data['scope'], factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/weyl-middle-coefficients.json')
    args = parser.parse_args()
    if args.generate:
        data = generate()
        args.evidence.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    else:
        data = json.loads(args.evidence.read_text())
    print(json.dumps(verify(data)))


if __name__ == '__main__':
    main()
