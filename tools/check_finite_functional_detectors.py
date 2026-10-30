#!/usr/bin/env python3
"""Exact finite-field verification of weighted functional rank detectors."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import time

from check_binary_two_piece_f4_boundary import Field, QUADS, cuts, rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    d = 5
    pairs = list(product(range(d), repeat=2))
    quads = list(product(range(d), repeat=4))
    cert2 = [(7, (1, 0, 0, 0)), (5, (1, 1, 0, 0)),
             (3, (1, 0, 1, 0)), (2, (1, 0, 0, 1)), (1, (1, 0, 1, 1))]
    cert3 = [(6, (1, 0, 0, 0)), (1, (1, 0, 0, 1)),
             (3, (1, 0, 1, 0)), (1, (1, 0, 2, 2)),
             (3, (1, 1, 0, 0)), (1, (1, 2, 0, 1))]
    field_results = []
    for q in (2, 3, 4):
        field = Field(q)
        cert = cert3 if q == 3 else cert2
        bases = []
        for x, y, u, v in quads:
            e, r, c = int(x == v and y == u), int(x == y and u == v), int(x == u and y == v)
            h = field.sub(field.sub(e, r), c) if max(x, y, u, v) < 3 else 0
            bases.append([e, r, c, h])
        assert rank(bases, field) == 4
        projective = [coefficients for coefficients in product(range(q), repeat=4)
                      if any(coefficients) and next(x for x in coefficients if x) == 1]
        assert len(projective) == (q**4-1)//(q-1)
        profiles = []
        for coefficients in projective:
            tensor = {quad: field.dot(coefficients, basis) for quad, basis in zip(quads, bases)}
            rr = rank([[tensor[x, y, u, v] for u, v in pairs] for x, y in pairs], field)
            cr = rank([[tensor[x, y, u, v] for y, v in pairs] for x, u in pairs], field)
            detector_rank = sum(weight for weight, ell in cert if field.dot(ell, coefficients))
            assert detector_rank**2 <= rr*cr
            profiles.append({'coefficients': coefficients, 'cut_ranks': [rr, cr],
                             'detector_rank': detector_rank})
        crossing_rank = sum(weight for weight, ell in cert if ell[0])
        assert crossing_rank == (15 if q == 3 else 18)
        field_results.append({'field_size': q, 'projective_cases': len(projective),
                              'basis_rank': 4, 'certificate': cert,
                              'verified_spectral_lower_bound': crossing_rank,
                              'profiles': profiles})

    scope_controls = []
    for q, coefficients in ((2, (1, 0)), (2, (0, 1)), (2, (1, 1)), (4, (1, 2))):
        field = Field(q)
        a, b = coefficients
        identity, companion = [[1, 0], [0, 1]], [[0, 1], [1, 1]]
        matrix = [[field.add(field.mul(a, identity[x][u]), field.mul(b, companion[x][u]))
                   for u in range(2)] for x in range(2)]
        tensor = {(x, y, u, v): matrix[x][u] if y == v == 0 else 0
                  for x, y, u, v in QUADS}
        rr, cr = cuts(tensor, field)
        weighted = 7*int(a != 0)
        valid = weighted**2 <= 25*rr*cr
        assert valid == (q == 2)
        scope_controls.append({'field_size': q, 'coefficients': coefficients,
                               'cut_ranks': [rr, cr], 'base_field_certificate_valid': valid})
    result = {'alphabet': d, 'subalphabet': 3,
              'basis': ['E', 'R', 'C', 'restricted(E-R-C)'],
              'field_certificates': field_results,
              'field_scope_negative_control': scope_controls,
              'scope': 'Each finite-field certificate is separately checked; no automatic extension-field claim.',
              'seconds': time.monotonic()-started}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({**{key: value for key, value in result.items() if key != 'field_certificates'},
                      'field_certificates': [{key: value for key, value in entry.items() if key != 'profiles'}
                                             for entry in field_results]}, indent=2))


if __name__ == '__main__':
    main()
