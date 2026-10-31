#!/usr/bin/env python3
"""Replay exact arbitrary-mode and parameter-dependent coded quartic controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from build_coded_quartic_crossing import build


ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [((0, 1), 2, ((1, 0),) * 4, (1, 0)),
             ((0, 1, 2, 3), 4, ((1, 0),) * 4, (1, 0)),
             ((0, 1, 2, 4), 8, ((1, 0),) * 4, (1, 0)),
             ((0, 1, 3, 7), 8, ((1, 0),) * 4, (1, 0)),
             ((0, 1, 2, 3, 4), 8, ((1, 0),) * 4, (1, 0)),
             ((0, 1, 2, 3), 4, ((2, 0), (4, 0), (7, 0), (1, -2)), (11, 0)),
             ((0, 1, 2, 3), 4, ((4099, -6), (4101, -7), (8201, -6), (8209, -7)), (1, 0)),
             ((0, 1, 3, 7), 8, ((7, -1), (11, 2), (13, -3), (1, 4)), (21, -2))]
    records = []
    for indices, ambient, coefficients, quartic in cases:
        data = build(indices, ambient, coefficients, quartic)
        values = data.pop('values')
        data['coefficient_sha256'] = sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()
        if indices == (0, 1, 2, 3) and coefficients == ((1, 0),) * 4 and quartic == (1, 0):
            baseline = json.loads((ROOT / 'evidence/direct-residual-code6-n4.json').read_text())
            assert values == baseline['values']
            data['previous_independent_exterior_control_matches'] = True
        records.append(data)
    report = {'status': 'exact polynomial and typed-endpoint controls; no pressure claim',
              'field': 2, 'cases': records,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                                for name in ('tools/build_coded_quartic_crossing.py',
                                             'tools/check_coded_quartic_mode_families.py',
                                             'tools/check_schur_downset_crossing.py',
                                             'tools/weighted_smith_oracle.py')}}
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
