#!/usr/bin/env python3
"""Finite source-recursion checks of the high-mode creation obstruction."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_modular_crossing import determinant_operator_two, modular_operator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = []
    for h in (8, 16, 32, 64):
        for n in range(3, h//2, 2):
            m = h-n
            for cutoff in (0, 7, 31):
                state = tuple(range(1, cutoff+1, 2))
                high_bit = 1 << max(h, cutoff+1).bit_length()
                p, q = n+high_bit, m+2*high_bit
                output = tuple(sorted(state+(p, q)))
                assert p > cutoff and q > cutoff
                source = modular_operator(2, (n, m), output, state)
                determinant = determinant_operator_two((n, m), output, state)
                assert source == determinant == 1
                records.append({'operator_parts': [n, m], 'frozen_state': state,
                                'created_modes': [p, q], 'weight_shift': p+q,
                                'coefficient': source})
    paths = [Path(__file__), Path(__file__).with_name('check_modular_crossing.py')]
    result = {'status': 'source operators leak above every fixed finite core; no general rank bound',
              'cases': records, 'source_recursion_comparisons': len(records),
              'sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'source_recursion_comparisons': len(records)}))


if __name__ == '__main__':
    main()
