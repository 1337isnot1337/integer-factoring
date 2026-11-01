#!/usr/bin/env python3
"""Complete 27-case F4 screen with coupled height-eight and height-twelve constraints.

The operators overlap in physical modes and change the protected low-height
core. Every coefficient triple (a,b,c) in (F4*)^3 is included. Coordinate
retractions always exist; orthogonal retractions are additionally checked
when available. Numerical pressures are discovery diagnostics only.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds_f4 import profiles
from orthogonal_differential_pairs_f4 import orthogonal_retraction
from search_coherent_hot_retractions import pressure
from tensor_f4_five_pair_catalogue_profiles import orthogonal_zero_profile


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--modes', type=int, default=6)
    ap.add_argument('--height', type=int, default=36)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    names = [Path(__file__).name, 'coherent_differential_kernels_f4.py',
             'coherent_hot_compounds_f4.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'orthogonal_differential_pairs_f4.py', 'search_coherent_hot_retractions.py',
             'tensor_f4_five_pair_catalogue_profiles.py']
    report = {'scope': __doc__, 'modes': args.modes, 'height_cap': args.height,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'records': [], 'completed': False, 'in_progress': None}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()

    def save():
        report['elapsed_seconds'] = time.monotonic() - start
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    for coefficients in product((1, 2, 3), repeat=3):
        a, b, c = coefficients
        operators = [{9: 1, 6: a}, {33: 1, 18: b, 12: c}]
        code = differential_kernel(operators, args.modes, args.height)
        orthogonal = orthogonal_retraction(code)
        for kind, left in [('coordinate', code['retraction'])] + ([('orthogonal', orthogonal)] if orthogonal else []):
            report['in_progress'] = {'coefficients': coefficients, 'retraction_kind': kind, 'q': code['q']}
            save()
            if kind == 'coordinate':
                profile = profiles(code['embedding'], left, args.modes)
            else:
                profile = orthogonal_zero_profile(code, left)
            fit = pressure(profile)
            report['records'].append({'coefficients': coefficients, 'retraction_kind': kind,
                                      'code': code, 'retraction': left, 'profile': profile, 'pressure': fit})
            report['in_progress'] = None
            save()
            print(json.dumps({'completed': coefficients, 'kind': kind, 'q': code['q'],
                              'rho': fit['rho'], 'seconds': profile['seconds']}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
