#!/usr/bin/env python3
"""Complete-profile F4 inverse-system discovery in one homogeneous top cell."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from coherent_hot_compounds import height
from coherent_hot_compounds_f4 import profiles
from coherent_inverse_systems_f4 import inverse_system
from search_coherent_hot_retractions import pressure


def projective_vectors(size):
    for pivot in range(size):
        for tail in product(range(4), repeat=size-pivot-1):
            yield (0,)*pivot+(1,)+tail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, default=8)
    parser.add_argument('--particles', type=int, default=4)
    parser.add_argument('--height', type=int, default=32)
    parser.add_argument('--max-trials', type=int, default=64)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--candidate-file', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.modes <= 16 or args.max_trials < 1 or args.seconds <= 0:
        parser.error('Use one through sixteen modes and positive trial/time budgets.')
    words = [word for word in range(1 << args.modes)
             if word.bit_count() == args.particles and height(word) == args.height]
    if not words:
        parser.error('The declared homogeneous cell is empty.')
    if args.candidate_file:
        candidates = json.loads(args.candidate_file.read_text())
        total = len(candidates)
        for values in candidates:
            if len(values) != len(words) or any(v not in (0, 1, 2, 3) for v in values) or not any(values):
                parser.error('Each candidate must be a nonzero F4 vector in the declared cell.')
    else:
        candidates = projective_vectors(len(words))
        total = (4**len(words)-1)//3
    names = (Path(__file__).name, 'coherent_inverse_systems_f4.py', 'coherent_hot_compounds_f4.py',
             'coherent_hot_compounds_f4.cpp', 'coherent_hot_compounds.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'search_coherent_hot_retractions.py')
    report = dict(modes=args.modes, particles=args.particles, height=args.height,
                  top_monomials=words, declared_candidates=total, records=[], unresolved=0,
                  scope='Complete exact two-cut profiles for the completed finite prefix; numerical pressure fits are discovery diagnostics, not rate certificates.',
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in names},
                  candidate_file_sha256=sha256(args.candidate_file.read_bytes()).hexdigest()
                  if args.candidate_file else None)
    started, exhausted, best = time.monotonic(), True, None
    for index, values in enumerate(candidates):
        if len(report['records']) >= args.max_trials or time.monotonic()-started >= args.seconds:
            exhausted = False
            break
        top = {word: value for word, value in zip(words, values) if value}
        code = inverse_system([top], args.modes)
        record = dict(index=index, coefficients=list(values), code=code, completed=False)
        try:
            result = profiles(code['embedding'], code['retraction'], args.modes)
            fit = pressure(result)
            record.update(completed=True, q=result['q'], zero_orders_t=result['zero_orders_t'],
                          infinity_orders_t=result['infinity_orders_t'], pressure=fit,
                          block_precisions=[[block['precision_attempts'] for block in cut['blocks']]
                                            for cut in result['cuts']],
                          endpoint_entries=sum(block['endpoint_entries'] for cut in result['cuts']
                                               for block in cut['blocks']),
                          construction_seconds=result['seconds'])
            if best is None or fit['rho'] < best:
                best = fit['rho']
                print(json.dumps(dict(index=index, q=result['q'], best_rho=best,
                                      coefficients=list(values), seconds=time.monotonic()-started)), flush=True)
        except (ValueError, ArithmeticError) as error:
            report['unresolved'] += 1
            record['error'] = f'{type(error).__name__}: {error}'
        report['records'].append(record)
        report['elapsed_seconds'] = time.monotonic()-started
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':'))+'\n')
        temporary.replace(args.output)
    report.update(completed=True, finite_list_exhausted=exhausted,
                  elapsed_seconds=time.monotonic()-started, best_rho=best)
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(dict(tested=len(report['records']), unresolved=report['unresolved'],
                          exhausted=exhausted, best_rho=best, seconds=report['elapsed_seconds'])))


if __name__ == '__main__':
    main()
