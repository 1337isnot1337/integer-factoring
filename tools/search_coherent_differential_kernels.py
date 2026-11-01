#!/usr/bin/env python3
"""Complete-profile screen of bounded-height Boolean differential kernels."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

from coherent_differential_kernels import differential_kernel
from coherent_hot_compounds_fast import profiles
from search_coherent_hot_retractions import _atomic_write, pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--caps', default='16,20,22')
    parser.add_argument('--families', default='low,high')
    parser.add_argument('--transpose', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    caps = list(map(int, args.caps.split(',')))
    families = args.families.split(',')
    low = (1 << 6) ^ (1 << 9)
    high = sum(1 << word for word in (24, 36, 66, 129))
    operators = {'low': [low], 'high': [high], 'joint': [low, high]}
    if any(name not in operators for name in families):
        raise ValueError('Families must be low, high, or joint.')
    names = (Path(__file__).name, 'coherent_differential_kernels.py', 'coherent_inverse_systems.py',
             'coherent_hot_compounds_fast.py', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds.py', 'search_coherent_hot_retractions.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': False,
              'parameters': {'caps': caps, 'families': families, 'transpose': args.transpose},
              'started_utc': datetime.now(timezone.utc).isoformat(), 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    for cap in caps:
        for name in families:
            modes = max(8, (cap + 1) // 2)
            code = differential_kernel(operators[name], modes, cap)
            embedding, retraction = code['embedding'], code['retraction']
            if args.transpose:
                embedding, retraction = retraction, embedding
            report['in_progress'] = {'family': name, 'height_cap': cap, 'code': code}
            _atomic_write(args.output, report)
            print(json.dumps({'starting': name, 'cap': cap, 'q': code['q']}), flush=True)
            profile = profiles(embedding, retraction, modes)
            fit = pressure(profile)
            record = {'family': name, 'height_cap': cap, 'code': code,
                      'profile': profile, 'pressure': fit}
            report['records'].append(record)
            report['in_progress'] = None
            report['elapsed_seconds'] = time.monotonic() - started
            _atomic_write(args.output, report)
            print(json.dumps({'family': name, 'cap': cap, 'q': code['q'], 'pressure': fit,
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, finite_list_exhausted=True,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
