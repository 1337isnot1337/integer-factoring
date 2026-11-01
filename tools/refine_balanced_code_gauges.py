#!/usr/bin/env python3
"""Bind a sustained existing gauge optimizer to a frozen selected code seed."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--family', choices=('symmetric', 'asymmetric'), required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--restarts', type=int, default=32)
    parser.add_argument('--maxfun', type=int, default=20000)
    parser.add_argument('--height-bound', type=float, default=80.)
    parser.add_argument('--random-seed', type=int, default=2026103221)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    assert source['completed_requested_schedule']
    selected = min((row for row in source['candidates'] if row['family'] == args.family),
                   key=lambda row: row['best']['rho'])
    seed_path = Path(str(args.prefix)+'-seed.json')
    start_path = Path(str(args.prefix)+'-start.json')
    output_path = Path(str(args.prefix)+'-search.json')
    provenance_path = Path(str(args.prefix)+'-provenance.json')
    seed = selected['literal_seed']
    seed_path.write_text(json.dumps(seed, indent=2)+'\n')
    start = {'dimensions': seed['dimensions'], 'parameters': selected['best']['params'],
             'rho_floating': selected['best']['rho'],
             'source_record': str(args.source), 'source_record_sha256': sha256(args.source.read_bytes()).hexdigest()}
    start_path.write_text(json.dumps(start, indent=2)+'\n')
    command = [sys.executable, str(Path(__file__).with_name('search_rectangular_f4_pressure.py')),
               str(seed_path), '--starts', str(start_path), '--restarts', str(args.restarts),
               '--maxfun', str(args.maxfun), '--seconds', str(args.seconds),
               '--height-bound', str(args.height_bound), '--random-seed', str(args.random_seed),
               '--output', str(output_path)]
    names = (Path(__file__).name, 'search_rectangular_f4_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Sustained exploratory all-height paired-Gibbs optimization; source tensor and starting state are frozen literally.',
              'family': args.family, 'source_record': str(args.source),
              'source_record_sha256': sha256(args.source.read_bytes()).hexdigest(),
              'seed': str(seed_path), 'seed_sha256': sha256(seed_path.read_bytes()).hexdigest(),
              'start': str(start_path), 'start_sha256': sha256(start_path.read_bytes()).hexdigest(),
              'output': str(output_path), 'initial_rho': selected['best']['rho'],
              'command': command,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'status': 'running'}
    provenance_path.write_text(json.dumps(report, indent=2)+'\n')
    started = time.monotonic()
    completed = subprocess.run(command, check=False)
    report.update(status='completed' if completed.returncode == 0 else 'failed',
                  returncode=completed.returncode, elapsed_seconds=time.monotonic()-started)
    if output_path.exists():
        report['output_sha256'] = sha256(output_path.read_bytes()).hexdigest()
    provenance_path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'family': args.family, 'status': report['status'],
                      'provenance': str(provenance_path)}), flush=True)
    if completed.returncode:
        raise SystemExit(completed.returncode)


if __name__ == '__main__':
    main()
