#!/usr/bin/env python3
"""Resume only the unspent wall budget of one frozen exploratory worker.

The original source and atomic prefix remain unchanged.  Only the main
loop's proposal index and carried fused winners are restored.  A declared
new random stream is used; this is not an exact RNG-state continuation.
Every new checkpoint records both source hashes and the prior elapsed time.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import inspect
import json
from pathlib import Path
import sys


ALLOWED = {'search_fused_physical_bases.py', 'search_reflection16_coupled_frames.py',
           'search_reflection16_square_zero_frames.py'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('worker', type=Path)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed-offset', type=int, default=10000)
    args = parser.parse_args()
    worker = args.worker.resolve()
    if worker.name not in ALLOWED:
        parser.error('This harness supports only its three named frozen workers.')
    if args.output.exists():
        parser.error('Use a new continuation output; never overwrite a saved prefix.')
    raw = args.prefix.read_bytes()
    previous = json.loads(raw)
    for name, expected in previous['source_sha256'].items():
        if sha256(worker.with_name(name).read_bytes()).hexdigest() != expected:
            raise ArithmeticError('A source differs from the interrupted prefix: ' + name)
    arguments = dict(previous['arguments'])
    spent = float(previous['elapsed_seconds'])
    remaining = float(arguments['seconds']) - spent
    if remaining <= 0:
        parser.error('The atomic prefix already spent its declared time budget.')
    records = previous['records']
    start = max((record['trial'] for record in records), default=-1) + 1
    metadata = {'prefix_path': str(args.prefix.resolve()), 'prefix_sha256': sha256(raw).hexdigest(),
                'original_worker': worker.name, 'original_worker_sha256': sha256(worker.read_bytes()).hexdigest(),
                'continuation_harness_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
                'previous_checkpoint_seconds': spent, 'previous_saved_attempts': len(records),
                'first_continuation_trial': start, 'requested_total_seconds': arguments['seconds'],
                'requested_continuation_seconds': remaining,
                'random_stream': 'New declared seed; no claim to restore the lost in-memory RNG state.'}
    arguments['seconds'], arguments['output'] = remaining, str(args.output)
    arguments['seed'] += args.seed_offset
    metadata['continuation_seed'] = arguments['seed']
    spec = importlib.util.spec_from_file_location('_frozen_continuation_worker', worker)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    module.RESUME_START = start
    module.RESUME_PREVIOUS = previous
    code = inspect.getsource(module.main)
    if worker.name == 'search_fused_physical_bases.py':
        changes = {'    winners, records = {}, []': "    winners, records = dict(RESUME_PREVIOUS['winners']), []",
                   '    for trial in range(args.trials):': '    for trial in range(RESUME_START, RESUME_START + args.trials):'}
        persistence_name = 'atomic_write'
    elif worker.name == 'search_reflection16_coupled_frames.py':
        changes = {'    trial = 0': '    trial = RESUME_START',
                   'trial < args.trials:': 'trial < RESUME_START + args.trials:'}
        persistence_name = 'save'
    else:
        changes = {'    trial, candidate = 0, None': "    trial, candidate = RESUME_START, dict(RESUME_PREVIOUS['records'][-1]['candidate'])",
                   'while trial < args.trials and': 'while trial < RESUME_START + args.trials and'}
        persistence_name = 'save'
    for old, new in changes.items():
        if code.count(old) != 1:
            raise ArithmeticError('The frozen loop no longer has the checked resume anchor: ' + old)
        code = code.replace(old, new)
    # Every saved continuation has its own attempt count, even though record
    # trial numbers continue the original proposal numbering.
    persist = getattr(module, persistence_name)

    def checkpoint(path, report):
        report['continuation'] = metadata
        report['attempted_trials'] = len(report['records'])
        report['combined_checkpoint_seconds'] = spent + float(report.get('elapsed_seconds', 0))
        persist(path, report)

    setattr(module, persistence_name, checkpoint)
    exec(compile(code, str(worker) + ':continuation-main', 'exec'), module.__dict__)
    positional = 'baseline' if worker.name == 'search_fused_physical_bases.py' else 'catalogue' if worker.name == 'search_reflection16_coupled_frames.py' else None
    argv = [str(worker)]
    if positional:
        argv.append(str(arguments.pop(positional)))
    for key, value in arguments.items():
        if value is None or value is False:
            continue
        argv.append('--' + key.replace('_', '-'))
        if value is not True:
            argv.extend(map(str, value if isinstance(value, list) else [value]))
    sys.argv = argv
    module.main()


if __name__ == '__main__':
    main()
