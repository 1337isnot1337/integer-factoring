#!/usr/bin/env python3
"""Continue interrupted physical-frame searches without changing their sources.

Saved prefixes are immutable. This harness restores candidate deduplication,
available elites, and the unused wall budget. It starts a declared new random
stream; no exact restoration of lost in-memory optimizer or RNG state is claimed.
Only initialization and queue bookkeeping are instrumented. Mathematical
constructors, exact controls, objectives, and optimizers remain unchanged.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import inspect
import json
from pathlib import Path
import sys


WORKERS = {
    'basis4': 'search_reflection_physical_bases.py',
    'asymmetric4': 'search_asymmetric_physical_bases.py',
    'basis16': 'search_reflection16_local_bases.py',
    'polynomial4': 'search_physical_polynomial_frames.py',
    'squarezero4': 'search_reflection4_square_zero_frames.py',
    'affine4': 'search_affine_coupled_reflection_frames.py',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=WORKERS)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--seed-offset', type=int, default=20000)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('The continuation output must be new.')
    raw = args.prefix.read_bytes()
    previous = json.loads(raw)
    if previous.get('completed'):
        parser.error('The supplied prefix already reports completion.')
    directory = Path(__file__).resolve().parent
    worker = directory / WORKERS[args.kind]
    for name, expected in previous['source_sha256'].items():
        if sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ArithmeticError('Changed frozen source: ' + name)
    spent = float(previous['elapsed_seconds'])
    remaining = float(previous['requested_seconds']) - spent
    if remaining <= 0:
        parser.error('There is no unused declared wall budget.')
    records = previous['cases']
    metadata = {
        'prefix_path': str(args.prefix.resolve()),
        'prefix_sha256': sha256(raw).hexdigest(),
        'prefix_was_completed': False,
        'previous_saved_cases': len(records),
        'previous_checkpoint_seconds': spent,
        'requested_total_seconds': previous['requested_seconds'],
        'requested_continuation_seconds': remaining,
        'original_worker': worker.name,
        'original_worker_sha256': sha256(worker.read_bytes()).hexdigest(),
        'harness_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'random_stream': 'New declared seed; lost RNG and optimizer states are not restored.',
        'saved_work_accounting': 'Only saved completed candidate records count; any unsaved interrupted candidate may be proposed again.',
    }
    spec = importlib.util.spec_from_file_location('_physical_continuation_worker', worker)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    module.RESUME_PREVIOUS = previous
    module.RESUME_COUNT = len(records)
    code = inspect.getsource(module.main)
    changes = []

    def replace(old, new):
        nonlocal code
        if code.count(old) != 1:
            raise ArithmeticError('Nonunique frozen initialization anchor: ' + old)
        code = code.replace(old, new)
        changes.append({'old': old, 'new': new})

    argv = [str(worker), '--output', str(args.output), '--seconds', str(remaining),
            '--seed', str(previous['seed'] + args.seed_offset),
            '--iterations', str(previous['iteration_limit'])]
    metadata['continuation_seed'] = previous['seed'] + args.seed_offset
    persistence = 'save'
    if args.kind in ('basis4', 'asymmetric4', 'basis16'):
        canonical = module.canonical_basis
        module.RESUME_SEEN = {
            (canonical(r['basis_x']), canonical(r['basis_y'])) +
            ((r['cut'],) if args.kind == 'asymmetric4' else ()) for r in records}
        score = 'support_value' if args.kind == 'asymmetric4' else 'rho'
        resolved = sorted((r for r in records if 'best' in r), key=lambda r: r['best'][score])
        module.RESUME_ELITES = resolved[:10 if args.kind == 'basis16' else 12]
        module.RESUME_BEST = resolved[0] if resolved else None
    if args.kind == 'basis4':
        argv += ['--interval', str(previous['interval']), '--starts', str(previous['starts_per_basis'])]
        module.RESUME_START = max(r['index'] for r in records) + 1
        replace('    seen, elites, best = set(), [], None',
                '    seen, elites, best = set(RESUME_SEEN), list(RESUME_ELITES), RESUME_BEST')
        replace('    started, index = time.monotonic(), 0',
                '    started, index = time.monotonic(), RESUME_START')
        persistence = 'write_atomic'
    elif args.kind == 'asymmetric4':
        if args.source is None:
            parser.error('The asymmetric continuation requires its current immutable basis source.')
        argv += ['--source', str(args.source)]
        source_raw = args.source.read_bytes()
        source_data = json.loads(source_raw)
        module.RESUME_START = max(r['index'] for r in source_data['cases']) + 1
        metadata['new_source_snapshot_sha256'] = sha256(source_raw).hexdigest()
        metadata['source_policy'] = 'Use the later saved basis prefix and skip all already-tested basis/cut pairs.'
        replace('    source_cursor = generation = 0', '    source_cursor, generation = 0, RESUME_START')
        replace('    seen, elites = set(), []', '    seen, elites = set(RESUME_SEEN), list(RESUME_ELITES)')
        replace('    started, best_case = time.monotonic(), None',
                '    started, best_case = time.monotonic(), RESUME_BEST')
    elif args.kind == 'basis16':
        argv += ['--trial-seconds', str(previous['trial_seconds'])]
        module.RESUME_START = max(r['generation'] for r in records) + 1
        replace('    seen, elites, best = set(), [], None',
                '    seen, elites, best = set(RESUME_SEEN), list(RESUME_ELITES), RESUME_BEST')
        replace('    index = 0', '    index = RESUME_START')
    elif args.kind == 'polynomial4':
        argv += ['--trial-seconds', str(previous['trial_seconds'])]
        balanced = sorted((r for r in records if 'balanced' in r), key=lambda r: r['balanced']['rho'])
        asymmetric = sorted((r for r in records if 'asymmetric' in r), key=lambda r: r['asymmetric']['support_value'])
        module.RESUME_BALANCED = balanced[0] if balanced else None
        module.RESUME_ASYMMETRIC = asymmetric[0] if asymmetric else None
        module.RESUME_ELITES = balanced[:12]
        module.RESUME_SEEN = {r['seed_sha256'] for r in records}
        module.RESUME_START = max(r['generation'] for r in records) + 1
        replace('    started, index = time.monotonic(), 0',
                '    started, index = time.monotonic(), RESUME_START')
        replace('    best_balanced = best_asymmetric = None',
                '    best_balanced, best_asymmetric = RESUME_BALANCED, RESUME_ASYMMETRIC')
        replace('    elites, seen = [], set()', '    elites, seen = list(RESUME_ELITES), set(RESUME_SEEN)')
    elif args.kind == 'squarezero4':
        argv += ['--trial-seconds', str(previous['trial_seconds'])]
        key = lambda g: tuple(x for row in g for x in row)
        done = {(key(r['generator']), r['pole']) for r in records}
        module.RESUME_DONE = done
        module.RESUME_SEEN = {
            key(r['generator']) for r in previous['admissibility']
            if all((key(r['generator']), p) in done for p in r['allowed_poles'])}
        pending = [r for r in previous['admissibility'] if key(r['generator']) not in module.RESUME_SEEN]
        original_generators = module.generators

        def generators(rng):
            for r in pending:
                yield r['family'] + '-interrupted-pole-queue', r['generator']
            yield from original_generators(rng)

        module.generators = generators
        balanced = sorted((r for r in records if r['phase'] == 'balanced' and 'best' in r),
                          key=lambda r: r['best']['value'])
        module.RESUME_WARM = (balanced[0]['best']['params'] if balanced and
                              balanced[0]['best']['value'] <= previous['best_rho'] + 1e-9 else None)
        replace('    best_rho, best_support = baseline, math.inf',
                "    best_rho, best_support = RESUME_PREVIOUS['best_rho'], (RESUME_PREVIOUS['best_support_value'] or math.inf)")
        replace('    seen = set()', '    seen = set(RESUME_SEEN)')
        replace('    warm = initial.copy()',
                '    warm = initial.copy() if RESUME_WARM is None else np.asarray(RESUME_WARM, float)')
        replace('        for pole, seed in candidates:',
                '        for pole, seed in candidates:\n            if (key, pole) in RESUME_DONE:\n                continue')
        replace("            phase = 'balanced' if len(report['cases']) % 3 != 2 else 'asymmetric'",
                "            phase = 'balanced' if (RESUME_COUNT + len(report['cases'])) % 3 != 2 else 'asymmetric'")
    else:
        if not previous.get('screen_phase_complete'):
            parser.error('This affine continuation requires the complete saved screen queue.')
        argv += ['--trial-seconds', str(previous['trial_seconds']), '--screen', str(previous['requested_screens']),
                 '--refine', str(previous['requested_refinements'])]
        module.RESUME_DONE = {r['seed_sha256'] for r in records}
        balanced = sorted((r for r in records if 'balanced' in r), key=lambda r: r['balanced']['value'])
        module.RESUME_WARM = (balanced[0]['balanced']['params'] if balanced and
                              balanced[0]['balanced']['value'] <= previous['best_rho'] + 1e-9 else None)
        begin = code.index('    seen = set()\n')
        end = code.index('    threshold = False\n', begin)
        replace(code[begin:end], """    started = time.monotonic()
    selected = sorted((r for r in RESUME_PREVIOUS['screen_records'] if r['natural_rho'] is not None),
                      key=lambda r: r['natural_rho'])[:args.refine]
    selected = [r for r in selected if r['seed_sha256'] not in RESUME_DONE]
    report['screen_phase_complete'] = True
    report['completed_screens'] = 0
    report['screen_records_reused_from_prefix'] = len(RESUME_PREVIOUS['screen_records'])
    best_rho = RESUME_PREVIOUS['best_rho']
    best_support = RESUME_PREVIOUS['best_support_value'] or math.inf
    warm = initial.copy() if RESUME_WARM is None else np.asarray(RESUME_WARM, float)
""")
    metadata['initialization_changes'] = changes
    metadata['instrumented_main_sha256'] = sha256(code.encode()).hexdigest()
    persist = getattr(module, persistence)

    def checkpoint(path, report):
        report['continuation'] = metadata
        report['combined_checkpoint_seconds'] = spent + float(report.get('elapsed_seconds', 0))
        report['combined_saved_cases'] = len(records) + len(report['cases'])
        persist(path, report)

    setattr(module, persistence, checkpoint)
    exec(compile(code, str(worker) + ':continuation-main', 'exec'), module.__dict__)
    print(json.dumps({'worker': worker.name, 'previous_cases': len(records), 'remaining_seconds': remaining,
                      'output': str(args.output), 'prefix_sha256': metadata['prefix_sha256']}), flush=True)
    sys.argv = argv
    module.main()


if __name__ == '__main__':
    main()
