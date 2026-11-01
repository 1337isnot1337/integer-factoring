#!/usr/bin/env python3
"""Compact a finished search without losing trial outcomes or retained winners."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=Path)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    source = json.loads(raw)
    baseline_bytes = args.baseline.read_bytes()
    assert sha256(baseline_bytes).hexdigest() == source['baseline_sha256']
    assert source['budget_completed'] or source['trial_limit_completed']
    messages, records, winners = [], [], {}
    for item in source['records']:
        outcomes = []
        for stage in item['outcomes']:
            message = stage['optimizer_message']
            if message not in messages:
                messages.append(message)
            outcomes.append([stage['smoothing'], stage['optimizer_success'], stage['iterations'],
                             stage['function_evaluations'], messages.index(message)])
        best = item['best']
        records.append([item['trial'], item['alphabet'], item['mode'], item['include_inactive'],
                        item['train_allocations'], item['completed_trial'], item['queries'],
                        item['elapsed_seconds'], outcomes,
                        best['rho'] if best else None, best['entropy_sum'] if best else None,
                        digest(best['parameters']) if best else None, item['interruption']])
        key = str(item['alphabet'])
        if best and (key not in winners or best['rho'] < winners[key]['best']['rho']):
            winners[key] = item
    baselines = []
    for entry in json.loads(baseline_bytes)['seeds']:
        baselines.append({'seed': {key: entry['seed'][key] for key in ('m', 'theta', 'alphabet', 'degree')},
                          'partitions': sorted(entry['completed_partitions'], key=lambda row: row['rho'])[:3]})
    result = {key: source[key] for key in ('scope', 'arguments', 'baseline_sha256', 'source_sha256',
                                         'budget_completed', 'trial_limit_completed', 'elapsed_seconds', 'completed_trials')}
    result.update({'schema_version': 1,
                   'compact_scope': 'Every trial outcome and best value is retained; all winning parameters and local records are retained. Redundant nonwinning local profiles remain identified by the full-run hash.',
                   'full_run_sha256': sha256(raw).hexdigest(), 'full_run_bytes': len(raw),
                   'trial_fields': ['trial', 'alphabet', 'mode', 'include_inactive', 'train_allocations',
                                    'completed_trial', 'queries', 'elapsed_seconds', 'optimizer_stages',
                                    'best_rho', 'best_entropy_sum', 'best_parameters_sha256', 'interruption'],
                   'optimizer_stage_fields': ['smoothing', 'success', 'iterations', 'function_evaluations', 'message_index'],
                   'optimizer_messages': messages, 'trials': records,
                   'attempted_trials': len(records), 'winners': winners,
                   'baseline_partitions': baselines,
                   'compactor_sha256': sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps({'attempted_trials': len(records), 'completed_trials': source['completed_trials'],
                      'raw_bytes': len(raw), 'compact_bytes': args.output.stat().st_size,
                      'winners': {q: row['best']['rho'] for q, row in winners.items()}}))


if __name__ == '__main__':
    main()
