#!/usr/bin/env python3
"""Run one declared bounded prepared selector and save its complete raw result."""

import argparse
from hashlib import sha256
import json
import lzma
from pathlib import Path
import resource
from time import process_time

from check_fresh_prime_cofactor_control import verify_input
from check_prepared_oriented_power_search import SOURCES as MATCHED_SOURCES, binding_controls
from prepared_oriented_power_search import factor_with_prepared_oriented_refinement

ROOT = Path(__file__).resolve().parents[1]
SOURCES = MATCHED_SOURCES + ('prepared_oriented_fresh_run.py',
                            'check_prepared_oriented_fresh_run.py')


def source_hashes():
    return {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCES}


def write_xz(path, data, *, preset=3):
    if path.suffix != '.xz':
        raise ValueError('The complete result path must end in .xz.')
    temporary = path.with_name(path.name + '.partial')
    with lzma.open(temporary, 'wt', encoding='utf-8', preset=preset) as stream:
        json.dump(data, stream, sort_keys=True, separators=(',', ':'))
        stream.write('\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--policy', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    policy_bytes = args.policy.read_bytes()
    policy_hash = sha256(policy_bytes).hexdigest()
    policy = json.loads(policy_bytes)
    hashes = source_hashes()
    assert hashes == policy['source_sha256']
    binding_controls()
    limit = policy['memory_limit_mib'] * 2**20
    assert limit >= 512 * 2**20
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    input_path = ROOT / policy['input_file']
    assert sha256(input_path.read_bytes()).hexdigest() == policy['input_sha256']
    source = json.loads(input_path.read_text())
    input_check = verify_input(source)
    assert source['integer'].bit_length() == policy['bits']
    print(json.dumps({'event': 'started', 'bits': policy['bits'],
        'integer': source['integer'], 'input_verification': input_check,
        'policy_sha256': policy_hash}), flush=True)

    def progress(row):
        print(json.dumps({key: row.get(key) for key in ('queries', 'keys', 'atoms',
            'old_base_rank', 'core_columns', 'core_rows', 'independent_rank',
            'extraction_skipped', 'refined_factor', 'elapsed_seconds')}), flush=True)

    started = process_time()
    result = factor_with_prepared_oriented_refinement(source['integer'], policy['cutoff'],
        policy['max_queries'], checkpoint_interval=policy['checkpoint_interval'],
        max_extraction_columns=policy['max_extraction_columns'], progress=progress)
    cpu = process_time() - started
    assert source_hashes() == hashes
    assert sha256(args.policy.read_bytes()).hexdigest() == policy_hash
    report = {'input': source, 'result': result, 'cpu_seconds': cpu,
              'policy': policy, 'policy_sha256': policy_hash,
              'policy_file': str(args.policy.resolve().relative_to(ROOT)),
              'search_source_sha256': hashes,
              'scope': 'Complete bounded finite search; independent saved-artifact replay is required.'}
    write_xz(args.output, report)
    print(json.dumps({'event': 'saved', 'factor': result['factor'],
        'queries': result['norm_queries'], 'method': result['method'],
        'elapsed_seconds': result['elapsed_seconds'], 'cpu_seconds': cpu,
        'complete_result_sha256': sha256(args.output.read_bytes()).hexdigest(),
        'maximum_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        'saved_artifact_replay_complete': False}), flush=True)


if __name__ == '__main__':
    main()
