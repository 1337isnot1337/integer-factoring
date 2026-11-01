#!/usr/bin/env python3
"""Independently replay every witness of a complete prepared-selector result."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource

from check_sparse_oriented_power_search import read_data, summarize
from prepared_oriented_fresh_run import source_hashes, write_xz
from prepared_oriented_power_search import verify_prepared_search

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--import-result', type=Path)
    parser.add_argument('--full-legacy-result', type=Path)
    parser.add_argument('--memory-limit-mib', type=int, default=2048)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--summary-output', type=Path)
    args = parser.parse_args()
    if bool(args.replay) == bool(args.import_result):
        parser.error('Choose exactly one complete saved result to replay or import.')
    if args.memory_limit_mib < 512:
        parser.error('The address-space limit must be at least 512 MiB.')
    limit = args.memory_limit_mib * 2**20
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    hashes = source_hashes()
    data = read_data(args.replay or args.import_result)
    for name, digest in data['search_source_sha256'].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    if args.replay:
        for name, digest in data['source_sha256'].items():
            assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    policy = data['policy']
    policy_path = ROOT / data['policy_file']
    assert sha256(policy_path.read_bytes()).hexdigest() == data['policy_sha256']
    assert json.loads(policy_path.read_text()) == policy
    input_path = ROOT / policy['input_file']
    assert sha256(input_path.read_bytes()).hexdigest() == policy['input_sha256']
    assert json.loads(input_path.read_text()) == data['input']
    assert policy['source_sha256'] == data['search_source_sha256'] == hashes
    for name in ('cutoff', 'max_queries', 'checkpoint_interval', 'max_extraction_columns'):
        assert data['result'][name] == policy[name], name
    assert data['input']['integer'].bit_length() == policy['bits']
    if args.full_legacy_result:
        data['full_legacy'] = read_data(args.full_legacy_result)
    verification = verify_prepared_search(data)
    assert verification.get('units', {}).get('complete_original_graph_checked', True)
    assert source_hashes() == hashes
    report = {**data, 'verification': verification, 'source_sha256': hashes,
              'scope': 'Independently verified complete finite prepared selector; no uniform exponent improvement.'}
    if args.output:
        write_xz(args.output, report, preset=6)
    summary = {**summarize(report), 'cpu_seconds': data['cpu_seconds'],
               'component_replacements': data['result']['component_replacements']}
    if args.summary_output:
        args.summary_output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
