#!/usr/bin/env python3
"""Run or independently replay a larger optimized oriented power prefix."""

import argparse
from hashlib import sha256
import json
import lzma
from pathlib import Path
import resource
from time import process_time

from check_fresh_prime_cofactor_control import generate_input, verify_input
from check_sparse_oriented_power_search import SOURCES as REFERENCE_SOURCES, canonical, read_data, summarize, verify_search
from optimized_sparse_oriented_power_search import factor_with_optimized_oriented_refinement

ROOT = Path(__file__).resolve().parents[1]
SOURCES = REFERENCE_SOURCES + ('optimized_sparse_oriented_power_search.py',
    'check_optimized_oriented_fresh_run.py', 'quotient_update_oriented_basis.py',
    'sparse_unit_presentation.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--import-result', type=Path)
    parser.add_argument('--full-legacy-result', type=Path)
    parser.add_argument('--bits', type=int, default=96)
    parser.add_argument('--cutoff', type=int, default=1021)
    parser.add_argument('--max-queries', type=int, default=65536)
    parser.add_argument('--checkpoint-interval', type=int, default=2048)
    parser.add_argument('--max-extraction-columns', type=int, default=256)
    parser.add_argument('--memory-limit-mib', type=int, default=1536)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--summary-output', type=Path)
    args = parser.parse_args()
    if args.memory_limit_mib < 512:
        parser.error('The address-space limit must be at least 512 MiB.')
    limit = args.memory_limit_mib * 2**20
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    hashes = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
              for name in SOURCES}
    if args.replay:
        data = read_data(args.replay)
        for name, digest in data['source_sha256'].items():
            assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    elif args.import_result:
        data = read_data(args.import_result)
        if 'unit_source_sha256_at_search' in data:
            assert data['unit_source_sha256_at_search'] == hashes['tools/sparse_unit_presentation.py']
        if 'search_source_sha256' in data:
            for name, digest in data['search_source_sha256'].items():
                assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    else:
        source = generate_input(args.bits)
        verify_input(source)
        started = process_time()
        result = factor_with_optimized_oriented_refinement(source['integer'], args.cutoff,
            args.max_queries, checkpoint_interval=args.checkpoint_interval,
            max_extraction_columns=args.max_extraction_columns,
            progress=lambda row: print(json.dumps({key: row.get(key) for key in
                ('queries', 'keys', 'atoms', 'core_columns', 'core_rows',
                 'refined_factor', 'elapsed_seconds')}), flush=True))
        cpu = process_time() - started
        data = canonical({'input': source, 'result': result, 'cpu_seconds': cpu,
                          'search_source_sha256': hashes})
    if args.full_legacy_result:
        data['full_legacy'] = read_data(args.full_legacy_result)
    verification = verify_search(data)
    for name, digest in hashes.items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    report = {**data, 'verification': verification, 'source_sha256': hashes,
              'scope': 'Certified finite optimized prefix with exact original-traversal comparison; no uniform exponent gain.'}
    if args.output:
        if args.output.suffix != '.xz':
            parser.error('The complete output must end in .xz.')
        temporary = args.output.with_name(args.output.name + '.partial')
        with lzma.open(temporary, 'wb') as stream:
            for piece in json.JSONEncoder(sort_keys=True, separators=(',', ':')).iterencode(report):
                stream.write(piece.encode())
            stream.write(b'\n')
        temporary.replace(args.output)
    summary = summarize(report)
    summary['cpu_seconds'] = data.get('cpu_seconds')
    summary['component_replacements'] = data['result'].get('component_replacements')
    if args.summary_output:
        args.summary_output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
