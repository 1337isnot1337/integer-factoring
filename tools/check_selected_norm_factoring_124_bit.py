#!/usr/bin/env python3
"""Replay the 124-bit fixture through the integer-only selector interface."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_selected_norm_factoring import public_summary, verify_certificates
from selected_norm_factoring import factor_from_selected_norms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/selected-norm-factoring-124-bit.json')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source_path = root / 'evidence/partial-norm-factoring-certificates.json'
    data = json.loads(source_path.read_text())
    stored = data['records'][0]
    # The consumer receives no former representation or relation data.
    selected = [(witness['generator_prime'], witness['represented_norm'])
                for witness in stored['witnesses']]
    result = factor_from_selected_norms(
        stored['integer'], stored['smooth_bound'], stored['word_budget'],
        selected, extract_each=False, keep_certificates=True,
    )
    verify_certificates(result)
    assert result['factor'] is not None
    assert result['stats']['accepted_relations'] == len(selected)
    report = public_summary(result)
    report.update({
        'status': 'pass',
        'scope': 'Only unsigned selected integers supplied; source search cost remains separate.',
        'source_fixture': source_path.name,
        'witnesses_supplied': len(selected),
        'independent_relation_replay': True,
        'engine_sha256': sha256((root / 'tools/selected_norm_factoring.py').read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_fixture_sha256': sha256(source_path.read_bytes()).hexdigest(),
    })
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in (
        'status', 'input_bits', 'factor', 'factor_base_size', 'accepted_target_count',
    )}, indent=2))


if __name__ == '__main__':
    main()
