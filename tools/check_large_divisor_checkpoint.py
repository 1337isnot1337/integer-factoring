#!/usr/bin/env python3
"""Larger exact divisor checkpoints and the intentional large-integer CLI fix."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import difflib
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

try:
    import certified_divisor_oracle as oracle
except ModuleNotFoundError:
    staging = Path(__file__).resolve().parents[2]
    sys.path[:0] = [str(staging / name / 'tools') for name in
                   ('algorithm-certified', 'algorithm-adaptive', 'algorithm-polynomial',
                    'integer-factoring')]
    import certified_divisor_oracle as oracle


GUARD = ('    # This CLI intentionally accepts and emits arbitrary-size exact integers.\n'
         '    if hasattr(sys, "set_int_max_str_digits"):\n'
         '        sys.set_int_max_str_digits(0)\n')


def cli_sources():
    source = Path(oracle.__file__).read_text()
    unguarded = source.replace(GUARD, '').replace('import sys\n', '')
    assert unguarded.count('def main():\n') == 1
    guarded = unguarded.replace('import json\n', 'import json\nimport sys\n', 1)
    guarded = guarded.replace('def main():\n', 'def main():\n' + GUARD, 1)
    return unguarded, guarded


def cli_control(result, record_path, unguarded, guarded):
    if not hasattr(sys, 'set_int_max_str_digits'):
        return {'applicable': False, 'reason': 'interpreter predates the decimal conversion guard'}
    bound_text = str(result.exclusive_bound)
    assert len(bound_text) > 4300
    environment = dict(os.environ)
    environment.pop('PYTHONINTMAXSTRDIGITS', None)
    environment['PYTHONPATH'] = os.pathsep.join(str(Path(path)) for path in sys.path if path)
    with tempfile.TemporaryDirectory(prefix='divisor-cli-control-') as temporary:
        root = Path(temporary)
        broken = root / 'unguarded_cli.py'
        fixed = root / 'guarded_cli.py'
        broken.write_text(unguarded)
        fixed.write_text(guarded)
        # Exercise the actual CLI parsing and serialization, while reusing the
        # independently verified solver result. Repeating a large search is not
        # necessary to test this decimal-conversion boundary.
        wrapper = root / 'cli_wrapper.py'
        wrapper.write_text('''import json, runpy, sys
source, checkpoint = sys.argv[1:3]
arguments = sys.argv[3:]
namespace = runpy.run_path(source, run_name="cli_fixture")
limit = sys.get_int_max_str_digits()
sys.set_int_max_str_digits(0)
with open(checkpoint) as stream:
    record = json.load(stream)
sys.set_int_max_str_digits(limit)
result = namespace["CertifiedMaximum"](**record)
def checked_result(bound):
    assert bound == result.exclusive_bound
    return result
namespace["main"].__globals__["maximum_below"] = checked_result
sys.argv = [source] + arguments
namespace["main"]()
''')
        before = subprocess.run([sys.executable, str(wrapper), str(broken), str(record_path.resolve()),
                                 '--below', bound_text],
                                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                text=True, env=environment)
        assert before.returncode != 0
        assert 'invalid int value' in before.stderr or 'Exceeds the limit' in before.stderr
        cli_output = root / 'cli-result.json'
        after = subprocess.run([sys.executable, str(wrapper), str(fixed), str(record_path.resolve()),
                                '--below', bound_text,
                                '--output', str(cli_output)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                               text=True, env=environment)
        assert after.returncode == 0, after.stderr
        expected = json.loads(json.dumps(asdict(result)))
        assert json.loads(cli_output.read_text()) == expected
        before_bits = subprocess.run([sys.executable, str(wrapper), str(broken), str(record_path.resolve()),
                                      '--bits', str(result.input_bits)],
                                     stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                     text=True, env=environment)
        assert before_bits.returncode != 0 and 'Exceeds the limit' in before_bits.stderr
        after_bits = subprocess.run([sys.executable, str(wrapper), str(fixed), str(record_path.resolve()),
                                     '--bits', str(result.input_bits), '--output', str(cli_output)],
                                    stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                    text=True, env=environment)
        assert after_bits.returncode == 0, after_bits.stderr
        assert json.loads(cli_output.read_text()) == expected
    return {'applicable': True, 'input_decimal_digits': len(bound_text),
            'unguarded_exit_code': before.returncode, 'guarded_exit_code': after.returncode,
            'unguarded_bits_exit_code': before_bits.returncode,
            'guarded_bits_exit_code': after_bits.returncode,
            'full_result_matches_verified_checkpoint': True,
            'solver_reused_verified_checkpoint': True}


def main():
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, nargs='+', default=[16384, 32768])
    parser.add_argument('--source-dir', type=Path, default=Path('evidence'),
                        help='completed result and independent-audit files; no search is repeated')
    parser.add_argument('--output-dir', type=Path, default=Path('evidence'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    unguarded, guarded = cli_sources()
    patch_dir = args.output_dir.parent / 'patches'
    patch_dir.mkdir(parents=True, exist_ok=True)
    patch = ''.join(difflib.unified_diff(unguarded.splitlines(True), guarded.splitlines(True),
                    fromfile='a/tools/certified_divisor_oracle.py',
                    tofile='b/tools/certified_divisor_oracle.py'))
    (patch_dir / 'certified-divisor-cli.patch').write_text(patch)
    cases = []
    for bits in args.bits:
        assert bits >= 16384
        source_result = args.source_dir / f'result-{bits}.json'
        if not source_result.exists():
            source_result = args.source_dir / f'divisor-maximum-{bits}.json'
        source_audit = args.source_dir / f'audit-{bits}.json'
        if not source_audit.exists():
            source_audit = args.source_dir / f'divisor-audit-{bits}.json'
        record_text = source_result.read_text()
        audit_text = source_audit.read_text()
        record = json.loads(record_text)
        audit = json.loads(audit_text)
        result = oracle.CertifiedMaximum(**record)
        assert result.input_bits == audit['bits'] == bits
        assert result.exclusive_bound == 1 << bits
        assert result.least_witness < result.exclusive_bound
        assert result.divisor_parameter_bits == (result.maximum_divisors - 1).bit_length()
        assert result.stages[-1]['certified'] is True
        assert audit['independent_profiles'] > 0
        record_path = args.output_dir / f'divisor-maximum-{bits}.json'
        record_path.write_text(record_text)
        audit_path = args.output_dir / f'divisor-audit-{bits}.json'
        audit_path.write_text(audit_text)
        cli = cli_control(result, record_path, unguarded, guarded)
        cases.append({'input_bits': bits, 'maximum_divisors': result.maximum_divisors,
                      'divisor_parameter_bits': result.divisor_parameter_bits,
                      'support_size': len(result.exponents),
                      'least_witness_decimal_digits': len(str(result.least_witness)),
                      'visited_feasible_profiles': result.visited_profiles,
                      'total_search_nodes': result.total_search_nodes,
                      'independent_terminal_profiles': audit['independent_profiles'],
                      'certificate_comparisons': audit['marginals'],
                      'independent_replay': 'Imported completed exact prime-order replay; not rerun by this CLI regression.',
                      'cli_control': cli,
                      'result_file': record_path.name,
                      'result_sha256': hashlib.sha256(record_text.encode()).hexdigest(),
                      'audit_file': audit_path.name,
                      'audit_sha256': hashlib.sha256(audit_text.encode()).hexdigest()})
    report = {'scope': 'Consolidated completed exact larger-bound checkpoints with a decimal-guard CLI regression; no large independent search is repeated, and timings are not proof evidence.',
              'cases': cases,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'guarded_oracle_source_sha256': hashlib.sha256(guarded.encode()).hexdigest(),
              'patch_sha256': hashlib.sha256(patch.encode()).hexdigest()}
    report_text = json.dumps(report, indent=2) + '\n'
    (args.output_dir / 'large-divisor-checkpoint.json').write_text(report_text)
    print(report_text, end='')


if __name__ == '__main__':
    main()
