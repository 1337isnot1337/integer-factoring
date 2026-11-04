#!/usr/bin/env python3
"""Replay parity-refined theta programs, certificates, and evaluator controls."""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path

from theta_coordinate_program import (
    compile_program as compile_original, literal_planes,
    theta_data as original_theta_data,
)
from theta_parity_program import (
    compile_program, normalized_factors, original_coefficient,
    require_odd_prime, sturm_bound, theta_data, validate_program,
)

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT/'evidence/theta-parity-programs.json.gz'
RECEIPT = ROOT/'evidence/theta-parity-programs.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-programs', action='store_true',
                        help='write the deterministic compressed program artifact')
    args = parser.parse_args()
    records = []
    programs = []
    evaluator_steps = controls = 0
    repaired = []
    for prime in (3, 5, 7, 11, 13, 17, 19):
        original = compile_original(prime)
        result = compile_program(prime)
        assert result['status'] == 'certified_by_ordinary_modular_Sturm_proof'
        assert not result['validation']['failures']
        assert result['Sturm_bound_inclusive'] == sturm_bound(prime)
        if original['missing_coordinates']:
            assert original['missing_coordinates'] == [0]
            repaired.append(prime)
        targets = literal_planes(prime, 511)
        local_steps = 0
        for integer in (1, 3, 7, 15, 31, 63, 127, 255, 511):
            value = original_coefficient(integer, result)
            expected = sum(((plane >> integer) & 1) << i
                           for i, plane in enumerate(targets))
            assert value['value'] == expected
            evaluator_steps += value['charged_steps']
            local_steps += value['charged_steps']
            controls += 1
        old_variables, old_trace = original_theta_data(prime, 511)
        new_variables, new_trace = theta_data(prime, 511)
        assert old_trace == new_trace
        for a in range((prime-1)//2):
            assert new_variables[2*a] ^ new_variables[2*a+1] == old_variables[a]
        summary = {k: v for k, v in result.items() if k != 'program'}
        summary['program_instructions'] = len(result['program'])
        summary['maximum_quadratic_variables'] = max(
            [len(normalized_factors(prime, x['monomial']))
             for x in result['program'] if x['monomial']] or [0])
        summary['original_module_missing_coordinates'] = original['missing_coordinates']
        summary['original_module_status'] = original['status']
        summary['original_evaluator_controls'] = 9
        summary['charged_original_evaluator_steps'] = local_steps
        programs.append(result)
        records.append(summary)
    assert repaired == [5, 13, 17]
    # A deliberately underdetermined fit is still refused.
    underfit = compile_program(17)
    underfit['program'] = underfit['program'][:1]
    underfit['status'] = 'certified_by_ordinary_modular_Sturm_proof'
    underfit_validation = validate_program(underfit, 1024)
    assert underfit_validation['failures']
    # The original compiler's partial output is still refused by the parity
    # evaluator, so the new decoder does not silently promote old programs.
    partial = compile_original(13)
    assert partial['status'] == 'bounded_module_excluded'
    try:
        original_coefficient(15, partial)
        raise AssertionError('partial original program accepted as a full engine')
    except ValueError:
        pass
    for invalid in (0, 1, 2, 4, 9, 21):
        try:
            require_odd_prime(invalid)
            raise AssertionError('invalid phase accepted')
        except ValueError:
            pass
    for invalid in (0, -1, 2):
        try:
            original_coefficient(invalid, programs[0])
            raise AssertionError('invalid coefficient index accepted')
        except ValueError:
            pass
    raw = json.dumps(programs, separators=(',', ':'), sort_keys=True).encode()
    if args.write_programs:
        ARTIFACT.write_bytes(gzip.compress(raw, mtime=0))
    assert gzip.decompress(ARTIFACT.read_bytes()) == raw, 'regenerated programs differ from artifact'
    receipt = {
        'scope': ('Parity-refined degree-four theta programs with conservative '
                  'ordinary modular/Sturm certificates; no uniform success theorem '
                  'and no faster factoring bound'),
        'phase_results': records,
        'certified_phases': [p['phase'] for p in programs],
        'repaired_phases': repaired,
        'original_module_excluded_phases': [r['phase'] for r in records
                                            if r['original_module_missing_coordinates']],
        'evaluator_controls': controls,
        'charged_evaluator_steps': evaluator_steps,
        'underfit_regression': {'phase': 17, 'validation': underfit_validation},
        'original_partial_refused': {'phase': 13, 'status': partial['status']},
        'program_artifact': {'path': str(ARTIFACT.relative_to(ROOT)),
                             'uncompressed_bytes': len(raw),
                             'compressed_bytes': ARTIFACT.stat().st_size,
                             'uncompressed_sha256': sha256(raw).hexdigest(),
                             'compressed_sha256': sha256(ARTIFACT.read_bytes()).hexdigest()},
        'source_sha256': {name: sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/theta_parity_program.py',
                                       'tools/check_theta_parity_program.py',
                                       'tools/theta_coordinate_program.py',
                                       'tools/tate_coordinate_compiler.py',
                                       'tools/cyclotomic_character_decoder.py')},
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
