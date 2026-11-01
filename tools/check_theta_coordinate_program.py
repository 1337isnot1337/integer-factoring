#!/usr/bin/env python3
"""Replay a partial theta compiler, finite modular certificates and exclusions."""

import argparse
import gzip
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from lambert_theta_sections import AdditivePhaseAlgebra
from tate_coordinate_compiler import compile_coordinate, embed_coefficient
from theta_coordinate_program import (
    compile_program, fit_program, literal_planes, native_phase_weights,
    normalized_factors, original_coefficient, require_odd_prime, validate_program,
)

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT/'evidence/theta-coordinate-programs.json.gz'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-programs', action='store_true',
                        help='write the deterministic compressed program artifact')
    args = parser.parse_args()
    records = []
    programs = []
    evaluator_steps = controls = 0
    for prime in (3, 5, 7, 11, 13, 17, 19, 23, 31):
        result = compile_program(prime)
        expected_success = prime % 4 == 3
        assert (result['status'] == 'certified_by_ordinary_modular_Sturm_proof') == expected_success
        if not expected_success:
            assert result['missing_coordinates'] == [0]
        # Every phase weight is independently replayed under the actual embedding.
        constant, weights = native_phase_weights(prime)
        plan = compile_coordinate(prime, intrinsic=True)
        phase = AdditivePhaseAlgebra(prime)
        assert embed_coefficient(constant, plan) == compile_coordinate(prime)['constant']
        for residue, value in enumerate(weights):
            assert embed_coefficient(value, plan) == phase.phase(residue) ^ phase.phase(-residue)
        summary = {k: v for k, v in result.items() if k != 'program'}
        summary['program_instructions'] = len(result['program'])
        summary['maximum_quadratic_variables'] = max(
            [len(normalized_factors(x['monomial'])) for x in result['program']
             if x['monomial']] or [0])
        if expected_success:
            targets = literal_planes(prime, 255)
            local_steps = 0
            for integer in (1, 3, 7, 15, 31, 63, 127, 255):
                value = original_coefficient(integer, result)
                expected = sum(((plane >> integer) & 1) << i for i, plane in enumerate(targets))
                assert value['value'] == expected
                evaluator_steps += value['charged_steps']
                local_steps += value['charged_steps']
                controls += 1
            summary['original_evaluator_controls'] = 8
            summary['charged_original_evaluator_steps'] = local_steps
        else:
            try:
                original_coefficient(15, result)
                raise AssertionError('partial-coordinate program accepted as a full engine')
            except ValueError:
                pass
        programs.append(result)
        records.append(summary)
    # An underdetermined training fit is deliberately forced to fail validation.
    underfit = fit_program(31, 4096)
    assert not underfit['missing_coordinates']
    underfit_validation = validate_program(underfit, 8192)
    assert underfit_validation['failures']
    try:
        original_coefficient(31, underfit)
        raise AssertionError('unvalidated fit accepted as a full engine')
    except ValueError:
        pass
    # Large ORIGINAL-index controls at phase three use only the proved unary trace.
    q3 = programs[0]
    large_controls = []
    for bits in (128, 256, 512):
        for exponent in (0, 1, 3):
            for multiplier in (1, 5):
                integer = 3**exponent * multiplier * ((1 << bits)+1)**2
                reduced = integer
                while reduced % 3 == 0:
                    reduced //= 3
                value = original_coefficient(integer, q3)
                assert value['value'] == int(isqrt(reduced)**2 == reduced)
                assert value['charged_steps'] <= 2
                large_controls.append({'bits': integer.bit_length(), 'value': value['value']})
    for invalid in (0, 1, 2, 4, 9, 21):
        try:
            require_odd_prime(invalid)
            raise AssertionError('invalid phase accepted')
        except ValueError:
            pass
    for invalid in (0, -1, 2):
        try:
            original_coefficient(invalid, q3)
            raise AssertionError('invalid coefficient index accepted')
        except ValueError:
            pass
    raw = json.dumps(programs, separators=(',', ':'), sort_keys=True).encode()
    if args.write_programs:
        ARTIFACT.write_bytes(gzip.compress(raw, mtime=0))
    assert gzip.decompress(ARTIFACT.read_bytes()) == raw, 'regenerated programs differ from artifact'
    receipt = {
        'scope': 'Partial compiler with ordinary modular/Sturm certificates at six listed phases; no uniform success theorem or faster factoring bound',
        'phase_results': records,
        'certified_phases': [p['phase'] for p in programs if not p['missing_coordinates']],
        'bounded_module_excluded_phases': [p['phase'] for p in programs if p['missing_coordinates']],
        'original_evaluator_controls': controls,
        'charged_original_evaluator_steps': evaluator_steps,
        'large_phase_three_original_controls': large_controls,
        'underfit_regression': {'phase': 31, 'fit_limit': 4096, 'binary_rank': underfit['binary_rank'],
                               'missing_training_coordinates': [], 'validation': underfit_validation,
                               'uncertified_evaluator_refused': True},
        'program_artifact': {'path': str(ARTIFACT.relative_to(ROOT)),
                             'uncompressed_bytes': len(raw), 'compressed_bytes': ARTIFACT.stat().st_size,
                             'uncompressed_sha256': sha256(raw).hexdigest(),
                             'compressed_sha256': sha256(ARTIFACT.read_bytes()).hexdigest()},
        'source_sha256': {name: sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/theta_coordinate_program.py',
                                       'tools/check_theta_coordinate_program.py',
                                       'research/theta-coordinate-programs.md',
                                       'tools/tate_coordinate_compiler.py',
                                       'tools/cyclotomic_character_decoder.py')},
    }
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
