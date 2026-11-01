#!/usr/bin/env python3
"""Replay every power, ideal, and integer presentation in a sparse search."""

import argparse
from hashlib import sha256
import json
import lzma
from math import gcd, isqrt, prod
from pathlib import Path
from time import perf_counter

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from check_fresh_prime_cofactor_control import generate_input, verify_input
from cofactor_collision_relations import RationalRowBasis, factor_with_cofactor_collisions, partial_leading_word
from incremental_oriented_coprime_basis import verify_incremental_oriented_basis
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from signed_residual_graph import verify_signed_graph_elimination
from sparse_oriented_power_search import factor_with_sparse_oriented_refinement, transport_sparse_witnesses
from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('sparse_oriented_power_search.py', 'check_sparse_oriented_power_search.py',
           'incremental_oriented_coprime_basis.py', 'oriented_coprime_basis.py',
           'signed_residual_graph.py', 'signed_residual_unit_elimination.py',
           'cofactor_collision_relations.py', 'leading_power_smooth_relations.py',
           'relation_two_saturation.py', 'quadratic_ideal_hnf.py',
           'check_fresh_prime_cofactor_control.py', 'certified_small_factor.py')


def canonical(value):
    return json.loads(json.dumps(value))


def independent_rank(rows, width):
    basis = RationalRowBasis(width)
    for row in rows:
        basis.add(row)
    return basis


def verify_search(data):
    input_check = verify_input(data['input'])
    record = data['result']
    integer = data['input']['integer']
    assert record['integer'] == integer and record['input_bits'] == integer.bit_length()
    delta = -integer if integer % 4 == 3 else -4 * integer
    assert record['discriminant'] == delta
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    base = build_prime_form_factor_base(delta, record['cutoff'])
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes, keys, indices, expected_witnesses = sorted(records), [], {}, []
    assert record['active_primes'] == primes and canonical(records) == record['prime_records']
    assert all(not 1 < gcd(p, integer) < integer for p in base['excluded_primes'])
    trial_divisions = 0
    direct = record['terminal_witness'] is not None
    assert len(record['witnesses']) + int(direct) == record['norm_queries']
    for query, witness in enumerate(record['witnesses']):
        offset, column = divmod(query, len(primes))
        prime = primes[column]
        exponent = least_power_above(prime, height) + offset
        assert witness['prime'] == prime and witness['exponent'] == exponent
        target = power(records[prime][1], exponent, identity, compose_forms_hnf)
        assert tuple(witness['target']) == target and gcd(target[0], integer) == 1
        if target != identity and target == inverse_form(target):
            assert ambiguous_factor(target, integer) is None
        split = partial_leading_word(integer, target, primes, records, prod(primes))
        trial_divisions += split['trial_divisions']
        assert canonical(split) == witness['split']
        key = split['key']
        if key is not None and key not in indices:
            indices[key] = len(keys)
            keys.append(key)
        word = {column: exponent}
        for p, value in split['word'].items():
            j = primes.index(p)
            word[j] = word.get(j, 0) - value
            if not word[j]:
                del word[j]
        index = indices[key] if key is not None else None
        assert witness['base_word'] == canonical(sorted(word.items()))
        assert witness['key_index'] == index
        expected_witnesses.append({'base_word': sorted(word.items()), 'key_index': index, 'split': split})
    assert canonical(keys) == record['keys'] and trial_divisions == record['trial_divisions']
    verification = {'input': input_check, 'literal_power_checks': len(expected_witnesses)}
    refined = record['refinement']
    if refined is not None:
        prefix_queries = len(record['graph_reduction']['old_row_combinations'])
        prefix = expected_witnesses[:prefix_queries]
        used = [w['key_index'] for w in prefix if w['key_index'] is not None]
        prefix_keys = keys[:max(used, default=-1) + 1]
        assert canonical(prefix_keys) == refined['keys']
        verification['ideals'] = verify_incremental_oriented_basis(refined)
        rows = transport_sparse_witnesses(prefix, len(primes), refined['sparse_words'])
        assert canonical(rows) == record['transported_rows']
        verification['graph'] = verify_signed_graph_elimination(rows, record['graph_reduction'])
        verification['units'] = verify_signed_graph_elimination(rows, record['unit_reduction'])
        compressed = record['unit_reduction']
        width = len(compressed['torsion_generator_indices'])
        if record['independent_rows'] is not None:
            independent = independent_rank(compressed['torsion_rows'], width)
            assert independent.rows == record['independent_rows']
            generators = [records[p][1] for p in primes] + [tuple(a['form']) for a in refined['atoms']]
            selected = [generators[i] for i in compressed['torsion_generator_indices']]
            top, factor, extraction = first_useful_relation_top(
                selected, independent.rows, identity, compose_forms_hnf, inverse_form,
                lambda form: ambiguous_factor(form, integer), verify_relations=True)
            assert factor == record['checkpoints'][-1]['refined_factor']
            if not direct:
                assert factor == record['factor']
            verification['independent_rank'] = len(independent.rows)
            verification['extraction'] = extraction
            verification['ambiguous_form'] = top
    if direct:
        query = record['norm_queries'] - 1
        offset, column = divmod(query, len(primes))
        prime = primes[column]
        exponent = least_power_above(prime, height) + offset
        target = power(records[prime][1], exponent, identity, compose_forms_hnf)
        terminal = record['terminal_witness']
        assert terminal['prime'] == prime and terminal['exponent'] == exponent
        assert tuple(terminal['target']) == target
        leading = gcd(target[0], integer)
        factor = leading if 1 < leading < integer else ambiguous_factor(target, integer)
        assert factor == record['factor']
    assert record['factor'] is None or record['factor'] in data['input']['input_primes']
    started = perf_counter()
    old = factor_with_cofactor_collisions(integer, record['cutoff'], record['norm_queries'],
        max_queries=record['norm_queries'], verify_candidates=True)
    comparison = {key: old[key] for key in ('factor', 'norm_queries', 'relation_rank',
        'retained_cofactor_keys', 'candidate_relations', 'dependent_relations')}
    comparison['elapsed_seconds'] = perf_counter() - started
    comparison['cap'] = record['norm_queries']
    if record['checkpoints'] and old['norm_queries'] == record['norm_queries'] and not direct:
        assert old['relation_rank'] == record['checkpoints'][-1]['old_base_rank']
        assert old['factor'] == record['checkpoints'][-1]['old_factor']
    verification['legacy_comparison'] = comparison
    if 'full_legacy' in data:
        stored = data['full_legacy']
        full = stored['result']
        assert full['integer'] == integer and full['cutoff'] == record['cutoff']
        tick = perf_counter()
        replayed = factor_with_cofactor_collisions(integer, full['cutoff'], full['window'],
            max_queries=full['max_queries'], verify_candidates=True)
        assert canonical(replayed) == full
        assert full['factor'] in data['input']['input_primes']
        verification['full_legacy_comparison'] = {
            **{key: full[key] for key in ('factor', 'norm_queries', 'relation_rank',
                 'retained_cofactor_keys', 'candidate_relations', 'dependent_relations')},
            'recorded_elapsed_seconds': stored['elapsed_seconds'],
            'recorded_cpu_seconds': stored['cpu_seconds'],
            'independent_replay_elapsed_seconds': perf_counter() - tick,
            'complete_original_traversal_reproduced': True}
    return verification


def read_data(path):
    return json.loads(lzma.decompress(path.read_bytes()) if path.suffix == '.xz' else path.read_bytes())


def summarize(report):
    record, check = report['result'], report['verification']
    final = record['checkpoints'][-1] if record['checkpoints'] else {}
    return {'integer': record['integer'], 'bits': record['input_bits'],
            'cutoff': record['cutoff'], 'max_queries': record['max_queries'],
            'checkpoint_interval': record['checkpoint_interval'],
            'factor': record['factor'], 'queries': record['norm_queries'],
            'elapsed_seconds': record['elapsed_seconds'], 'phase_seconds': record['phase_seconds'],
            'keys': len(record['keys']), 'old_base_rank': final.get('old_base_rank'),
            'core_columns': final.get('core_columns'), 'core_rows': final.get('core_rows'),
            'independent_rank': final.get('independent_rank'),
            'graph_counts': final.get('graph_counts'),
            'basis_statistics': record['refinement']['statistics'] if record['refinement'] else None,
            'verification': check, 'source_sha256': report['source_sha256'],
            'scope': report['scope']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--import-result', type=Path)
    parser.add_argument('--full-legacy-result', type=Path)
    parser.add_argument('--bits', type=int, default=64)
    parser.add_argument('--cutoff', type=int, default=127)
    parser.add_argument('--max-queries', type=int, default=8192)
    parser.add_argument('--checkpoint-interval', type=int, default=256)
    parser.add_argument('--max-extraction-columns', type=int, default=256)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--summary-output', type=Path)
    args = parser.parse_args()
    if args.replay:
        data = read_data(args.replay)
        for path, digest in data['source_sha256'].items():
            assert sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    elif args.import_result:
        data = read_data(args.import_result)
    else:
        source = generate_input(args.bits)
        verify_input(source)
        result = factor_with_sparse_oriented_refinement(source['integer'], args.cutoff,
            args.max_queries, checkpoint_interval=args.checkpoint_interval,
            max_extraction_columns=args.max_extraction_columns)
        data = canonical({'input': source, 'result': result})
    if args.full_legacy_result:
        data['full_legacy'] = read_data(args.full_legacy_result)
    verification = verify_search(data)
    report = {**data, 'verification': verification,
              'source_sha256': {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in SOURCES},
              'scope': 'Finite matched-prefix factor gain on a certified input; no uniform exponent claim.'}
    if args.output:
        if args.output.suffix != '.xz':
            parser.error('The complete certificate output must end in .xz.')
        temporary = args.output.with_name(args.output.name + '.partial')
        with lzma.open(temporary, 'wb') as stream:
            for piece in json.JSONEncoder(sort_keys=True, separators=(',', ':')).iterencode(report):
                stream.write(piece.encode())
            stream.write(b'\n')
        temporary.replace(args.output)
    summary = summarize(report)
    if args.summary_output:
        args.summary_output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
