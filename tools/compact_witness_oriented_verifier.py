"""Independently reconstruct the omitted fields of compact power witnesses.

The literal-power and complete-presentation checks follow
check_sparse_oriented_power_search.py at SHA256
06c58bc8af09fe2d76af3b0d847aee0f95fdec3b8606fa9e644caa131b2cba99.
No full witness expansion is retained during replay.
"""

from math import gcd, isqrt, prod
from time import perf_counter

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from check_fresh_prime_cofactor_control import verify_input
from check_sparse_oriented_power_search import canonical, independent_rank
from cofactor_collision_relations import factor_with_cofactor_collisions, partial_leading_word
from compact_oriented_power_search import _bind
from compact_witness_oriented_power_search import transport_compact_witnesses
from incremental_oriented_coprime_basis import verify_incremental_oriented_basis
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from sparse_projected_signed_graph import verify_sparse_projected_signed_graph
from unit_presentation_trace import verify_unit_presentation_trace
from weighted_norm_relations import least_power_above

_unit_verifier = _bind(verify_unit_presentation_trace,
    {'verify_signed_graph_elimination': verify_sparse_projected_signed_graph})


def verify_compact_witness_search(data):
    input_check = verify_input(data['input'])
    record = data['result']
    assert record['witness_encoding'] == 'base-word-key-orientation-v1'
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
        assert len(witness) == 3
        stored_word, stored_index, stored_sign = witness
        assert type(stored_sign) is int
        offset, column = divmod(query, len(primes))
        prime = primes[column]
        exponent = least_power_above(prime, height) + offset
        target = power(records[prime][1], exponent, identity, compose_forms_hnf)
        assert gcd(target[0], integer) == 1
        if target != identity and target == inverse_form(target):
            assert ambiguous_factor(target, integer) is None
        split = partial_leading_word(integer, target, primes, records, prod(primes))
        trial_divisions += split['trial_divisions']
        assert stored_sign == split['orientation']
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
        assert stored_word == canonical(sorted(word.items()))
        assert stored_index is None or type(stored_index) is int
        assert stored_index == index
        expected_witnesses.append(witness)
    assert canonical(keys) == record['keys'] and trial_divisions == record['trial_divisions']
    verification = {'input': input_check, 'literal_power_checks': len(expected_witnesses)}
    refined = record['refinement']
    if refined is not None:
        prefix_queries = len(record['graph_reduction']['old_row_combinations'])
        prefix = expected_witnesses[:prefix_queries]
        used = [w[1] for w in prefix if w[1] is not None]
        prefix_keys = keys[:max(used, default=-1) + 1]
        assert canonical(prefix_keys) == refined['keys']
        verification['ideals'] = verify_incremental_oriented_basis(refined)
        rows = transport_compact_witnesses(prefix, len(primes), refined['sparse_words'])
        stored_rows = record['transported_rows']
        assert len(rows) == len(stored_rows)
        for rebuilt, stored in zip(rows, stored_rows):
            assert rebuilt == {int(column): value for column, value in stored.items()}
        verification['graph'] = verify_sparse_projected_signed_graph(rows, record['graph_reduction'])
        verification['units'] = _unit_verifier(rows, record['graph_reduction'],
                                                record['unit_reduction'], verify_graph=True)
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
