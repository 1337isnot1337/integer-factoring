"""Exact power prefixes with incremental atoms and integer graph elimination."""

from math import gcd, isqrt, prod
from time import perf_counter

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from cofactor_collision_relations import RationalRowBasis, partial_leading_word, verify_partial_identity
from incremental_oriented_coprime_basis import IncrementalOrientedBasis
from leading_power_smooth_relations import compose_with_prime_form
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from signed_residual_graph import eliminate_signed_residual_graph
from signed_residual_unit_elimination import eliminate_unit_residual_columns
from weighted_norm_relations import least_power_above


def transport_sparse_witnesses(witnesses, base_size, sparse_words):
    """Build exact sparse rows from stored base words and residual key indices."""
    rows = []
    for witness in witnesses:
        row = dict(witness['base_word'])
        index = witness['key_index']
        if index is not None:
            sign = witness['split']['orientation']
            for column, exponent in sparse_words[index]:
                row[base_size + column] = -sign * exponent
        rows.append(row)
    return rows


def factor_with_sparse_oriented_refinement(integer, cutoff, max_queries, *,
                                           checkpoint_interval=256,
                                           max_extraction_columns=256,
                                           verify_candidates=True, progress=None):
    """Retain every queried witness and charge every checkpoint from the start.

The search knows only the integer and policy. It does not test residual
primality or factor residual norms. The column cap skips only extraction
at an oversized checkpoint; it never drops supplied rows or columns.
    """
    if (not isinstance(integer, int) or integer < 3 or integer % 2 == 0
            or not isinstance(cutoff, int) or cutoff < 2
            or not isinstance(max_queries, int) or max_queries < 1
            or not isinstance(checkpoint_interval, int) or checkpoint_interval < 1
            or not isinstance(max_extraction_columns, int) or max_extraction_columns < 1):
        raise ValueError('Require an odd integer and positive bounded search parameters.')
    started = perf_counter()
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    base = build_prime_form_factor_base(delta, cutoff)
    if any(1 < gcd(p, integer) < integer for p in base['excluded_primes']):
        raise ValueError('The input is already split by the initial factor base.')
    records = {record.prime: (record.root, record.positive) for record in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes = sorted(records)
    if not primes:
        raise ValueError('Require a nonempty initial prime-form base.')
    base_size, prime_product = len(primes), prod(primes)
    columns = {p: i for i, p in enumerate(primes)}
    first = {p: least_power_above(p, height) for p in primes}
    targets = {p: power(records[p][1], first[p], identity, compose_forms_hnf) for p in primes}
    state = IncrementalOrientedBasis(delta)
    key_indices, old_residual_words, witnesses, checkpoints = {}, {}, [], []
    old_basis = RationalRowBasis(base_size)
    trial_divisions = queries = candidate_collisions = 0
    phase_times = {'incremental_basis_seconds': 0.0, 'transport_seconds': 0.0,
                   'graph_seconds': 0.0, 'unit_seconds': 0.0, 'extraction_seconds': 0.0}
    final_refinement = final_rows = final_graph = final_units = final_independent = None
    inspect = lambda form: ambiguous_factor(form, integer)

    def finish(factor, method, terminal=None):
        return {'integer': integer, 'input_bits': integer.bit_length(), 'discriminant': delta,
                'cutoff': cutoff, 'max_queries': max_queries,
                'checkpoint_interval': checkpoint_interval,
                'max_extraction_columns': max_extraction_columns,
                'active_primes': primes, 'prime_records': records,
                'norm_queries': queries, 'trial_divisions': trial_divisions,
                'keys': state.keys, 'witnesses': witnesses, 'checkpoints': checkpoints,
                'old_base_rows': old_basis.rows, 'candidate_collisions': candidate_collisions,
                'refinement': final_refinement, 'transported_rows': final_rows,
                'graph_reduction': final_graph, 'unit_reduction': final_units,
                'independent_rows': final_independent,
                'factor': factor, 'method': method, 'terminal_witness': terminal,
                'elapsed_seconds': perf_counter() - started, 'phase_seconds': phase_times,
                'factoring_oracle_calls': 0, 'residual_primality_tests': 0,
                'scope': 'Finite deterministic prefix and exact supplied-presentation reduction.'}

    for query in range(max_queries):
        offset, column = divmod(query, base_size)
        prime = primes[column]
        target, exponent = targets[prime], first[prime] + offset
        queries += 1
        divisor = inspect(target) if target != identity and target == inverse_form(target) else None
        leading_gcd = gcd(target[0], integer)
        if 1 < leading_gcd < integer:
            divisor = leading_gcd
        if divisor is not None:
            return finish(divisor, 'direct target factor',
                          {'prime': prime, 'exponent': exponent, 'target': target})
        split = partial_leading_word(integer, target, primes, records, prime_product)
        trial_divisions += split['trial_divisions']
        if verify_candidates:
            assert verify_partial_identity(integer, target, split, records)
        key = split['key']
        if key is not None and key not in key_indices:
            tick = perf_counter()
            key_indices[key] = state.add(*key)
            phase_times['incremental_basis_seconds'] += perf_counter() - tick
        word = {column: exponent}
        for p, value in split['word'].items():
            j = columns[p]
            value = word.get(j, 0) - value
            if value:
                word[j] = value
            else:
                word.pop(j, None)
        witnesses.append({'prime': prime, 'exponent': exponent, 'target': target,
                          'split': split, 'base_word': sorted(word.items()),
                          'key_index': key_indices[key] if key is not None else None})
        candidate = None
        if key is None:
            candidate = word
        else:
            signed = {j: split['orientation'] * value for j, value in word.items()}
            if key not in old_residual_words:
                old_residual_words[key] = signed
            else:
                candidate_collisions += 1
                candidate = dict(signed)
                for j, value in old_residual_words[key].items():
                    candidate[j] = candidate.get(j, 0) - value
        if candidate is not None:
            old_basis.add([candidate.get(j, 0) for j in range(base_size)])
        if query + 1 < max_queries:
            targets[prime] = compose_with_prime_form(target, prime, *records[prime])
        if (query + 1) % checkpoint_interval and query + 1 != max_queries:
            continue
        tick = perf_counter()
        final_refinement = state.result(verify=False)
        final_rows = transport_sparse_witnesses(witnesses, base_size,
                                                final_refinement['sparse_words'])
        phase_times['transport_seconds'] += perf_counter() - tick
        tick = perf_counter()
        final_graph = eliminate_signed_residual_graph(final_rows, base_size,
            generator_count=base_size + len(final_refinement['atoms']),
            allow_extra_rows=True, verify=True)
        phase_times['graph_seconds'] += perf_counter() - tick
        tick = perf_counter()
        final_units = eliminate_unit_residual_columns(final_rows, final_graph, verify=True)
        phase_times['unit_seconds'] += perf_counter() - tick
        count = len(final_units['torsion_generator_indices'])
        checkpoint = {'queries': queries, 'keys': len(state.keys),
                      'atoms': len(final_refinement['atoms']),
                      'old_base_rank': len(old_basis.rows),
                      'core_columns': count, 'core_rows': len(final_units['torsion_rows']),
                      'graph_counts': final_units['counts'],
                      'basis_gcd_decisions': final_refinement['decision_gcd_calls'],
                      'extraction_skipped': count > max_extraction_columns}
        tick = perf_counter()
        _, old_factor, old_stats = first_useful_relation_top(
            [records[p][1] for p in primes], old_basis.rows, identity,
            compose_forms_hnf, inverse_form, inspect, verify_relations=True)
        checkpoint.update(old_factor=old_factor, old_stats=old_stats)
        factor = top = None
        final_independent = None
        if count <= max_extraction_columns:
            generators = [records[p][1] for p in primes] + [a['form'] for a in final_refinement['atoms']]
            selected = [generators[i] for i in final_units['torsion_generator_indices']]
            independent = RationalRowBasis(count)
            for row in final_units['torsion_rows']:
                independent.add(row)
            final_independent = independent.rows
            top, factor, stats = first_useful_relation_top(
                selected, final_independent, identity, compose_forms_hnf,
                inverse_form, inspect, verify_relations=True)
            checkpoint.update(refined_factor=factor, refined_stats=stats,
                              independent_rank=len(final_independent), ambiguous_form=top)
        phase_times['extraction_seconds'] += perf_counter() - tick
        checkpoint['elapsed_seconds'] = perf_counter() - started
        checkpoints.append(checkpoint)
        if progress is not None:
            progress(checkpoint)
        if factor is not None:
            return finish(factor, 'oriented residual presentation')
    return finish(None, 'query cap exhausted')
