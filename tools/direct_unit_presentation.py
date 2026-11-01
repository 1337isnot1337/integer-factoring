"""Exact unit traces directly on a supplied sparse integer presentation.

The producer reuses the existing compact unit algorithm. Its independent
verifier indexes nonzero residual columns, then rebuilds every affected row
from the defining integer substitution. No graph generator maps are formed.
"""

from hashlib import sha256
import json

from compact_oriented_power_search import _bind
from unit_presentation_trace import (
    _finish as _reference_finish,
    eliminate_unit_residual_columns_compact as _reference_producer,
)


def _normalize_rows(rows, initial_size, generator_count):
    if (type(initial_size) is not int or type(generator_count) is not int
            or not 0 <= initial_size <= generator_count):
        raise ValueError('Require valid initial and total generator counts.')
    normalized = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError('Supply sparse integer-column dictionaries.')
        integer_keys = all(type(j) is int for j in row)
        clean = row if integer_keys and all(row.values()) else {}
        seen = set()
        for column, value in row.items():
            if not (type(column) is int or isinstance(column, str)):
                raise ValueError('A column index must be an integer or its decimal JSON key.')
            try:
                j = int(column)
            except ValueError as error:
                raise ValueError('A column index must be an integer.') from error
            if (j in seen or not 0 <= j < generator_count or type(value) is not int):
                raise ValueError('A sparse row has a repeated, out-of-range or noninteger entry.')
            seen.add(j)
            if clean is not row and value:
                clean[j] = value
        normalized.append(clean)
    return normalized


def _stream_fingerprint(reduction):
    """The reference fingerprint, serialized one row at a time."""
    digest = sha256()

    def emit(value):
        digest.update(json.dumps(value, separators=(',', ':')).encode())

    digest.update(b'[')
    emit(reduction['initial_size'])
    digest.update(b',')
    emit(reduction['new_generator_indices'])
    digest.update(b',[')
    for index, row in enumerate(reduction['output_rows']):
        if index:
            digest.update(b',')
        emit(sorted((int(j), value) for j, value in row.items() if value))
    digest.update(b']]')
    return digest.hexdigest()


def _identity_presentation(rows, initial_size, generator_count):
    return {'initial_size': initial_size,
            'new_generator_indices': list(range(generator_count)),
            'output_rows': rows,
            'counts': {'input_rows': len(rows), 'input_generator_count': generator_count}}


_finish = _bind(_reference_finish, {'_fingerprint': _stream_fingerprint})
_producer = _bind(_reference_producer, {'_finish': _finish})


def eliminate_direct_unit_presentation(rows, initial_size, generator_count, *, verify=True):
    """Eliminate residual unit columns while preserving the entire presentation."""
    rows = _normalize_rows(rows, initial_size, generator_count)
    identity = _identity_presentation(rows, initial_size, generator_count)
    result = _producer(rows, identity, verify=False, verify_graph=False)
    result['input_kind'] = 'identity-original-sparse-presentation-v1'
    if verify:
        result['verification'] = verify_direct_unit_presentation(
            rows, initial_size, generator_count, result)
    return result


def verify_direct_unit_presentation(original_rows, initial_size, generator_count, result):
    """Replay all local integer isomorphisms with a separately rebuilt index.

For each touched row this verifier discards all its previous incidence
memberships, computes the entire substituted row from the union of supports,
and inserts the new memberships. It does not reuse the producer's incremental
coefficient updates, heap, version numbers or pivot-selection decisions.
    """
    source = _normalize_rows(original_rows, initial_size, generator_count)
    identity = _identity_presentation(source, initial_size, generator_count)
    assert result['input_kind'] == 'identity-original-sparse-presentation-v1'
    assert result['schema'] == 'integer-unit-substitution-trace-v1'
    assert result['input_presentation_sha256'] == _stream_fingerprint(identity)
    assert result['initial_size'] == initial_size
    assert result['input_row_count'] == len(source)
    assert result['input_generator_count'] == generator_count
    rows = [dict(row) for row in source]
    active = set(range(len(rows)))
    removed = set()
    incidence = {}
    for index, row in enumerate(rows):
        for column in row:
            if column >= initial_size:
                incidence.setdefault(column, set()).add(index)
    substitutions = support_visits = 0
    maximum_row_support = max(map(len, rows), default=0)

    def remove_memberships(index, row):
        for column in row:
            if column >= initial_size:
                incident = incidence[column]
                assert index in incident
                incident.remove(index)
                if not incident:
                    del incidence[column]

    for step, entry in enumerate(result['pivots']):
        assert len(entry) == 3
        index, column, sign = entry
        assert type(index) is int and type(column) is int and type(sign) is int
        assert index in active and initial_size <= column < generator_count
        assert column not in removed and sign in (-1, 1)
        pivot = rows[index]
        assert pivot.get(column) == sign
        assert index in incidence[column]
        if step == 0:
            active = {i for i in active if rows[i]}
        affected = sorted(incidence[column] - {index})
        remove_memberships(index, pivot)
        active.remove(index)
        removed.add(column)
        for other in affected:
            assert other in active
            previous = rows[other]
            coefficient = previous[column]
            assert coefficient
            remove_memberships(other, previous)
            support = previous.keys() | pivot.keys()
            support_visits += len(support)
            updated = {}
            for j in support:
                value = previous.get(j, 0) - coefficient * sign * pivot.get(j, 0)
                if value:
                    updated[j] = value
            assert column not in updated
            rows[other] = updated
            substitutions += 1
            maximum_row_support = max(maximum_row_support, len(updated))
            if updated:
                for j in updated:
                    if j >= initial_size:
                        assert j not in removed
                        incidence.setdefault(j, set()).add(other)
            else:
                active.remove(other)
        rows[index] = {}
        assert column not in incidence
    assert all(not (j >= initial_size and abs(value) == 1)
               for i in active for j, value in rows[i].items())
    assert all(i in active and rows[i].get(j)
               for j, indices in incidence.items() for i in indices)
    rebuilt_incidence = {}
    for i in active:
        for j in rows[i]:
            if j >= initial_size:
                rebuilt_incidence.setdefault(j, set()).add(i)
    assert incidence == rebuilt_incidence

    row_ids = sorted(active)
    remaining = [j for j in range(generator_count) if j not in removed]
    positions = {j: i for i, j in enumerate(remaining)}
    output = [{positions[j]: value for j, value in rows[i].items()} for i in row_ids]
    present = {j for row in output for j in row}
    free = [j for j in range(initial_size, len(remaining)) if j not in present]
    free_set = set(free)
    torsion = [j for j in range(len(remaining)) if j not in free_set]
    expected = {
        'output_row_indices': row_ids,
        'new_generator_indices': remaining,
        'free_root_columns': free,
        'torsion_columns': torsion,
        'torsion_generator_indices': [remaining[j] for j in torsion],
        'torsion_rows': [[row.get(j, 0) for j in torsion] for row in output],
        'counts': {'input_rows': len(source), 'input_generator_count': generator_count,
                   'output_rows': len(output), 'free_components': len(free),
                   'torsion_generator_count': len(torsion),
                   'later_unit_columns_eliminated': len(result['pivots']),
                   'remaining_nonfree_residual_generators': len(torsion) - initial_size},
    }
    assert len(result['output_rows']) == len(output)
    for stored, rebuilt in zip(result['output_rows'], output):
        assert {int(j): value for j, value in stored.items()} == rebuilt
    for name, value in expected.items():
        assert result[name] == value, name
    return {'unit_isomorphisms': len(result['pivots']),
            'row_substitutions': substitutions,
            'substitution_support_visits': support_visits,
            'maximum_replayed_row_support': maximum_row_support,
            'full_integer_presentation_preserved': True,
            'original_sparse_matrix_checked': True,
            'expanded_composite_maps_materialized': False,
            'graph_maps_materialized': False}
