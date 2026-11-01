"""Compact certificates for exact residual unit substitutions.

The supplied graph certificate is shared. A pivot trace composes elementary
integer presentation isomorphisms without expanding every old source map.
"""
from hashlib import sha256
from heapq import heappop, heappush
import json

from signed_residual_graph import verify_signed_graph_elimination


def _rows(reduction):
    return [{int(j): value for j, value in row.items() if value}
            for row in reduction['output_rows']]


def _fingerprint(reduction):
    data = [reduction['initial_size'], reduction['new_generator_indices'],
            [sorted(row.items()) for row in _rows(reduction)]]
    return sha256(json.dumps(data, separators=(',', ':')).encode()).hexdigest()


def _finish(reduction, rows, active, removed, pivots):
    base = reduction['initial_size']
    width = len(reduction['new_generator_indices'])
    column_map = {old: new for new, old in enumerate(j for j in range(width) if j not in removed)}
    indices = [reduction['new_generator_indices'][j] for j in column_map]
    row_ids = sorted(active)
    output = [{column_map[j]: value for j, value in rows[i].items()} for i in row_ids]
    present = {j for row in output for j in row}
    free = [j for j in range(base, len(indices)) if j not in present]
    free_set = set(free)
    nonfree = [j for j in range(len(indices)) if j not in free_set]
    counts = dict(reduction['counts'])
    counts.update(output_rows=len(output), free_components=len(free),
                  torsion_generator_count=len(nonfree),
                  later_unit_columns_eliminated=len(pivots),
                  remaining_nonfree_residual_generators=len(nonfree)-base)
    return {'schema': 'integer-unit-substitution-trace-v1',
            'input_presentation_sha256': _fingerprint(reduction),
            'initial_size': base, 'input_row_count': len(rows),
            'input_generator_count': width, 'pivots': pivots,
            'output_row_indices': row_ids, 'new_generator_indices': indices,
            'output_rows': output, 'free_root_columns': free,
            'torsion_columns': nonfree,
            'torsion_generator_indices': [indices[j] for j in nonfree],
            'torsion_rows': [[row.get(j, 0) for j in nonfree] for row in output],
            'counts': counts}


def eliminate_unit_residual_columns_compact(original_rows, reduction, *, verify=True, verify_graph=True):
    """Produce a compact exact certificate with the reference pivot ordering."""
    rows = _rows(reduction)
    base = reduction['initial_size']
    width = len(reduction['new_generator_indices'])
    active, removed = set(range(len(rows))), set()
    incidence = [set() for _ in range(width)]
    versions = [0] * len(rows)
    heap, pivots = [], []

    def offer(index):
        for column, value in rows[index].items():
            if column >= base and abs(value) == 1:
                heappush(heap, (len(rows[index]), column, index, versions[index]))

    for index, row in enumerate(rows):
        for column in row:
            if column >= base:
                incidence[column].add(index)
        offer(index)

    while heap:
        support, column, pivot_index, version = heappop(heap)
        if pivot_index not in active or version != versions[pivot_index]:
            continue
        pivot = rows[pivot_index]
        if support != len(pivot) or abs(pivot.get(column, 0)) != 1:
            continue
        sign = pivot[column]
        expression = {j: -sign * value for j, value in pivot.items() if j != column}
        pivots.append([pivot_index, column, sign])
        if len(pivots) == 1:
            active.difference_update(i for i, row in enumerate(rows) if not row)
        active.remove(pivot_index)
        removed.add(column)
        for j in pivot:
            if j >= base:
                incidence[j].discard(pivot_index)
        for index in tuple(incidence[column]):
            row = rows[index]
            coefficient = row.pop(column)
            incidence[column].discard(index)
            for j, value in expression.items():
                previous = row.get(j, 0)
                updated = previous + coefficient * value
                if updated:
                    row[j] = updated
                    if j >= base and not previous:
                        incidence[j].add(index)
                else:
                    row.pop(j, None)
                    if j >= base:
                        incidence[j].discard(index)
            versions[index] += 1
            if row:
                offer(index)
            else:
                active.remove(index)
    result = _finish(reduction, rows, active, removed, pivots)
    if verify:
        result['verification'] = verify_unit_presentation_trace(
            original_rows, reduction, result, verify_graph=verify_graph)
    return result


def verify_unit_presentation_trace(original_rows, reduction, result, *, verify_graph=True):
    """Check each elementary isomorphism by an independent full row replay.

This verifier uses neither the incidence index nor the pivot heap. It never
needs to materialize the composite maps to the original generator set.
    """
    graph_check = verify_signed_graph_elimination(original_rows, reduction) if verify_graph else None
    assert result['schema'] == 'integer-unit-substitution-trace-v1'
    assert result['input_presentation_sha256'] == _fingerprint(reduction)
    rows = _rows(reduction)
    base = reduction['initial_size']
    width = len(reduction['new_generator_indices'])
    assert result['initial_size'] == base
    assert result['input_row_count'] == len(rows) and result['input_generator_count'] == width
    active, removed = set(range(len(rows))), set()
    substitutions = 0
    for step, (index, column, sign) in enumerate(result['pivots']):
        assert isinstance(index, int) and isinstance(column, int)
        assert index in active and base <= column < width and column not in removed
        assert sign in (-1, 1) and rows[index].get(column, 0) == sign
        if step == 0:
            active = {i for i in active if rows[i]}
        pivot = dict(rows[index])
        active.remove(index)
        removed.add(column)
        for i in tuple(active):
            coefficient = rows[i].get(column, 0)
            if coefficient:
                updated = dict(rows[i])
                for j, value in pivot.items():
                    total = updated.get(j, 0) - coefficient * sign * value
                    if total:
                        updated[j] = total
                    else:
                        updated.pop(j, None)
                rows[i] = updated
                substitutions += 1
                if not updated:
                    active.remove(i)
    assert all(not (j >= base and abs(value) == 1)
               for i in active for j, value in rows[i].items())
    expected = _finish(reduction, rows, active, removed, result['pivots'])
    for key, value in expected.items():
        assert json.loads(json.dumps(result[key])) == json.loads(json.dumps(value)), key
    return {'unit_isomorphisms': len(result['pivots']), 'row_substitutions': substitutions,
            'complete_original_graph_checked': graph_check is not None,
            'full_integer_presentation_preserved': True,
            'expanded_composite_maps_materialized': False}
