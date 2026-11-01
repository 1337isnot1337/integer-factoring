"""Sparse unit elimination with the reference integer presentation certificates."""
from copy import deepcopy
from heapq import heappop, heappush

from signed_residual_graph import verify_signed_graph_elimination


def _add(target, source, coefficient):
    for column, value in source.items():
        updated = target.get(column, 0) + coefficient * value
        if updated:
            target[column] = updated
        else:
            target.pop(column, None)


def eliminate_sparse_unit_residual_columns(original_rows, reduction, *, verify=True):
    """Preserve the reference pivot rule, using incidence and stable identifiers.

Rows and columns retain their initial identifiers throughout elimination.
Only vectors containing a pivot column receive its substitution. A lazy
heap preserves the (support, column, row) pivot ordering. The final pass
renumbers survivors once. Expanded certificates remain fully charged.
    """
    result = deepcopy(reduction)
    names = ('generator_words', 'generator_definition_sources', 'output_rows',
             'output_sources', 'old_row_combinations')
    for name in names:
        result[name] = [{int(i): x for i, x in row.items() if x}
                        for row in result[name]]
    base = result['initial_size']
    width = len(result['new_generator_indices'])
    rows, sources = result['output_rows'], result['output_sources']
    words, definitions = result['generator_words'], result['generator_definition_sources']
    row_incidence = [set() for _ in range(width)]
    word_incidence = [set() for _ in range(width)]
    active_rows = set(range(len(rows)))
    removed_columns = set()
    versions = [0] * len(rows)
    heap, steps = [], []

    def offer(index):
        row = rows[index]
        for column, value in row.items():
            if column >= base and abs(value) == 1:
                heappush(heap, (len(row), column, index, versions[index]))

    for index, row in enumerate(rows):
        for column in row:
            if column >= base:
                row_incidence[column].add(index)
        offer(index)
    for index, word in enumerate(words):
        for column in word:
            if column >= base:
                word_incidence[column].add(index)

    def substitute(vector, incidence, index, column, expression):
        coefficient = vector.pop(column)
        incidence[column].discard(index)
        for j, value in expression.items():
            previous = vector.get(j, 0)
            updated = previous + coefficient * value
            if updated:
                vector[j] = updated
                if j >= base and not previous:
                    incidence[j].add(index)
            else:
                vector.pop(j, None)
                if j >= base:
                    incidence[j].discard(index)
        return coefficient

    while heap:
        support, column, pivot_index, version = heappop(heap)
        if pivot_index not in active_rows or version != versions[pivot_index]:
            continue
        pivot = rows[pivot_index]
        if support != len(pivot) or abs(pivot.get(column, 0)) != 1:
            continue
        source = dict(sources[pivot_index])
        sign = pivot[column]
        expression = {j: -sign * x for j, x in pivot.items() if j != column}
        steps.append({'old_generator_index': result['new_generator_indices'][column],
                      'pivot_coefficient': sign, 'source': source})
        removed_columns.add(column)

        # The reference discards pre-existing zero rows at its first pivot.
        if len(steps) == 1:
            active_rows.difference_update(i for i, row in enumerate(rows) if not row)
        active_rows.remove(pivot_index)
        for j in pivot:
            if j >= base:
                row_incidence[j].discard(pivot_index)

        for index in tuple(word_incidence[column]):
            coefficient = substitute(words[index], word_incidence, index, column, expression)
            _add(definitions[index], source, coefficient * sign)

        for index in tuple(row_incidence[column]):
            coefficient = substitute(rows[index], row_incidence, index, column, expression)
            _add(sources[index], source, -coefficient * sign)
            versions[index] += 1
            if rows[index]:
                offer(index)
            else:
                active_rows.remove(index)

    column_map = {old: new for new, old in enumerate(
        j for j in range(width) if j not in removed_columns)}
    kept_rows = sorted(active_rows)
    row_map = {old: new for new, old in enumerate(kept_rows)}
    result['new_generator_indices'] = [result['new_generator_indices'][j] for j in column_map]
    result['generator_words'] = [{column_map[j]: value for j, value in word.items()}
                                 for word in words]
    result['output_rows'] = [{column_map[j]: value for j, value in rows[i].items()}
                            for i in kept_rows]
    result['output_sources'] = [sources[i] for i in kept_rows]
    result['old_row_combinations'] = [
        {row_map[i]: value for i, value in combination.items() if i in row_map}
        for combination in result['old_row_combinations']]
    for component in result['components']:
        current = component['root_column']
        if current is not None:
            if current in removed_columns:
                component['root_column'] = None
                component['root_removed_by_later_unit'] = True
            else:
                component['root_column'] = column_map[current]

    new_count = len(result['new_generator_indices'])
    present = {j for row in result['output_rows'] for j in row}
    free = [j for j in range(base, new_count) if j not in present]
    free_set = set(free)
    columns = [j for j in range(new_count) if j not in free_set]
    result['free_root_columns'] = free
    result['torsion_columns'] = columns
    result['torsion_generator_indices'] = [result['new_generator_indices'][j] for j in columns]
    result['torsion_rows'] = [[row.get(j, 0) for j in columns] for row in result['output_rows']]
    result['unit_elimination_steps'] = steps
    result['counts'].update(output_rows=len(result['output_rows']), free_components=len(free),
                            torsion_generator_count=len(columns),
                            later_unit_columns_eliminated=len(steps),
                            remaining_nonfree_residual_generators=len(columns)-base)
    if verify:
        result['verification'] = verify_signed_graph_elimination(original_rows, result)
    return result
