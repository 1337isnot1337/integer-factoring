"""Remove remaining unit residual pivots from a certified graph presentation."""
from copy import deepcopy

from signed_residual_graph import verify_signed_graph_elimination


def _add(target, source, coefficient):
    for column, value in source.items():
        result = target.get(column, 0) + coefficient*value
        if result:
            target[column] = result
        else:
            target.pop(column, None)


def eliminate_unit_residual_columns(original_rows, reduction, *, verify=True):
    """Preserve all integer presentation maps while eliminating +/-1 pivots.

Only residual columns can be removed. Base columns remain. Double-root
relations are retained unless another genuine unit row eliminates their
root, in which case substitution transfers the relation exactly.
    """
    result = deepcopy(reduction)
    names = ('generator_words', 'generator_definition_sources', 'output_rows',
             'output_sources', 'old_row_combinations')
    for name in names:
        result[name] = [{int(i): x for i, x in row.items() if x} for row in result[name]]
    base = result['initial_size']
    steps = []
    while True:
        candidates = [(len(row), column, index) for index, row in enumerate(result['output_rows'])
                      for column, value in row.items() if column >= base and abs(value) == 1]
        if not candidates:
            break
        _, column, pivot_index = min(candidates)
        pivot = dict(result['output_rows'][pivot_index])
        source = dict(result['output_sources'][pivot_index])
        sign = pivot[column]
        expression = {j: -sign*x for j, x in pivot.items() if j != column}
        steps.append({'old_generator_index': result['new_generator_indices'][column],
                      'pivot_coefficient': sign, 'source': source})

        for word, certificate in zip(result['generator_words'], result['generator_definition_sources']):
            coefficient = word.pop(column, 0)
            if coefficient:
                _add(word, expression, coefficient)
                _add(certificate, source, coefficient*sign)

        retained_rows, retained_sources, row_map = [], [], {}
        for index, (row, certificate) in enumerate(zip(result['output_rows'], result['output_sources'])):
            if index == pivot_index:
                continue
            row, certificate = dict(row), dict(certificate)
            coefficient = row.pop(column, 0)
            if coefficient:
                _add(row, expression, coefficient)
                _add(certificate, source, -coefficient*sign)
            if row:
                row_map[index] = len(retained_rows)
                retained_rows.append(row)
                retained_sources.append(certificate)
        result['old_row_combinations'] = [
            {row_map[i]: value for i, value in combination.items() if i in row_map}
            for combination in result['old_row_combinations']]
        result['output_rows'], result['output_sources'] = retained_rows, retained_sources
        result['new_generator_indices'].pop(column)
        for name in ('generator_words', 'output_rows'):
            result[name] = [{j-(j > column): value for j, value in row.items()}
                            for row in result[name]]
        for component in result['components']:
            current = component['root_column']
            if current is not None:
                if current == column:
                    component['root_column'] = None
                    component['root_removed_by_later_unit'] = True
                elif current > column:
                    component['root_column'] -= 1

    new_count = len(result['new_generator_indices'])
    free = [j for j in range(base, new_count) if all(row.get(j, 0) == 0 for row in result['output_rows'])]
    columns = [j for j in range(new_count) if j not in free]
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
