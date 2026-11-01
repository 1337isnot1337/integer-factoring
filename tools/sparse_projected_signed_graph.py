"""The complete signed-graph presentation with a sparse projected view.

Derived from signed_residual_graph.py at SHA256
f7374233f0b6d2284316b5ddd7a5d9c32e24f7d2f9e81b2121ce77d8f142d025.
All four integer presentation maps and the graph traversal are unchanged.
The projected view stays sparse, support discovery uses sets, and completed
construction scratch is released before verification. Read-only input rows
and normalized dictionaries are reused; projection checks process one row.
"""
from collections import deque


def _add(target, source, multiple=1):
    for key, value in source.items():
        result = target.get(key, 0) + multiple * value
        if result:
            target[key] = result
        else:
            target.pop(key, None)
    return target


def _sparse(row):
    if isinstance(row, dict):
        return {i: value for i, value in row.items() if value}
    return {i: value for i, value in enumerate(row) if value}


def eliminate_sparse_projected_signed_graph(rows, initial_size, *, generator_count=None,
                                    allow_extra_rows=False, verify=True):
    """Eliminate tree/unit relations and preserve necessary double roots.

    Rows are dense integer lists or sparse integer-column dictionaries.
    Supply generator_count to retain trailing isolated generators in a
    sparse input. Rows have fixed-base columns followed by residual columns. The latter
have absolute weight <=2, allowing a single coefficient +/-2. All
    unanchored roots occur in the full output; free root columns are omitted
    from sparse_torsion_rows because they occur in no relation. With allow_extra_rows,
    larger residual words are transported through the same exact substitution.
    """
    # Neither construction nor verification mutates these supplied rows.
    rows = [row if isinstance(row, (dict, list, tuple)) else list(row) for row in rows]
    if not isinstance(initial_size, int) or initial_size < 0:
        raise ValueError('Require a nonnegative fixed-base size.')
    dense_widths = [len(row) for row in rows if not isinstance(row, dict)]
    sparse_indices = [i for row in rows if isinstance(row, dict) for i in row]
    if any(not isinstance(i, int) or i < 0 for i in sparse_indices):
        raise ValueError('Sparse column indices must be nonnegative integers.')
    inferred = dense_widths[0] if dense_widths else max([initial_size] + [i+1 for i in sparse_indices])
    count = inferred if generator_count is None else generator_count
    if (not isinstance(count, int) or count < initial_size
            or any(width != count for width in dense_widths)
            or any(i >= count for i in sparse_indices)
            or any(any(not isinstance(x, int) for x in
                       (row.values() if isinstance(row, dict) else row)) for row in rows)
            or (not allow_extra_rows and
                any(sum(abs(x) for i, x in _sparse(row).items() if i >= initial_size) > 2
                    for row in rows))):
        raise ValueError('Require rectangular integer rows with residual absolute weight <=2.')
    residual_count = count - initial_size
    sparse_rows = [_sparse(row) for row in rows]
    base_rows = [{j: x for j, x in row.items() if j < initial_size} for row in sparse_rows]
    residual_rows = [{j-initial_size: x for j, x in row.items() if j >= initial_size}
                     for row in sparse_rows]
    eligible = {i for i, word in enumerate(residual_rows) if sum(map(abs, word.values())) <= 2}
    adjacency = [[] for _ in range(residual_count)]
    for index, word in enumerate(residual_rows):
        if index in eligible and len(word) == 2:
            u, v = word
            adjacency[u].append((v, index))
            adjacency[v].append((u, index))

    components, visited = [], set()
    constants = [{} for _ in range(residual_count)]
    signs = [0]*residual_count
    definitions = [{} for _ in range(residual_count)]
    component_of = [-1]*residual_count
    for root in range(residual_count):
        if root in visited:
            continue
        identifier = len(components)
        vertices, tree = [], []
        signs[root] = 1
        visited.add(root)
        queue = deque([root])
        while queue:
            u = queue.popleft()
            vertices.append(u)
            component_of[u] = identifier
            for v, edge in adjacency[u]:
                if v in visited:
                    continue
                visited.add(v)
                queue.append(v)
                tree.append(edge)
                su, sv = residual_rows[edge][u], residual_rows[edge][v]
                multiplier = -su*sv
                signs[v] = multiplier*signs[u]
                constants[v] = _add({j: -sv*x for j, x in base_rows[edge].items()}, constants[u], multiplier)
                definitions[v] = _add({edge: sv}, definitions[u], multiplier)
        components.append({'root': root, 'vertices': vertices, 'tree_rows': tree,
                           'anchor_row': None, 'root_column': None})

    for index, word in enumerate(residual_rows):
        if index in eligible and len(word) == 1:
            vertex, coefficient = next(iter(word.items()))
            if abs(coefficient) == 1:
                component = components[component_of[vertex]]
                if component['anchor_row'] is None:
                    component['anchor_row'] = index

    new_generators = list(range(initial_size))
    for component in components:
        anchor = component['anchor_row']
        if anchor is None:
            component['root_column'] = len(new_generators)
            new_generators.append(initial_size + component['root'])
            continue
        vertex, coefficient = next(iter(residual_rows[anchor].items()))
        root_coefficient = coefficient * signs[vertex]
        assert abs(root_coefficient) == 1
        q = _add(dict(base_rows[anchor]), constants[vertex], coefficient)
        source = _add({anchor: 1}, definitions[vertex], -coefficient)
        for v in component['vertices']:
            _add(constants[v], q, -signs[v]*root_coefficient)
            _add(definitions[v], source, signs[v]*root_coefficient)
            signs[v] = 0

    generator_words = [{i: 1} for i in range(initial_size)]
    definition_sources = [{} for _ in range(initial_size)]
    for vertex in range(residual_count):
        word = dict(constants[vertex])
        if signs[vertex]:
            word[components[component_of[vertex]]['root_column']] = signs[vertex]
        generator_words.append(word)
        definition_sources.append(definitions[vertex])

    mapped_rows, mapped_sources = [], []
    for index, row in enumerate(sparse_rows):
        mapped, source = {}, {index: 1}
        for column, value in row.items():
            _add(mapped, generator_words[column], value)
            _add(source, definition_sources[column], -value)
        mapped_rows.append(mapped)
        mapped_sources.append(source)

    pivots = {}
    for index, row in enumerate(mapped_rows):
        if index not in eligible:
            continue
        roots = [(column, value) for column, value in row.items() if column >= initial_size]
        assert len(roots) <= 1
        if roots:
            column, value = roots[0]
            assert abs(value) == 2
            pivots.setdefault(column, index)
    output_rows, output_sources, pivot_outputs = [], [], {}
    for column, index in sorted(pivots.items()):
        sign = mapped_rows[index][column] // 2
        pivot_outputs[column] = len(output_rows)
        output_rows.append({j: sign*x for j, x in mapped_rows[index].items()})
        output_sources.append({j: sign*x for j, x in mapped_sources[index].items()})

    old_row_combinations = []
    for index, original in enumerate(mapped_rows):
        reduced, source, combination = dict(original), dict(mapped_sources[index]), {}
        for column, pivot_index in pivot_outputs.items():
            coefficient = reduced.get(column, 0)
            if coefficient:
                multiple = coefficient // 2
                _add(reduced, output_rows[pivot_index], -multiple)
                _add(source, output_sources[pivot_index], -multiple)
                combination[pivot_index] = multiple
        if reduced:
            if index in eligible:
                assert all(j < initial_size for j in reduced)
            combination[len(output_rows)] = 1
            output_rows.append(reduced)
            output_sources.append(source)
        old_row_combinations.append(combination)

    present = {j for row in output_rows for j, value in row.items() if value}
    free_columns = [i for i in range(initial_size, len(new_generators)) if i not in present]
    free_set = set(free_columns)
    torsion_columns = [i for i in range(len(new_generators)) if i not in free_set]
    projected_positions = {old: new for new, old in enumerate(torsion_columns)}
    result = {
        'schema': 'sparse-projected-signed-graph-v1',
        'initial_size': initial_size, 'old_generator_count': count,
        'new_generator_indices': new_generators, 'generator_words': generator_words,
        'generator_definition_sources': definition_sources,
        'output_rows': output_rows, 'output_sources': output_sources,
        'old_row_combinations': old_row_combinations,
        'free_root_columns': free_columns, 'torsion_columns': torsion_columns,
        'torsion_generator_indices': [new_generators[i] for i in torsion_columns],
        'sparse_torsion_rows': [{projected_positions[j]: value for j, value in row.items()}
                                for row in output_rows],
        'components': components,
        'counts': {'input_rows': len(rows), 'input_residual_generators': residual_count,
                   'graph_rows': len(eligible), 'extra_rows': len(rows)-len(eligible),
                   'tree_rows_eliminated': sum(len(c['tree_rows']) for c in components),
                   'anchored_components': sum(c['anchor_row'] is not None for c in components),
                   'free_components': len(free_columns), 'unbalanced_components': len(pivots),
                   'output_rows': len(output_rows), 'torsion_generator_count': len(torsion_columns)},
    }
    # These construction arrays are not part of the certificate. Retaining
    # them during independent map replay needlessly doubles the live state.
    del sparse_rows, base_rows, residual_rows, adjacency, constants, signs
    del definitions, component_of, mapped_rows, mapped_sources, visited, eligible
    if verify:
        result['verification'] = verify_sparse_projected_signed_graph(rows, result)
    return result


def verify_sparse_projected_signed_graph(rows, result):
    """Replay both integer presentation maps without graph traversal."""
    def normalize(row):
        if all(isinstance(i, int) and value for i, value in row.items()):
            return row
        return {int(i): value for i, value in row.items() if value}

    old = [normalize(row) if isinstance(row, dict) else _sparse(row) for row in rows]
    words = [normalize(row) for row in result['generator_words']]
    output = [normalize(row) for row in result['output_rows']]
    sources = [normalize(row) for row in result['output_sources']]
    definitions = [normalize(row) for row in result['generator_definition_sources']]
    combinations = [normalize(row) for row in result['old_row_combinations']]
    chosen = result['new_generator_indices']
    old_count = result['old_generator_count']
    assert all((all(isinstance(i, int) and 0 <= i < old_count for i in row)
                if isinstance(row, dict) else len(row) == old_count) for row in rows)
    assert len(words) == len(definitions) == old_count
    assert len(combinations) == len(old) and len(sources) == len(output)
    assert len(set(chosen)) == len(chosen) and all(0 <= i < old_count for i in chosen)

    def combine(coefficients, vectors):
        answer = {}
        for i, coefficient in coefficients.items():
            assert 0 <= i < len(vectors) and isinstance(coefficient, int)
            for j, value in vectors[i].items():
                answer[j] = answer.get(j, 0) + coefficient*value
        return {j: value for j, value in answer.items() if value}

    def lift(word):
        assert all(0 <= i < len(chosen) and isinstance(value, int) for i, value in word.items())
        return {chosen[i]: value for i, value in word.items()}

    for index, old_index in enumerate(chosen):
        assert words[old_index] == {index: 1}
    for index, word in enumerate(words):
        difference = {j: -value for j, value in lift(word).items()}
        difference[index] = difference.get(index, 0) + 1
        difference = {j: value for j, value in difference.items() if value}
        assert difference == combine(definitions[index], old)
    for row, certificate in zip(output, sources):
        assert lift(row) == combine(certificate, old)
    for row, certificate in zip(old, combinations):
        assert combine(row, words) == combine(certificate, output)
    present = {column for row in output for column in row}
    for column in result['free_root_columns']:
        assert column not in present
    columns = result['torsion_columns']
    free_set = set(result['free_root_columns'])
    assert columns == [i for i in range(len(chosen)) if i not in free_set]
    assert result['torsion_generator_indices'] == [chosen[i] for i in columns]
    assert result['schema'] == 'sparse-projected-signed-graph-v1'
    assert len(result['sparse_torsion_rows']) == len(output)
    projected_positions = {old: new for new, old in enumerate(columns)}
    for stored, row in zip(result['sparse_torsion_rows'], output):
        assert normalize(stored) == {projected_positions[j]: value for j, value in row.items()}
    return {'inverse_generator_identities': len(chosen),
            'old_generator_identities': old_count,
            'output_relation_identities': len(output),
            'input_relation_identities': len(old), 'integer_presentation_isomorphism': True}


def materialize_dense_projection(result):
    """Compatibility expansion for bounded comparisons, never used in search."""
    assert result['schema'] == 'sparse-projected-signed-graph-v1'
    expanded = dict(result)
    del expanded['schema']
    sparse = expanded.pop('sparse_torsion_rows')
    width = len(expanded['torsion_columns'])
    dense = []
    for row in sparse:
        normalized = {int(j): value for j, value in row.items()}
        dense.append([normalized.get(i, 0) for i in range(width)])
    expanded['torsion_rows'] = dense
    return expanded
