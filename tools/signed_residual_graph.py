"""Integer presentation elimination for residual words of absolute weight <=2.

Keep one double-root relation per unbalanced component. Return mutually
inverse presentation maps with exact integer certificates. No factoring
or rational row selection is used by this transformation.
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


def eliminate_signed_residual_graph(rows, initial_size, *, generator_count=None,
                                    allow_extra_rows=False, verify=True):
    """Eliminate tree/unit relations and preserve necessary double roots.

    Rows are dense integer lists or sparse integer-column dictionaries.
    Supply generator_count to retain trailing isolated generators in a
    sparse input. Rows have fixed-base columns followed by residual columns. The latter
have absolute weight <=2, allowing a single coefficient +/-2. All
    unanchored roots occur in the full output; free root columns are omitted
    from torsion_rows because they occur in no relation. With allow_extra_rows,
    larger residual words are transported through the same exact substitution.
    """
    rows = [dict(row) if isinstance(row, dict) else list(row) for row in rows]
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

    free_columns = [i for i in range(initial_size, len(new_generators))
                    if all(row.get(i, 0) == 0 for row in output_rows)]
    torsion_columns = [i for i in range(len(new_generators)) if i not in free_columns]
    result = {
        'initial_size': initial_size, 'old_generator_count': count,
        'new_generator_indices': new_generators, 'generator_words': generator_words,
        'generator_definition_sources': definition_sources,
        'output_rows': output_rows, 'output_sources': output_sources,
        'old_row_combinations': old_row_combinations,
        'free_root_columns': free_columns, 'torsion_columns': torsion_columns,
        'torsion_generator_indices': [new_generators[i] for i in torsion_columns],
        'torsion_rows': [[row.get(i, 0) for i in torsion_columns] for row in output_rows],
        'components': components,
        'counts': {'input_rows': len(rows), 'input_residual_generators': residual_count,
                   'graph_rows': len(eligible), 'extra_rows': len(rows)-len(eligible),
                   'tree_rows_eliminated': sum(len(c['tree_rows']) for c in components),
                   'anchored_components': sum(c['anchor_row'] is not None for c in components),
                   'free_components': len(free_columns), 'unbalanced_components': len(pivots),
                   'output_rows': len(output_rows), 'torsion_generator_count': len(torsion_columns)},
    }
    if verify:
        result['verification'] = verify_signed_graph_elimination(rows, result)
    return result


def verify_signed_graph_elimination(rows, result):
    """Replay both integer presentation maps without graph traversal."""
    old = [_sparse(row) for row in rows]
    normalize = lambda row: {int(i): value for i, value in row.items() if value}
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
    for column in result['free_root_columns']:
        assert all(row.get(column, 0) == 0 for row in output)
    columns = result['torsion_columns']
    assert columns == [i for i in range(len(chosen)) if i not in result['free_root_columns']]
    assert result['torsion_generator_indices'] == [chosen[i] for i in columns]
    assert result['torsion_rows'] == [[row.get(i, 0) for i in columns] for row in output]
    return {'inverse_generator_identities': len(chosen),
            'old_generator_identities': old_count,
            'output_relation_identities': len(output),
            'input_relation_identities': len(old), 'integer_presentation_isomorphism': True}
