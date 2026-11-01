#!/usr/bin/env python3
"""Binary polynomial row reduction and complete bounded Laurent spaces.

The 0/1 saturation is provisional until a polynomial right inverse is
checked. This generalizes the algorithm in joint_hot_retraction_module.py.
"""
from joint_hot_retraction_module import degree, dependency, divide_linear


def row_reduced_basis(source):
    rows = [list(row) for row in source]
    width = len(rows[0])
    steps = []
    for root in (0, 1):
        while True:
            evaluations = [sum(((p & 1) if root == 0 else p.bit_count() % 2) << j
                               for j, p in enumerate(row)) for row in rows]
            relation = dependency(evaluations)
            if not relation:
                break
            selected = [i for i in range(len(rows)) if relation >> i & 1]
            pivot = max(selected, key=lambda i: degree(rows[i]))
            new = [0]*width
            for i in selected:
                new = [a ^ b for a, b in zip(new, rows[i])]
            if any(new):
                rows[pivot] = [divide_linear(p, root) for p in new]
            else:
                rows.pop(pivot)
            steps.append(dict(root=root, rows=selected, pivot=pivot))
    reductions = 0
    while True:
        ds = [degree(row) for row in rows]
        relation = dependency([sum(((p >> d) & 1) << j for j, p in enumerate(row))
                               for row, d in zip(rows, ds)])
        if not relation:
            break
        selected = [i for i in range(len(rows)) if relation >> i & 1]
        pivot = max(selected, key=lambda i: ds[i])
        new = [0]*width
        for i in selected:
            new = [a ^ (p << (ds[pivot]-ds[i])) for a, p in zip(new, rows[i])]
        if any(new):
            assert degree(new) < ds[pivot]
            rows[pivot] = new
        else:
            rows.pop(pivot)
        reductions += 1
    return rows, dict(saturation_steps=steps, leading_reductions=reductions)


def bounded_laurent_space(basis, lower, upper):
    """All strict rows in span(basis), with common lower and entrywise caps.

    A row is returned as z^-clearing times its polynomial numerator.
    Independent constant and leading rows give complete coefficient bounds.
    """
    width = len(basis[0])
    if len(upper) != width or not isinstance(lower, int):
        raise ValueError('Supply one integer lower bound and a degree cap per column.')
    assert not dependency([sum((p & 1) << j for j, p in enumerate(row)) for row in basis])
    assert not dependency([sum(((p >> degree(row)) & 1) << j for j, p in enumerate(row)) for row in basis])
    strict = [[p ^ (p << 1) for p in row] for row in basis]
    clearing = max(0, -lower)
    packing = max(upper)+clearing+1
    allowed = [((1 << (u+clearing+1))-1) ^ ((1 << (lower+clearing))-1)
               if u >= lower else 0 for u in upper]
    generators = [[p << (shift+clearing) for p in row] for row in strict
                  for shift in range(lower, max(upper)-degree(row)+1)]
    pivots, output = {}, []
    for i, generator in enumerate(generators):
        image = sum((p & ~mask) << (j*packing) for j, (p, mask) in enumerate(zip(generator, allowed)))
        combo = 1 << i
        while image:
            k = image.bit_length()-1
            if k in pivots:
                image ^= pivots[k][0]
                combo ^= pivots[k][1]
            else:
                pivots[k] = image, combo
                break
        if not image:
            row = [0]*width
            for j, g in enumerate(generators):
                if combo >> j & 1:
                    row = [a ^ b for a, b in zip(row, g)]
            assert any(row) and all(not (p & ~m) for p, m in zip(row, allowed))
            assert all(p.bit_count() % 2 == 0 for p in row)
            output.append(row)
    return dict(rows=output, clearing=clearing, dimension=len(output),
                coefficient_variables=len(generators), lower=lower, upper=list(upper))
