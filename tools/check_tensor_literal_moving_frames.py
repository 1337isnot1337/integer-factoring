#!/usr/bin/env python3
"""Independent raw-Laurent-kernel checks for larger literal moving frames."""
import argparse
import ast
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import check_tensor_hot_compounds_f16 as controls
import tensor_hot_compounds_f16 as field
import tensor_literal_hot_pfaffians_f16 as unweighted
import tensor_literal_weighted_maps_f16 as weighted
from tensor_literal_moving_frames import cut_blocks, frame_code, serialized_code


@lru_cache(None)
def raw_pfaffian(modes):
    source = Path(__file__).with_name('reflection_dual_rail_gaussian.py')
    tree = ast.parse(source.read_text())
    wanted = {'validate', 'laurent_kernel'}
    selected = ast.Module(body=[node for node in tree.body
                               if isinstance(node, ast.FunctionDef) and node.name in wanted], type_ignores=[])
    namespace = {}
    exec(compile(selected, str(source), 'exec'), namespace)
    interval = 4
    while interval < 4*modes:
        interval *= 2
    size = interval//2
    full = namespace['laurent_kernel'](interval)
    ids = [leg*size+size-1-hot for leg in range(4) for hot in range(modes)]
    assert all(e % 2 == 0 for i in ids for j in ids for e in full[i][j])
    kernel = [[tuple(e//2 for e in full[i][j]) for j in ids] for i in ids]

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return frozenset([0])
        if mask.bit_count() % 2:
            return frozenset()
        first = mask & -mask
        i = first.bit_length()-1
        rest, remaining, result = mask ^ first, mask ^ first, set()
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length()-1
            for shift in kernel[i][j]:
                for e in pfaffian(rest ^ chosen):
                    if e+shift in result:
                        result.remove(e+shift)
                    else:
                        result.add(e+shift)
        return frozenset(result)
    return pfaffian


def reference_entry(code, x, y, u, v, cut, clearing):
    p, left, modes, heights = (code[k] for k in ('embedding', 'retraction', 'modes', 'heights'))
    pfaffian = raw_pfaffian(modes)
    planes = [0]*4
    for tx, ty, tu, tv in product(left[x].items(), left[y].items(), p[u].items(), p[v].items()):
        terms = (tx, ty, tu, tv)
        mask = sum(w << (leg*modes) for leg, ((w, _), _) in enumerate(terms))
        coefficient = 1
        exponent = 0
        for (_, e), c in terms:
            coefficient = field.mul(coefficient, c)
            exponent += e
        for raw in pfaffian(mask):
            if cut == 'zero':
                power = raw+exponent+heights[x]+heights[y]+clearing
            else:
                power = -raw-exponent+heights[u]+heights[v]+clearing
            assert power >= 0
            for bit in range(4):
                if coefficient >> bit & 1:
                    planes[bit] ^= 1 << power
    return tuple(planes)


def reverse(polynomial, degree):
    assert polynomial.bit_length()-1 <= degree
    return sum(((polynomial >> i) & 1) << (degree-i) for i in range(degree+1))


def check(code, homogeneous_control=False):
    q, modes = code['q'], code['modes']
    matrices, clearings, profiles = [], [], []
    entry_count, toeplitz_count, capped_count = 0, 0, 0
    for cut in ('zero', 'infinity'):
        blocks, clearing = cut_blocks(code, cut)
        clearings.append(clearing)
        dense = [[(0, 0, 0, 0) for _ in range(q*q)] for _ in range(q*q)]
        records = []
        for block in blocks:
            native = weighted.block(modes, block['rows'], block['columns'],
                                    clearing=clearing, coefficients=True)
            tests = controls.toeplitz_controls(native['matrix'], native['cleared_orders_t'])
            toeplitz_count += len(tests)
            for cap in (8, 24, 48):
                capped = weighted.block(modes, block['rows'], block['columns'],
                                        clearing=clearing, order_cap=cap)
                assert capped['orders_t'] == [min(cap, order) for order in native['orders_t']]
                capped_count += len(block['pairs'])
            for a, (i, j) in enumerate(block['pairs']):
                for b, (k, l) in enumerate(block['pairs']):
                    dense[q*i+j][q*k+l] = native['matrix'][a][b]
            if homogeneous_control:
                assert clearing == 0
                row_plain = [{w: c for (w, shift), c in vector.items() if shift == 0} for vector in block['rows']]
                column_plain = [{w: c for (w, shift), c in vector.items() if shift == 0} for vector in block['columns']]
                assert all(len(a) == len(b) for a, b in zip(row_plain, block['rows']))
                assert all(len(a) == len(b) for a, b in zip(column_plain, block['columns']))
                old = unweighted.compound_block(modes, row_plain, column_plain, coefficients=True)
                assert old['matrix'] == native['matrix'] and old['orders_t'] == native['orders_t']
            records.append({'parity': block['parity'], 'orders_t': native['orders_t'],
                            'cleared_orders_t': native['cleared_orders_t'], 'toeplitz_checks': tests,
                            'degree_bound': native['cleared_polynomial_degree_bound'],
                            'precision_attempts': native['precision_attempts']})
        matrices.append(dense)
        profiles.append({'cut': cut, 'clearing_power': clearing, 'blocks': records})
    endpoint_count, reciprocal_count = 0, 0
    for x, y, u, v in product(range(q), repeat=4):
        zero = matrices[0][q*x+y][q*v+u]
        infinity = matrices[1][q*x+u][q*v+y]
        assert zero == reference_entry(code, x, y, u, v, 'zero', clearings[0])
        assert infinity == reference_entry(code, x, y, u, v, 'infinity', clearings[1])
        entry_count += 2
        assert sum((plane.bit_count() % 2) << bit for bit, plane in enumerate(zero)) == int(x == v and y == u)
        endpoint_count += 1
        degree = sum(code['heights'][i] for i in (x, y, u, v))+sum(clearings)
        assert infinity == tuple(reverse(plane, degree) for plane in zero)
        reciprocal_count += 1
    return {'code': serialized_code(code), 'homogeneous_control': homogeneous_control,
            'profiles': profiles, 'raw_laurent_reference_entries': entry_count,
            'independent_toeplitz_rank_checks': toeplitz_count, 'capped_order_comparisons': capped_count,
            'endpoint_checks': endpoint_count, 'infinity_reciprocity_checks': reciprocal_count,
            'raw_pfaffian_cache_size': raw_pfaffian(modes).cache_info().currsize}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-moving-frame-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    cases = []
    specifications = [(0, 0, 1, 1, 1), (1, 3, 1, 1, 1),
                      (3, 7, 1, 2, 3), (15, 15, 0, 1, 1), (7, 15, -1, 1, 1)]
    for a_mask, b_mask, flow, a_scalar, b_scalar in specifications:
        code = frame_code(10, 9, [0, 3, 5, 6], a_mask, b_mask, flow, a_scalar, b_scalar)
        cases.append(check(code, homogeneous_control=(a_mask == b_mask == 0)))
    names = [Path(__file__).name, 'tensor_literal_moving_frames.py',
             'tensor_literal_weighted_maps_f16.py', 'tensor_literal_weighted_maps_f16.cpp',
             'tensor_literal_hot_pfaffians_f16.py', 'tensor_literal_hot_pfaffians_f16.cpp',
             'reflection_dual_rail_gaussian.py', 'check_tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'cases': cases,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'seconds': time.monotonic()-started}
    for key in ('raw_laurent_reference_entries', 'independent_toeplitz_rank_checks',
                'capped_order_comparisons', 'endpoint_checks', 'infinity_reciprocity_checks'):
        report[key] = sum(case[key] for case in cases)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('cases', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
