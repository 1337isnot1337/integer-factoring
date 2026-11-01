#!/usr/bin/env python3
"""Literal mixed-grade radical endpoints and two exact crossing repairs."""

import argparse
import json
import time
import hashlib
import subprocess
from pathlib import Path
from collections import defaultdict
from itertools import product
from functools import lru_cache
from explore_symplectic_boundary_codes import boundary_code, scalar_coefficient_polynomials, common_linear_factor_power
from check_double_endpoint_retracts import bit_indices, native_text, optimize
from explore_complement_edge_pencils import evaluate_local, independent_coefficients
from complement_shell_tensor import multiply
from check_modular_crossing import rank_two

def mixed_code(P, grades, radical=True):
    out = {'pairs': P, 'length': max(grades), 'labels': 0, 'ambient_occupations': [], 'input_vectors': [], 'output_indices': [], 'syndromes': [], 'input_Gram_rows': [], 'induced_mode_complement_rows': [], 'grades': grades}
    for L in grades:
        c = boundary_code(P, L, radical)
        offset = len(out['ambient_occupations'])
        n = out['labels']
        out['ambient_occupations'] += c['ambient_occupations']
        out['input_vectors'] += [v << offset for v in c['input_vectors']]
        out['output_indices'] += [x + offset for x in c['output_indices']]
        out['syndromes'] += c['syndromes']
        out['labels'] += c['labels']
        out['input_Gram_rows'] += [r << n for r in c['input_Gram_rows']]
        out['induced_mode_complement_rows'] += [r << n for r in c['induced_mode_complement_rows']]
    return out

def coefficients(code, mode='normalize'):
    labels, inputs, outputs = (code['ambient_occupations'], code['input_vectors'], code['output_indices'])
    d = code['labels']
    Lmax = code['length']
    base = Lmax + 1
    local = evaluate_local((2, 1 << base, 1, 1, 1, 1))
    supports = [bit_indices(v) for v in inputs]
    grades = [sum((m.bit_count() for m in labels[x])) for x in outputs]
    groups = defaultdict(list)
    for x, y in product(range(d), repeat=2):
        groups[code['syndromes'][x] ^ code['syndromes'][y]].append((x, y))
    keys = sorted(groups)
    pos = {pair: (b, i) for b, k in enumerate(keys) for i, pair in enumerate(groups[k])}
    size = Lmax * base + 1
    rows = [[[0] * len(groups[k]) for k in keys] for _ in range(size)]
    cols = [[[0] * len(groups[k]) for k in keys] for _ in range(size)]

    @lru_cache(None)
    def value(pattern):
        z = 1
        for index in pattern:
            z = multiply(z, local[index])
            if not z:
                break
        return z
    failures = []
    count = 0
    failure_count = 0
    for b, key in enumerate(keys):
        group = groups[key]
        for ri, (x, y) in enumerate(group):
            for ci, (v, u) in enumerate(group):
                z = 0
                for yy, vv in product(supports[y], supports[v]):
                    pattern = tuple(sorted((((a * 4 + bb) * 4 + c) * 4 + e for a, bb, c, e in zip(labels[outputs[x]], labels[yy], labels[outputs[u]], labels[vv]))))
                    z ^= value(pattern)
                if mode == 'balanced' and (grades[x] != grades[v] or grades[y] != grades[u]):
                    z = 0
                if z & 1 != int(x == v and y == u):
                    failure_count += 1
                    if len(failures) < 12:
                        failures.append((x, y, u, v, grades[x], grades[y], grades[u], grades[v]))
                cb, cr = pos[x, u]
                other, cc = pos[v, y]
                assert cb == other
                for degree in bit_indices(z):
                    rows[degree][b][ri] ^= 1 << ci
                    cols[degree][cb][cr] ^= 1 << cc
                count += 1
    if mode == 'normalize':
        inverses = []
        for M in rows[0]:
            n = len(M)
            work = [r | 1 << n + i for i, r in enumerate(M)]
            for i in range(n):
                pivot = next((j for j in range(i, n) if work[j] >> i & 1))
                work[i], work[pivot] = (work[pivot], work[i])
                for j in range(n):
                    if j != i and work[j] >> i & 1:
                        work[j] ^= work[i]
            inverses.append([v >> n for v in work])
        for M in rows:
            for b, inv in enumerate(inverses):
                source = M[b]
                M[b] = []
                for ir in inv:
                    value = 0
                    for i in bit_indices(ir):
                        value ^= source[i]
                    M[b].append(value)
        assert all((M == [1 << i for i in range(len(M))] for M in rows[0]))
        cols = [[[0] * len(groups[k]) for k in keys] for _ in range(size)]
        for a, M in enumerate(rows):
            for b, key in enumerate(keys):
                group = groups[key]
                for ri, value in enumerate(M[b]):
                    x, y = group[ri]
                    for ci in bit_indices(value):
                        v, u = group[ci]
                        cb, cr = pos[x, u]
                        other, cc = pos[v, y]
                        assert cb == other
                        cols[a][cb][cr] ^= 1 << cc
    elif failures:
        return ({'r_endpoint_failures': failures, 'r_endpoint_failure_count': failure_count, 'entries': count}, None, None)
    active = [i for i, M in enumerate(rows) if any((any(b) for b in M))]
    rows = [rows[i] for i in active]
    cols = [cols[i] for i in active]
    selected = independent_coefficients(rows)
    assert active[selected[0]] == 0
    polynomials, coords = scalar_coefficient_polynomials(rows, active, selected)
    data = {'labels': d, 'length': Lmax, 'basis_dimension': len(selected), 'block_sizes': [len(groups[k]) for k in keys], 'basis_monomials_p_q': [[active[i] % base, active[i] // base] for i in selected], 'endpoint_coordinates': [1] + [0] * (len(selected) - 1), 'universal_coefficient_profile': [[a % base, a // base, sum((rank_two(b) for b in r)), sum((rank_two(b) for b in c))] for a, r, c in zip(active, rows, cols)], 'basis_scalar_polynomials_p_q': [[[a % base, a // base] for a in bit_indices(v)] for v in polynomials], 'entries': count, 'grades_by_label': grades, 'raw_endpoint_failure_count': failure_count, 'raw_endpoint_failure_examples': failures}
    return (data, [rows[i] for i in selected], [cols[i] for i in selected])


def main():
    from analyze_pfaffian_binary_screen import analyze, unpack_binary
    from check_pfaffian_span_f4 import verify
    from check_symplectic_seed_tensor_fields import integer_gate
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=4)
    parser.add_argument("--grades", default="2,4,6")
    parser.add_argument("--modes", default="raw,normalize,balanced")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--field4", action="store_true")
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    code = mixed_code(args.pairs, list(map(int, args.grades.split(","))))
    args.export_dir.mkdir(parents=True, exist_ok=True)
    result = {"code": code, "cases": [],
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    for mode in args.modes.split(","):
        assert mode in ("raw", "normalize", "balanced")
        began = time.monotonic()
        data, rows, columns = coefficients(code, mode)
        record = dict(data, mode=mode)
        if rows is not None:
            name = f"omega_mixed_P{args.pairs}_{mode}"
            basis_path = args.export_dir / (name + ".json")
            exported = dict(data, name=name, code=code,
                            basis_row_blocks=[[[hex(value) for value in block] for block in matrix]
                                              for matrix in rows],
                            basis_column_blocks=[[[hex(value) for value in block] for block in matrix]
                                                 for matrix in columns])
            basis_path.write_text(json.dumps(exported, indent=2, sort_keys=True) + "\n")
            native = native_text(data, rows, columns, integer_gate(code["labels"]))
            (args.export_dir / (name + ".txt")).write_text(native)
            screen_path = args.export_dir / (name + "-F2.jsonl")
            completed = subprocess.run([str(args.native), "--field2"], input=native,
                                       text=True, stdout=subprocess.PIPE, check=True)
            screen_path.write_text(completed.stdout)
            record["F2_screen"] = analyze(basis_path, screen_path)
            audited = 0
            for event in map(json.loads, completed.stdout.splitlines()):
                if event["event"] != "eligible":
                    continue
                mask = unpack_binary(event["code"], len(rows))
                ranks = []
                for basis in (rows, columns):
                    total = 0
                    for block, size in enumerate(data["block_sizes"]):
                        work = [0] * size
                        for i in bit_indices(mask):
                            work = [a ^ b for a, b in zip(work, basis[i][block])]
                        total += rank_two(work)
                    ranks.append(total)
                assert ranks == [event["row_rank"], event["column_rank"]]
                audited += 1
            record["F2_eligible_vectors_independently_audited"] = audited
            if len(rows) <= 10:
                record["F2_optimization"] = optimize(data, rows, columns)
            if args.field4 and len(rows) <= 10:
                field_path = args.export_dir / (name + "-F4.jsonl")
                completed = subprocess.run([str(args.native)], input=native,
                                           text=True, stdout=subprocess.PIPE, check=True)
                field_path.write_text(completed.stdout)
                record["F4_screen"] = verify(basis_path, field_path, 40)
        record["seconds"] = time.monotonic() - began
        result["cases"].append(record)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"mode": mode, "labels": code["labels"],
                          "span": data.get("basis_dimension"),
                          "raw_endpoint_failures": data.get("r_endpoint_failure_count",
                                                             data.get("raw_endpoint_failure_count")),
                          "F2": record.get("F2_screen"), "F4": record.get("F4_screen"),
                          "seconds": record["seconds"]}), flush=True)


if __name__ == "__main__":
    main()

