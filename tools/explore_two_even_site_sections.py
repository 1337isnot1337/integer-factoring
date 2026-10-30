#!/usr/bin/env python3
"""Exhaust graph-supported two-even-site Pfaffian retracts and their sections.

For each edge of K4, choose omission, its lower endpoint, or its upper
endpoint as output.  All four single-mode orientations are retained on
each selected edge.  The resulting 3^6-1 finite spaces are distinct
retract choices, not tensor powers of one fixed 24-label coefficient space.
"""

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import subprocess

from analyze_pfaffian_binary_screen import analyze, unpack_binary
from check_double_endpoint_retracts import bit_indices, native_text
from check_modular_crossing import rank_two
from check_symplectic_seed_tensor_fields import integer_gate
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local, independent_coefficients
from explore_symplectic_boundary_codes import boundary_code, coefficient_basis


def two_even_code(pairs):
    original = boundary_code(pairs, pairs)
    selected = [i for i, sigma in enumerate(original["syndromes"])
                if sigma.bit_count() == pairs - 2]
    code = dict(original)
    code["labels"] = len(selected)
    for key in ("input_vectors", "output_indices", "syndromes"):
        code[key] = [original[key][i] for i in selected]
    for key in ("input_Gram_rows", "induced_mode_complement_rows"):
        code[key] = [sum(((original[key][i] >> j) & 1) << k
                        for k, j in enumerate(selected)) for i in selected]
    assert code["labels"] == pairs * (pairs - 1) // 2 * (1 << (pairs - 2))
    assert all(value.bit_count() == 2 for value in code["input_vectors"])
    return code


def parent_coefficients():
    """Rectangular ambient blocks containing every uniform edge section."""
    code = two_even_code(4)
    words, inputs, sigmas = (code[key] for key in
                            ("ambient_occupations", "input_vectors", "syndromes"))
    edges = list(combinations(range(4), 2))
    edge_index = {sum(1 << i for i in edge): j for j, edge in enumerate(edges)}
    edge_of = [edge_index[15 ^ sigma] for sigma in sigmas]
    supports = [bit_indices(value) for value in inputs]
    outputs = []
    for support in supports:
        ordered = sorted(support, key=lambda i: next(j for j, x in enumerate(words[i]) if x == 3))
        outputs.extend(ordered)
    groups = []
    positions = []
    # Each tuple below is a row or column index in one of the two cuts.
    for first, second in ((48, 24), (24, 48), (48, 48), (24, 24)):
        bucket = defaultdict(list)
        for a, b in product(range(first), range(second)):
            aa, bb = a // 2 if first == 48 else a, b // 2 if second == 48 else b
            bucket[sigmas[aa] ^ sigmas[bb]].append((a, b))
        grouped = [bucket[key] for key in sorted(bucket)]
        groups.append(grouped)
        positions.append({pair: (block, i) for block, group in enumerate(grouped)
                          for i, pair in enumerate(group)})
    rows = [[[0] * len(group) for group in groups[0]] for _ in range(21)]
    cols = [[[0] * len(group) for group in groups[2]] for _ in range(21)]
    local = evaluate_local((2, 1 << 5, 1, 1, 1, 1))

    @lru_cache(None)
    def evaluate(pattern):
        value = 1
        for entry in pattern:
            value = multiply(value, local[entry])
            if not value:
                break
        return value

    count = 0
    for block, group in enumerate(groups[0]):
        for ri, (x, y) in enumerate(group):
            for ci, (v, u) in enumerate(groups[1][block]):
                value = 0
                for yy, vv in product(supports[y], supports[v]):
                    pattern = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + e
                                           for a, b, c, e in zip(words[outputs[x]], words[yy],
                                                                 words[outputs[u]], words[vv])))
                    value ^= evaluate(pattern)
                assert value & 1 == int(x // 2 == v and y == u // 2)
                assert not (value >> 4) & 1
                cb, cr = positions[2][x, u]
                other, cc = positions[3][v, y]
                assert cb == other
                for degree in bit_indices(value):
                    assert degree % 5 + degree // 5 <= 4
                    rows[degree][block][ri] ^= 1 << ci
                    cols[degree][cb][cr] ^= 1 << cc
                count += 1
    active = [degree for degree, matrix in enumerate(rows) if any(any(b) for b in matrix)]
    selected = independent_coefficients([rows[degree] for degree in active])
    assert active[selected[0]] == 0
    return {"code": code, "edges": edges, "edge_of": edge_of, "positions": positions,
            "rows": [rows[active[i]] for i in selected],
            "columns": [cols[active[i]] for i in selected],
            "ambient_coefficient_monomials": [[active[i] % 5, active[i] // 5] for i in selected],
            "ambient_entries_checked": count, "distinct_local_patterns": evaluate.cache_info().currsize}


def restrict(parent, choices):
    selected = [i for i, edge in enumerate(parent["edge_of"]) if choices[edge]]
    variants = {i: 2 * i + choices[parent["edge_of"][i]] - 1 for i in selected}
    sigmas = parent["code"]["syndromes"]
    buckets = defaultdict(list)
    for i, j in product(selected, repeat=2):
        buckets[sigmas[i] ^ sigmas[j]].append((i, j))
    groups = [buckets[key] for key in sorted(buckets)]
    cuts = []
    for cut, key in enumerate(("rows", "columns")):
        selectors = []
        for group in groups:
            if cut == 0:
                rowpairs = [(variants[i], j) for i, j in group]
                colpairs = [(i, variants[j]) for i, j in group]
            else:
                rowpairs = [(variants[i], variants[j]) for i, j in group]
                colpairs = list(group)
            rowpos = [parent["positions"][2 * cut][pair] for pair in rowpairs]
            colpos = [parent["positions"][2 * cut + 1][pair] for pair in colpairs]
            blocks = {b for b, i in rowpos + colpos}
            assert len(blocks) == 1
            selectors.append((blocks.pop(), [i for b, i in rowpos], [i for b, i in colpos]))
        matrices = []
        for matrix in parent[key]:
            blocks = []
            for block, ri, ci in selectors:
                blocks.append([sum(((matrix[block][r] >> c) & 1) << j for j, c in enumerate(ci))
                               for r in ri])
            matrices.append(blocks)
        cuts.append(matrices)
    row, col = cuts
    basis = independent_coefficients(row)
    assert basis[0] == 0
    row, col = [[cut[i] for i in basis] for cut in (row, col)]
    for matrix in (row[0], col[0]):
        assert all(block == [1 << i for i in range(len(block))] for block in matrix)
    data = {"name": "two_even_sections_" + "".join(map(str, choices)),
            "labels": len(selected), "length": 4, "pairs": 4,
            "section_choices": list(choices), "selected_original_labels": selected,
            "selected_output_variants": [variants[i] for i in selected],
            "basis_dimension": len(basis), "block_sizes": list(map(len, groups)),
            "basis_parent_indices": basis, "endpoint_coordinates": [1] + [0] * (len(basis) - 1)}
    return data, row, col


def save_basis(path, data, rows, columns):
    result = dict(data, basis_row_blocks=[[[hex(x) for x in b] for b in m] for m in rows],
                  basis_column_blocks=[[[hex(x) for x in b] for b in m] for m in columns])
    path.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fixed-upper-only", action="store_true")
    parser.add_argument("--stop-after", type=int)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    parent = parent_coefficients()
    result = {"construction": "two-even-site edge subcodes with uniform binary output sections",
              "choices_per_edge": ["omit", "lower", "upper"],
              "edges": parent["edges"],
              "ambient_entries_checked": parent["ambient_entries_checked"],
              "ambient_coefficient_monomials": parent["ambient_coefficient_monomials"],
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in (Path(__file__).name, "explore_symplectic_boundary_codes.py",
                               "explore_complement_edge_pencils.py", "screen_pfaffian_span_f4.cpp",
                               "analyze_pfaffian_binary_screen.py")},
              "cases": [], "complete": False}
    alphabet = (0, 2) if args.fixed_upper_only else (0, 1, 2)
    for choices in product(alphabet, repeat=6):
        if not any(choices):
            continue
        data, rows, cols = restrict(parent, choices)
        name = data["name"]
        basis_path = args.work_dir / (name + ".json")
        screen_path = args.work_dir / (name + "-F2.jsonl")
        save_basis(basis_path, data, rows, cols)
        native = native_text(data, rows, cols, integer_gate(data["labels"]))
        run = subprocess.run([str(args.native), "--field2"], input=native, text=True,
                             stdout=subprocess.PIPE, check=True)
        screen_path.write_text(run.stdout)
        record = analyze(basis_path, screen_path)
        record.update(section_choices=list(choices), name=name)
        events = [json.loads(line) for line in run.stdout.splitlines()]
        checked = 0
        for event in events:
            if event["event"] != "eligible":
                continue
            mask = unpack_binary(event["code"], len(rows))
            ranks = []
            for cut in (rows, cols):
                rank = 0
                for bi, size in enumerate(data["block_sizes"]):
                    work = [0] * size
                    for i in bit_indices(mask):
                        work = [a ^ b for a, b in zip(work, cut[i][bi])]
                    rank += rank_two(work)
                ranks.append(rank)
            assert ranks == [event["row_rank"], event["column_rank"]]
            checked += 1
        record["eligible_rank_pairs_independently_checked"] = checked
        result["cases"].append(record)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(record), flush=True)
        if args.stop_after and len(result["cases"]) >= args.stop_after:
            break
    result["complete"] = len(result["cases"]) == len(alphabet) ** 6 - 1
    result["all_excluded_over_F2"] = result["complete"] and all(
        x["status"].startswith("excluded:") for x in result["cases"])
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
