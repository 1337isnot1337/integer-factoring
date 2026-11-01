#!/usr/bin/env python3
"""Exact Pfaffian retracts cancelling one or both noncrossing endpoints."""

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import random
import time

from complement_shell_tensor import multiply, pair_occupations
from check_complement_shell_crossings import shell
from check_modular_crossing import rank_two
from explore_complement_edge_pencils import (
    all_length_labels, evaluate_local, independent_coefficients,
    minimum_decomposition, syndrome,
)


def mode_complement(word):
    return tuple((0, 2, 1, 3)[mask] for mask in word)


def sparse_retract(labels, seed=0, within_syndrome=True):
    labels = tuple(labels)
    lookup = {word: i for i, word in enumerate(labels)}
    J = [lookup[mode_complement(word)] for word in labels]
    groups = defaultdict(list)
    fixed = 0
    rng = random.Random(seed)
    for i, j in enumerate(J):
        if i == j:
            fixed += 1
        elif i < j:
            representative = i if seed == 0 else (i, j)[rng.randrange(2)]
            key = syndrome(labels[i]) if within_syndrome else 0
            groups[key].append(representative)
    pairs = []
    discarded = 0
    for key in sorted(groups):
        values = groups[key]
        if seed:
            rng.shuffle(values)
        discarded += len(values) % 2
        pairs.extend(zip(values[::2], values[1::2]))
    outputs = [a for a, b in pairs]
    assert len(set(outputs)) == len(outputs)
    assert len({i for pair in pairs for i in pair}) == 2 * len(pairs)
    for x, (a, b) in enumerate(pairs):
        for y, (c, d) in enumerate(pairs):
            assert int(a in (c, d)) == int(x == y)
            assert len({a, b} & {c, d}) % 2 == 0
            assert a not in (J[c], J[d])
    return {
        "ambient_occupations": labels, "mode_complement_permutation": J,
        "input_pairs": pairs, "output_indices": outputs,
        "fixed_complement_labels": fixed, "discarded_representatives": discarded,
        "within_syndrome": within_syndrome, "seed": seed,
    }


def invariant_orbit_retract(labels, seed=0, pair_fixed=True):
    labels = tuple(labels)
    lookup = {word: i for i, word in enumerate(labels)}
    J = [lookup[mode_complement(word)] for word in labels]
    rng = random.Random(seed)
    pairs = [(i, j) if seed == 0 or rng.randrange(2) == 0 else (j, i)
             for i, j in enumerate(J) if i < j]
    fixed = [i for i, j in enumerate(J) if i == j]
    if pair_fixed:
        if seed:
            rng.shuffle(fixed)
        pairs.extend(zip(fixed[::2], fixed[1::2]))
    outputs = [a for a, b in pairs]
    for x, (a, b) in enumerate(pairs):
        for y, (c, e) in enumerate(pairs):
            assert int(a in (c, e)) == int(x == y)
            assert len({a, b} & {c, e}) % 2 == 0
            assert int(a in (J[c], J[e])) == int(x == y)
    return {
        "ambient_occupations": labels, "mode_complement_permutation": J,
        "input_pairs": pairs, "output_indices": outputs,
        "fixed_complement_labels": len(fixed),
        "discarded_representatives": len(fixed) % 2 if pair_fixed else len(fixed),
        "within_syndrome": True, "seed": seed, "q_endpoint": "noncrossing",
        "construction": "invariant_orbit_with_fixed_pairs" if pair_fixed else "invariant_orbit",
    }


def bit_indices(value):
    result = []
    while value:
        bit = value & -value
        result.append(bit.bit_length() - 1)
        value ^= bit
    return result


def check_linear_maps(retract):
    inputs = retract["input_vectors"]
    outputs = retract["output_vectors"]
    J = retract["mode_complement_permutation"]
    for i, first in enumerate(inputs):
        transformed = sum(1 << J[j] for j in bit_indices(first))
        for j, second in enumerate(inputs):
            assert (first & second).bit_count() % 2 == 0
            assert (outputs[j] & first).bit_count() % 2 == int(i == j)
            assert (outputs[j] & transformed).bit_count() % 2 == 0
    return len(inputs) ** 2 * 3


def triangular_retract(labels, use_fixed=True, group_mode="one", seed=0):
    labels = tuple(labels)
    lookup = {word: i for i, word in enumerate(labels)}
    J = [lookup[mode_complement(word)] for word in labels]
    orbits = [(i, j) for i, j in enumerate(J) if i < j]
    fixed = [i for i, j in enumerate(J) if i == j]
    m = len(orbits)
    rng = random.Random(seed)
    if seed:
        orbits = [(a, b) if rng.randrange(2) == 0 else (b, a) for a, b in orbits]
        rng.shuffle(orbits)
    if group_mode != "one":
        if group_mode == "syndrome":
            buckets = defaultdict(list)
            for orbit in orbits:
                buckets[syndrome(labels[orbit[0]])].append(orbit)
            groups = [buckets[key] for key in sorted(buckets)]
        else:
            number = len(fixed) if group_mode == "fixed" else int(group_mode)
            if number < 1:
                raise ValueError("at least one triangular group is required")
            groups = [orbits[i::number] for i in range(number)]
        inputs, outputs = [], []
        for group_index, group in enumerate(groups):
            if not group:
                continue
            anchor = fixed[group_index] if use_fixed and group_index < len(fixed) else None
            for i, (a, b) in enumerate(group if anchor is not None else group[:-1]):
                if anchor is not None:
                    first = (1 << a) | (1 << anchor)
                    second = 1 << a
                    for c, e in group[:i]:
                        first ^= (1 << c) | (1 << e)
                    for c, e in group[i + 1:]:
                        second ^= (1 << c) | (1 << e)
                else:
                    first = (1 << b) | (1 << group[-1][1])
                    second = 1 << b
                    for c, e in group[i + 1:-1]:
                        first ^= (1 << c) | (1 << e)
                    for c, e in group[:i]:
                        second ^= (1 << c) | (1 << e)
                inputs.append(first)
                outputs.append(second)
        construction = "grouped_triangular_" + group_mode
    elif fixed and use_fixed:
        inputs = []
        outputs = []
        for i, (a, b) in enumerate(orbits):
            first = (1 << a) | (1 << fixed[0])
            second = 1 << a
            for c, e in orbits[:i]:
                first ^= (1 << c) | (1 << e)
            for c, e in orbits[i + 1:]:
                second ^= (1 << c) | (1 << e)
            inputs.append(first)
            outputs.append(second)
        construction = "triangular_with_fixed_coordinate"
    else:
        inputs = []
        outputs = []
        for i, (a, b) in enumerate(orbits[:-1]):
            first = (1 << b) | (1 << orbits[-1][1])
            second = 1 << b
            for c, e in orbits[i + 1:-1]:
                first ^= (1 << c) | (1 << e)
            for c, e in orbits[:i]:
                second ^= (1 << c) | (1 << e)
            inputs.append(first)
            outputs.append(second)
        construction = "triangular_without_fixed_coordinate"
    result = {
        "ambient_occupations": labels, "mode_complement_permutation": J,
        "input_vectors": inputs, "output_vectors": outputs,
        "fixed_complement_labels": len(fixed), "construction": construction,
        "within_syndrome": False, "group_mode": group_mode, "seed": seed,
    }
    result["linear_identity_entries_checked"] = check_linear_maps(result)
    return result


@lru_cache(maxsize=4)
def ambient_coefficients(labels):
    ambient = len(labels)
    length = sum(mask.bit_count() for mask in labels[0])
    base = length + 1
    local = evaluate_local((2, 1 << base, 1, 1, 1, 1))
    groups = defaultdict(list)
    for x, y in product(range(ambient), repeat=2):
        groups[syndrome(labels[x]) ^ syndrome(labels[y])].append((x, y))
    ambient_rows = [[0] * (ambient * ambient) for _ in range(length * base + 1)]

    @lru_cache(None)
    def pattern_value(pattern):
        result = 1
        for index in pattern:
            result = multiply(result, local[index])
            if not result:
                break
        return result

    evaluated = 0
    for group in groups.values():
        for x, y in group:
            row = x * ambient + y
            for v, u in group:
                pattern = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + e
                                       for a, b, c, e in zip(labels[x], labels[y], labels[u], labels[v])))
                value = pattern_value(pattern)
                column = 1 << (v * ambient + u)
                while value:
                    bit = value & -value
                    ambient_rows[bit.bit_length() - 1][row] ^= column
                    value ^= bit
                evaluated += 1
    return length, ambient_rows, evaluated, pattern_value.cache_info().currsize


def dense_coefficient_blocks(retract):
    """Contract all four legs by exact binary linear maps, then form both cuts."""
    labels = tuple(tuple(word) for word in retract["ambient_occupations"])
    ambient = len(labels)
    inputs = retract["input_vectors"]
    outputs = retract["output_vectors"]
    d = len(inputs)
    length, ambient_rows, evaluated, local_patterns = ambient_coefficients(labels)
    base = length + 1
    input_supports = [bit_indices(value) for value in inputs]
    output_supports = [bit_indices(value) for value in outputs]
    left_supports = [[a * ambient + b for a in output_supports[x]
                     for b in input_supports[y]] for x, y in product(range(d), repeat=2)]
    input_transpose = [sum(1 << i for i, value in enumerate(inputs) if value >> a & 1)
                       for a in range(ambient)]
    output_transpose = [sum(1 << i for i, value in enumerate(outputs) if value >> a & 1)
                        for a in range(ambient)]
    right_columns = [sum(output_transpose[b] << (v * d)
                         for v in bit_indices(input_transpose[a]))
                     for a, b in product(range(ambient), repeat=2)]
    tables = []
    for first in range(0, ambient * ambient, 8):
        columns = right_columns[first:first + 8]
        table = [0] * (1 << len(columns))
        for mask in range(1, len(table)):
            bit = mask & -mask
            table[mask] = table[mask ^ bit] ^ columns[bit.bit_length() - 1]
        tables.append(table)

    def transform_right(row):
        result = 0
        index = 0
        while row:
            result ^= tables[index][row & 255]
            row >>= 8
            index += 1
        return result

    degrees = []
    rows = []
    cols = []
    digest = hashlib.sha256()
    for degree, matrix in enumerate(ambient_rows):
        if not any(matrix):
            continue
        intermediate = [transform_right(row) for row in matrix]
        transformed = []
        for support in left_supports:
            value = 0
            for index in support:
                value ^= intermediate[index]
            transformed.append(value)
        if degree in (length, length * base):
            assert not any(transformed)
        if degree == 0:
            assert transformed == [1 << i for i in range(d * d)]
        if not any(transformed):
            continue
        columns = [0] * (d * d)
        for index, value in enumerate(transformed):
            x, y = divmod(index, d)
            digest.update(value.to_bytes((d * d + 7) // 8, "little"))
            while value:
                bit = value & -value
                v, u = divmod(bit.bit_length() - 1, d)
                columns[x * d + u] ^= 1 << (v * d + y)
                value ^= bit
        degrees.append(degree)
        rows.append([transformed])
        cols.append([columns])
    selected = independent_coefficients(rows)
    assert degrees[selected[0]] == 0
    data = {
        "labels": d, "length": length,
        "nonzero_universal_monomials": [[degree % base, degree // base,
                                          length - degree % base - degree // base]
                                         for degree in degrees],
        "basis_coefficient_indices": selected, "basis_dimension": len(selected),
        "universal_coefficient_profile": [
            [degree % base, degree // base, length - degree % base - degree // base,
             sum(rank_two(block) for block in row), sum(rank_two(block) for block in column)]
            for degree, row, column in zip(degrees, rows, cols)],
        "block_sizes": [d * d], "endpoint_coordinates": [1] + [0] * (len(selected) - 1),
        "complete_entries_checked": d ** 4, "ambient_entries_evaluated": evaluated,
        "ambient_distinct_local_patterns": local_patterns,
        "dense_entry_count": d ** 4, "universal_profile_sha256": digest.hexdigest(),
        "both_noncrossing_endpoints_zero": True,
    }
    return data, [rows[i] for i in selected], [cols[i] for i in selected]


def coefficient_blocks(retract):
    labels = retract["ambient_occupations"]
    pairs = retract["input_pairs"]
    outputs = retract["output_indices"]
    d = len(pairs)
    if d < 2:
        raise ValueError("at least two retracted labels are required")
    length = sum(mask.bit_count() for mask in labels[0])
    assert all(sum(mask.bit_count() for mask in word) == length for word in labels)
    base = length + 1
    local = evaluate_local((2, 1 << base, 1, 1, 1, 1))
    groups = defaultdict(list)
    sigmas = [syndrome(labels[i]) for i in outputs]
    for x, y in product(range(d), repeat=2):
        key = sigmas[x] ^ sigmas[y] if retract["within_syndrome"] else 0
        groups[key].append((x, y))
    keys = sorted(groups)
    positions = {pair: (block, i) for block, key in enumerate(keys)
                 for i, pair in enumerate(groups[key])}
    size = length * base + 1
    rows = [[[0] * len(groups[key]) for key in keys] for _ in range(size)]
    cols = [[[0] * len(groups[key]) for key in keys] for _ in range(size)]

    @lru_cache(None)
    def entry(x, y, u, v):
        value = 1
        for a, b, c, e in zip(labels[x], labels[y], labels[u], labels[v]):
            value = multiply(value, local[((a * 4 + b) * 4 + c) * 4 + e])
            if not value:
                break
        return value

    evaluated = 0
    digest = hashlib.sha256()
    for block, key in enumerate(keys):
        group = groups[key]
        for i, (x, y) in enumerate(group):
            for j, (v, u) in enumerate(group):
                value = 0
                for yy, vv in product(pairs[y], pairs[v]):
                    value ^= entry(outputs[x], yy, outputs[u], vv)
                assert value & 1 == int(x == v and y == u)
                assert not (value >> length) & 1
                expected_q = int(x == y and u == v) if retract.get("q_endpoint") == "noncrossing" else 0
                assert (value >> (length * base)) & 1 == expected_q
                digest.update(value.to_bytes((size + 7) // 8, "little"))
                cblock, ci = positions[x, u]
                cblock_check, cj = positions[v, y]
                assert cblock == cblock_check
                while value:
                    bit = value & -value
                    degree = bit.bit_length() - 1
                    assert degree % base + degree // base <= length
                    rows[degree][block][i] ^= 1 << j
                    cols[degree][cblock][ci] ^= 1 << cj
                    value ^= bit
                evaluated += 1
    degrees = [degree for degree in range(size)
               if any(any(block) for block in rows[degree])]
    rows = [rows[degree] for degree in degrees]
    cols = [cols[degree] for degree in degrees]
    selected = independent_coefficients(rows)
    assert degrees[selected[0]] == 0
    data = {
        "labels": d, "length": length,
        "nonzero_universal_monomials": [
            [degree % base, degree // base, length - degree % base - degree // base]
            for degree in degrees],
        "basis_coefficient_indices": selected, "basis_dimension": len(selected),
        "universal_coefficient_profile": [
            [degree % base, degree // base, length - degree % base - degree // base,
             sum(rank_two(block) for block in row), sum(rank_two(block) for block in column)]
            for degree, row, column in zip(degrees, rows, cols)],
        "block_sizes": [len(groups[key]) for key in keys],
        "endpoint_coordinates": [1] + [0] * (len(selected) - 1),
        "complete_entries_checked": evaluated, "dense_entry_count": d ** 4,
        "universal_profile_sha256": digest.hexdigest(),
        "both_noncrossing_endpoints_zero": retract.get("q_endpoint") != "noncrossing",
        "p_endpoint_zero": True, "q_endpoint": retract.get("q_endpoint", "zero"),
    }
    return data, [rows[i] for i in selected], [cols[i] for i in selected]


def optimize(data, rows, cols):
    width = len(rows)
    if width > 16:
        raise ValueError("use the native exported basis when the binary span exceeds 16")
    rwork = [[0] * size for size in data["block_sizes"]]
    cwork = [[0] * size for size in data["block_sizes"]]
    both = [[0, 0] for _ in range(1 << width)]
    previous = 0
    for index in range(1, 1 << width):
        gray = index ^ (index >> 1)
        changed = (gray ^ previous).bit_length() - 1
        for target, source in ((rwork, rows[changed]), (cwork, cols[changed])):
            for block, delta in zip(target, source):
                for i, value in enumerate(delta):
                    block[i] ^= value
        both[gray] = [sum(rank_two(block) for block in rwork),
                      sum(rank_two(block) for block in cwork)]
        previous = gray
    upper, masks = minimum_decomposition([math.sqrt(r * c) for r, c in both], target=1)
    lower, _ = minimum_decomposition([math.isqrt(r * c) for r, c in both], target=1)
    d = data["labels"]
    return {
        "combination_ranks": both, "found_decomposition_cost": upper,
        "decomposition_masks": masks, "integer_shortest_lower_bound": lower,
        "finite_seed_threshold": d ** math.sqrt(2),
        "all_F2_decompositions_fail_seed_criterion": lower ** 70 > d ** 99,
    }


def orientation_sweep(labels):
    """Exhaust the nonfixed orbit representatives, modulo simultaneous J."""
    initial = invariant_orbit_retract(labels)
    J = initial["mode_complement_permutation"]
    moving = sum(a != J[a] for a, b in initial["input_pairs"])
    if moving > 12:
        raise ValueError("the complete orientation sweep is limited to twelve moving orbits")
    records = []
    best = None
    lower = None
    digest = hashlib.sha256()
    for mask in range(1 << max(0, moving - 1)):
        retract = dict(initial)
        retract["input_pairs"] = [(b, a) if mask >> i & 1 else (a, b)
                                  for i, (a, b) in enumerate(initial["input_pairs"])]
        retract["output_indices"] = [a for a, b in retract["input_pairs"]]
        data, rows, cols = coefficient_blocks(retract)
        result = optimize(data, rows, cols)
        digest.update(json.dumps(result["combination_ranks"], separators=(",", ":")).encode())
        item = {
            "orientation_mask": mask, "basis_dimension": data["basis_dimension"],
            "integer_shortest_lower_bound": result["integer_shortest_lower_bound"],
            "found_decomposition_cost": result["found_decomposition_cost"],
        }
        records.append(item)
        lower = item["integer_shortest_lower_bound"] if lower is None else min(lower, item["integer_shortest_lower_bound"])
        if best is None or item["found_decomposition_cost"] < best["found_decomposition_cost"]:
            best = dict(item, decomposition_masks=result["decomposition_masks"],
                        decomposition_ranks=[result["combination_ranks"][i] for i in result["decomposition_masks"]])
    d = len(initial["input_pairs"])
    return {
        "labels": d, "nonfixed_orbits": moving, "orientations_checked": len(records),
        "symmetry": "Simultaneous J on both output legs leaves the tensor unchanged because both input legs are J invariant.",
        "fixed_monomial_pairing": initial["input_pairs"][moving:],
        "minimum_integer_shortest_lower_bound": lower,
        "all_orientations_fail_F2_seed_criterion": lower ** 70 > d ** 99,
        "best_found_decomposition": best, "rank_tables_sha256": digest.hexdigest(),
        "per_orientation": records,
    }


def export(name, retract, data, rows, cols):
    result = dict(data)
    result.update(
        name=name, retract=retract,
        basis_monomials_p_q_r=[data["nonzero_universal_monomials"][i]
                              for i in data["basis_coefficient_indices"]],
        basis_row_blocks=[[[hex(row) for row in block] for block in matrix] for matrix in rows],
        basis_column_blocks=[[[hex(row) for row in block] for block in matrix] for matrix in cols],
        exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        format="Direct-sum cut blocks, with each binary row encoded as a hexadecimal bit mask.",
    )
    return result


def verify_dense_contraction():
    """Compare the four-leg bitset contraction with literal four-term entries."""
    sparse = sparse_retract(all_length_labels(3, 3))
    data, rows, cols = coefficient_blocks(sparse)
    dense = dict(sparse)
    dense["input_vectors"] = [sum(1 << i for i in pair) for pair in sparse["input_pairs"]]
    dense["output_vectors"] = [1 << i for i in sparse["output_indices"]]
    check_linear_maps(dense)
    other, dense_rows, dense_cols = dense_coefficient_blocks(dense)
    assert data["nonzero_universal_monomials"] == other["nonzero_universal_monomials"]
    assert data["basis_coefficient_indices"] == other["basis_coefficient_indices"]
    sigmas = [syndrome(sparse["ambient_occupations"][i]) for i in sparse["output_indices"]]
    groups = defaultdict(list)
    d = data["labels"]
    for a, b in product(range(d), repeat=2):
        groups[sigmas[a] ^ sigmas[b]].append(a * d + b)
    indices = [groups[key] for key in sorted(groups)]
    checks = 0
    for blocked_cut, dense_cut in ((rows, dense_rows), (cols, dense_cols)):
        for matrix, expected in zip(blocked_cut, dense_cut):
            expanded = [0] * (d * d)
            for block, positions in zip(matrix, indices):
                for row, position in zip(block, positions):
                    expanded[position] = sum(1 << positions[i] for i in bit_indices(row))
            assert expanded == expected[0]
            checks += d ** 4
    return {"scalar_matrix_entries_compared": checks, "labels": d,
            "independent_coefficients": data["basis_dimension"],
            "both_cuts_match_literal_four_term_expansion": True}


def native_text(data, rows, cols, gate):
    lines = [f'{data["labels"]} {len(rows)} {len(data["block_sizes"])} {gate}',
             " ".join(map(str, data["block_sizes"]))]
    for cut in (rows, cols):
        for matrix in cut:
            for block, size in zip(matrix, data["block_sizes"]):
                stride = (size + 63) // 64
                for row in block:
                    lines.append(" ".join(f'{row >> (64 * k) & ((1 << 64) - 1):x}'
                                          for k in range(stride)))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="full:3:3,full:4:2,full:4:3,full:4:4")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--cross-syndrome", action="store_true")
    parser.add_argument("--construction", choices=("sparse", "triangular", "invariant"), default="sparse")
    parser.add_argument("--ignore-fixed-coordinate", action="store_true")
    parser.add_argument("--groups", default="one", help="one, syndrome, fixed, or a positive integer")
    parser.add_argument("--export-dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-optimization", action="store_true")
    parser.add_argument("--orientation-sweep", action="store_true")
    args = parser.parse_args()
    report = {
        "scope": "exact linear retracts, full two-cut ranks, and finite coefficient-space costs",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "cases": [],
        "dense_contraction_audit": verify_dense_contraction(),
    }
    for description in args.cases.split(","):
        start = time.monotonic()
        kind, *numbers = description.split(":")
        if kind == "full":
            pairs, length = map(int, numbers)
            labels = all_length_labels(pairs, length)
            name = f"double_endpoint_P{pairs}_L{length}_seed{args.seed}"
        elif kind == "shell":
            modulus, length, weight = map(int, numbers)
            labels = tuple(pair_occupations(parts, modulus)
                           for parts in shell(modulus, length, weight))
            if 2 * weight != length * modulus:
                raise ValueError("the shell must have central energy for complement stability")
            name = f"double_endpoint_M{modulus}_L{length}_h{weight}_seed{args.seed}"
        else:
            raise ValueError("use full:P:L or shell:M:L:h")
        if args.construction == "triangular":
            retract = triangular_retract(labels, not args.ignore_fixed_coordinate, args.groups, args.seed)
            name += "_triangular" + ("_nofix" if args.ignore_fixed_coordinate else "")
            if args.groups != "one":
                name += "_groups" + args.groups
            data, rows, cols = dense_coefficient_blocks(retract)
        elif args.construction == "invariant":
            retract = invariant_orbit_retract(labels, args.seed, not args.ignore_fixed_coordinate)
            name += "_invariant" + ("_nofix" if args.ignore_fixed_coordinate else "")
            data, rows, cols = coefficient_blocks(retract)
        else:
            retract = sparse_retract(labels, args.seed, not args.cross_syndrome)
            data, rows, cols = coefficient_blocks(retract)
        result = dict(data, name=name, retract=retract)
        if args.export_dir:
            args.export_dir.mkdir(parents=True, exist_ok=True)
            (args.export_dir / (name + ".json")).write_text(
                json.dumps(export(name, retract, data, rows, cols), indent=2, sort_keys=True) + "\n")
            gate = math.ceil(data["labels"] ** (99 / 70))
            (args.export_dir / (name + ".txt")).write_text(native_text(data, rows, cols, gate))
        if not args.skip_optimization:
            result["F2_optimization"] = optimize(data, rows, cols)
        if args.orientation_sweep:
            if args.construction != "invariant":
                raise ValueError("the orientation sweep applies to invariant-orbit retracts")
            result["orientation_sweep"] = orientation_sweep(labels)
        result["seconds"] = time.monotonic() - start
        report["cases"].append(result)
        print(json.dumps({
            "name": name, "labels": data["labels"], "span": data["basis_dimension"],
            "largest_block": max(data["block_sizes"]), "seconds": result["seconds"],
            **({key: result["F2_optimization"][key] for key in (
                "found_decomposition_cost", "integer_shortest_lower_bound", "finite_seed_threshold")}
               if "F2_optimization" in result else {}),
        }), flush=True)
        if args.output:
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
