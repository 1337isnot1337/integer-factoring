#!/usr/bin/env python3
"""Exact Pfaffian profiles on full-grade symplectic boundary codes."""

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import subprocess
import time

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import all_length_labels, evaluate_local, independent_coefficients, syndrome
from check_double_endpoint_retracts import bit_indices, native_text, optimize
from analyze_pfaffian_binary_screen import analyze, unpack_binary


def reduced_basis(vectors):
    pivots = {}
    for value in vectors:
        for pivot in sorted(pivots, reverse=True):
            if value >> pivot & 1:
                value ^= pivots[pivot]
        if not value:
            continue
        pivot = value.bit_length() - 1
        for old in pivots:
            if pivots[old] >> pivot & 1:
                pivots[old] ^= value
        pivots[pivot] = value
    return [(pivot, pivots[pivot]) for pivot in sorted(pivots, reverse=True)]


def nullspace(rows, width):
    reduced = reduced_basis(rows)
    pivots = {pivot for pivot, row in reduced}
    result = []
    for free in range(width):
        if free in pivots:
            continue
        value = 1 << free
        for pivot, row in reduced:
            if row >> free & 1:
                value ^= 1 << pivot
        assert all((row & value).bit_count() % 2 == 0 for row in rows)
        result.append(value)
    return result


def predicted_dimensions(pairs, length):
    def choose(n, k):
        return math.comb(n, k) if 0 <= k <= n else 0
    boundary = gram = 0
    for singles in range(min(pairs, length) + 1):
        if (length - singles) & 1:
            continue
        even = pairs - singles
        if not even:
            continue
        degree = (length - singles) // 2
        copies = math.comb(pairs, singles) * (1 << singles)
        dimension = choose(even - 1, degree - 1)
        boundary += copies * dimension
        gram += copies * (dimension if even & 1 else choose(even - 2, degree - 2))
    return {"boundary_dimension": boundary, "boundary_Gram_rank": gram,
            "radical_dimension": boundary - gram}


def boundary_code(pairs, length, radical=True):
    labels = all_length_labels(pairs, length)
    physical_masks = [sum(mask << (2 * i) for i, mask in enumerate(word)) for word in labels]
    positions = {mask: i for i, mask in enumerate(physical_masks)}
    groups = defaultdict(list)
    for subset in combinations(range(2 * pairs), length - 2):
        source = sum(1 << i for i in subset)
        value = 0
        for i in range(pairs):
            pair = 3 << (2 * i)
            if not source & pair:
                value ^= 1 << positions[source | pair]
        if value:
            first = (value & -value).bit_length() - 1
            groups[syndrome(labels[first])].append(value)
    inputs, outputs, sigmas = [], [], []
    original_dimension = gram_rank = 0
    for sigma in sorted(groups):
        basis = reduced_basis(groups[sigma])
        vectors = [value for _, value in basis]
        original_dimension += len(vectors)
        gram = [sum(((first & second).bit_count() % 2) << j
                    for j, second in enumerate(vectors)) for first in vectors]
        gram_rank += rank_two(gram)
        if radical:
            kernel = nullspace(gram, len(vectors))
            new = []
            for coefficients in kernel:
                value = 0
                for index in bit_indices(coefficients):
                    value ^= vectors[index]
                new.append(value)
            basis = reduced_basis(new)
        for pivot, value in basis:
            outputs.append(pivot)
            inputs.append(value)
            sigmas.append(sigma)
    d = len(inputs)
    lookup = {word: i for i, word in enumerate(labels)}
    J = [lookup[tuple((0, 2, 1, 3)[mask] for mask in word)] for word in labels]
    induced = [0] * d
    gram = [0] * d
    for y, value in enumerate(inputs):
        shifted = sum(1 << J[i] for i in bit_indices(value))
        reconstructed = 0
        for x, pivot in enumerate(outputs):
            assert ((value >> pivot) & 1) == int(x == y)
            if shifted >> pivot & 1:
                induced[x] ^= 1 << y
                reconstructed ^= inputs[x]
            gram[x] |= ((inputs[x] & value).bit_count() % 2) << y
        assert reconstructed == shifted
    assert rank_two(induced) == d
    if radical:
        assert not any(gram)
    predicted = predicted_dimensions(pairs, length)
    assert original_dimension == predicted["boundary_dimension"]
    assert gram_rank == predicted["boundary_Gram_rank"]
    assert d == predicted["radical_dimension" if radical else "boundary_dimension"]
    return {
        "pairs": pairs, "length": length, "labels": d,
        "ambient_occupations": labels, "input_vectors": inputs,
        "output_indices": outputs, "syndromes": sigmas,
        "induced_mode_complement_rows": induced, "input_Gram_rows": gram,
        "original_boundary_dimension": original_dimension,
        "original_boundary_Gram_rank": gram_rank,
        "radical": radical,
        "predicted_dimensions": predicted,
        "maximum_input_terms": max(map(int.bit_count, inputs), default=0),
    }


def scalar_coefficient_polynomials(rows, active, selected):
    """Recover scalar coefficients from an invertible entry minor, then audit all entries."""
    basis = [rows[i] for i in selected]
    width = len(basis)
    samples = []
    pivots = {}
    for block in range(len(basis[0])):
        for ri in range(len(basis[0][block])):
            union = 0
            for matrix in basis:
                union |= matrix[block][ri]
            for ci in bit_indices(union):
                signature = sum(((matrix[block][ri] >> ci) & 1) << i
                                for i, matrix in enumerate(basis))
                reduced = signature
                while reduced:
                    pivot = reduced.bit_length() - 1
                    if pivot not in pivots:
                        pivots[pivot] = reduced
                        value = sum(((matrix[block][ri] >> ci) & 1) << degree
                                    for degree, matrix in zip(active, rows))
                        samples.append([signature, value])
                        break
                    reduced ^= pivots[pivot]
                if len(samples) == width:
                    break
            if len(samples) == width:
                break
        if len(samples) == width:
            break
    assert len(samples) == width
    for col in range(width):
        pivot = next(i for i in range(col, width) if samples[i][0] >> col & 1)
        samples[col], samples[pivot] = samples[pivot], samples[col]
        for i in range(width):
            if i != col and samples[i][0] >> col & 1:
                samples[i][0] ^= samples[col][0]
                samples[i][1] ^= samples[col][1]
    polynomials = [value for _, value in samples]
    coordinates = []
    for degree, matrix in zip(active, rows):
        mask = sum(((value >> degree) & 1) << i for i, value in enumerate(polynomials))
        coordinates.append(mask)
        for block, expected in enumerate(matrix):
            for ri, value in enumerate(expected):
                actual = 0
                for i in bit_indices(mask):
                    actual ^= basis[i][block][ri]
                assert actual == value
    return polynomials, coordinates


def common_linear_factor_power(polynomials, length, variable="q"):
    base = length + 1
    answer = length
    for value in polynomials:
        slices = defaultdict(int)
        for degree in bit_indices(value):
            a, b = degree % base, degree // base
            fixed, moving = (a, b) if variable == "q" else (b, a)
            slices[fixed] ^= 1 << moving
        for polynomial in slices.values():
            shifted = 0
            for exponent in bit_indices(polynomial):
                subset = exponent
                while True:
                    shifted ^= 1 << subset
                    if not subset:
                        break
                    subset = (subset - 1) & exponent
            if shifted:
                answer = min(answer, (shifted & -shifted).bit_length() - 1)
    return answer


def coefficient_basis(code):
    labels, inputs, outputs = code["ambient_occupations"], code["input_vectors"], code["output_indices"]
    d, length = code["labels"], code["length"]
    if d < 2:
        raise ValueError("at least two labels are needed")
    supports = [bit_indices(value) for value in inputs]
    base = length + 1
    local = evaluate_local((2, 1 << base, 1, 1, 1, 1))
    groups = defaultdict(list)
    for x, y in product(range(d), repeat=2):
        groups[code["syndromes"][x] ^ code["syndromes"][y]].append((x, y))
    keys = sorted(groups)
    positions = {pair: (block, index) for block, key in enumerate(keys)
                 for index, pair in enumerate(groups[key])}
    size = length * base + 1
    rows = [[[0] * len(groups[key]) for key in keys] for _ in range(size)]
    columns = [[[0] * len(groups[key]) for key in keys] for _ in range(size)]

    @lru_cache(None)
    def pattern_value(pattern):
        value = 1
        for index in pattern:
            value = multiply(value, local[index])
            if not value:
                break
        return value

    evaluated = literal = 0
    for block, key in enumerate(keys):
        group = groups[key]
        for row, (x, y) in enumerate(group):
            for col, (v, u) in enumerate(group):
                value = 0
                for yy, vv in product(supports[y], supports[v]):
                    pattern = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + e
                                           for a, b, c, e in zip(labels[outputs[x]], labels[yy],
                                                                 labels[outputs[u]], labels[vv])))
                    value ^= pattern_value(pattern)
                    literal += 1
                assert (value & 1) == int(x == v and y == u)
                expected_p = int(x == u) * ((code["input_Gram_rows"][y] >> v) & 1)
                expected_q = ((code["induced_mode_complement_rows"][x] >> y) & 1) * (
                    (code["induced_mode_complement_rows"][u] >> v) & 1)
                assert ((value >> length) & 1) == expected_p
                assert ((value >> (length * base)) & 1) == expected_q
                cblock, crow = positions[x, u]
                other, ccol = positions[v, y]
                assert cblock == other
                while value:
                    bit = value & -value
                    degree = bit.bit_length() - 1
                    rows[degree][block][row] ^= 1 << col
                    columns[degree][cblock][crow] ^= 1 << ccol
                    value ^= bit
                evaluated += 1
    active = [i for i, matrix in enumerate(rows) if any(any(block) for block in matrix)]
    rows = [rows[i] for i in active]
    columns = [columns[i] for i in active]
    selected = independent_coefficients(rows)
    assert active[selected[0]] == 0
    scalar_polynomials, coordinates = scalar_coefficient_polynomials(rows, active, selected)
    profile = [[degree % base, degree // base, length - degree % base - degree // base,
                sum(rank_two(block) for block in row),
                sum(rank_two(block) for block in col)]
               for degree, row, col in zip(active, rows, columns)]
    data = {
        "labels": d, "length": length, "basis_dimension": len(selected),
        "block_sizes": [len(groups[key]) for key in keys],
        "basis_monomials_p_q_r": [profile[i][:3] for i in selected],
        "basis_scalar_polynomials_p_q_r": [
            [[degree % base, degree // base, length - degree % base - degree // base]
             for degree in bit_indices(value)] for value in scalar_polynomials],
        "universal_coefficient_coordinates": coordinates,
        "common_q_plus_r_power": common_linear_factor_power(scalar_polynomials, length, "q"),
        "common_p_plus_r_power": common_linear_factor_power(scalar_polynomials, length, "p"),
        "universal_coefficient_profile": profile,
        "endpoint_coordinates": [1] + [0] * (len(selected) - 1),
        "contracted_entries_checked": evaluated, "literal_entries_summed": literal,
        "distinct_local_patterns": pattern_value.cache_info().currsize,
        "p_endpoint_zero": not any(code["input_Gram_rows"]),
        "q_endpoint_ranks": [int(any(code["induced_mode_complement_rows"])),
                             rank_two(code["induced_mode_complement_rows"]) ** 2],
    }
    return data, [rows[i] for i in selected], [columns[i] for i in selected]


def half_orientation_code(code):
    assert code["radical"]
    labels = code["ambient_occupations"]
    selected = []
    for i, (value, sigma) in enumerate(zip(code["input_vectors"], code["syndromes"])):
        if not sigma:
            continue
        first_site = (sigma & -sigma).bit_length() - 1
        support = bit_indices(value)
        orientations = {labels[j][first_site] for j in support}
        assert len(orientations) == 1 and orientations <= {1, 2}
        if orientations == {1}:
            selected.append(i)
    result = dict(code)
    result["labels"] = len(selected)
    result["input_vectors"] = [code["input_vectors"][i] for i in selected]
    result["output_indices"] = [code["output_indices"][i] for i in selected]
    result["syndromes"] = [code["syndromes"][i] for i in selected]
    result["input_Gram_rows"] = [0] * len(selected)
    result["induced_mode_complement_rows"] = [0] * len(selected)
    result["half_orientation_selector"] = "least single-occupied pair uses its low mode; all zero-syndrome vectors are omitted"
    lookup = {tuple(word): i for i, word in enumerate(labels)}
    for y, value in enumerate(result["input_vectors"]):
        shifted = sum(1 << lookup[tuple((0, 2, 1, 3)[mask] for mask in labels[i])]
                      for i in bit_indices(value))
        for x, pivot in enumerate(result["output_indices"]):
            assert ((value >> pivot) & 1) == int(x == y)
            assert not ((shifted >> pivot) & 1)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="4:4")
    parser.add_argument("--full-boundary", action="store_true")
    parser.add_argument("--half-orientations", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--export-dir", type=Path, default=Path("/tmp/symplectic-boundary-bases"))
    parser.add_argument("--native", type=Path)
    parser.add_argument("--optimize", action="store_true")
    args = parser.parse_args()
    records = []
    args.export_dir.mkdir(parents=True, exist_ok=True)
    for spec in args.cases.split(","):
        started = time.monotonic()
        pairs, length = map(int, spec.split(":"))
        code = boundary_code(pairs, length, not args.full_boundary)
        if args.half_orientations:
            code = half_orientation_code(code)
        data, rows, columns = coefficient_basis(code)
        name = f'omega_{"boundary" if args.full_boundary else "radical"}_P{pairs}_L{length}'
        if args.half_orientations:
            name += "_half"
        exported = dict(data, name=name, code=code,
                        basis_row_blocks=[[[hex(v) for v in block] for block in matrix] for matrix in rows],
                        basis_column_blocks=[[[hex(v) for v in block] for block in matrix] for matrix in columns],
                        exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        basis_path = args.export_dir / (name + ".json")
        basis_path.write_text(json.dumps(exported, indent=2, sort_keys=True) + "\n")
        gate = math.ceil(code["labels"] ** (99 / 70))
        native = native_text(data, rows, columns, gate)
        native_path = args.export_dir / (name + ".txt")
        native_path.write_text(native)
        result = dict(data, name=name, code=code)
        if args.optimize:
            result["F2_optimization"] = optimize(data, rows, columns)
        if args.native:
            completed = subprocess.run([str(args.native), "--field2"], input=native, text=True,
                                       stdout=subprocess.PIPE, check=True)
            screen = args.export_dir / (name + ".jsonl")
            screen.write_text(completed.stdout)
            result["F2_screen"] = analyze(basis_path, screen)
            events = [json.loads(line) for line in completed.stdout.splitlines()]
            audited = 0
            for event in events:
                if event["event"] != "eligible":
                    continue
                mask = unpack_binary(event["code"], len(rows))
                ranks = []
                for basis in (rows, columns):
                    work = [[0] * size for size in data["block_sizes"]]
                    for i, matrix in enumerate(basis):
                        if mask >> i & 1:
                            for target, source in zip(work, matrix):
                                for j, value in enumerate(source):
                                    target[j] ^= value
                    ranks.append(sum(rank_two(block) for block in work))
                assert ranks == [event["row_rank"], event["column_rank"]]
                audited += 1
            result["eligible_vectors_independently_audited"] = audited
        result["seconds"] = time.monotonic() - started
        records.append(result)
        print(json.dumps({"name": name, "labels": data["labels"], "span": data["basis_dimension"],
                          "max_input_terms": code["maximum_input_terms"],
                          "largest_block": max(data["block_sizes"]),
                          "seconds": result["seconds"],
                          "screen": result.get("F2_screen"),
                          "optimum": result.get("F2_optimization", {}).get("found_decomposition_cost")}), flush=True)
        if args.output:
            args.output.write_text(json.dumps({"cases": records,
                                               "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                                              indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

