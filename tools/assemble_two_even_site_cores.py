#!/usr/bin/env python3
"""Assemble two-even-site seeds from at most four exceptional sites."""

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import json
import hashlib
from pathlib import Path
import random
import time

from explore_symplectic_boundary_codes import boundary_code, coefficient_basis


ROOT = Path(__file__).resolve().parents[1]


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def add(first, second):
    return [[a ^ b for a, b in zip(left, right)] for left, right in zip(first, second)]


def normalize(data, matrices):
    pairs = data["code"]["pairs"]
    result = matrices[:pairs]
    for i in range(pairs - 2):
        value = matrices[pairs + i]
        if (pairs - i) % 2:
            value = add(value, matrices[i])
        if (i + 1) % 2:
            value = add(value, matrices[i + 1])
        result.append(value)
    return result


def monomials(pairs):
    return [(i, 0) for i in range(pairs)] + [(i, 1) for i in range(pairs - 2)]


def pair_axes(pairs, first, second):
    return [(site, role) for site in range(pairs)
            for role, edge in enumerate((first, second)) if site not in edge]


def pair_words(pairs, first, second, index):
    words = [[0] * pairs for _ in range(2)]
    for role, edge in enumerate((first, second)):
        words[role][max(edge)] = 3
    for bit, (site, role) in enumerate(pair_axes(pairs, first, second)):
        words[role][site] = 1 + ((index >> bit) & 1)
    return tuple(map(tuple, words))


@lru_cache(None)
def core(pairs):
    if pairs == 2:
        return {"monomials": [(0, 0)], "cuts": [[[[1]]], [[[1]]]],
                "labels": {(0, 3): 0}, "positions": {(0, 0): (0, 0)}}
    if pairs == 3:
        code = boundary_code(3, 3)
        data, row, column = coefficient_basis(code)
        data["code"] = code
    else:
        assert pairs == 4
        data = json.loads((ROOT / "evidence/symplectic-radical-sigma2-basis.json").read_text())
        code = data["code"]
        row, column = [[[[int(value, 16) for value in block] for block in matrix]
                        for matrix in data[name]]
                       for name in ("basis_row_blocks", "basis_column_blocks")]
    labels = {tuple(code["ambient_occupations"][output]): i
              for i, output in enumerate(code["output_indices"])}
    grouped = defaultdict(list)
    for i, first in enumerate(code["syndromes"]):
        for j, second in enumerate(code["syndromes"]):
            grouped[first ^ second].append((i, j))
    positions = {pair: (block, i) for block, key in enumerate(sorted(grouped))
                 for i, pair in enumerate(grouped[key])}
    return {"monomials": monomials(pairs), "cuts": [normalize(data, row), normalize(data, column)],
            "labels": labels, "positions": positions}


@lru_cache(None)
def core_subblock(pairs, cut, first, second, third, fourth):
    data = core(pairs)
    size = 1 << (2 * pairs - 4)
    rowpos = [data["positions"][tuple(data["labels"][word] for word in pair_words(
        pairs, first, second, i))] for i in range(size)]
    colpos = [data["positions"][tuple(data["labels"][word] for word in pair_words(
        pairs, third, fourth, i))] for i in range(size)]
    assert len({block for block, _ in rowpos + colpos}) == 1
    block = rowpos[0][0]
    return [[sum(((matrix[block][oldrow] >> oldcol) & 1) << j
                 for j, (_, oldcol) in enumerate(colpos)) for _, oldrow in rowpos]
            for matrix in data["cuts"][cut]]


@lru_cache(None)
def free_rows(cut, copies):
    if copies == 0:
        return [[1]]
    previous = free_rows(cut, copies - 1)
    oldsize = len(previous)
    result = []
    first_order = [1, 4, 2, 8] if cut == 0 else [9, 0, 0, 9]
    for local in range(4):
        for old in range(oldsize):
            coefficients = [0] * (copies + 1)
            for degree, row in enumerate(previous[old]):
                coefficients[degree] ^= row << (local * oldsize)
                for column in bits(first_order[local]):
                    coefficients[degree + 1] ^= row << (column * oldsize)
            result.append(coefficients)
    return result


def kron_row(core_row, free_row, core_size):
    result = 0
    for column in bits(free_row):
        result ^= core_row << (column * core_size)
    return result


@lru_cache(None)
def extended_subblock(pairs, core_pairs, cut, first, second, third, fourth):
    small = core_subblock(core_pairs, cut, first, second, third, fourth)
    free = free_rows(cut, pairs - core_pairs)
    core_size = len(small[0])
    size = core_size * len(free)
    target = {term: i for i, term in enumerate(monomials(pairs))}
    result = [[0] * size for _ in target]
    for index, (degree, charge) in enumerate(core(core_pairs)["monomials"]):
        for core_row, row in enumerate(small[index]):
            if not row:
                continue
            for free_row, coefficients in enumerate(free):
                position = core_row + core_size * free_row
                for extra, other in enumerate(coefficients):
                    if other:
                        result[target[degree + extra, charge]][position] ^= kron_row(row, other, core_size)
    return result


@lru_cache(None)
def axis_permutation(positions):
    count = len(positions)
    size = 1 << count
    layout = list(range(count))
    target = [positions.index(i) for i in range(count)]
    swaps = []
    for i in range(count):
        if layout[i] != target[i]:
            j = layout.index(target[i])
            assert j > i
            delta = (1 << j) - (1 << i)
            mask = sum(1 << value for value in range(size)
                       if value >> i & 1 and not value >> j & 1)
            swaps.append((mask, delta))
            layout[i], layout[j] = layout[j], layout[i]
    return swaps


def permute_row(row, swaps):
    for mask, delta in swaps:
        difference = (row ^ (row >> delta)) & mask
        row ^= difference ^ (difference << delta)
    return row


def union_axis_positions(pairs, first, second, exceptional):
    axes = pair_axes(pairs, first, second)
    return tuple([i for i, (site, _) in enumerate(axes) if site in exceptional] +
                 [i for i, (site, _) in enumerate(axes) if site not in exceptional])


@lru_cache(None)
def source_rows(positions):
    return [sum(((value >> old) & 1) << new for new, old in enumerate(positions))
            for value in range(1 << len(positions))]


def assemble(pairs):
    edges = list(combinations(range(pairs), 2))
    buckets = defaultdict(list)
    for first in edges:
        for second in edges:
            syndrome = sum(1 << i for i in set(first) ^ set(second))
            buckets[syndrome].append((first, second))
    syndromes = sorted(buckets)
    orientations = 1 << (2 * pairs - 4)
    sizes = [len(buckets[key]) * orientations for key in syndromes]
    cuts = []
    for cut in range(2):
        matrices = [[[0] * size for size in sizes] for _ in monomials(pairs)]
        for block, syndrome in enumerate(syndromes):
            edgepairs = buckets[syndrome]
            for ri, (first, second) in enumerate(edgepairs):
                for ci, (third, fourth) in enumerate(edgepairs):
                    exceptional = sorted(set(first + second + third + fourth))
                    assert 2 <= len(exceptional) <= 4
                    rel = {site: i for i, site in enumerate(exceptional)}
                    local_edges = tuple(tuple(rel[i] for i in edge)
                                        for edge in (first, second, third, fourth))
                    little = extended_subblock(pairs, len(exceptional), cut, *local_edges)
                    rowmap = source_rows(union_axis_positions(pairs, first, second, exceptional))
                    permutation = axis_permutation(union_axis_positions(pairs, third, fourth, exceptional))
                    row_offset, col_offset = ri * orientations, ci * orientations
                    for i, matrix in enumerate(little):
                        target = matrices[i][block]
                        for row, old in enumerate(rowmap):
                            value = matrix[old]
                            if value:
                                target[row_offset + row] ^= permute_row(value, permutation) << col_offset
            print(json.dumps({"cut": cut, "syndrome": syndrome, "size": sizes[block]}), flush=True)
        assert all(rows == [1 << i for i in range(len(rows))] for rows in matrices[0])
        cuts.append(matrices)
    return {"pairs": pairs, "labels": len(edges) * (1 << (pairs - 2)),
            "basis_dimension": len(monomials(pairs)), "normalized_monomials_t_z": monomials(pairs),
            "syndromes": syndromes, "block_sizes": sizes, "edgepairs": [buckets[s] for s in syndromes]}, cuts


def compare_literal(metadata, cuts, path):
    data = json.loads(path.read_text())
    code = data["code"]
    assert metadata["labels"] == code["labels"]
    labels = {tuple(code["ambient_occupations"][output]): i
              for i, output in enumerate(code["output_indices"])}
    groups = defaultdict(list)
    for i, first in enumerate(code["syndromes"]):
        for j, second in enumerate(code["syndromes"]):
            groups[first ^ second].append((i, j))
    positions = [{pair: i for i, pair in enumerate(groups[s])} for s in metadata["syndromes"]]
    maps = []
    orientations = 1 << (2 * metadata["pairs"] - 4)
    for block, edgepairs in enumerate(metadata["edgepairs"]):
        maps.append([positions[block][tuple(labels[word] for word in pair_words(
            metadata["pairs"], first, second, i))]
                     for first, second in edgepairs for i in range(orientations)])
    count = 0
    for cut, name in enumerate(("basis_row_blocks", "basis_column_blocks")):
        literal = normalize(data, [[[int(value, 16) for value in block] for block in matrix]
                                   for matrix in data[name]])
        for expected, actual in zip(literal, cuts[cut]):
            for old, new, indices in zip(expected, actual, maps):
                for row, source in zip(new, indices):
                    assert row == sum(((old[source] >> column) & 1) << j
                                      for j, column in enumerate(indices))
                    count += 1
    return count


def write_native(path, metadata, cuts, gate):
    with path.open("w") as stream:
        stream.write(f"{metadata['labels']} {metadata['basis_dimension']} {gate}\n")
        for matrices in cuts:
            stream.write(f"{len(metadata['block_sizes'])}\n")
            for block, size in enumerate(metadata["block_sizes"]):
                stream.write(f"{size} 1\n")
                for matrix in matrices:
                    for row in matrix[block]:
                        stream.write(" ".join(f"{(row >> shift) & ((1 << 64) - 1):x}"
                                              for shift in range(0, size, 64)) + "\n")


def independent_coordinates(matrices):
    pivots, witnesses = {}, []
    dimension = len(matrices)
    for block in range(len(matrices[0])):
        for row in range(len(matrices[0][block])):
            values = [matrix[block][row] for matrix in matrices]
            union = 0
            for value in values:
                union |= value
            for column in bits(union):
                original = sum(((value >> column) & 1) << i for i, value in enumerate(values))
                vector = original
                while vector:
                    pivot = vector.bit_length() - 1
                    if pivot not in pivots:
                        pivots[pivot] = vector
                        witnesses.append([block, row, column, original])
                        break
                    vector ^= pivots[pivot]
                if len(pivots) == dimension:
                    return witnesses
    raise AssertionError("the proposed generators are not independent")


def literal_controls(metadata, cuts, trials=128):
    from complement_shell_tensor import multiply
    from explore_complement_edge_pencils import evaluate_local
    pairs = metadata["pairs"]
    variable_p, variable_q = 2, 1 << (pairs + 1)
    local = evaluate_local((variable_p, variable_q, 1, 1, 1, 1))

    def power(value, exponent):
        result = 1
        for _ in range(exponent):
            result = multiply(result, value)
        return result

    weights = [multiply(multiply(power(1 ^ variable_q, pairs - degree - 2 * charge),
                                power(variable_p ^ variable_q, degree)),
                        power(variable_q, charge))
               for degree, charge in metadata["normalized_monomials_t_z"]]
    orientations = 1 << (2 * pairs - 4)

    def input_words(output, edge):
        result = []
        for full in edge:
            word = list(output)
            for site in edge:
                word[site] = 3 * int(site == full)
            result.append(word)
        return result

    rng = random.Random(771503)
    reports = []
    for cut in range(2):
        checked, union_sizes = 0, defaultdict(int)
        samples = []
        for i, edgepair in enumerate(metadata["edgepairs"][0]):
            for j in range(min(3, len(metadata["edgepairs"][0]))):
                samples.append((0, i * orientations + rng.randrange(orientations),
                                j * orientations + rng.randrange(orientations)))
        samples += [(block, rng.randrange(metadata["block_sizes"][block]),
                     rng.randrange(metadata["block_sizes"][block]))
                    for block in [rng.randrange(len(metadata["block_sizes"])) for _ in range(trials)]]
        for block, row, column in samples:
            first, second = metadata["edgepairs"][block][row // orientations]
            third, fourth = metadata["edgepairs"][block][column // orientations]
            rword = pair_words(pairs, first, second, row % orientations)
            cword = pair_words(pairs, third, fourth, column % orientations)
            if cut == 0:
                x, y, v, u = rword[0], rword[1], cword[0], cword[1]
                yedge, vedge = second, third
            else:
                x, u, v, y = rword[0], rword[1], cword[0], cword[1]
                yedge, vedge = fourth, third
            actual = 0
            for yy in input_words(y, yedge):
                for vv in input_words(v, vedge):
                    value = 1
                    for a, b, c, d in zip(x, yy, u, vv):
                        value = multiply(value, local[((a * 4 + b) * 4 + c) * 4 + d])
                    actual ^= value
            expected = 0
            for matrix, weight in zip(cuts[cut], weights):
                if matrix[block][row] >> column & 1:
                    expected ^= weight
            assert actual == expected
            union_sizes[len(set(first + second + third + fourth))] += 1
            checked += 1
        assert set(union_sizes) == {2, 3, 4}
        reports.append({"cut": cut, "complete_polynomial_entries_checked": checked,
                        "exceptional_union_size_counts": dict(union_sizes)})
    return reports


def audit_binary_screen(metadata, cuts, path, samples=12):
    from check_modular_crossing import rank_two
    events = [json.loads(line) for line in path.read_text().splitlines()]
    first, last = events[0], events[-1]
    dimension, labels = metadata["basis_dimension"], metadata["labels"]
    assert first["event"] == "start" and last["event"] == "complete"
    assert first["field_size"] == 2 and first["labels"] == labels
    assert first["dimension"] == dimension and first["code_start"] == 1
    assert first["code_count"] == last["projective_tested"] == (1 << dimension) - 1
    gate = first["integer_cost_gate"]
    assert gate ** 70 > labels ** 99 and 99 ** 2 > 2 * 70 ** 2
    eligible = {entry["index"]: entry for entry in events if entry["event"] == "eligible"}
    assert len(eligible) == last["eligible"]
    selected = set(eligible)
    rng = random.Random(775098)
    while len(selected) < min((1 << dimension) - 1, len(eligible) + samples):
        selected.add(rng.randrange(1, 1 << dimension))
    checks = []
    for mask in sorted(selected):
        ranks = []
        for matrices in cuts:
            total = 0
            for block, size in enumerate(metadata["block_sizes"]):
                rows = [0] * size
                for index in bits(mask):
                    rows = [a ^ b for a, b in zip(rows, matrices[index][block])]
                total += rank_two(rows)
            ranks.append(total)
        affordable = ranks[0] * ranks[1] < gate * gate
        assert affordable == (mask in eligible)
        if affordable:
            assert ranks == [eligible[mask]["row_rank"], eligible[mask]["column_rank"]]
        checks.append({"mask": mask, "ranks": ranks, "affordable": affordable})
    pivots = {}
    for mask in eligible:
        while mask:
            pivot = mask.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = mask
                break
            mask ^= pivots[pivot]
    remainder = 1
    while remainder and remainder.bit_length() - 1 in pivots:
        remainder ^= pivots[remainder.bit_length() - 1]
    return {"field_size": 2, "vectors_tested": last["projective_tested"], "integer_cost_gate": gate,
            "eligible_vectors": list(eligible.values()), "eligible_span_dimension": len(pivots),
            "endpoint_remainder": remainder, "excluded_by_affordable_span": bool(remainder),
            "independent_rank_controls": checks,
            "screen_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "native_source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                      for name in ("screen_weighted_pfaffian_span.cpp", "screen_pfaffian_span_f4.cpp")}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--compare-basis", type=Path)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--screen", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    began = time.monotonic()
    metadata, cuts = assemble(args.pairs)
    if args.compare_basis:
        metadata["literal_normalized_bit_rows_compared"] = compare_literal(metadata, cuts, args.compare_basis)
    metadata["independent_coordinate_witnesses"] = independent_coordinates(cuts[0])
    metadata["literal_polynomial_controls"] = literal_controls(metadata, cuts)
    if args.screen:
        metadata["binary_screen_audit"] = audit_binary_screen(metadata, cuts, args.screen)
    if args.native:
        from check_symplectic_seed_tensor_fields import integer_gate
        write_native(args.native, metadata, cuts, integer_gate(metadata["labels"]))
        with args.native.open("rb") as stream:
            metadata["native_input_sha256"] = hashlib.file_digest(stream, "sha256").hexdigest()
    metadata.pop("edgepairs")
    metadata["elapsed_seconds"] = time.monotonic() - began
    metadata["source_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                  for name in (Path(__file__).name, "explore_symplectic_boundary_codes.py",
                                               "explore_complement_edge_pencils.py")}
    path = ROOT / "evidence/symplectic-radical-sigma2-basis.json"
    metadata["four_site_basis_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    encoded = json.dumps(metadata, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
