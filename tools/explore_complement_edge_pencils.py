#!/usr/bin/env python3
"""Search exact endpoint-preserving Pfaffian pencils on fixed-length codes.

The six edge polynomials are UA, AY, AV, UV, VY, UY. At t=1 the first
two vanish and AV=UY=1. Equal global A/V lengths then force the exact
crossing, even without an energy-shell restriction.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import heapq
import json
import math
from pathlib import Path
import time

from check_complement_shell_crossings import shell
from check_modular_crossing import rank_two
from complement_shell_tensor import LOCAL_TENSOR, multiply, pair_occupations


@lru_cache(None)
def local_monomials(x, y, u, v):
    vertices = ([(0, bit) for bit in range(2) if u >> bit & 1]
                + [(1, bit) for bit in range(2) if x >> bit & 1]
                + [(2, bit) for bit in range(2) if v >> bit & 1]
                + [(3, bit) for bit in range(2) if y >> bit & 1])
    n = len(vertices)
    if n & 1:
        return frozenset()
    edges = [[None] * n for _ in range(n)]
    for i, (first, a) in enumerate(vertices):
        for j in range(i + 1, n):
            second, b = vertices[j]
            if (first, second) == (0, 1) and a == b:
                edges[i][j] = 0
            elif (first, second) == (1, 3) and a != b:
                edges[i][j] = 1
            elif (first, second) == (1, 2) and a == b:
                edges[i][j] = 2
            elif (first, second) == (0, 2) and a != b:
                edges[i][j] = 3
            elif (first, second) == (2, 3) and a == b:
                edges[i][j] = 4
            elif (first, second) == (0, 3) and a == b:
                edges[i][j] = 5

    @lru_cache(None)
    def recurse(mask):
        if not mask:
            return frozenset({(0,) * 6})
        first_bit = mask & -mask
        first = first_bit.bit_length() - 1
        rest = mask ^ first_bit
        result = set()
        candidates = rest
        while candidates:
            second_bit = candidates & -candidates
            second = second_bit.bit_length() - 1
            edge = edges[first][second]
            if edge is not None:
                for monomial in recurse(rest ^ second_bit):
                    term = list(monomial)
                    term[edge] += 1
                    term = tuple(term)
                    if term in result:
                        result.remove(term)
                    else:
                        result.add(term)
            candidates ^= second_bit
        return frozenset(result)
    return recurse((1 << n) - 1)


MONOMIALS = tuple(local_monomials(x, y, u, v)
                  for x, y, u, v in product(range(4), repeat=4))


def evaluate_local(edges):
    powers = []
    for edge in edges:
        row = [1]
        for _ in range(4):
            row.append(multiply(row[-1], edge))
        powers.append(row)
    values = []
    for terms in MONOMIALS:
        value = 0
        for term in terms:
            product_ = 1
            for i, exponent in enumerate(term):
                product_ = multiply(product_, powers[i][exponent])
            value ^= product_
        values.append(value)
    return values


def verify_local_family():
    assert tuple(evaluate_local((6, 3, 2, 1, 1, 1))) == LOCAL_TENSOR
    return {"original_local_entries_compared": len(LOCAL_TENSOR)}


def syndrome(label):
    return sum((mask.bit_count() & 1) << i for i, mask in enumerate(label))


def all_length_labels(pairs, length):
    return tuple(tuple(sum(1 << (k % 2) for k in subset if k // 2 == i)
                       for i in range(pairs))
                 for subset in combinations(range(2 * pairs), length))


def compile_code(parameters, states=None):
    parameters = tuple(parameters)
    states = parameters if states is None else tuple(states)
    d = len(parameters)
    assert d == len(states) and d > 0
    length = sum(mask.bit_count() for mask in parameters[0])
    assert all(sum(mask.bit_count() for mask in word) == length for word in parameters + states)
    assert len(set(parameters)) == d and len(set(states)) == d
    assert len({len(word) for word in parameters + states}) == 1
    blocks = defaultdict(list)
    for x, y in product(range(d), repeat=2):
        blocks[syndrome(parameters[x]) ^ syndrome(states[y])].append((x, y))
    keys = sorted(blocks)
    positions = {pair: (block, i) for block, key in enumerate(keys)
                 for i, pair in enumerate(blocks[key])}
    entries_by_type = defaultdict(list)
    diagonal_entries = 0
    for block, key in enumerate(keys):
        pairs = blocks[key]
        for row, (x, y) in enumerate(pairs):
            for col, (v, u) in enumerate(pairs):
                pattern = tuple(sorted(
                    ((a * 4 + b) * 4 + c) * 4 + e
                    for a, b, c, e in zip(parameters[x], states[y], states[u], parameters[v])
                    if a or b or c or e))
                other_block, other_row = positions[x, u]
                other_block_check, other_col = positions[v, y]
                assert other_block == other_block_check
                crossing = int(x == v and y == u)
                diagonal_entries += crossing
                entries_by_type[pattern].append(
                    (block, row, col, other_block, other_row, other_col, crossing))
    assert diagonal_entries == d * d
    return {
        "labels": d, "length": length, "pairs": len(parameters[0]),
        "block_sizes": [len(blocks[key]) for key in keys],
        "types": list(entries_by_type.items()),
        "evaluated_entries": sum(len(entries) for entries in entries_by_type.values()),
        "distinct_local_types": len(entries_by_type),
    }


def coefficient_blocks(code, edges, constant_crossing=False):
    if not constant_crossing:
        assert not edges[0].bit_count() % 2 and not edges[1].bit_count() % 2
        assert edges[2].bit_count() % 2 == edges[5].bit_count() % 2 == 1
    local = evaluate_local(edges)
    values = []
    maximum_degree = 0
    for pattern, entries in code["types"]:
        value = 1
        for local_index in pattern:
            value = multiply(value, local[local_index])
            if not value:
                break
        specialization = (value & 1) if constant_crossing else value.bit_count() % 2
        assert all(specialization == entry[-1] for entry in entries)
        maximum_degree = max(maximum_degree, value.bit_length() - 1)
        values.append(value)
    sizes = code["block_sizes"]
    rows = [[[0] * size for size in sizes] for _ in range(maximum_degree + 1)]
    cols = [[[0] * size for size in sizes] for _ in range(maximum_degree + 1)]
    for (_, entries), value in zip(code["types"], values):
        while value:
            bit = value & -value
            degree = bit.bit_length() - 1
            for block, row, col, cblock, crow, ccol, _ in entries:
                rows[degree][block][row] ^= 1 << col
                cols[degree][cblock][crow] ^= 1 << ccol
            value ^= bit
    active = [degree for degree in range(maximum_degree + 1)
              if any(any(block) for block in rows[degree])]
    return active, [rows[degree] for degree in active], [cols[degree] for degree in active]


def rank_profile(code, edges, return_blocks=False):
    degrees, rows, cols = coefficient_blocks(code, edges)
    profile = [[degree, sum(rank_two(block) for block in row),
                sum(rank_two(block) for block in col)]
               for degree, row, col in zip(degrees, rows, cols)]
    score = sum(math.sqrt(r * c) for _, r, c in profile)
    result = {
        "edges_UA_AY_AV_UV_VY_UY": list(edges), "profile": profile,
        "rank_sum": score, "threshold": code["labels"] ** math.sqrt(2),
        "margin": code["labels"] ** math.sqrt(2) / score,
    }
    return (result, (degrees, rows, cols)) if return_blocks else result


def minimum_decomposition(costs, target=None):
    size = len(costs)
    target = size - 1 if target is None else target
    distances = [float("inf")] * size
    distances[0] = 0
    parents = [None] * size
    upper = costs[target]
    upper_masks = [target]
    for mask in range(1, size):
        if costs[mask] + costs[mask ^ target] < upper:
            upper = costs[mask] + costs[mask ^ target]
            upper_masks = [mask, mask ^ target]
    queue = [(0, 0)]
    edges = sorted((costs[mask], mask) for mask in range(1, size))
    while queue:
        distance, state = heapq.heappop(queue)
        if distance != distances[state]:
            continue
        for cost, mask in edges:
            candidate = distance + cost
            if candidate >= upper - 1e-10:
                break
            destination = state ^ mask
            if candidate < distances[destination] - 1e-10:
                distances[destination] = candidate
                parents[destination] = state, mask
                if destination == target:
                    upper = candidate
                    upper_masks = None
                heapq.heappush(queue, (candidate, destination))
    if upper_masks is not None:
        return upper, [mask for mask in upper_masks if mask]
    masks = []
    state = target
    while state:
        state, mask = parents[state]
        masks.append(mask)
    return upper, masks


def independent_coefficients(rows):
    pivots = {}
    selected = []
    for degree, blocks in enumerate(rows):
        vector = 0
        offset = 0
        for block in blocks:
            for row in block:
                vector |= row << offset
                offset += len(block)
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                selected.append(degree)
                break
            vector ^= pivots[pivot]
    return selected


def universal_span(code):
    """Optimize the full F_2 span of the homogeneous ternary tensor."""
    started = time.monotonic()
    base = code["length"] + 1
    degrees, rows, cols = coefficient_blocks(
        code, (2, 1 << base, 1, 1, 1, 1), constant_crossing=True)
    assert degrees[0] == 0
    for degree in degrees:
        assert degree % base + degree // base <= code["length"]
    selected = independent_coefficients(rows)
    assert selected[0] == 0
    basis_rows = [rows[i] for i in selected]
    basis_cols = [cols[i] for i in selected]
    width = len(selected)
    if width > 18:
        raise ValueError("the full universal span exceeds the exhaustive binary bound")
    rwork = [[0] * size for size in code["block_sizes"]]
    cwork = [[0] * size for size in code["block_sizes"]]
    ranks = [[0, 0] for _ in range(1 << width)]
    previous = 0
    for index in range(1, 1 << width):
        gray = index ^ (index >> 1)
        changed = (gray ^ previous).bit_length() - 1
        for target, source in ((rwork, basis_rows[changed]), (cwork, basis_cols[changed])):
            for block, delta in zip(target, source):
                for row, value in enumerate(delta):
                    block[row] ^= value
        ranks[gray] = [sum(rank_two(block) for block in rwork),
                       sum(rank_two(block) for block in cwork)]
        previous = gray
    cost, masks = minimum_decomposition([math.sqrt(r * c) for r, c in ranks], target=1)
    lower, _ = minimum_decomposition([math.isqrt(r * c) for r, c in ranks], target=1)
    excluded = lower ** 70 > code["labels"] ** 99
    return {
        "nonzero_universal_coefficients": [[degree % base, degree // base,
                                            code["length"] - degree % base - degree // base]
                                           for degree in degrees],
        "basis_coefficient_indices": selected, "span_dimension": width,
        "combination_ranks": ranks, "optimal_cost": cost,
        "optimal_combination_masks": masks, "integer_lower_bound": lower,
        "threshold": code["labels"] ** math.sqrt(2),
        "margin": code["labels"] ** math.sqrt(2) / cost,
        "all_F2_decompositions_fail_seed_criterion": excluded,
        "integer_threshold_certificate": {
            "sqrt_two_upper_numerator": 99, "sqrt_two_upper_denominator": 70,
            "lower_bound_power_70_exceeds_label_power_99": excluded,
        },
        "seconds": time.monotonic() - started,
    }


def export_basis(name, labels):
    """Canonical independent binary generators for extension-field screens."""
    code = compile_code(labels)
    base = code["length"] + 1
    degrees, rows, cols = coefficient_blocks(
        code, (2, 1 << base, 1, 1, 1, 1), constant_crossing=True)
    selected = independent_coefficients(rows)
    assert selected[0] == 0
    return {
        "name": name, "labels": len(labels), "length": code["length"],
        "occupations": labels, "block_sizes": code["block_sizes"],
        "basis_dimension": len(selected),
        "basis_monomials_p_q_r": [
            [degrees[i] % base, degrees[i] // base,
             code["length"] - degrees[i] % base - degrees[i] // base]
            for i in selected],
        "endpoint_coordinates": [1] + [0] * (len(selected) - 1),
        "basis_row_blocks": [[[hex(row) for row in block] for block in rows[i]]
                             for i in selected],
        "basis_column_blocks": [[[hex(row) for row in block] for block in cols[i]]
                                for i in selected],
        "format": (
            "Each basis entry is a list of direct-sum blocks; each block is a list "
            "of binary row bit masks in hexadecimal. Both cut ranks are sums of "
            "block ranks over the coefficient field. Labels are identical on all "
            "four legs. The crossing endpoint is basis vector zero."),
        "exporter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def recombine(code, edges):
    result, (degrees, rows, cols) = rank_profile(code, edges, return_blocks=True)
    width = len(degrees)
    if width > 12:
        raise ValueError("exact recombination is restricted to at most 12 coefficients")
    rwork = [[0] * size for size in code["block_sizes"]]
    cwork = [[0] * size for size in code["block_sizes"]]
    ranks = [[0, 0] for _ in range(1 << width)]
    previous = 0
    for index in range(1, 1 << width):
        gray = index ^ (index >> 1)
        changed = (gray ^ previous).bit_length() - 1
        for target, source in ((rwork, rows[changed]), (cwork, cols[changed])):
            for block, delta in zip(target, source):
                for row, value in enumerate(delta):
                    block[row] ^= value
        ranks[gray] = [sum(rank_two(block) for block in rwork),
                       sum(rank_two(block) for block in cwork)]
        previous = gray
    cost, masks = minimum_decomposition([math.sqrt(r * c) for r, c in ranks])
    lower, _ = minimum_decomposition([math.isqrt(r * c) for r, c in ranks])
    result.update(
        coefficient_degrees=degrees, combination_ranks=ranks,
        optimal_cost=cost, optimal_combination_masks=masks,
        integer_lower_bound=lower, optimized_margin=result["threshold"] / cost,
    )
    return result


def search_case(name, labels, keep=12):
    started = time.monotonic()
    code = compile_code(labels)
    results = []
    for f, g, h, k, ell in product((0, 3, 5, 6), (0, 3, 5, 6), (1, 2, 4, 7), (0, 1), (0, 1)):
        edges = (f, g, h, k, ell, 1)
        result = rank_profile(code, edges)
        results.append(result)
    results.sort(key=lambda r: (r["rank_sum"], r["edges_UA_AY_AV_UV_VY_UY"]))
    finalists = [recombine(code, tuple(result["edges_UA_AY_AV_UV_VY_UY"]))
                 for result in results[:keep]]
    finalists.sort(key=lambda r: (r["optimal_cost"], r["edges_UA_AY_AV_UV_VY_UY"]))
    return {
        "name": name, "labels": code["labels"], "length": code["length"],
        "complementary_pairs": code["pairs"],
        "candidate_count": len(results), "evaluated_entries_per_candidate": code["evaluated_entries"],
        "distinct_local_types": code["distinct_local_types"],
        "best_monomial_profiles": results[:keep], "exact_recombination_finalists": finalists,
        "seconds": time.monotonic() - started,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="shell16:4:32,length2:2,length3:2,length3:3,length4:2")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--keep", type=int, default=12)
    parser.add_argument("--universal-span", action="store_true")
    parser.add_argument("--export-basis-dir", type=Path)
    args = parser.parse_args()
    report = {"scope": "exact endpoint-preserving Pfaffian seed screen; no stronger factoring exponent",
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "local_verification": verify_local_family(), "cases": []}
    for description in args.cases.split(","):
        if description.startswith("shell"):
            modulus, length, weight = map(int, description[5:].split(":"))
            labels = tuple(pair_occupations(parts, modulus) for parts in shell(modulus, length, weight))
            name = f"shell{modulus}_L{length}_j{weight}"
        elif description.startswith("length"):
            pairs, length = map(int, description[6:].split(":"))
            labels = all_length_labels(pairs, length)
            name = f"all_P{pairs}_L{length}"
        else:
            raise ValueError("use shellM:L:j or lengthP:L")
        if args.export_basis_dir:
            args.export_basis_dir.mkdir(parents=True, exist_ok=True)
            path = args.export_basis_dir / ("pfaffian-basis-" + name + ".json")
            data = export_basis(name, labels)
            path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
            print(json.dumps({"basis": str(path), "labels": len(labels),
                              "dimension": data["basis_dimension"]}), flush=True)
            continue
        if args.universal_span:
            code = compile_code(labels)
            result = {"name": description, "labels": len(labels),
                      "universal_span": universal_span(code)}
            best = result["universal_span"]
            print(json.dumps({"case": description, "labels": len(labels),
                              "span_dimension": best["span_dimension"],
                              "cost": best["optimal_cost"], "threshold": best["threshold"],
                              "seconds": best["seconds"]}), flush=True)
            report["cases"].append(result)
            if args.output:
                args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            continue
        result = search_case(description, labels, args.keep)
        report["cases"].append(result)
        best = result["exact_recombination_finalists"][0]
        print(json.dumps({"case": description, "labels": len(labels),
                          "cost": best["optimal_cost"], "threshold": best["threshold"],
                          "edges": best["edges_UA_AY_AV_UV_VY_UY"],
                          "seconds": result["seconds"]}), flush=True)
        if args.output:
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
