#!/usr/bin/env python3
"""Literal checks for the two-cut Euler refinement of the Ward factorization."""

from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json

from check_modular_crossing import rank_two
from check_noncommutative_crossing_ward_rank_descent import (
    entry, labels, rectangle_profile, supports,
)
from check_noncommutative_crossing_symplectic_radicals import omega, contraction
from check_noncommutative_crossing_radical_ward import local_clifford
from explore_complement_edge_pencils import evaluate_local
from explore_overlapping_symplectic_codes import overlapping_code
from explore_symplectic_boundary_codes import boundary_code, half_orientation_code


def mask(word):
    return sum(value << (2 * i) for i, value in enumerate(word))


def occupations(value, pairs):
    return tuple((value >> (2 * i)) & 3 for i in range(pairs))


def homogeneous_preimages(code):
    pairs, length = code["pairs"], code["length"]
    source = labels(pairs, length - 2)
    pivots = {}
    for index, word in enumerate(source):
        value = mask(word)
        image = omega({value}, pairs) ^ contraction({value}, pairs)
        column = sum(1 << item for item in image)
        coordinates = 1 << index
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column, coordinates
                break
            old, old_coordinates = pivots[pivot]
            column ^= old
            coordinates ^= old_coordinates
    result = []
    ambient = code["ambient_occupations"]
    for vector in code["input_vectors"]:
        target = sum(1 << mask(ambient[i]) for i in supports(vector))
        coordinates = 0
        while target:
            pivot = target.bit_length() - 1
            assert pivot in pivots, "Homogeneous Ward preimage missing."
            old, old_coordinates = pivots[pivot]
            target ^= old
            coordinates ^= old_coordinates
        words = tuple(source[i] for i in supports(coordinates))
        reconstructed = set()
        for word in words:
            value = mask(word)
            reconstructed ^= omega({value}, pairs) ^ contraction({value}, pairs)
        assert reconstructed == {mask(ambient[i]) for i in supports(vector)}
        result.append(words)
    return result


def parity_filter(value, degree, stride):
    result = 0
    for exponent in supports(value):
        a, b = exponent % stride, exponent // stride
        if (degree + a + b) & 1:
            result ^= 1 << exponent
    return result


def local_derivative_checks():
    stride = 8
    p, q, r = 2, 1 << stride, 1 << (stride * stride)
    table = evaluate_local((p, q, r, 1, 1, 1))
    clifford = local_clifford()

    def derivative(value):
        result = 0
        for exponent in supports(value):
            for place in (1, stride, stride * stride):
                if exponent // place % stride & 1:
                    result ^= 1 << (exponent - place)
        return result

    def centered(value):
        result = 0
        for exponent in supports(value):
            a, b = exponent % stride, exponent // stride % stride
            for first in range(a + 1):
                if first & ~a:
                    continue
                for second in range(b + 1):
                    if not second & ~b:
                        result ^= 1 << (first + 3 * second)
        return result

    for x, y, u, v in product(range(4), repeat=4):
        value = table[((x * 4 + y) * 4 + u) * 4 + v]
        acted = 0
        for mode in range(2):
            external = mode ^ 1
            if x >> external & 1:
                source = x ^ (1 << external)
                for state in supports(clifford[mode][y]):
                    acted ^= table[((source * 4 + state) * 4 + u) * 4 + v]
        assert acted == derivative(value)
        assert centered(acted) == parity_filter(centered(value), x.bit_count(), 3)
    return {"generic_incidence_derivative_equations": 256,
            "centered_Euler_equations": 256}


def gram_rows(columns, size):
    return [sum((sum(((value >> a) & 1) * ((value >> b) & 1)
                         for value in columns) % 2) << b for b in range(size))
            for a in range(size)]


def principal_rank(rows, selected):
    indices = supports(selected)
    return rank_two([sum(((rows[a] >> b) & 1) << i for i, b in enumerate(indices))
                     for a in indices])


def typed_scalar_checks():
    records = []
    for size, columns in ((1, (0, 1)), (2, (0, 1, 2, 3)),
                          (3, (1, 2, 3, 4, 5, 6, 7))):
        count = 0
        for active in range(1 << len(columns)):
            active_columns = [columns[i] for i in supports(active)]
            gram = gram_rows(active_columns, size)
            for rows in range(1 << size):
                chosen = supports(rows)
                degree = len(chosen)
                scalar = 0
                for subset in combinations(active_columns, degree):
                    minor = [sum(((value >> a) & 1) << b for b, value in enumerate(subset))
                             for a in chosen]
                    scalar ^= int(rank_two(minor) == degree)
                assert scalar == int(principal_rank(gram, rows) == degree)
                if degree > rank_two(gram):
                    assert scalar == 0
                count += 1
        records.append({"forms": size, "types": list(columns),
                        "all_typed_jet_parity_and_row_subset_checks": count})
    return records


def coefficient_ranks(polynomials, d, stride):
    rows, columns = {}, {}
    for (x, y, u, v), value in polynomials.items():
        for exponent in supports(value):
            key = exponent % stride, exponent // stride
            if key not in rows:
                rows[key] = [0] * (d * d)
                columns[key] = [0] * (d * d)
            rows[key][x * d + y] ^= 1 << (v * d + u)
            columns[key][x * d + u] ^= 1 << (v * d + y)
    return {key: (rank_two(rows[key]), rank_two(columns[key])) for key in rows}


def audit_case(code, name, perturb=False, full_rectangle_bounds=False):
    pairs, length, d = code["pairs"], code["length"], code["labels"]
    ambient = code["ambient_occupations"]
    injection = [tuple(ambient[i] for i in supports(vector)) for vector in code["input_vectors"]]
    outputs = [(ambient[i],) for i in code["output_indices"]]
    if perturb:
        # A completely single-occupied monomial annihilates every radical input.
        extra = next(word for word in ambient if all(value in (0, 1, 2) for value in word))
        outputs = [tuple(dict.fromkeys((*word, extra))) for word in outputs]
        assert all(not (vector >> ambient.index(extra) & 1) for vector in code["input_vectors"])
    preimages = homogeneous_preimages(code)
    stride = 2 * pairs + 1
    polynomials = {}
    equations = 0
    for x, y, u, v in product(range(d), repeat=4):
        actual = first = second = 0
        for xx, uu, vv in product(outputs[x], outputs[u], injection[v]):
            for yy in injection[y]:
                actual ^= entry(pairs, xx, yy, uu, vv)
            for yy in preimages[y]:
                first ^= entry(pairs, xx, yy, uu, vv)
                for site, value in enumerate(xx):
                    if value == 3:
                        reduced = xx[:site] + (0,) + xx[site + 1:]
                        second ^= entry(pairs, reduced, yy, uu, vv)
        expected = (parity_filter(first, length, stride) << stride) ^ (second << (2 * stride))
        assert actual == expected, (name, x, y, u, v, actual, expected)
        polynomials[x, y, u, v] = actual
        equations += 1
    ranks = coefficient_ranks(polynomials, d, stride)
    records = []
    if full_rectangle_bounds:
        for (a, b), (row, column) in sorted(ranks.items()):
            flag = (length + a + b - 1) & 1
            r1 = rectangle_profile(pairs, length, length - 2, length, length).get((a, b - 1), 0)
            c1 = rectangle_profile(pairs, length, length, length, length - 2).get((b - 1, a), 0)
            r2 = rectangle_profile(pairs, length - 2, length - 2, length, length).get((a, b - 2), 0)
            c2 = rectangle_profile(pairs, length - 2, length, length, length - 2).get((b - 2, a), 0)
            row_bound, column_bound = flag * r1 + r2, flag * c1 + c2
            assert row <= row_bound and column <= column_bound, (name, a, b, row, column, row_bound, column_bound)
            records.append([a, b, row, column, row_bound, column_bound])
    return {"name": name, "pairs": pairs, "grade": length, "labels": d,
            "noncoordinate_outputs": perturb, "complete_tensor_polynomial_equations": equations,
            "maximum_preimage_terms": max(map(len, preimages), default=0),
            "jets_actual_row_actual_column_bound_row_bound_column": records}


def restricted_case(pairs, length, count):
    code = boundary_code(pairs, length, True)
    result = dict(code)
    for key in ("input_vectors", "output_indices", "syndromes"):
        result[key] = code[key][:count]
    result["labels"] = count
    return result


def overlapping_multiplicity(columns, length, limit=None):
    code = overlapping_code(columns, length)
    original = code["labels"]
    if limit is not None:
        for key in ("input_vectors", "output_indices", "syndromes"):
            code[key] = code[key][:limit]
        code["labels"] = min(limit, original)
    pairs, d, forms = code["pairs"], code["labels"], code["form_count"]
    ambient = code["ambient_occupations"]
    inputs = [tuple(ambient[i] for i in supports(vector)) for vector in code["input_vectors"]]
    outputs = [ambient[i] for i in code["output_indices"]]
    stride = 2 * pairs + 1
    total_rank = rank_two(gram_rows(columns, forms))
    required = 2 * forms - total_rank
    actual_minimum = None
    equations = 0
    for x, y, u, v in product(range(d), repeat=4):
        value = 0
        for yy, vv in product(inputs[y], inputs[v]):
            value ^= entry(pairs, outputs[x], yy, outputs[u], vv)
        for exponent in supports(value):
            a, b = exponent % stride, exponent // stride
            assert b >= forms and a + b >= required
            actual_minimum = a + b if actual_minimum is None else min(actual_minimum, a + b)
        equations += 1
    return {"columns": columns, "pairs": pairs, "grade": length,
            "forms": forms, "full_labels": original, "checked_labels": d,
            "total_Gram_rank": total_rank, "proved_joint_multiplicity": required,
            "literal_minimum_total_degree": actual_minimum,
            "complete_tensor_polynomial_entries": equations}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    p3 = boundary_code(3, 3, True)
    cases = [audit_case(boundary_code(2, 2, True), "P2_L2", full_rectangle_bounds=True),
             audit_case(p3, "P3_L3", full_rectangle_bounds=True),
             audit_case(half_orientation_code(p3), "P3_L3_half", full_rectangle_bounds=True),
             audit_case(p3, "P3_L3_noncoordinate", perturb=True, full_rectangle_bounds=True),
             audit_case(restricted_case(4, 4, 3), "P4_L4_three_labels", perturb=True)]
    result = {"status": "pass", "scope": "Complete polynomial identities and both-cut rank bounds for the global-form Euler refinement; no improved exponent asserted.",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "local_derivative_and_Euler": local_derivative_checks(),
              "typed_Gram_scalar_controls": typed_scalar_checks(),
              "overlapping_joint_multiplicity_controls": [
                  overlapping_multiplicity([1, 1, 2, 2, 3, 3], 6),
                  overlapping_multiplicity([1, 1, 1, 2, 2, 3, 3], 7, 3)],
              "cases": cases, "literal_rectangles_constructed": rectangle_profile.cache_info().currsize}
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
