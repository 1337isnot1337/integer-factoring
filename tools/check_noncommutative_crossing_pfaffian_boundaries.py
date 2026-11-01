#!/usr/bin/env python3
"""Exact audit and finite basis export for Pfaffian boundary sections.

The proof is in research/noncommutative-crossing-pfaffian-boundaries.md.
No search outcome or factoring improvement is assumed by this checker.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import json
import math
from pathlib import Path

from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local
from check_modular_crossing import rank_two


def omega_boundary(source, pairs):
    """Images of a squarefree monomial under wedge with sum e_2i e_(2i+1)."""
    return tuple(source | (3 << (2 * i)) for i in range(pairs)
                 if not source & (3 << (2 * i)))


def occupations(mask, pairs):
    return tuple((mask >> (2 * i)) & 3 for i in range(pairs))


def parity_syndrome(mask, pairs):
    return sum((((mask >> (2 * i)) & 3).bit_count() & 1) << i
               for i in range(pairs))


def powers(value, maximum):
    result = [1]
    for _ in range(maximum):
        result.append(multiply(result[-1], value))
    return result


def even_boundary_gram(pairs):
    """Independent multiplication of all middle-grade even-sector boundaries."""
    assert pairs > 0 and pairs % 2 == 0
    half = pairs // 2
    base = pairs + 1
    p, q, r = 2, 1 << base, 1 << (base * base)
    s = p ^ q ^ r
    local = ((p, q), (r, s))
    a, b = multiply(p, s), multiply(q, r)
    apow, bpow = powers(a, half), powers(b, half)
    sources = [sum(1 << i for i in subset)
               for subset in combinations(range(pairs), half - 1)]
    checks = 0
    for first, second in product(sources, repeat=2):
        actual = 0
        for i in range(pairs):
            if first >> i & 1:
                continue
            target_first = first | (1 << i)
            for j in range(pairs):
                if second >> j & 1:
                    continue
                target_second = second | (1 << j)
                value = 1
                for k in range(pairs):
                    value = multiply(
                        value, local[(target_first >> k) & 1]
                                    [(target_second >> k) & 1])
                actual ^= value
        intersection = (first & second).bit_count()
        different = half - 1 - intersection
        expected = 0
        if intersection & 1:
            expected ^= multiply(apow[intersection + 1],
                                 bpow[half - intersection - 1])
        if different & 1:
            expected ^= multiply(apow[intersection + 2],
                                 bpow[half - intersection - 2])
        assert actual == expected, (pairs, first, second)
        checks += 1
    return {"even_pairs": pairs, "source_monomials": len(sources),
            "exact_gram_entries": checks}


def literal_boundary_audit(pairs):
    """Contract literal local Pfaffians, without the proposed Gram formula."""
    base = pairs + 1
    p, q, r = 2, 1 << base, 1 << (base * base)
    local = evaluate_local((p, q, r, 1, 1, 1))
    groups = defaultdict(list)
    for subset in combinations(range(2 * pairs), pairs - 2):
        source = sum(1 << i for i in subset)
        boundary = omega_boundary(source, pairs)
        if boundary:
            groups[parity_syndrome(source, pairs)].append(boundary)
    singles = tuple(product((1, 2), repeat=pairs))
    expected_scalar = powers(multiply(p, p ^ q ^ r), 2)[2]

    @lru_cache(None)
    def entry(x, u, v, y):
        value = 1
        for a, b, c, d in zip(x, occupations(y, pairs), u,
                              occupations(v, pairs)):
            value = multiply(value, local[((a * 4 + b) * 4 + c) * 4 + d])
            if not value:
                break
        return value

    checks = 0
    for syndrome, boundaries in groups.items():
        for x, u in product(singles, repeat=2):
            for first, second in product(boundaries, repeat=2):
                actual = 0
                for v, y in product(first, second):
                    actual ^= entry(x, u, v, y)
                expected = (expected_scalar if pairs == 4 and syndrome == 0
                            and x == u else 0)
                assert actual == expected, (pairs, syndrome, x, u,
                                            first, second)
                checks += 1
    return {"pairs": pairs, "labels": 1 << pairs,
            "boundary_generators": sum(map(len, groups.values())),
            "syndrome_groups": len(groups), "literal_contractions": checks}


def binary_basis(pairs=4):
    """Coefficient tensors of product_i (a E_i + b C_i).

    Label i has occupation 1+bit_k(i) on complementary pair k.
    Row cut: (x,y) | (u,v). Column cut: (x,u) | (v,y).
    Coefficient j is a**(pairs-j) b**j; j=0 is the crossing.
    """
    d = 1 << pairs
    rows = [[0] * (d * d) for _ in range(pairs + 1)]
    cols = [[0] * (d * d) for _ in range(pairs + 1)]
    for x, y, u, v in product(range(d), repeat=4):
        value = 1
        for i in range(pairs):
            a, b, c, e = ((word >> i) & 1 for word in (x, y, u, v))
            local = int(a == e and b == c) ^ (int(a == c and b == e) << 1)
            value = multiply(value, local)
            if not value:
                break
        while value:
            bit = value & -value
            degree = bit.bit_length() - 1
            rows[degree][d * x + y] ^= 1 << (d * u + v)
            cols[degree][d * x + u] ^= 1 << (d * v + y)
            value ^= bit
    return rows, cols


def correction_basis(f, pairs=4):
    d = 1 << pairs
    if not 0 <= f < 1 << d:
        raise ValueError("f must be a binary label functional")
    rows, cols = [0] * (d * d), [0] * (d * d)
    input_pair = sum(1 << (d * v + y) for v, y in product(range(d), repeat=2)
                     if f >> v & 1 and f >> y & 1)
    for x in range(d):
        cols[d * x + x] = input_pair
        for y in range(d):
            if f >> y & 1:
                rows[d * x + y] = f << (d * x)
    return rows, cols


def all_field_charge_certificate(binary=None):
    """Integer certificate plus a binary audit of the proved rank law."""
    v = (256, 120, 46, 10, 1)
    weights = (15, 30, 36, 23, 16)
    intervals = []
    for low in range(5):
        for high in range(low, 5):
            charge = sum(weights[low:high + 1])
            baseline = v[low] * v[4 - high]
            perturbed = (max(1, abs(v[low] - 16)) *
                         max(1, abs(v[4 - high] - 1)))
            assert charge * charge <= min(baseline, perturbed)
            intervals.append({"least": low, "greatest": high,
                              "charge": charge,
                              "baseline_product": baseline,
                              "perturbed_product_lower": perturbed})
    row, col = binary if binary is not None else binary_basis()
    n_cuts = []
    for cut in (row, col):
        transformed = []
        for degree in range(5):
            matrix = [0] * 256
            for j in range(5):
                exponent = degree - 4 + j
                if 0 <= exponent <= j and math.comb(j, exponent) & 1:
                    for index, value in enumerate(cut[j]):
                        matrix[index] ^= value
            transformed.append(matrix)
        n_cuts.append(transformed)
    for mask in range(1, 32):
        low = (mask & -mask).bit_length() - 1
        high = mask.bit_length() - 1
        ranks = []
        for cut in n_cuts:
            matrix = [0] * 256
            for j in range(5):
                if mask >> j & 1:
                    for index, value in enumerate(cut[j]):
                        matrix[index] ^= value
            ranks.append(rank_two(matrix))
        assert ranks == [v[low], v[4 - high]], (mask, ranks)
    return {"rank_sequence": list(v), "coordinate_weights": list(weights),
            "certified_rank_sum_lower": sum(weights),
            "binary_rank_law_combinations_checked": 31,
            "intervals": intervals,
            "field_scope": "Every characteristic-two field by the ordinary proof; "
                           "binary tests only audit the matrix convention."}


def export_basis(prefix, f, binary=None):
    pairs, d = 4, 16
    row, col = binary if binary is not None else binary_basis(pairs)
    rows, cols = list(row), list(col)
    if f:
        correction_row, correction_col = correction_basis(f, pairs)
        assert rank_two(correction_row) == d
        assert rank_two(correction_col) == 1
        rows.append(correction_row)
        cols.append(correction_col)
    dimension = len(rows)
    payload = {
        "name": "pfaffian_boundary_P4_f_%04x" % f,
        "labels": d, "length": pairs, "pairs": pairs, "functional_f": f,
        "occupations": [[1 + ((label >> i) & 1) for i in range(pairs)]
                        for label in range(d)],
        "block_sizes": [d * d], "basis_dimension": dimension,
        "basis_description": [
            "coefficient of a^%d b^%d in product(a E_i+b C_i)" % (pairs-i, i)
            for i in range(pairs + 1)] + (["Z_f=delta(x,u) f(y) f(v)"] if f else []),
        "endpoint_coordinates": [1] + [0] * (dimension - 1),
        "binary_variables": {"a": "r+q", "b": "p+q"},
        "correction_polynomial_p_q_r": [[4, 0, 0], [2, 2, 0], [2, 0, 2]],
        "correction_polynomial_a_b_q": [[2, 2, 0], [0, 4, 0],
                                       [2, 0, 2], [0, 0, 4]],
        "basis_row_blocks": [[[hex(value) for value in matrix]]
                             for matrix in rows],
        "basis_column_blocks": [[[hex(value) for value in matrix]]
                                for matrix in cols],
        "format": "Complete dense binary cut matrices as hexadecimal row masks; "
                  "row cut (x,y)|(u,v), column cut (x,u)|(v,y).",
        "exporter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    prefix = Path(prefix)
    prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    gate = math.ceil(d ** math.sqrt(2))
    words = (d * d + 63) // 64
    with prefix.with_suffix(".txt").open("w") as handle:
        handle.write("%d %d 1 %d\n%d\n" %
                     (d, dimension, gate, d * d))
        for cut in (rows, cols):
            for matrix in cut:
                for value in matrix:
                    handle.write(" ".join(
                        format((value >> (64 * i)) & ((1 << 64) - 1), "x")
                        for i in range(words)) + "\n")
    return {"prefix": str(prefix), "functional_f": f,
            "basis_dimension": dimension, "integer_cost_gate": gate,
            "basis_ranks": [[rank_two(a), rank_two(b)]
                            for a, b in zip(rows, cols)]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    parser.add_argument("--export-prefix")
    parser.add_argument("--functional", type=lambda value: int(value, 0), default=1)
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()
    report = {"claim": "Exact boundary-section formula; no seed success asserted.",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    binary = binary_basis() if (not args.skip_audit or args.export_prefix) else None
    if not args.skip_audit:
        report["literal_audits"] = [literal_boundary_audit(pairs)
                                    for pairs in (2, 3, 4)]
        report["all_size_gram_identity_checks"] = [even_boundary_gram(pairs)
                                                   for pairs in (2, 4, 6, 8)]
        report["all_field_charge_certificate"] = all_field_charge_certificate(binary)
    if args.export_prefix:
        report["export"] = export_basis(args.export_prefix, args.functional, binary)
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
