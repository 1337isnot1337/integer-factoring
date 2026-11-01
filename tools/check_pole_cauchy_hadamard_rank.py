#!/usr/bin/env python3
"""Exact Smith and separation-rank controls for entrywise Cauchy powers.

This checks a representation theorem, not a compressed factor-count algorithm.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random

from check_pole_smith_stabilization import Quad, matrix_at, smith_orders
from check_pole_cauchy_conditioning import inverse


def valuation(value: int) -> int | None:
    return None if value == 0 else (abs(value) & -abs(value)).bit_length() - 1


def power(value: Quad, exponent: int) -> Quad:
    result = Quad.of(1)
    while exponent:
        if exponent & 1:
            result = result * value
        value = value * value
        exponent >>= 1
    return result


def ceil_div(n: int, d: int) -> int:
    return (n + d - 1) // d


def matching_data(n: int, rows: tuple[int, ...], cols: tuple[int, ...], k: int):
    m = len(rows)
    assert m == len(cols) and m & (m - 1) == 0
    s = m.bit_length() - 1
    modulus = 2 * m
    column_at = {b % modulus: i for i, b in enumerate(cols)}
    assert len(column_at) == m
    assert len({a % modulus for a in rows}) == m
    assert all(a % 2 for a in rows + cols)
    pairing, orders = [], []
    for a in rows:
        index = column_at[n * pow(a, -1, modulus) % modulus]
        value = valuation(n - a * cols[index])
        pairing.append(index)
        orders.append(k - min(k, value if value is not None else k))
    assert len(set(pairing)) == m
    assert all(0 <= value <= k - s - 1 for value in orders)
    return s, pairing, orders


def check_matrix(n: int, rows: tuple[int, ...], cols: tuple[int, ...], p: int):
    m = len(rows)
    h = max(1, max(abs(n - a * b).bit_length() for a in rows for b in cols))
    h = max(h, m.bit_length())
    k = h + p
    s, pairing, orders = matching_data(n, rows, cols, k)
    matrix = matrix_at(n, rows, cols, k)
    profile_checks = approximation_checks = ratio_checks = 0
    records = []
    for exponent in (1, 2, 3, 5):
        powered = [[power(value, exponent) for value in row] for row in matrix]
        expected = sorted(exponent * value for value in orders)
        assert smith_orders(powered) == expected
        profile_checks += 1
        for i, row in enumerate(powered):
            pivot = row[pairing[i]]
            assert pivot.valuation() == exponent * orders[i]
            for j, value in enumerate(row):
                if j != pairing[i]:
                    assert (value / pivot).valuation() >= exponent
                    ratio_checks += 1
        limit = k - s - 1
        for q in sorted({1, p + limit, exponent * limit, exponent * limit + 1}):
            q = max(1, q)
            retained = [i for i, value in enumerate(orders) if exponent * value < q]
            level = max(s + 1, k - (q - 1) // exponent)
            modulus = 1 << level
            overlap = sum((n - a * cols[pairing[i]]) % modulus == 0
                          for i, a in enumerate(rows))
            assert overlap == len(retained)
            pruned = [powered[i] for i in retained]
            assert smith_orders(pruned) == [value for value in expected if value < q]
            for i, row in enumerate(powered):
                if i not in retained:
                    assert all(value.valuation() >= q for value in row)
            approximation_checks += 1
        records.append({"entrywise_power": exponent, "smith_orders": expected})
    q = p + k - s - 1
    level = max(s + 1, k - (q - 1) // 2)
    assert level == max(s + 1, ceil_div(h + s + 2, 2))
    rank = sum(2 * value < q for value in orders)
    padded_k = max(h + p, s + 1, 2 * h + p - s - 3)
    padded_q = p + padded_k - s - 1
    padded_matrix = matrix_at(n, rows, cols, padded_k)
    padded_square = [[value * value for value in row] for row in padded_matrix]
    factor_count = sum(n == a * b for a in rows for b in cols)
    for i, a in enumerate(rows):
        for j, b in enumerate(cols):
            if n == a * b:
                assert padded_square[i][j] == Quad.of(1)
            else:
                assert padded_square[i][j].valuation() >= padded_q
    padded_orders = smith_orders(padded_square)
    assert sum(value < padded_q for value in padded_orders) == factor_count
    inv = inverse(padded_matrix)
    trace = sum((inv[j][i] * padded_square[i][j]
                 for i in range(m) for j in range(m)), Quad())
    incidence_trace = sum((inv[j][i] for i, a in enumerate(rows)
                           for j, b in enumerate(cols) if a * b == n), Quad())
    assert not trace - factor_count or (trace - factor_count).valuation() >= p
    assert not incidence_trace - factor_count or (incidence_trace - factor_count).valuation() >= p
    return {"n": n, "rows": rows, "columns": cols, "p": p, "h": h, "k": k,
            "matching_columns": pairing, "base_row_orders": orders,
            "power_profiles": records, "trace_precision": q,
            "trace_overlap_level": level, "minimum_trace_separation_rank": rank,
            "higher_padding": {"k": padded_k, "kernel_precision": padded_q,
                               "minimum_separation_rank": factor_count,
                               "factor_count": factor_count,
                               "hadamard_square_orders": padded_orders}}, {
                "smith_profiles": profile_checks,
                "optimal_approximation_controls": approximation_checks,
                "off_diagonal_ratio_controls": ratio_checks,
                "higher_padding_rank_controls": 1,
                "higher_padding_entry_controls": m * m,
                "higher_padding_trace_controls": 2,
            }


def witness_family(s: int, root_free: bool, literal: bool) -> dict:
    m = 1 << s
    n = (16 if root_free else 1) * m * m + (m + 1) ** 2
    h = 2 * s + (5 if root_free else 2)
    assert n.bit_length() == h
    assert max(abs(n - 1), abs(n - (2 * m - 1) ** 2)).bit_length() == h
    if root_free:
        assert n > (2 * m - 1) ** 2
    else:
        assert 1 < n < (2 * m - 1) ** 2
    level = ceil_div(h + s + 2, 2)
    u = ceil_div(level, 2)
    assert level >= s + 1 and level <= 2 * s
    assert u <= s - 1
    radius = (1 << (s - u)) - 1
    witness_count = 2 * radius + 1
    if not root_free:
        assert u == s + 1 - s // 4
        assert witness_count == (1 << (s // 4)) - 1
    samples = {-radius, -1, 0, 1, radius}
    samples.update(1 << i for i in range(radius.bit_length()) if (1 << i) <= radius)
    samples.update(-j for j in tuple(samples))
    for j in samples:
        a, b = m + 1 + (1 << u) * j, m + 1 - (1 << u) * j
        assert 1 <= a <= 2 * m - 1 and 1 <= b <= 2 * m - 1
        assert a % 2 and b % 2
        difference = (16 if root_free else 1) * m * m + (1 << (2 * u)) * j * j
        assert n - a * b == difference > 0
        assert difference % (1 << level) == 0
        assert n * pow(a, -1, 2 * m) % (2 * m) == b
        assert not (n % a == 0 and 1 <= n // a <= 2 * m - 1)
    rank_controls = []
    for p in (1, 8, s + 1, 3 * s + 7):
        k = h + p
        scale = k - s - 1
        precision = p + scale
        assert k - (precision - 1) // 2 == level
        rank_controls.append({"p": p, "k": k, "inverse_scale": scale,
                              "kernel_precision": precision, "overlap_level": level})
    record = {"s": s, "dimension": m, "n": n, "input_bits": h,
              "zero_factor_incidence_from_corners": root_free,
              "overlap_level": level, "witness_step_bits": u,
              "certified_nonfactor_rows": witness_count,
              "pointwise_witness_samples": len(samples),
              "precision_controls": rank_controls,
              "literal_overlap_enumerated": literal}
    if literal:
        modulus = 1 << level
        overlap = factors = 0
        for a in range(1, 2 * m, 2):
            b = n * pow(a, -1, 2 * m) % (2 * m)
            overlap += (n - a * b) % modulus == 0
            factors += n == a * b
        assert overlap >= factors + witness_count
        if root_free:
            assert factors == 0
        record.update({"literal_overlap": overlap, "literal_factor_count": factors})
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(2026110101)
    totals = {"smith_profiles": 0, "optimal_approximation_controls": 0,
              "off_diagonal_ratio_controls": 0,
              "higher_padding_rank_controls": 0,
              "higher_padding_entry_controls": 0,
              "higher_padding_trace_controls": 0}
    cases = []
    for s in range(4):
        m = 1 << s
        for trial in range(5 if s < 3 else 3):
            residues = tuple(range(1, 2 * m, 2))
            rows = tuple(a + 2 * m * rng.randrange(4) for a in residues)
            cols = tuple(b + 2 * m * rng.randrange(4) for b in residues)
            n = rows[trial % m] * cols[(trial + 1) % m] if trial % 2 else rng.randrange(3, 300, 2)
            result, counts = check_matrix(n, rows, cols, 1 + trial % 3)
            cases.append(result)
            for key, value in counts.items():
                totals[key] += value
    families = [witness_family(s, root_free, s <= 12)
                for root_free in (False, True) for s in (12, 16, 32, 80, 200, 400)]
    result = {
        "theorem": "Exact finite-precision separation ranks of entrywise Cauchy powers",
        "arithmetic": "Exact Q[w]/(w^2+w+1), integer modular inverse controls",
        "matrix_cases": len(cases), **totals, "cases": cases,
        "explicit_exponential_families": families,
        "scope": "Whole-kernel separated approximation lower bound. Large families use proved symmetric witnesses, not global enumeration. No factoring lower bound or improved score.",
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
