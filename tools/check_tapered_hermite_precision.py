#!/usr/bin/env python3
"""Check a polynomial witness for the tapered Hermite precision scale."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import monotonic


def valuation2(value):
    if not value:
        raise ValueError("zero has no finite valuation")
    value = abs(value)
    return (value & -value).bit_length() - 1


def multiply_linear(coefficients, constant, linear):
    answer = [0] * (len(coefficients) + 1)
    for j, coefficient in enumerate(coefficients):
        answer[j] += constant * coefficient
        answer[j + 1] += linear * coefficient
    return answer


def literal_shifted_numerator(a, u, M, t):
    """Expand the whole numerator using only ordinary integer arithmetic."""
    L, alpha, h = 1 << a, 1 << (a + 1), M // u
    coefficients = [1]
    denominator = 1
    for q in range(-L, L):
        if not q:
            continue
        multiplicity = h * max(0, u - valuation2(q))
        for _ in range(multiplicity):
            coefficients = multiply_linear(coefficients, q - t, -alpha)
            denominator *= q
    return coefficients, denominator


def predicted_content(a, u, M, t):
    if not t or valuation2(t) >= u:
        return 0
    i = valuation2(t)
    h, D = M // u, 1 << (a - i)
    m = h * (u - i)
    return m * (a + 1 - i) + h * sum(
        (D >> j) * j * min(j, u - i) for j in range(1, a - i + 1)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    expanded_polynomials = coefficient_integrality_checks = active_checks = 0
    rows = []
    for a in range(1, 5):
        L = 1 << a
        for u in range(1, a + 1):
            for M in sorted({u, u + 1, 2 * u + 1}):
                h = M // u
                expected_degree = sum((L >> i) * h * (u - i) for i in range(u))
                assert expected_degree <= M * (2 * L - 1)
                contents = {}
                for t in range(-L, L):
                    coefficients, denominator = literal_shifted_numerator(a, u, M, t)
                    assert len(coefficients) - 1 == expected_degree
                    denominator_v = valuation2(denominator)
                    actual_content = min(valuation2(c) for c in coefficients if c) - denominator_v
                    expected_content = predicted_content(a, u, M, t)
                    assert actual_content == expected_content >= 0
                    for coefficient in coefficients:
                        if coefficient:
                            assert valuation2(coefficient) >= denominator_v
                        coefficient_integrality_checks += 1
                    if t == 0:
                        assert coefficients[0] == denominator
                    elif valuation2(t) < u:
                        i = valuation2(t)
                        D = 1 << (a - i)
                        assert 2 * actual_content >= h * D
                        assert 4 * u * actual_content >= M * D
                        contents[i] = actual_content
                        active_checks += 1
                    expanded_polynomials += 1
                rows.append({"a": a, "u": u, "M": M, "degree": expected_degree,
                             "maximum_allowed_degree": M * (2 * L - 1),
                             "content_by_active_valuation": contents})

    # Broader integer checks of the closed formula and degree estimate.
    all_size_parameter_checks = 0
    for a in range(1, 65):
        L = 1 << a
        for u in range(1, a + 1):
            for M in [u, u + 1, 3 * u - 1, 100 * u + 17]:
                h = M // u
                degree = sum((L >> i) * h * (u - i) for i in range(u))
                assert degree <= M * (2 * L - 1)
                for i in range(u):
                    value = predicted_content(a, u, M, 1 << i)
                    D = 1 << (a - i)
                    assert 2 * value >= h * D
                    assert 4 * u * value >= M * D
                    all_size_parameter_checks += 1

    report = {
        "status": "PASS",
        "literal_shifted_polynomial_expansions": expanded_polynomials,
        "literal_coefficient_integrality_checks": coefficient_integrality_checks,
        "literal_active_center_content_checks": active_checks,
        "closed_formula_parameter_checks": all_size_parameter_checks,
        "small_witnesses": rows,
        "seconds": round(monotonic() - started, 3),
        "scope": "An indistinguishability witness for the source polynomial class and residue grid, not a factoring lower bound.",
    }
    output = json.dumps(report, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
