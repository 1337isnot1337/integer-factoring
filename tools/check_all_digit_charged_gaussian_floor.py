#!/usr/bin/env python3
"""Certify positive Gaussian digit derivatives on all required compact bases.

The research note supplies the unbounded-parameter proof. This checker
rebuilds complete compact covers using only outward integer arithmetic,
checks their tree coverage independently, and verifies the rational
inequalities used in the analytic tails and scalar kernels.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from all_digit_charged_gaussian_floor import DigitDerivativeBounds
from charged_gaussian_floor import split_rectangle
from check_small_charged_gaussian_bases import cancelled_pressure, exact_profile_controls


def interval_arguments(arithmetic, rectangle):
    endpoints = [arithmetic.rational(x.numerator, x.denominator) for x in rectangle]
    return (endpoints[0][0], endpoints[1][1]), (endpoints[2][0], endpoints[3][1])


def construct_cover(arithmetic, function, rectangle):
    pending = [("", rectangle)]
    leaves = []
    evaluated = 0
    while pending:
        path, current = pending.pop()
        lower = function(*interval_arguments(arithmetic, current))
        evaluated += 1
        if lower > 0:
            leaves.append([path, str(lower)])
        else:
            assert len(path) < 32, (path, current, lower)
            pending.extend((path + symbol, child)
                           for symbol, child in reversed(split_rectangle(current)))
        assert evaluated < 200000
    return {
        "rectangle": [str(value) for value in rectangle],
        "evaluated_rectangles": evaluated,
        "positive_cover_rectangles": len(leaves),
        "minimum_lower_numerator": str(min(int(value) for path, value in leaves)),
        "cover": leaves,
    }


def audit_cover(arithmetic, function, record):
    """Reconstruct the full partition from the saved leaves, independently."""
    root = tuple(Fraction(value) for value in record["rectangle"])
    entries = dict(record["cover"])
    assert len(entries) == len(record["cover"]), "Duplicate leaf path"
    prefixes = {path[:length] for path in entries for length in range(len(path) + 1)}
    visited = set()
    nodes = 0

    def visit(path, rectangle):
        nonlocal nodes
        nodes += 1
        assert path in prefixes, "Missing portion of the claimed cover"
        if path in entries:
            assert not any(len(candidate) > len(path) and candidate.startswith(path)
                           for candidate in entries), "Overlapping leaf paths"
            lower = function(*interval_arguments(arithmetic, rectangle))
            assert int(entries[path]) == lower > 0, "Incorrect outward lower bound"
            visited.add(path)
            return
        for symbol, child in split_rectangle(rectangle):
            visit(path + symbol, child)

    visit("", root)
    assert visited == set(entries)
    assert nodes == record["evaluated_rectangles"]
    assert len(visited) == record["positive_cover_rectangles"]
    return nodes


def exact_analytic_controls():
    """Finite rational certificates named in the ordinary analytic proof."""
    F = Fraction
    checks = []

    def require(name, condition):
        assert condition, name
        checks.append(name)

    require("sqrt2 between 7/5 and 10/7", F(7, 5) ** 2 < 2 < F(10, 7) ** 2)
    require("sqrt2 below 3/2", 2 < F(3, 2) ** 2)
    require("e^3 above 20", sum(F(3) ** k / factorial(k) for k in range(9)) > 20)
    require("e^4 above 50", sum(F(4) ** k / factorial(k) for k in range(10)) > 50)
    eta = F(1, 20 ** 14)
    require("large-scale geometric coefficient below 21/20",
            1 / (1 - eta) + 44 * eta / (1 - eta) ** 2 < F(21, 20))
    require("small-scale Fourier error below t/900",
            F(160, 160000) / (1 - F(1, 160000)) ** 3 < F(1, 900))
    require("small-scale positive continuation", F(1, 6) - F(1, 900) == F(149, 900))
    require("two-mode exponential kernel", F(79, 7) ** 41 > F(159, 100) ** 200)
    require("edge-kernel lower radical", F(106, 75) ** 2 < 2)
    require("edge-kernel upper radical", 2 < F(1061, 750) ** 2)
    require("edge-kernel base ratio", F(1317, 122) > F(107, 10))
    require("edge exponential kernel", F(107, 10) ** 3 > F(33, 25) ** 25)
    require("periodic tail prefactor", F(40, 7) * F(21, 20) == 6)
    require("periodic low-charge large-scale ratio", F(150, 400) < 1)
    require("base zero-ground large-scale ratio", F(225, 400) < 1)
    require("base positive-ground large-scale ratio", F(711, 15 * 50) < 1)
    require("base ground covers physical terms", F(3, 5) * 16 > 9)
    require("residual range at length four", 2 < 6 - 2 * F(10, 7))
    for n in (1, 2, 3):
        require(f"small-scale finite base {n}", F(8 * n, 7) - F(4 * n * n + 8, 16) > 0)
    return checks


def exact_numerical_identity_controls(arithmetic):
    """Outward overlap checks for two independent complete expressions."""
    count = 0
    for n, left, right, side in ((1, 0, 2, 0), (2, 0, 2, 0),
                                (2, 2, 4, 1), (3, 0, 2, 0)):
        for quarter in range(4 * left, 4 * right + 1):
            charge = arithmetic.rational(quarter, 4)
            for scale in (Fraction(1, 4), Fraction(1), Fraction(4), Fraction(16)):
                time = arithmetic.rational(scale.numerator, scale.denominator)
                original = arithmetic.finite_pressure(n, time, charge)
                cancelled = cancelled_pressure(arithmetic, n, time, charge, side)
                assert max(original[0], cancelled[0]) <= min(original[1], cancelled[1])
                count += 1
    assert arithmetic.softplus_point(0) == arithmetic.one
    for numerator in (1 << 60, 1 << 80, 7 << 80):
        positive = arithmetic.softplus_point(numerator)
        negative = arithmetic.softplus_point(-numerator)
        assert negative == (positive[0] + numerator, positive[1] + numerator)
    return count


def source_hashes():
    directory = Path(__file__).resolve().parent
    names = ("charged_gaussian_floor.py", "check_small_charged_gaussian_bases.py",
             "all_digit_charged_gaussian_floor.py", Path(__file__).name)
    return {name: sha256((directory / name).read_bytes()).hexdigest() for name in names}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=root / "evidence" / "all-digit-charged-gaussian-floor.json")
    args = parser.parse_args()
    arithmetic = DigitDerivativeBounds(80)
    bases = []
    audit_nodes = 0
    for n, left, right, side in ((1, 0, 2, 0), (2, 0, 2, 0),
                                (2, 2, 4, 1), (3, 0, 2, 0)):
        function = lambda time, charge, n=n, side=side: cancelled_pressure(
            arithmetic, n, time, charge, side)[0]
        row = construct_cover(arithmetic, function,
                              (Fraction(1, 16), Fraction(16), Fraction(left), Fraction(right)))
        audit_nodes += audit_cover(arithmetic, function, row)
        row.update(n=n, physical_branch=side)
        bases.append(row)
    function = lambda time, charge: arithmetic.continuation(time, charge, terms=12)
    continuation = construct_cover(arithmetic, function,
                                    (Fraction(1), Fraction(16), Fraction(0), Fraction(2)))
    audit_nodes += audit_cover(arithmetic, function, continuation)

    # The independent partition audit must reject a missing leaf and an
    # inflated lower endpoint. Neither mutation is used in the proof record.
    malformed = dict(continuation, cover=continuation["cover"][1:])
    try:
        audit_cover(arithmetic, function, malformed)
    except AssertionError:
        pass
    else:
        raise AssertionError("Missing-leaf negative control was accepted")
    malformed_leaves = [list(row) for row in continuation["cover"]]
    malformed_leaves[0][1] = str(int(malformed_leaves[0][1]) + 1)
    try:
        audit_cover(arithmetic, function, dict(continuation, cover=malformed_leaves))
    except AssertionError:
        pass
    else:
        raise AssertionError("Inflated-bound negative control was accepted")

    record = {
        "status": "All digit-derivative compact bases have positive rigorous lower bounds.",
        "scope": "Compact certificates for the all-real-charge, all-real-digit-count ordinary proof; no factoring exponent improvement.",
        "source_base": "c2e951fa1f391edd97ba587d15ec13d3d6c8892d",
        "dyadic_precision": 80,
        "common_denominator": str(arithmetic.scale),
        "continuation_terms": 12,
        "positive_cover_rectangles": sum(row["positive_cover_rectangles"] for row in bases)
                                     + continuation["positive_cover_rectangles"],
        "independently_audited_tree_nodes": audit_nodes,
        "rejected_certificate_negative_controls": 2,
        "analytic_rational_controls": exact_analytic_controls(),
        "exact_signed_energy_profile_controls": exact_profile_controls(),
        "complete_expression_overlap_controls": exact_numerical_identity_controls(arithmetic),
        "source_sha256": source_hashes(),
        "bases": bases,
        "continuation": continuation,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, separators=(",", ":")) + "\n")
    summary = {key: value for key, value in record.items() if key not in ("bases", "continuation")}
    summary["bases"] = [{key: value for key, value in row.items() if key != "cover"} for row in bases]
    summary["continuation"] = {key: value for key, value in continuation.items() if key != "cover"}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
