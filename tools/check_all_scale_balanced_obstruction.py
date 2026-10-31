#!/usr/bin/env python3
"""Replay the algebraic inputs and analytic controls for all-scale exclusion.

The exact rational duals are certificates. The Decimal Gaussian and tail
checks are independent finite controls, not proofs of an asymptotic bound
or of an effective eventual interval threshold.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def exact_dual_replay(root: Path) -> dict:
    from check_rectangular_uniform_entry_bounds import verify_supported_duals
    from check_full_alphabet_leading_obstruction import dual_checks

    rectangular_path = root / "evidence/rectangular-full-leading-face-lp.json"
    witness_path = root / "evidence/rectangular-uniform-entry-inputs.json"
    full_path = root / "evidence/full-alphabet-leading-inputs.json"
    rectangular = json.loads(rectangular_path.read_text())
    witnesses = json.loads(witness_path.read_text())["witnesses"]
    full = json.loads(full_path.read_text())
    rect_result = verify_supported_duals(rectangular, witnesses)
    full_result = dual_checks(full)

    groups = [
        (rectangular["cases"], rectangular["C0_physical_gaps"], 1, [4, 3]),
        (full["cases"], full["gaps"], 2, [4, 4]),
    ]
    gap_records = []
    for cases, gaps, multiplier_bound, cardinalities in groups:
        for case in cases:
            chosen = (case["selected_x"], case["selected_y"])
            assert Fraction(case["optimum"]) >= (0 if chosen == (1, 1) else 2)
        for case in gaps:
            equality = dict(case["dual_equalities"])
            multiplier = Fraction(equality[2])
            assert -multiplier_bound <= multiplier <= -1
            assert Fraction(case["optimum"]) >= 2
        # If gap >= 2 tau - m*T - 4 ell and
        # T <= (2B+8)ell+b, its ell coefficient is (2m)B+(8m+4).
        polynomial_in_B = (2 * multiplier_bound, 8 * multiplier_bound + 4)
        assert polynomial_in_B == ((2, 12) if multiplier_bound == 1 else (4, 20))
        gap_records.append({
            "physical_cardinalities": cardinalities,
            "nonlow_states": sum(cardinalities) - 2,
            "maximum_off_face_multiplier_magnitude": multiplier_bound,
            "gap_loss_B_coefficient": polynomial_in_B[0],
            "gap_loss_constant_coefficient": polynomial_in_B[1],
        })

    return {
        "rectangular_exact_duals": rect_result["exact_dual_certificates"],
        "full_alphabet_exact_duals": full_result["exact_primal_dual_certificates"],
        "gap_reductions": gap_records,
        "input_sha256": {
            path.name: sha256(path.read_bytes()).hexdigest()
            for path in (rectangular_path, witness_path, full_path)
        },
        "scope": "Replays rational duals; the existing algebraic entry proofs are reused.",
    }


def log1p_positive(value: Decimal) -> Decimal:
    """Avoid losing tiny positive summands in log(1+value)."""
    assert value >= 0
    if value >= Decimal("0.1"):
        return (1 + value).ln()
    if not value:
        return value
    total = term = value
    power = value
    epsilon = Decimal("1e-70")
    for k in range(2, 10000):
        power *= value
        term = power / k
        total += -term if k % 2 == 0 else term
        if term <= value * epsilon:
            return total
    raise AssertionError("log1p series did not converge within its explicit cap")


def h(value: Decimal) -> Decimal:
    if value < 0:
        return -value + h(-value)
    return log1p_positive((-value).exp())


def gaussian_paired(length: int, ell: Decimal, root_two: Decimal) -> Decimal:
    assert length >= 2 and length % 2 == 0 and ell > 0
    first = sum((h((2 * j - 1) * ell) for j in range(1, length + 1)), Decimal())
    second = sum((h(root_two * (2 * j - 1) * ell)
                  for j in range(1, length // 2 + 1)), Decimal())
    return 4 * first - 4 * root_two * second


def gaussian_original(length: int, ell: Decimal, root_two: Decimal) -> Decimal:
    log_profile = length * length * ell + 4 * sum(
        (h((2 * j - 1) * ell) for j in range(1, length + 1)), Decimal())
    log_physical = 2 * sum(
        (h((4 * j - 2 * length + 2) * ell / root_two) for j in range(length)),
        Decimal())
    return log_profile - root_two * log_physical


def geometric_log_tail(first: Decimal, step: Decimal) -> tuple[Decimal, int]:
    """Sum log(1+exp(-first-k*step)) through its alternating geometric series."""
    assert first > 0 and step > 0
    initial = (-first).exp() / (1 - (-step).exp())
    total = Decimal()
    for k in range(1, 100000):
        term = (-k * first).exp() / (k * (1 - (-k * step).exp()))
        total += term if k % 2 else -term
        if term <= initial * Decimal("1e-65"):
            return total, k
    raise AssertionError("tail series did not converge within its explicit cap")


def close(left: Decimal, right: Decimal) -> None:
    assert abs(left - right) <= Decimal("1e-57") * max(Decimal(1), abs(left), abs(right))


def analytic_controls() -> dict:
    with localcontext() as context:
        context.prec = 90
        root_two = Decimal(2).sqrt()
        log_two = Decimal(2).ln()
        finite_cases = []
        tail_cases = []
        negative_controls = 0

        for length in (4, 8, 16, 64, 256):
            for tau_text in ("0.01", "0.1", "1", "4", "16", "64"):
                tau = Decimal(tau_text)
                ell = tau / length
                paired = gaussian_paired(length, ell, root_two)
                original = gaussian_original(length, ell, root_two)
                close(paired, original)
                lower = 2 * (2 - root_two) * length * h(2 * tau) - 4 * (1 + root_two) * log_two
                assert paired >= lower
                assert paired > 0
                finite_cases.append({
                    "M": length, "tau": tau_text,
                    "critical_pressure": str(paired),
                    "bounded_scale_lower_bound": str(lower),
                })

        for length in (6, 14, 62, 254):
            M = length + 2
            for tau_text in ("1", "4", "16", "64"):
                tau = Decimal(tau_text)
                ell = tau / M
                p_sum, p_terms = geometric_log_tail(root_two * (length + 1) * ell,
                                                  2 * root_two * ell)
                n_sum, n_terms = geometric_log_tail((2 * length + 1) * ell, 2 * ell)
                positive, negative = 4 * root_two * p_sum, 4 * n_sum
                p_lower = (4 * root_two * (-root_two * (length + 1) * ell).exp()
                           / ((1 - (-2 * root_two * ell).exp())
                              * (1 + (-root_two * (length + 1) * ell).exp())))
                n_upper = 4 * (-(2 * length + 1) * ell).exp() / (1 - (-2 * ell).exp())
                ratio_lower = ((((2 - root_two) * length + 1 - root_two) * ell).exp()
                               / (1 + (-root_two * (length + 1) * ell).exp()))
                scale_lower = (-root_two * (M - 1) * ell).exp() / ell
                assert positive >= p_lower >= scale_lower
                assert negative <= n_upper
                assert positive / negative >= ratio_lower
                coarse_exponent = (2 - root_two) * tau / 2
                assert ((2 - root_two) * length + 1 - root_two) * ell >= coarse_exponent
                if coarse_exponent.exp() > 4:
                    assert positive > 2 * negative
                # The 1/ell factor is substantial, not an innocuous fixed constant.
                if ell < Decimal("0.1"):
                    assert positive > (-root_two * tau).exp() / (2 * ell)
                    negative_controls += 1
                # Check the exact cancellation of the common tau exponent in (16).
                for C, D in ((Decimal(1), Decimal(0)), (Decimal(7), Decimal(9))):
                    error = C * (-root_two * tau + D * ell).exp()
                    ratio_upper = 2 * C * ell * ((D - root_two) * ell).exp()
                    assert error / (positive / 2) <= ratio_upper
                tail_cases.append({
                    "M": M, "tau": tau_text,
                    "positive_to_negative_tail_ratio": str(positive / negative),
                    "formula_lower_ratio": str(ratio_lower),
                    "positive_tail_to_scale_lower_bound": str(positive / scale_lower),
                    "series_terms": [p_terms, n_terms],
                })

        # Independently reconnect the finite pressure to its two omitted tails.
        identity_cases = 0
        for length, ell in ((4, Decimal("0.5")), (14, Decimal("0.25")),
                            (62, Decimal("0.125"))):
            first, unused = geometric_log_tail(ell, 2 * ell)
            second, unused = geometric_log_tail(root_two * ell, 2 * root_two * ell)
            infinity = 4 * first - 4 * root_two * second
            positive, unused = geometric_log_tail(root_two * (length + 1) * ell,
                                                  2 * root_two * ell)
            negative, unused = geometric_log_tail((2 * length + 1) * ell, 2 * ell)
            close(gaussian_paired(length, ell, root_two),
                  infinity + 4 * root_two * positive - 4 * negative)
            assert infinity > 0
            identity_cases += 1

        return {
            "arithmetic": "90-digit Decimal; finite diagnostic controls, not interval certificates",
            "finite_pressure_controls": len(finite_cases),
            "tail_controls": len(tail_cases),
            "finite_infinite_tail_identities": identity_cases,
            "small_scale_amplification_controls": negative_controls,
            "finite_cases": finite_cases,
            "tail_cases": tail_cases,
        }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--analytic-only", action="store_true",
                        help="Run analytic controls while skipping the separate exact dual replay.")
    args = parser.parse_args()
    result = {
        "status": "All-scale balanced-gauge obstruction: ordinary proof with exact inherited dual inputs.",
        "scope": "All positive scales for the two stated partial repairs, once M is sufficiently large.",
        "exact_duals": None if args.analytic_only else exact_dual_replay(root),
        "analytic_controls": analytic_controls(),
        "limitation": "No effective eventual M threshold or new factoring exponent is certified.",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256((root / "research/all-scale-balanced-obstruction.md").read_bytes()).hexdigest(),
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == "__main__":
    main()
