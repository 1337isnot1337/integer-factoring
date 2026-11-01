#!/usr/bin/env python3
"""Exact finite recurrence certificates for the relaxed piece interface.

These checks concern the declared rounded volume inequalities, not an
implementation of new crossing pieces or a faster factoring algorithm.
"""

from fractions import Fraction
import json


def ceil_fraction(value: Fraction) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def schedule_logs(maximum_level: int, x0: int = 18) -> list[int]:
    values = [x0]
    for level in range(maximum_level):
        values.append(2 * values[-1] + 6 * level + 2 * x0)
    for level, value in enumerate(values):
        assert value == (3 * x0 + 6) * (1 << level) - 6 * level - 2 * x0 - 6
        assert value >= x0 * (1 << level)
        if level:
            assert values[level - 1] * 2 <= value <= values[level - 1] * 10
    return values


def verify() -> dict:
    logs = schedule_logs(6)
    rho_values = (Fraction(1, 2), Fraction(2, 3), Fraction(3, 4), Fraction(7, 8))
    paths = rounded_steps = lower_bounds = early_cutoffs = entropy_sums = 0
    product_bounds = []

    for level in range(1, 7):
        root_log = logs[level]
        input_log = 40 * root_log
        n = 1 << input_log
        cutoff = 1 << (20 * root_log)
        label_lower = 1 << (18 * root_log)
        model_slack = Fraction(1, level + 1)
        rounding_allowance = Fraction(1, 100 * (level + 1))
        uniform_slack = model_slack + rounding_allowance
        assert Fraction(1, 1 << (19 * root_log)) <= rounding_allowance

        for rho in rho_values:
            envelope = Fraction(1)
            for previous in range(level):
                zeta = Fraction(1, 1 << (previous + 5))
                envelope *= 1 + (3 * zeta + uniform_slack) / rho
            assert envelope < 27
            product_bounds.append({
                "level": level,
                "rho": str(rho),
                "envelope_numerator": envelope.numerator,
                "envelope_denominator": envelope.denominator,
            })

            for forced_cut in sorted({0, level, max(1, level // 2)}):
                current_level, volume = level, n
                visited_sum = 0
                query_charge = 0
                search_charge = 0
                while current_level > 0 and volume > cutoff:
                    parent_log = logs[current_level]
                    child_log = logs[current_level - 1]
                    parent_squared = 1 << (2 * parent_log)
                    child_squared = 1 << (2 * child_log)
                    assert volume % parent_squared == 0
                    label = volume // parent_squared
                    assert label > label_lower
                    lower_bounds += 1

                    # A uniform O(lambda/loglog lambda) slack covers this model.
                    loglog_upper = (label.bit_length() + 1).bit_length()
                    assert 4 * model_slack * loglog_upper <= 64

                    zeta = Fraction(1, 1 << (current_level + 4))
                    raw = (rho + 3 * zeta + model_slack) * volume
                    child = child_squared * ceil_fraction(raw / child_squared)
                    assert child <= (rho + 3 * zeta + uniform_slack) * volume
                    assert child_squared <= 1 << parent_log
                    assert Fraction(child_squared, volume) <= rounding_allowance
                    rounded_steps += 1

                    if current_level == forced_cut:
                        child = child_squared * max(
                            1, (cutoff + 2 * child_squared - 1) // (2 * child_squared))
                        assert child <= cutoff
                        assert child <= raw + child_squared

                    assert volume <= envelope * n * rho ** (level - current_level)
                    visited_sum += volume
                    query_charge += volume // (1 << child_log)
                    search_charge += (1 << (parent_log // 2)) * child_log
                    volume, current_level = child, current_level - 1

                assert volume <= envelope * n * rho ** (level - current_level)
                assert visited_sum <= envelope * n / (1 - rho)
                assert query_charge <= envelope * n * rho ** level / ((1 << 17) - 1)
                assert search_charge <= level * (1 << (root_log // 2)) * root_log
                assert search_charge <= cutoff
                assert volume <= envelope * n * rho ** level + cutoff
                if current_level > 0:
                    assert volume <= cutoff
                    early_cutoffs += 1
                for delta in (1, 2):
                    branching_charge = Fraction(visited_sum, 4 * input_log ** delta)
                    assert branching_charge <= (
                        envelope * n / (4 * (1 - rho) * input_log ** delta))
                    entropy_sums += 1
                paths += 1

    counterexamples = 0
    for inverse_target in range(1, 9):
        previous = Fraction(0)
        for power in range(4, 33):
            t = inverse_target * power
            ratio = Fraction(1 << power, t)
            assert ratio > previous
            assert (1 << power) >= power * power
            assert ratio >= Fraction(power, inverse_target)
            previous = ratio
            counterexamples += 1
    assert Fraction(1, 256) < Fraction(1, 100)
    assert Fraction(1 << 32, 256) > 10 ** 6

    return {
        "scope": "exact declared-volume recurrence and cost-envelope checks; no new crossing family or factoring exponent",
        "rational_schedule_paths": paths,
        "rounded_grid_steps": rounded_steps,
        "explicit_label_lower_bound_checks": lower_bounds,
        "early_cutoff_paths": early_cutoffs,
        "geometric_entropy_sum_checks": entropy_sums,
        "subexponential_entropy_counterexample_cases": counterexamples,
        "uniform_finite_envelope_bound": "27",
        "label_lower_bound": "lambda > n^(9/20) at every grid step",
        "local_query_series_bound": "1/(2^17-1)",
        "counterexample": {
            "target_exponent": "1/8",
            "entropy_fraction": "1/256",
            "ratio_to_target": str(Fraction(1 << 32, 256)),
        },
        "product_envelopes": product_bounds,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
