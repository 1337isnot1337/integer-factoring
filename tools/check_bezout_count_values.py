#!/usr/bin/env python3
"""Check fixed-field count values and the actual boundary parity provider."""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path

from balanced_bezout_modular import _unit_product
from bezout_count_values import (
    boundary_parity_provider, count_value_contract,
    smallest_prime_from_count_values, smallest_prime_from_count_value_windows,
    value_from_reversed_jet,
)
from check_leading_bezout_factoring import factorial_odd_unit, literal_binomial_mod
from check_pole_cauchy_conditioning import mod_pair
from check_pole_smith_stabilization import Quad
from check_regularized_bezout_factor_moments import dense_regularized_jet
from interval_unit_jets import QuadraticJets


def known_factor_value(contract, proper_divisors):
    n, first, count = contract["n"], contract["first_a"], contract["rows"]
    precision, exponent = contract["precision"], contract["size_exponent"]
    modulus = 1 << precision
    common_square = (factorial_odd_unit(exponent, precision) ** 2
                     if precision >= 3 else 1)
    ring = QuadraticJets(precision, 0)
    pp0 = _unit_product(ring, (0, count), ring.constant(-first), ring.constant(-2))[0]
    assert pp0 == (1, 0)
    normalized_count = 0
    visits = 0
    for a in proper_divisors:
        if a < first:
            continue
        b = n // a
        i, j = (a - first) // 2, (b - 3) // 2
        relative = (pow(a, -count, modulus) * (-1) ** (i + j)
                    * literal_binomial_mod(exponent, i, precision)
                    * literal_binomial_mod(exponent, j, precision)) % modulus
        assert relative == 1
        normalized_count = (normalized_count + relative) % modulus
        visits += 1
    assert normalized_count == visits % modulus
    return -pp0[0] * pow(common_square, -1, modulus) * normalized_count % modulus, visits


def wrapper_controls():
    cases = calls = weights = window_cases = window_calls = 0
    examples = []
    for n in list(range(3, 258, 2)) + [17 * 19, 101 * 103, 1009 * 1013]:
        proper = [a for a in range(3, isqrt(n) + 1, 2) if n % a == 0]
        proper = sorted(set(proper + [n // a for a in proper]))
        for trial_bits in (1, 2, 4):
            def oracle(contract):
                nonlocal calls, weights
                calls += 1
                value, support = known_factor_value(contract, proper)
                weights += support
                return value
            result = smallest_prime_from_count_values(n, oracle, trial_bits=trial_bits)
            expected = min(proper) if proper else n
            assert result["smallest_prime_factor"] == expected
            if n in (45, 169, 323, 10403, 1022117) and trial_bits in (1, 4):
                examples.append({"n": n, "trial_bits": trial_bits, **result})
            cases += 1
            def window_oracle(contract):
                nonlocal window_calls
                window_calls += 1
                return known_factor_value(contract, proper)[0]
            window_result = smallest_prime_from_count_value_windows(n, window_oracle, trial_bits=trial_bits)
            assert window_result["smallest_prime_factor"] == expected
            window_cases += 1
    return {"roughness_wrapper_controls": cases, "known_factor_value_queries": calls,
            "prime_window_wrapper_controls": window_cases,
            "known_factor_window_value_queries": window_calls,
            "exact_flat_inverse_weight_controls": weights, "selected_controls": examples}


def dense_controls():
    values = parities = 0
    for index in range(48):
        n = (9, 15, 25, 35, 45, 77)[index % 6]
        rows, columns = 1 << (index % 4), 1 << ((index // 4) % 4)
        a0, b0 = 1 + 2 * (index % 3), 1 + 2 * ((index // 3) % 3)
        precision = (1, 2, 4)[index % 3]
        dense = dense_regularized_jet(n, a0, rows, b0, columns, precision)
        value = sum(((-1) ** j * coefficient for j, coefficient in enumerate(dense["v"])),
                    start=Quad())
        residue = mod_pair(value, precision)
        assert residue[1] == 0
        assert value_from_reversed_jet(rows, precision, dense["v_jet"]) == residue[0]
        parity = boundary_parity_provider(n, a0, rows, b0, columns, boundary_limit=rows + columns)
        assert parity["leading_coefficient_mod2"] == dense["v_jet"][0] % 2
        assert parity["value_at_zero_mod2"] == residue[0] % 2
        assert parity["interior_nodes_enumerated"] == 0
        expected_count = sum(a * b == n for a in dense["rows"] for b in dense["columns"])
        assert parity["leading_coefficient_mod2"] == expected_count % 2
        values += 1
        parities += 1
    return {"dense_exact_value_and_jet_controls": values,
            "dense_boundary_provider_controls": parities}


def implicit_provider_controls():
    root = (1 << 150) + 3
    square = boundary_parity_provider(root * root, 3, 1 << 200, 3, 1 << 200, boundary_limit=0)
    nonsquare = boundary_parity_provider(root * root + 2, 3, 1 << 200, 3, 1 << 200, boundary_limit=0)
    count = 1 << 300
    shifted = boundary_parity_provider(3 * (2 * count + 7), 3, count, 9, count, boundary_limit=6)
    assert square["leading_coefficient_mod2"] == 1
    assert nonsquare["leading_coefficient_mod2"] == 0
    assert shifted["leading_coefficient_mod2"] == 1
    assert shifted["raw_boundary_nodes"] == 6
    assert shifted["visited_boundary_nodes"] == 3
    try:
        boundary_parity_provider((1 << 400) + 1, 3, 1 << 200,
                                 (1 << 201) + 3, 1 << 200, boundary_limit=10)
    except ValueError as error:
        assert "boundary" in str(error)
    else:
        raise AssertionError("The excessive boundary was not rejected")
    def brief(data):
        return {"input_bits": data["contract"]["n"].bit_length(),
                "row_count_bits": data["contract"]["rows"].bit_length(),
                "column_count_bits": data["contract"]["columns"].bit_length(),
                "value_at_zero_mod2": data["value_at_zero_mod2"],
                "leading_coefficient_mod2": data["leading_coefficient_mod2"],
                "raw_boundary_nodes": data["raw_boundary_nodes"],
                "visited_boundary_nodes": data["visited_boundary_nodes"],
                "interior_nodes_enumerated": data["interior_nodes_enumerated"],
                "scope": data["scope"]}
    return {"actual_implicit_parity_providers": [brief(x) for x in (square, nonsquare, shifted)],
            "rejected_excessive_boundary_controls": 1}


def zero_suffix_control():
    n = 5 * 3 ** 16
    support = sorted(3 ** i * 5 ** j for i in range(17) for j in range(2)
                     if 1 < 3 ** i * 5 ** j < n)
    queried = {}
    def explicit_known_factor_oracle(contract):
        precision = contract["precision"]
        common = factorial_odd_unit(contract["size_exponent"], precision)
        count = sum(a >= contract["first_a"] for a in support)
        value = -pow(common * common, -1, 1 << precision) * count % (1 << precision)
        queried[contract["first_a"]] = value
        return value
    result = smallest_prime_from_count_value_windows(n, explicit_known_factor_oracle)
    assert result["smallest_prime_factor"] == 3
    assert queried[3] == 0 and queried[5] != 0
    return {"n": n, "precision": result["precision"], "proper_divisor_count": len(support),
            "initial_nonempty_suffix_value": queried[3], "next_suffix_value": queried[5],
            "smallest_prime_factor": result["smallest_prime_factor"],
            "scope": "Known-factor oracle control; nonempty suffix wraps to zero, but its prime-window difference detects the factor"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = {"theorem": "Fixed-field count values and a boundary parity provider",
               **wrapper_controls(), **dense_controls(), **implicit_provider_controls(),
               "wrapped_suffix_prime_window_control": zero_suffix_control(),
               "scope": "General count values are supplied by finite known-factor test oracles. The bounded-symmetric-difference precision-one provider is implemented and exercised on implicit large domains. No improved complete factoring score."}
    large = count_value_contract((1 << 401) + 45, 3, 22)
    receipt["large_general_value_contract"] = {"input_bits": large["n"].bit_length(),
        "precision": large["precision"], "domain_size_exponent": large["size_exponent"],
        "padding_exponent": large["padding_exponent"],
        "scope": "Contract construction only; no general large value computed"}
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
