#!/usr/bin/env python3
"""Finite oracle controls; the supplied divisor lists are explicitly charged."""

from math import isqrt
import json

from divisor_moment_recovery import factor_from_odd_moments, factor_with_odd_moments


def divisors(integer):
    small = [d for d in range(1, isqrt(integer)+1) if integer % d == 0]
    return sorted(set(small+[integer//d for d in small]))


def expected_polynomial(values):
    polynomial = [1]
    for divisor in values:
        result = [0]*(len(polynomial)+1)
        for i, coefficient in enumerate(polynomial):
            result[i] -= divisor*divisor*coefficient
            result[i+1] += coefficient
        polynomial = result
    return polynomial


def check_one(integer, values, bound):
    queries = []

    def moment(index):
        queries.append(index)
        return sum(divisor**(2*index+3) for divisor in values)

    result = factor_from_odd_moments(integer, bound, moment)
    assert result["smallest_prime_factor"] == values[1]
    assert result["divisor_count"] == len(values)
    assert result["node_polynomial"] == expected_polynomial(values)
    assert queries == list(range(2*bound))
    return result


def main():
    cases = queries = maximum_degree = 0
    for integer in list(range(2, 257))+[360, 900, 2310]:
        values = divisors(integer)
        for slack in (0, 2):
            result = check_one(integer, values, len(values)+slack)
            queries += result["moment_queries"]
            maximum_degree = max(maximum_degree, result["divisor_count"])
            cases += 1
    wrapper_cases = wrapper_queries = trial_cases = 0
    for integer in range(2, 257):
        values = divisors(integer)
        result = factor_with_odd_moments(integer, lambda j: sum(d**(2*j+3) for d in values))
        assert result["smallest_prime_factor"] == values[1]
        wrapper_cases += 1
        wrapper_queries += result["moment_queries"]
        trial_cases += result["stage"] == "trial_division"
    large = []
    for p, q in ((1000000007, 1000000009), (824633708473, 1099511627689)):
        # This expensive finite-provider work is charged explicitly. The
        # recovery routine never receives p or q and never performs this scan.
        assert all(p % divisor for divisor in range(2, isqrt(p)+1))
        assert all(q % divisor for divisor in range(2, isqrt(q)+1))
        integer = p*q
        result = check_one(integer, [1, p, q, integer], 8)
        result.pop("node_polynomial")
        large.append({"integer": integer, **result})
    print(json.dumps({"scope": "Deterministic recovery from supplied exact moments; no fast coefficient oracle",
                      "finite_oracle": "divisors enumerated for small N; supplied large factors independently trial-proved prime",
                      "core_recoveries": cases, "core_moment_queries": queries,
                      "largest_hankel_rank": maximum_degree,
                      "uniform_wrapper_cases": wrapper_cases,
                      "uniform_wrapper_trial_returns": trial_cases,
                      "uniform_wrapper_moment_queries": wrapper_queries,
                      "large_supplied_divisor_controls": large}, indent=2))


if __name__ == "__main__":
    main()
