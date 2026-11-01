#!/usr/bin/env python3
"""Finite supplied-residue controls; no fast weighted evaluator is used."""

import json

from check_divisor_moment_recovery import divisors
from sparse_divisor_residues import recover_from_global_residues


def polynomial_from_divisors(values):
    polynomial = [1]
    for value in values:
        result = [0]*(len(polynomial)+1)
        for index, coefficient in enumerate(polynomial):
            result[index] -= value*coefficient
            result[index+1] += coefficient
        polynomial = result
    return polynomial


def check_one(integer, values, characteristic, slack):
    queries = []

    def oracle(prime):
        queries.append(prime)
        output = [0]*prime
        for divisor in values:
            residue = divisor % prime
            output[residue] = (output[residue]+1) % characteristic
        return output

    result = recover_from_global_residues(integer, characteristic, len(values)+slack, oracle)
    assert result["smallest_prime_factor"] == values[1]
    assert result["divisor_count"] == len(values)
    assert result["divisor_polynomial"] == polynomial_from_divisors(values)
    assert len(queries) == result["query_count"] and max(queries) == result["largest_cyclic_degree"]
    return result


def main():
    cases = queries = cancellations = 0
    maximum_degree = 0
    for integer in (2, 3, 4, 6, 8, 9, 12, 15, 21, 25, 30, 36, 60, 105):
        values = divisors(integer)
        for characteristic in (2, 3, 5):
            for slack in (0, 2):
                result = check_one(integer, values, characteristic, slack)
                cases += 1
                queries += result["query_count"]
                cancellations += result["residue_support_histogram"].get(0, 0)
                maximum_degree = max(maximum_degree, result["largest_cyclic_degree"])
    large = []
    for p, q in ((1000000007, 1000000009), (824633708473, 1099511627689)):
        result = check_one(p*q, [1, p, q, p*q], 2, 0)
        result.pop("divisor_polynomial")
        large.append({"integer": p*q, **result})
    print(json.dumps({"scope": "Garg--Schost recovery from supplied global residues, not a fast evaluator",
                      "finite_provider": "explicit small divisor enumeration and already trial-proved large primes",
                      "small_recoveries": cases, "small_residue_queries": queries,
                      "zero_residues_from_characteristic_cancellation": cancellations,
                      "maximum_small_cyclic_degree": maximum_degree,
                      "large_supplied_divisor_controls": large}, indent=2))


if __name__ == "__main__":
    main()
