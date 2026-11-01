#!/usr/bin/env python3
"""Exact cyclotomic conversion and independent signed-support recovery."""

import json

from check_divisor_moment_recovery import divisors
from odd_character_residues import (
    character_bit_decoder, cyclotomic_polynomial, multiply_cyclotomic, primitive_generator,
    root_of_unity_powers, signed_residue_from_character_bits, signed_residue_from_odd_characters,
)
from sparse_divisor_residues import recover_from_global_residues


def check_character_conversion():
    cases = zero = queries = 0
    modular_cases = modular_queries = 0
    binary_cases = binary_queries = 0
    precision_bits = set()
    inputs = [(n, divisors(n)) for n in range(2, 65)]
    for p, q in ((1000000007, 1000000009), (824633708473, 1099511627689)):
        inputs.append((p*q, [1, p, q, p*q]))
    for prime in (3, 5, 7, 11, 13, 17, 19, 23):
        generator = primitive_generator(prime)
        order = prime-1
        modulus, powers = root_of_unity_powers(order)
        assert modulus == cyclotomic_polynomial(order)
        logs = {pow(generator, exponent, prime): exponent for exponent in range(order)}
        assert len(logs) == order
        for integer, values in inputs:

            def coefficient(supplied_generator, index):
                nonlocal queries
                assert supplied_generator == generator and index % 2
                queries += 1
                return tuple(sum(powers[index*logs[d % prime] % order][i]
                                 for d in values if d % prime)
                             for i in range(len(modulus)-1))

            # A separate residue-count calculation avoids all Fourier algebra.
            direct = [0]*prime
            for divisor in values:
                direct[divisor % prime] += 1
                direct[-divisor % prime] -= 1
            for characteristic in (2, 3):
                recovered = signed_residue_from_odd_characters(prime, coefficient, characteristic)
                assert recovered == tuple(value % characteristic for value in direct)
                cases += 1
                zero += not any(recovered)
            coefficient_modulus = order & -order
            precision_bits.add(coefficient_modulus.bit_length()-1)
            before = queries

            def modular_coefficient(supplied_generator, index):
                return tuple(value % coefficient_modulus
                             for value in coefficient(supplied_generator, index))

            low_precision = signed_residue_from_odd_characters(
                prime, modular_coefficient, coefficient_modulus=coefficient_modulus)
            assert low_precision == tuple(value % 2 for value in direct)
            modular_cases += 1
            modular_queries += queries-before
            before = queries

            def coefficient_bits(supplied_generator, index):
                return tuple(value & 1 for value in coefficient(supplied_generator, index))

            binary = signed_residue_from_character_bits(prime, coefficient_bits)
            assert binary == tuple(value % 2 for value in direct)
            binary_cases += 1
            binary_queries += queries-before
        # Independent character multiplication and odd parity controls.
        for left in range(1, prime):
            for right in range(1, prime):
                product = multiply_cyclotomic(powers[logs[left]], powers[logs[right]], modulus)
                assert product == powers[logs[left*right % prime]]
        assert powers[order//2][0] == -1 and not any(powers[order//2][1:])
    return {"prime_levels": [3, 5, 7, 11, 13, 17, 19, 23],
            "conversion_cases": cases, "supplied_weight_one_coefficients": queries,
            "exact_coefficient_queries": queries-modular_queries-binary_queries,
            "completely_zero_signed_residues": zero,
            "finite_precision_cases": modular_cases,
            "finite_precision_coefficient_queries": modular_queries,
            "two_adic_precisions_tested": sorted(precision_bits),
            "single_bit_cases": binary_cases,
            "single_bit_coefficient_queries": binary_queries,
            "single_bit_decoder_ranks": [[q, len(character_bit_decoder(q)[3])]
                                         for q in (3, 5, 7, 11, 13, 17, 19, 23)]}


def check_signed_recovery():
    cases = queries = zero = 0
    records = []
    inputs = [(n, divisors(n)) for n in (2, 4, 6, 9, 15, 25, 30, 36, 105)]
    inputs.append((1000000016000000063, [1, 1000000007, 1000000009, 1000000016000000063]))
    for integer, values in inputs:

        def oracle(prime):
            output = [0]*prime
            for divisor in values:
                output[divisor % prime] ^= 1
                output[-divisor % prime] ^= 1
            return output

        result = recover_from_global_residues(integer, 2, len(values), oracle, signed=True)
        expected = [1]
        for divisor in values:
            product = [0]*(len(expected)+1)
            for i, value in enumerate(expected):
                product[i] -= divisor*divisor*value
                product[i+1] += value
            expected = product
        assert result["divisor_polynomial"] == expected
        assert result["divisor_count"] == len(values)
        assert result["smallest_prime_factor"] == values[1]
        cases += 1
        queries += result["query_count"]
        zero += result["residue_support_histogram"].get(0, 0)
        if integer.bit_length() > 50:
            result.pop("divisor_polynomial")
            records.append({"integer": integer, **result})
    return {"recoveries": cases, "supplied_signed_residue_queries": queries,
            "completely_zero_residues": zero, "large_control": records}


def main():
    print(json.dumps({"scope": "Exact weight-one character conversion and supplied signed residues; no fast coefficient engine",
                      "finite_provider": "explicit divisor lists; character coefficients computed from them",
                      "character_conversion": check_character_conversion(),
                      "signed_recovery": check_signed_recovery()}, indent=2))


if __name__ == "__main__":
    main()
