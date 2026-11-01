#!/usr/bin/env python3
"""Finite checks for direct exact counts and their shorter residue modulus.

This harness uses explicit small divisor sets as a reference. It does not
implement the asymptotic fast contraction algorithm or claim its memory
bound from finite experiments.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import comb, gcd, lcm, prod
from pathlib import Path
import random

from integral_schur_queries import schur_mode, translated_angle_mode
from check_partition_channels import partitions


def divisor_bound_bits(n):
    if n < 4:
        return n
    logarithm = n.bit_length()-1
    return min(n, (9*n+logarithm-1)//logarithm)


def prefix_matches_endpoint(value, endpoint, changed_bit, bits):
    """The endpoint and value are read-only; only positions are retained."""
    for position in range(changed_bit, bits):
        expected = 0 if position == changed_bit else (endpoint >> position) & 1
        if (value >> position) & 1 != expected:
            return False
    return True


def direct_below(divisors, endpoint, bits, modulus):
    accumulator = queries = 0
    for changed_bit in range(bits-1, -1, -1):
        if not (endpoint >> changed_bit) & 1:
            continue
        count = sum(prefix_matches_endpoint(value, endpoint, changed_bit, bits)
                    for value in divisors)
        residue = count % modulus
        assert residue == count
        accumulator += residue
        queries += 1
    assert accumulator == sum(value < endpoint for value in divisors)
    assert queries <= bits
    return accumulator, queries


def count_controls():
    rng = random.Random(130210)
    prefix_queries = interval_queries = intervals = short_modulus_cases = 0
    for number in range(2, 258):
        bits = number.bit_length()
        q = divisor_bound_bits(bits)
        modulus = 1 << (q+1)
        divisors = [value for value in range(1, number+1) if number % value == 0]
        assert len(divisors) <= 1 << q < modulus
        for remaining in range(bits+1):
            counts = {}
            for value in divisors:
                prefix = value >> remaining
                counts[prefix] = counts.get(prefix, 0)+1
            for prefix in range(1 << (bits-remaining)):
                exact = counts.get(prefix, 0)
                assert exact % modulus == exact
                prefix_queries += 1
        endpoints = list(range(1 << bits)) if number <= 32 else [
            0, 1, number-1, number, (1 << bits)-1,
            *(rng.randrange(1 << bits) for _ in range(4))]
        for endpoint in endpoints:
            _, queries = direct_below(divisors, endpoint, bits, modulus)
            interval_queries += queries
        for _ in range(4):
            left, right = sorted(rng.randrange(1 << bits) for _ in range(2))
            below_left, first_queries = direct_below(divisors, left, bits, modulus)
            below_right, second_queries = direct_below(divisors, right, bits, modulus)
            singleton = sum(value == right for value in divisors) % modulus
            exact = below_right+singleton-below_left
            assert exact == sum(left <= value <= right for value in divisors)
            assert first_queries+second_queries+1 <= 2*bits+1
            intervals += 1
            interval_queries += first_queries+second_queries+1

    long_cases = []
    for bits in (1024, 2048, 4096):
        # N=2^(bits-1) has exactly bits divisors. No factor search is needed.
        q = divisor_bound_bits(bits)
        modulus = 1 << (q+1)
        assert q < bits and bits <= 1 << q < modulus
        short_modulus_cases += 1
        tested = 0
        for remaining in (0, 1, bits//3, bits//2, bits-1, bits):
            prefixes = (0,) if remaining == bits else (0, 1, 2, 3,
                                                       1 << (bits-remaining-1))
            for prefix in prefixes:
                if prefix >= 1 << (bits-remaining):
                    continue
                exact = sum((1 << exponent) >> remaining == prefix
                            for exponent in range(bits))
                expected = (remaining if prefix == 0 else
                            int(prefix & (prefix-1) == 0
                                and prefix.bit_length()-1+remaining < bits))
                assert exact == expected and exact % modulus == exact
                tested += 1
        long_cases.append({"input_bits": bits, "count_bound_bits": q,
                           "ring_element_bits": q+1, "known_divisor_count": bits,
                           "prefix_queries": tested})
    assert 6 % 2 == 0 and 6 % (1 << (divisor_bound_bits(4)+1)) == 6
    return {"small_prefix_queries": prefix_queries,
            "read_only_interval_queries": interval_queries,
            "exact_interval_counts": intervals,
            "strictly_shorter_modulus_cases": short_modulus_cases,
            "long_known_factorization_cases": long_cases,
            "negative_small_modulus_erases_a_count": True}


def integral_crossing_controls():
    # These caches are confined to the finite harness.
    mode = lru_cache(None)(schur_mode)
    angle = lru_cache(None)(translated_angle_mode)
    weight = 2
    basis = partitions(weight)
    moduli = (4, 8, 16, 64)
    pieces = reductions = crossings = 0
    for x, y, v, u in product(basis, repeat=4):
        totals = [0]*len(moduli)
        for shift in range(-2*weight, 2*weight+1):
            row_value = column_value = 0
            for index in range(2*weight+1):
                for sigma in partitions(shift+index):
                    row_value += ((-1)**index * angle(x, index, sigma, y)
                                  * mode(v, u, sigma))
                for sigma in partitions(2*weight-shift-index):
                    column_value += ((-1)**index * angle(x, index, u, sigma)
                                     * mode(v, sigma, y))
            assert row_value == column_value
            for position, modulus in enumerate(moduli):
                assert row_value % modulus == column_value % modulus
                totals[position] = (totals[position]+row_value) % modulus
                reductions += 1
            pieces += 1
        assert totals == [int(x == v and y == u)]*len(moduli)
        crossings += len(moduli)
    assert 2*2 % 4 == 0
    return {"both_order_integer_pieces": pieces,
            "power_of_two_piece_reductions": reductions,
            "power_of_two_crossing_entries": crossings,
            "moduli_with_zero_divisors": list(moduli)}


def rational_controls():
    checks = vanished = 0
    for denominators in ((2,), (2, 3), (2, 3, 5), (4, 9, 25), (1,)):
        denominator_product = prod(denominators)
        for q in (1, 2, 5, 17, 512):
            modulus = (1 << (q+1))*denominator_product+1
            assert modulus > 1 << q
            for denominator in denominators:
                assert gcd(denominator, modulus) == 1
                value = Fraction(modulus, denominator)
                complement = 1-value
                reduced = value.numerator*pow(value.denominator, -1, modulus) % modulus
                other = complement.numerator*pow(complement.denominator, -1, modulus) % modulus
                assert value != 0 and reduced == 0 and (reduced+other) % modulus == 1
                vanished += 1
                for exact_count in (0, 1, 1 << q):
                    assert exact_count % modulus == exact_count
                    checks += 1
    return {"rational_unit_and_count_checks": checks,
            "vanishing_rational_terms_preserving_one": vanished}


def generalized_crt(residue, modulus, incoming, incoming_modulus):
    common = gcd(modulus, incoming_modulus)
    difference = (incoming-residue) % incoming_modulus
    if difference % common:
        raise ValueError("incompatible residues")
    quotient = incoming_modulus//common
    if quotient == 1:
        return residue, modulus
    inverse = pow((modulus//common) % quotient, -1, quotient)
    step = ((difference//common)*inverse) % quotient
    result = residue+modulus*step
    combined = modulus*quotient
    assert 0 <= result < combined
    assert result % modulus == residue
    assert result % incoming_modulus == incoming % incoming_modulus
    return result, combined


def lcm_supply_controls():
    binomial_divisibilities = supplies = 0
    full_lcm = 1
    for limit in range(1, 97):
        full_lcm = lcm(full_lcm, limit)
        for index in range(limit):
            assert full_lcm % (limit*comb(limit-1, index)) == 0
            binomial_divisibilities += 1
        assert full_lcm >= 1 << (limit-1)
    for denominator_product in range(1, 65):
        a = denominator_product.bit_length()+1
        for q in range(1, 13):
            limit = 4*(q+a*a+4)
            full_lcm = good_lcm = 1
            for modulus in range(1, limit+1):
                full_lcm = lcm(full_lcm, modulus)
                if gcd(modulus, 2*denominator_product) == 1:
                    good_lcm = lcm(good_lcm, modulus)
            assert full_lcm % good_lcm == 0
            assert full_lcm//good_lcm <= limit**a
            assert (1 << (limit-1)) > (1 << q)*limit**a
            assert good_lcm > 1 << q
            supplies += 1
    return {"binomial_lcm_divisibilities": binomial_divisibilities,
            "exact_small_modulus_supply_bounds": supplies}


def reconstruct_reference_count(exact, q, denominator_product):
    assert 0 <= exact <= 1 << q
    a = denominator_product.bit_length()+1
    limit = 4*(q+a*a+4)
    residue, combined = 0, 1
    queries = overlaps = skipped = 0
    for modulus in range(3, limit+1):
        if gcd(modulus, 2*denominator_product) != 1:
            skipped += 1
            continue
        old = combined
        common = gcd(old, modulus)
        overlaps += int(common > 1 and modulus//common > 1)
        # This is the exact residue returned by the proved count evaluator.
        incoming = exact % modulus
        residue, combined = generalized_crt(residue, combined, incoming, modulus)
        # A nonzero rational term may vanish without affecting the identity.
        zero_term = Fraction(modulus, denominator_product)
        complement = exact-zero_term
        first = zero_term.numerator*pow(zero_term.denominator, -1, modulus) % modulus
        second = complement.numerator*pow(complement.denominator, -1, modulus) % modulus
        assert first == 0 and (first+second) % modulus == incoming
        queries += 1
        assert combined <= (1 << q)*limit
        if combined > 1 << q:
            break
    assert combined > 1 << q and residue == exact
    assert combined.bit_length() <= q+limit.bit_length()+1
    return queries, overlaps, skipped


def crt_controls():
    reconstructions = queries = overlaps = skipped = 0
    for denominator_product in (1, 2, 3, 6, 30, 210, 2310):
        for q in range(1, 9):
            for exact in range((1 << q)+1):
                counts = reconstruct_reference_count(exact, q, denominator_product)
                queries += counts[0]
                overlaps += counts[1]
                skipped += counts[2]
                reconstructions += 1
    for q in (64, 256, 1024):
        counts = reconstruct_reference_count(1 << q, q, 30)
        queries += counts[0]
        overlaps += counts[1]
        skipped += counts[2]
        reconstructions += 1
    composite_merges = 0
    for exact in range(513):
        residue, combined = 0, 1
        for modulus in (9, 15, 21, 25):
            residue, combined = generalized_crt(residue, combined, exact % modulus, modulus)
            composite_merges += 1
        assert combined == 1575 and residue == exact
    try:
        generalized_crt(0, 9, 1, 15)
    except ValueError:
        pass
    else:
        raise AssertionError("incompatible residues were accepted")
    assert overlaps > 0 and skipped > 0
    number_field_checks = 0
    for modulus in range(3, 102, 2):
        if gcd(modulus, 10) != 1:
            continue
        # In (Z/jZ)[X]/(X^2-2), (X/5)*(5X/2)=1.
        first = pow(5, -1, modulus)
        second = 5*pow(2, -1, modulus) % modulus
        constant = 2*first*second % modulus
        assert constant == 1
        for exact in (0, 1, 7, 31):
            assert (exact*constant % modulus, 0) == (exact % modulus, 0)
            number_field_checks += 1
    return {"exhaustive_and_large_exact_reconstructions": reconstructions,
            "direct_residue_queries": queries,
            "nontrivial_overlapping_merges": overlaps,
            "nonunit_moduli_skipped": skipped,
            "explicit_overlapping_composite_merges": composite_merges,
            "incompatible_residue_negative_control": True,
            "fixed_quadratic_ring_constant_residues": number_field_checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"scope": "Direct count interface controls; not a full asymptotic evaluator implementation",
              "counts_and_read_only_prefixes": count_controls(),
              "integral_schur_ring_controls": integral_crossing_controls(),
              "rational_modulus_controls": rational_controls(),
              "elementary_lcm_supply": lcm_supply_controls(),
              "small_modulus_crt": crt_controls(),
              "least_factor_prefix_self_reduction_not_claimed_sublinear": True}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass", **result["counts_and_read_only_prefixes"],
                      **result["integral_schur_ring_controls"],
                      **result["rational_modulus_controls"],
                      **result["elementary_lcm_supply"],
                      **result["small_modulus_crt"]}))


if __name__ == "__main__":
    main()
