#!/usr/bin/env python3
"""Independent exact controls for the polynomial character decoder."""

from hashlib import sha256
import json
from pathlib import Path
from random import Random

from cyclotomic_character_decoder import (
    binary_cyclotomic, decode_character_polynomial, parameters, polynomial_bits,
    polynomial_divmod, signed_residue_from_algebraic_character_bits,
)
from odd_character_residues import (
    primitive_generator, root_of_unity_powers, signed_residue_from_character_bits,
)

ROOT = Path(__file__).resolve().parents[1]


def supplied_observations(order, mask):
    """Use the older exact integer root powers, independently of new ring code."""
    modulus, powers = root_of_unity_powers(order)
    u = order//(order & -order)
    observations = {}
    for j in range(1, u+1):
        if u % j:
            continue
        values = [0]*(len(modulus)-1)
        for exponent in range(order//2):
            if mask >> exponent & 1:
                values = [a ^ (b & 1) for a, b in zip(values, powers[j*exponent % order])]
        observations[j] = polynomial_bits(values)
    return observations


def coefficient_provider(prime, observations, calls):
    generator = primitive_generator(prime)
    dimension = len(root_of_unity_powers(prime-1)[0])-1

    def coefficient(g, index):
        assert g == generator and index in observations
        calls.append(index)
        return tuple(observations[index] >> position & 1 for position in range(dimension))

    return coefficient


def expected_residue(prime, mask):
    generator, residue = primitive_generator(prime), 1
    output = [0]*prime
    for exponent in range(prime-1):
        output[residue] = mask >> (exponent % ((prime-1)//2)) & 1
        residue = residue*generator % prime
    return tuple(output)


def check_roundtrips():
    rng = Random(20261101)
    levels = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 73, 97, 193, 241, 257, 337)
    records = []
    for prime in levels:
        width = (prime-1)//2
        if prime <= 19:
            masks = list(range(1 << width))
            exhaustive = True
        else:
            basis = range(width) if prime <= 41 else (0, width//2, width-1)
            masks = sorted(set([0, (1 << width)-1, *(1 << i for i in basis),
                                *(rng.getrandbits(width) for _ in range(6))]))
            exhaustive = False
        queried = 0
        for mask in masks:
            observations = supplied_observations(prime-1, mask)
            assert decode_character_polynomial(prime-1, observations) == mask
            calls = []
            provider = coefficient_provider(prime, observations, calls)
            algebraic = signed_residue_from_algebraic_character_bits(prime, provider)
            assert calls == sorted(observations)
            queried += len(calls)
            gaussian = signed_residue_from_character_bits(prime, provider)
            assert algebraic == gaussian == expected_residue(prime, mask)
        e, u, phi, modulus = parameters(prime-1)
        assert modulus == polynomial_bits(root_of_unity_powers(prime-1)[0])
        assert phi == binary_cyclotomic(u)
        records.append({"prime": prime, "nilpotent_length": e, "odd_order": u,
                        "paired_width": width, "vectors": len(masks),
                        "exhaustive": exhaustive, "coefficient_queries": queried})
    # The polynomial inverse also works at arbitrary even character orders.
    general_cases = 0
    for order in range(2, 82, 2):
        for _ in range(3):
            mask = rng.getrandbits(order//2)
            assert decode_character_polynomial(order, supplied_observations(order, mask)) == mask
            general_cases += 1
    return {"prime_levels": records, "total_prime_vectors": sum(x["vectors"] for x in records),
            "other_even_order_vectors": general_cases}


def check_essential_characters():
    witnesses = []
    for order in range(2, 98, 2):
        e, u, _, _ = parameters(order)
        for j in range(1, u+1):
            if u % j:
                continue
            d = u//j
            factor = polynomial_bits(root_of_unity_powers(d)[0])
            powered = sum(1 << (e*i) for i in range(factor.bit_length()) if factor >> i & 1)
            mask, remainder = polynomial_divmod((1 << (order//2)) ^ 1, powered)
            assert not remainder and mask and mask.bit_length() <= order//2
            observations = supplied_observations(order, mask)
            assert observations[j] and all(not value for k, value in observations.items() if k != j)
            witnesses.append({"order": order, "removed_character": j, "kernel_witness": hex(mask)})
    return witnesses


def check_divisor_vectors():
    cases = 0
    for integer in range(2, 65):
        divisors = [d for d in range(1, integer+1) if integer % d == 0]
        for prime in (3, 5, 7, 11, 13, 17, 19, 23):
            generator = primitive_generator(prime)
            order = prime-1
            modulus, powers = root_of_unity_powers(order)
            logs = {pow(generator, k, prime): k for k in range(order)}

            def coefficient(g, j):
                assert g == generator
                return tuple(sum(powers[j*logs[d % prime] % order][i]
                                 for d in divisors if d % prime) & 1
                             for i in range(len(modulus)-1))

            expected = [0]*prime
            for d in divisors:
                expected[d % prime] ^= 1
                expected[-d % prime] ^= 1
            assert signed_residue_from_algebraic_character_bits(prime, coefficient) == tuple(expected)
            cases += 1
    return cases


def check_rejections():
    invalid = [(6, {1: 0, 3: 2}), (10, {1: 0, 5: 2}), (6, {1: 0}),
               (6, {1: -1, 3: 0}), (6, {1: 4, 3: 0}), (3, {})]
    for order, observations in invalid:
        try:
            decode_character_polynomial(order, observations)
        except ValueError:
            pass
        else:
            raise AssertionError((order, observations))
    for prime in (0, 1, 4, 9, 15):
        try:
            signed_residue_from_algebraic_character_bits(prime, lambda *_: ())
        except ValueError:
            pass
        else:
            raise AssertionError(prime)
    assert signed_residue_from_algebraic_character_bits(2, lambda *_: ()) == (0, 0)
    return len(invalid)+5


def main():
    result = {"scope": "Exact polynomial inverse of supplied binary character data; no fast coefficient engine",
              "roundtrips": check_roundtrips(), "divisor_vectors": check_divisor_vectors(),
              "essential_character_witnesses": check_essential_characters(),
              "invalid_inputs_rejected": check_rejections()}
    result["source_sha256"] = {
        name: sha256((ROOT/name).read_bytes()).hexdigest()
        for name in ("tools/cyclotomic_character_decoder.py", "tools/check_cyclotomic_character_decoder.py",
                     "tools/odd_character_residues.py")}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
