#!/usr/bin/env python3
"""Literal ideal and independent prime-expansion controls for two-pass lifting."""

import argparse
from hashlib import sha256
import json
import lzma
from math import gcd
from pathlib import Path
from random import Random

from circuit_box_oriented_batch import two_pass_oriented_basis
from classical_coprime_basis import classical_oriented_coprime_basis


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATHS = [
    "tools/circuit_box_oriented_batch.py",
    "tools/check_circuit_box_oriented_batch.py",
    "tools/classical_coprime_basis.py",
    "tools/oriented_norm_basis_lift.py",
    "tools/oriented_coprime_basis.py",
    "tools/quadratic_ideal_hnf.py",
]


def _factor_small(value):
    result = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result[divisor] = result.get(divisor, 0) + 1
            value //= divisor
        divisor += 1
    if value > 1:
        result[value] = result.get(value, 0) + 1
    return result


def _prime_words(keys, record):
    factorings = [_factor_small(atom["norm"]) for atom in record["atoms"]]
    for (norm, root), word in zip(keys, record["sparse_words"]):
        expected = {}
        for prime, exponent in _factor_small(norm).items():
            residue = root % (2 * prime)
            canonical = min(residue, (-residue) % (2 * prime))
            expected[prime] = exponent if residue == canonical else -exponent
        actual = {}
        for column, exponent in word:
            atom = record["atoms"][column]
            for prime, multiplicity in factorings[column].items():
                residue = atom["root"] % (2 * prime)
                canonical = min(residue, (-residue) % (2 * prime))
                orientation = 1 if residue == canonical else -1
                actual[prime] = actual.get(prime, 0) + orientation * multiplicity * exponent
        assert actual == expected


def _mixed_keys():
    primes = [5, 13, 17, 29, 37, 41, 53, 61]
    keys = []
    for mask in range(256):
        exponent, root, norm = 1 + mask % 5, 0, 1
        for index, prime in enumerate(primes):
            local = next(value for value in range(prime) if (value * value + 4) % prime == 0)
            modulus = prime
            for _ in range(1, exponent):
                local += (((-4 - local * local) // modulus) * pow(2 * local, -1, prime) % prime) * modulus
                modulus *= prime
            if mask >> index & 1:
                local = -local % modulus
            root += ((local - root) * pow(norm, -1, modulus) % modulus) * norm
            norm *= modulus
        if root % 2:
            root += norm
        keys.append((norm, root))
    return primes, keys


def run_checks():
    source_hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCE_PATHS}
    rng = Random(2026110190731)
    counts = {"ordered_ideal_batches": 0, "literal_ideal_identities": 0, "independent_prime_expansions": 0, "single_pass_comparisons": 0, "raw_delta_variant_comparisons": 0, "invalid_inputs_rejected": 0}
    small_records = []
    for discriminant in (-4, -23, -84, -260, -575):
        keys = [(norm, root) for norm in range(1, 140, 2) if gcd(norm, discriminant) == 1 for root in range(2 * norm) if (root * root - discriminant) % (4 * norm) == 0]
        shuffled = list(keys)
        rng.shuffle(shuffled)
        for order in (keys, list(reversed(keys)), shuffled + shuffled[:20]):
            result = two_pass_oriented_basis(discriminant, order)
            _prime_words(order, result)
            raw = two_pass_oriented_basis(discriminant, order, unitary_preconditioning=False)
            _prime_words(order, raw)
            original = classical_oriented_coprime_basis(discriminant, order)
            _prime_words(order, original)
            counts["ordered_ideal_batches"] += 1
            counts["literal_ideal_identities"] += 2 * len(order)
            counts["independent_prime_expansions"] += 3 * len(order)
            counts["single_pass_comparisons"] += 1
            counts["raw_delta_variant_comparisons"] += 1
            statistics = result["statistics"]
            assert statistics["orientation_gcd_calls_after_basis_passes"] == 0
            assert statistics["second_pass_input_bits"] <= result["second_pass_data_bit_bound"]
            assert statistics["discarded_auxiliary_atoms"] == 0
            small_records.append({"discriminant": discriminant, "keys": len(order), "input_bits": result["input_bit_bound"], "atoms": len(result["atoms"]), "second_pass_input_bits": statistics["second_pass_input_bits"], "normalization_remainders": statistics["normalization_remainders"], "discarded_auxiliary_atoms": statistics["discarded_auxiliary_atoms"], "raw_delta_atoms": len(raw["atoms"]), "raw_delta_discarded_atoms": raw["statistics"]["discarded_auxiliary_atoms"]})

    edge_records = []
    for keys in ([], [(1, 0)] * 7, [(65, 16), (65, 36)], [(325, 114), (325, 614)], [(169, 140), (13, 10), (169, 198)], [(29, 24), (29, 34), (65, 16)]):
        result = two_pass_oriented_basis(-4, keys)
        _prime_words(keys, result)
        raw = two_pass_oriented_basis(-4, keys, unitary_preconditioning=False)
        _prime_words(keys, raw)
        counts["literal_ideal_identities"] += 2 * len(keys)
        counts["independent_prime_expansions"] += 2 * len(keys)
        counts["raw_delta_variant_comparisons"] += 1
        edge_records.append({"keys": keys, "atoms": result["atoms"], "words": result["sparse_words"], "statistics": result["statistics"], "raw_delta_atoms": raw["atoms"], "raw_delta_words": raw["sparse_words"], "raw_delta_statistics": raw["statistics"]})
    assert [atom["norm"] for atom in result["atoms"]] == [29, 65]
    assert [atom["norm"] for atom in raw["atoms"]] == [5, 13, 29]

    primes, keys = _mixed_keys()
    mixed = two_pass_oriented_basis(-4, keys)
    mixed_raw = two_pass_oriented_basis(-4, keys, unitary_preconditioning=False)
    assert [atom["norm"] for atom in mixed["atoms"]] == primes
    assert mixed["atoms"] == mixed_raw["atoms"] and mixed["sparse_words"] == mixed_raw["sparse_words"]
    for (norm, root), word in zip(keys, mixed["sparse_words"]):
        expected = []
        remainder = norm
        for column, prime in enumerate(primes):
            exponent = 0
            while remainder % prime == 0:
                remainder //= prime
                exponent += 1
            if exponent:
                residue = root % (2 * prime)
                canonical = mixed["atoms"][column]["root"]
                expected.append((column, exponent if residue == canonical else -exponent))
        assert remainder == 1 and word == expected
    counts["literal_ideal_identities"] += 2 * len(keys)
    counts["independent_prime_expansions"] += 2 * len(keys)
    counts["raw_delta_variant_comparisons"] += 1

    fixtures = []
    for bits in (80, 96):
        relative = f"evidence/sparse-oriented-power-{bits}.json.xz" if bits == 80 else f"evidence/optimized-oriented-power-{bits}.json.xz"
        raw = (ROOT / relative).read_bytes()
        saved = json.loads(lzma.decompress(raw))["result"]["refinement"]
        keys = saved["keys"][:128]
        result = two_pass_oriented_basis(saved["discriminant"], keys)
        raw_result = two_pass_oriented_basis(saved["discriminant"], keys, unitary_preconditioning=False)
        counts["literal_ideal_identities"] += 2 * len(keys)
        counts["raw_delta_variant_comparisons"] += 1
        fixtures.append({"input_path": relative, "input_sha256": sha256(raw).hexdigest(), "full_retained_key_count": len(saved["keys"]), "checked_prefix_keys": len(keys), "input_bits": result["input_bit_bound"], "atoms": len(result["atoms"]), "second_pass_data_bit_bound": result["second_pass_data_bit_bound"], "statistics": result["statistics"], "raw_delta_atoms": len(raw_result["atoms"]), "raw_delta_statistics": raw_result["statistics"]})

    for discriminant, keys in ((1, []), (-3, [(2, 1)]), (-4, [(5, 1)]), (-20, [(5, 0)]), (-4, [(0, 0)])):
        try:
            two_pass_oriented_basis(discriminant, keys)
        except ValueError:
            counts["invalid_inputs_rejected"] += 1
        else:
            raise AssertionError((discriminant, keys))
    assert source_hashes == {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCE_PATHS}, "a dependency changed during the checks"
    return {
        "status": "passed", "seed": 2026110190731, "counts": counts,
        "small_batches": small_records, "edge_cases": edge_records,
        "mixed_prime_power_family": {"keys": len(mixed["keys"]), "maximum_norm_bits": max(norm.bit_length() for norm, _ in mixed["keys"]), "first_pass_basis_size": mixed["statistics"]["first_pass_basis_size"], "final_atoms": len(mixed["atoms"]), "input_bits": mixed["input_bit_bound"], "second_pass_data_bit_bound": mixed["second_pass_data_bit_bound"], "statistics": mixed["statistics"], "raw_delta_statistics": mixed_raw["statistics"]},
        "retained_fixture_prefixes": fixtures,
        "source_sha256": source_hashes,
        "scope": "Complete two-pass construction on supplied keys. Classical fast-arithmetic complexity is separate from Python timings; no relation-supply or factoring-score claim.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    if args.output:
        args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": result["status"], **result["counts"], "mixed_family": {key: result["mixed_prime_power_family"][key] for key in ("keys", "maximum_norm_bits", "first_pass_basis_size", "final_atoms")}, "fixture_prefixes": [entry["checked_prefix_keys"] for entry in result["retained_fixture_prefixes"]]}))


if __name__ == "__main__":
    main()
