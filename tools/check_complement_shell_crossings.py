#!/usr/bin/env python3
"""Exact homogeneous duals and complete crossings on fixed-length mode shells.

The raw vacuum Gram is an identity after partwise complementation below a
power of two. Complete ranks are measured separately; no exponent is inferred
from the identity or finite Taylor vanishing alone.
"""

from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess

from check_homogeneous_central_profiles import compile_checker, encode_input, python_profile
from check_modular_crossing import determinant_operator_two
from complement_shell_tensor import (
    folded_profile, multiply, pair_occupations, unfold_profile,
)


@lru_cache(None)
def shell(modulus: int, length: int, weight: int):
    if modulus < 4 or modulus & (modulus - 1):
        raise ValueError("the mode cutoff must be a power of two at least four")
    return tuple(parts for parts in combinations(range(1, modulus, 2), length)
                 if sum(parts) == weight)


def complement(parts, modulus):
    return tuple(sorted(modulus - part for part in parts))


def mode_mask(parts):
    return sum(1 << ((part - 1) // 2) for part in parts)


def verify_gram(modulus, length, weight):
    partners = shell(modulus, length, weight)
    parameters = tuple(complement(part, modulus) for part in partners)
    for x, first in enumerate(parameters):
        for v, second in enumerate(partners):
            assert determinant_operator_two(first, (), second) == int(x == v)
    return parameters, partners


def canonical_shell(label_bits):
    """Deterministic all-width dual labels without a Gram-minor search."""
    if label_bits < 1:
        raise ValueError("the binary label width must be positive")
    modes = 2 * label_bits + 16
    length = modes // 2
    modulus = 1 << (4 * modes - 1).bit_length()

    @lru_cache(None)
    def count(start, remaining, weight):
        if remaining == 0:
            return int(weight == 0)
        if (modes - start < remaining
                or weight < remaining * (2 * start + remaining)
                or weight > remaining * (2 * modes - remaining)):
            return 0
        return count(start + 1, remaining - 1, weight - (2 * start + 1)) + count(
            start + 1, remaining, weight)

    weights = range(length * length, length * (2 * modes - length) + 1, 2)
    weight = max(weights, key=lambda j: (count(0, length, j), -j))
    available = count(0, length, weight)
    wanted = 1 << label_bits
    assert available >= wanted

    def unrank(rank):
        result = []
        start, remaining, target = 0, length, weight
        while remaining:
            size = count(start + 1, remaining - 1, target - (2 * start + 1))
            if rank < size:
                result.append(2 * start + 1)
                target -= 2 * start + 1
                remaining -= 1
            else:
                rank -= size
            start += 1
        assert target == 0 and rank == 0
        return tuple(result)

    parameters = tuple(unrank(rank) for rank in range(wanted))
    partners = tuple(complement(parts, modulus) for parts in parameters)
    assert len(set(parameters)) == wanted
    assert all(sum(a) == weight and len(a) == length for a in parameters)
    assert all(sum(b) == length * modulus - weight for b in partners)
    assert weight < length * modulus - weight < 8 * modes * modes
    return modulus, length, available, parameters, partners


def verify_canonical_shells(max_bits=5):
    records = []
    for label_bits in range(1, max_bits + 1):
        modulus, length, available, parameters, partners = canonical_shell(label_bits)
        for x, first in enumerate(parameters):
            for v, second in enumerate(partners):
                assert determinant_operator_two(first, (), second) == int(x == v)
        records.append({
            "label_bits": label_bits, "labels": len(parameters),
            "mode_cutoff": modulus, "length": length,
            "operator_weight": sum(parameters[0]), "partner_weight": sum(partners[0]),
            "available_in_selected_shell": available,
            "vacuum_gram_entries_checked": len(parameters) ** 2,
        })
    return records


def verify_paired_reduction():
    """Compare the local tensor with the primitive reciprocal entry formula."""
    checked = 0
    for pairs, selected in ((2, 1), (4, 1), (4, 2), (8, 1)):
        labels = tuple(tuple(int(i in subset) for i in range(pairs))
                       for subset in combinations(range(pairs), selected))
        for x, y, u, v in product(labels, repeat=4):
            counts = {}
            for pattern in zip(x, y, u, v):
                counts[pattern] = counts.get(pattern, 0) + 1
            if any(sum(pattern) & 1 for pattern in counts):
                expected = 0
            else:
                b = counts.get((1, 1, 1, 1), 0)
                r = counts.get((1, 1, 0, 0), 0)
                c = counts.get((1, 0, 1, 0), 0)
                e = counts.get((1, 0, 0, 1), 0)
                assert r == counts.get((0, 0, 1, 1), 0)
                assert c == counts.get((0, 1, 0, 1), 0)
                assert e == counts.get((0, 1, 1, 0), 0)
                primitive = 1
                for _ in range(b):
                    primitive = multiply(primitive, 0b111)
                for _ in range(r + c):
                    primitive = multiply(primitive, 0b11)
                primitive <<= c + e
                expected = multiply(primitive, primitive)
            actual = folded_profile(*(tuple(3 * bit for bit in word) for word in (x, y, u, v)))
            assert actual == expected
            checked += 1
    return checked


def run_case(binary, modulus, length, weight, emit=True):
    parameters, partners = verify_gram(modulus, length, weight)
    j, h = length * modulus - weight, weight
    if not 0 < j <= h or not 1 <= len(partners) <= 64:
        raise ValueError("the finite full-profile checker requires 1<=j<=h and at most 64 labels")
    # Both sides have the same dimension. The lower-energy complementary
    # domain minimizes the state-energy choice within this particular pair.
    states = parameters
    command = [str(binary)] + (["--emit-profiles"] if emit else [])
    completed = subprocess.run(command, input=encode_input(j, h, j, parameters, partners, states),
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    record = json.loads(completed.stdout)
    record.update({"mode_cutoff": modulus, "length": length,
                   "operator_partitions": [list(x) for x in parameters],
                   "partner_partitions": [list(x) for x in partners],
                   "state_partitions": [list(x) for x in states],
                   "vacuum_gram_entries_checked": len(partners) ** 2})
    assert all(shift % modulus == 0 for shift, _, _ in record["nonzero_profile"])
    assert len(record["nonzero_profile"]) <= 2 * length + 1
    if emit:
        words = record.pop("profiles_words")
        d = len(partners)
        profiles = [sum(word << (64 * k) for k, word in enumerate(parts)) for parts in words]
        rng = random.Random(20261030133 + modulus + length + weight)
        selected = sorted(set([0, len(profiles) - 1] + [rng.randrange(len(profiles)) for _ in range(62)]))
        for index in selected:
            x, rem = divmod(index, d ** 3)
            y, rem = divmod(rem, d ** 2)
            v, u = divmod(rem, d)
            assert profiles[index] == python_profile(parameters[x], states[u], partners[v], states[y])
        record["independent_Python_Pfaffian_entries"] = len(selected)
        occupations = [pair_occupations(parts, modulus) for parts in parameters]
        for index, value in enumerate(profiles):
            x, rem = divmod(index, d ** 3)
            y, rem = divmod(rem, d ** 2)
            v, u = divmod(rem, d)
            folded = folded_profile(occupations[x], occupations[y], occupations[u], occupations[v])
            assert folded.bit_length() <= 2 * length + 1
            assert unfold_profile(folded, modulus, length, j) == value
        record["independent_universal_local_tensor_entries"] = len(profiles)
        # The operator and state labels need not be the same monomials. The
        # conserved-bit test therefore uses the abstract label partners for
        # all four legs, respecting the exact crossing identification.
        masks = [mode_mask(parts) for parts in partners]
        witness = None
        for index, value in enumerate(profiles):
            if not value:
                continue
            x, rem = divmod(index, d ** 3)
            y, rem = divmod(rem, d ** 2)
            v, u = divmod(rem, d)
            if masks[x] ^ masks[y] ^ masks[u] ^ masks[v]:
                witness = {"indices_x_y_u_v": [x, y, u, v],
                           "labels": [list(partners[k]) for k in (x, y, u, v)],
                           "shift_support": [s - 2 * j for s in range(value.bit_length()) if value >> s & 1]}
                break
        record["abstract_mode_XOR_violation"] = witness
        digest = hashlib.sha256()
        for value in profiles:
            size = (value.bit_length() + 7) // 8
            digest.update(size.to_bytes(4, "little"))
            digest.update(value.to_bytes(size, "little"))
        record["complete_profile_sha256"] = digest.hexdigest()
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="16:4:32,16:3:25,32:3:49")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    binary, native_hash = compile_checker()
    canonical = verify_canonical_shells()
    paired_entries = verify_paired_reduction()
    records = []
    for description in args.cases.split(','):
        modulus, length, weight = map(int, description.split(':'))
        record = run_case(binary, modulus, length, weight)
        records.append(record)
        print(json.dumps({"case": description, "labels": record["labels"],
                          "rank_sum": record["rank_sum"], "threshold": record["finite_seed_threshold"],
                          "valuation": record["minimum_residual_taylor_valuation"],
                          "XOR_violation": record["abstract_mode_XOR_violation"] is not None}), flush=True)
        if args.output:
            report = {"scope": "exact complement-shell vacuum duality and finite full-profile ranks; no stronger factoring exponent",
                      "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "native_profile_sha256": native_hash,
                      "local_tensor_sha256": hashlib.sha256(
                          Path(__file__).with_name("complement_shell_tensor.py").read_bytes()).hexdigest(),
                      "canonical_all_width_construction_checks": canonical,
                      "paired_reciprocal_reduction_entries": paired_entries, "cases": records}
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
