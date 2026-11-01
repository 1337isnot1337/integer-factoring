#!/usr/bin/env python3
"""Exact finite checks for the one-scalar quadratic recovery interface.

The ordinary arithmetic helpers are shared with the unitary checker.
All roots here are found by literal residue enumeration, independently
of the orbit count and independent of a factorization of the input.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
from random import Random

from check_unitary_orbit_recovery import (
    add_one_trace,
    factor_reference,
    perfect_power,
    recover_nonzero,
)


def factor_by_roots(value: int, roots: dict[int, tuple[int, ...]]) -> tuple[list[int], int]:
    output: list[int] = []
    while value % 2 == 0:
        output.append(2)
        value //= 2
    pending = [] if value == 1 else [(value, 1)]
    queries = 0
    while pending:
        current, multiplicity = pending.pop()
        assert current % 2 == 1
        root = perfect_power(current)
        if root is not None:
            base, exponent = root
            pending.append((base, multiplicity * exponent))
            continue
        witnesses = tuple((z,) for z in roots[current])
        queries += 1
        if len(witnesses) % 2 == 0:
            output.extend([current] * multiplicity)
            continue
        (z,), used = recover_nonzero(witnesses, current.bit_length(), 1, 2)
        queries += used
        divisor = gcd(current, z - 1)
        assert 1 < divisor < current and current % divisor == 0
        pending.extend(((divisor, multiplicity), (current // divisor, multiplicity)))
    return sorted(output), queries


def check_all(limit: int) -> dict[str, int]:
    roots: dict[int, tuple[int, ...]] = {}
    stats: Counter[str] = Counter()
    for value in range(2, limit + 1):
        omega = len(set(factor_reference(value)))
        expected = 2 ** (omega - 1) - 1
        idempotents = tuple(x for x in range(value) if x * (x - 1) % value == 0)
        assert len(idempotents) == 2**omega
        assert all((1 - x) % value in idempotents and (1 - x) % value != x for x in idempotents)
        representatives = tuple(x for x in idempotents if 2 <= x <= value // 2)
        assert len(representatives) == expected
        for x in representatives:
            y = x * (x - 1) // value
            assert 0 <= y < value and x * (x - 1) == value * y
            assert 1 < gcd(value, x) < value
            stats["idempotent_quotients"] += 1
        stats["idempotent_inputs"] += 1
        stats["idempotent_residues_examined"] += value
        if value % 2:
            all_roots = tuple(z for z in range(value) if (z * z - 1) % value == 0)
            assert len(all_roots) == 2**omega
            assert all((-z) % value in all_roots and (-z) % value != z for z in all_roots)
            admitted = tuple(z for z in all_roots if 2 <= z <= value // 2)
            assert len(admitted) == expected
            roots[value] = admitted
            for z in admitted:
                y = (z * z - 1) // value
                assert 0 <= y < value and z * z == value * y + 1
                result, carry, trace = add_one_trace(value * y, 1, 2 * value.bit_length())
                assert result == z * z and carry == 0
                assert len(trace) == 2 * value.bit_length() + 1
                assert gcd(value, z) == 1
                assert 1 < gcd(value, z - 1) < value
                stats["square_root_quotients"] += 1
            stats["odd_square_root_inputs"] += 1
            stats["square_root_residues_examined"] += value
        found, queries = factor_by_roots(value, roots)
        assert found == factor_reference(value)
        stats["factorization_inputs"] += 1
        stats["maximum_factorization_queries"] = max(stats["maximum_factorization_queries"], queries)
    return dict(stats)


def check_negative_controls() -> dict[str, object]:
    roots8 = [z for z in range(8) if (z * z - 1) % 8 == 0]
    assert roots8 == [1, 3, 5, 7]
    assert len(roots8) != 2 ** len(set(factor_reference(8)))
    roots15 = [z for z in range(15) if (z * z - 1) % 15 == 0]
    admitted15 = [z for z in roots15 if 2 <= z <= 15 // 2]
    assert roots15 == [1, 4, 11, 14] and admitted15 == [4]
    assert len(roots15) % 2 == 0 and len(admitted15) % 2 == 1
    assert gcd(15, 4) == 1 and gcd(15, 4 - 1) == 3
    roots7_half = [z for z in range(1, 7 // 2 + 1) if (z * z - 1) % 7 == 0]
    assert roots7_half == [1]
    assert perfect_power(9) is not None
    assert [z for z in range(2, 9 // 2 + 1) if (z * z - 1) % 9 == 0] == []
    return {
        "odd_input_assumption_required": {"N": 8, "square_roots": roots8},
        "full_root_parity_loses_the_signal": {"N": 15, "full_roots": roots15, "canonical_nontrivial_roots": admitted15},
        "gcd_requires_z_minus_one": {"N": 15, "z": 4, "gcd_N_z": 1, "gcd_N_z_minus_one": 3},
        "trivial_root_must_be_removed": {"N": 7, "half_interval_roots": roots7_half},
        "perfect_power_preprocessing_required": {"N": 9, "canonical_nontrivial_root_count": 0},
    }


def fused_source_transform(z: int, u: int, fixed: int, scale: int) -> dict[str, int]:
    """Literal pinned source transforms with the changed pointwise operation.

    N*u is never used to produce a transform coefficient. Ordinary integer
    products and schoolbook convolution are independent final comparisons.
    The source's explicit transform inverse works over the composite ring.
    """
    length = 1 << (2 * scale + 1)
    a = 1 << (scale + 1)
    half = 1 << scale
    count = a * a
    base = 1 << length
    modulus = (1 << (4 * length)) + 1
    capacity = length * length
    qword = (1 << capacity) - 1
    assert all(0 <= value <= qword for value in (z, u, fixed))
    digit_lists = [[(value >> (length * i)) & (base - 1) for i in range(length)]
                   for value in (z, u, fixed)]
    inverse_root = pow(16, count - 1, modulus)
    inverse_count = pow(2, (8 * length - 1) * (2 * scale + 2), modulus)
    assert inverse_count * count % modulus == 1

    def first_layer(digits):
        return [[pow(16, r * k, modulus) * sum(
            digits[r + a * i] * pow(16, a * i * k, modulus)
            for i in range(half)) % modulus for k in range(a)] for r in range(a)]

    fa_z, fa_u, fa_fixed = map(first_layer, digit_lists)
    mid = [[0] * a for _ in range(a)]
    for k1 in range(a):
        pointwise = []
        for k2 in range(a):
            hats = [sum(family[r][k1] * pow(16, a * r * k2, modulus)
                        for r in range(a)) % modulus
                    for family in (fa_z, fa_u, fa_fixed)]
            pointwise.append((hats[0] * hats[0] + hats[2] * hats[1]) % modulus)
        for j1 in range(a):
            mid[k1][j1] = pow(inverse_root, k1 * j1, modulus) * sum(
                pointwise[k2] * pow(inverse_root, a * k2 * j1, modulus)
                for k2 in range(a)) % modulus
    residues = [0] * count
    for j1 in range(a):
        for j2 in range(a):
            residues[j1 + a * j2] = inverse_count * sum(
                mid[k1][j1] * pow(inverse_root, a * k1 * j2, modulus)
                for k1 in range(a)) % modulus

    expected = [0] * count
    zd, ud, nd = digit_lists
    for i in range(length):
        for j in range(length):
            expected[i + j] += zd[i] * zd[j] + nd[i] * ud[j]
    assert residues == expected
    assert 2 * length * (base - 1) ** 2 < modulus - 1
    assert all(0 <= coefficient < modulus for coefficient in residues)
    carry = 0
    output = 0
    for j, coefficient in enumerate(residues):
        total = coefficient + carry
        output |= (total % base) << (length * j)
        carry = total // base
        assert carry <= 2 * length * (base - 1)
        assert carry < 1 << (4 * length + 1)
    assert output + (carry << (length * count)) == z * z + fixed * u
    target = fixed * qword + 1
    accepted = int(carry == 0 and output == target)
    assert bool(accepted) == (z * z + fixed * u == target)
    return {"coefficients": count, "carry": carry, "accepted": accepted,
            "modulus": modulus, "low_output": output}


def check_fused_source() -> dict[str, object]:
    rng = Random(20261030)
    totals: Counter[str] = Counter()
    accepted_by_fixed: Counter[int] = Counter()
    for fixed in range(1, 16):
        for z in range(16):
            for u in range(16):
                record = fused_source_transform(z, u, fixed, 0)
                totals["literal_scale_zero_inputs"] += 1
                totals["coefficients"] += record["coefficients"]
                totals["overflow_cases"] += int(record["carry"] != 0)
                if 2 <= z <= fixed // 2:
                    accepted_by_fixed[fixed] += record["accepted"]
    for fixed in range(3, 16, 2):
        assert accepted_by_fixed[fixed] == 2 ** (len(set(factor_reference(fixed))) - 1) - 1
    for scale, repeats in ((1, 60), (2, 12), (3, 2)):
        capacity = (1 << (2 * scale + 1)) ** 2
        for _ in range(repeats):
            values = [rng.getrandbits(capacity) for _ in range(3)]
            record = fused_source_transform(*values, scale)
            totals["larger_scale_inputs"] += 1
            totals["coefficients"] += record["coefficients"]
    for fixed in (15, 21, 35, 45, 77, 105, 315, 1729):
        qword = (1 << 64) - 1
        for z in range(2, fixed // 2 + 1):
            if (z * z - 1) % fixed == 0:
                y = (z * z - 1) // fixed
                record = fused_source_transform(z, qword - y, fixed, 1)
                assert record["accepted"] == 1 and record["carry"] == 0
                totals["valid_root_inputs"] += 1
                totals["coefficients"] += record["coefficients"]
    overflow = fused_source_transform(15, 15, 15, 0)
    assert overflow["carry"] == 1 and overflow["low_output"] == 194
    assert overflow["accepted"] == 0
    composite_modulus = (1 << 32) + 1
    assert composite_modulus == 641 * 6700417
    return {**dict(totals), "composite_modulus_control": composite_modulus,
            "overflow_control": {"z": 15, "u": 15, "N": 15, "low_output": 194, "final_carry": 1},
            "source_commit": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=8192)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "PASS", "scope": "Exact finite transcription checks, not an asymptotic runtime measurement",
        "checks": check_all(args.limit), "negative_controls": check_negative_controls(),
        "fused_source": check_fused_source(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.output is not None:
        print(f"PASS: wrote {args.output}")


if __name__ == "__main__":
    main()
