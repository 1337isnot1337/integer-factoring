#!/usr/bin/env python3
"""Replay the pinned multiplication transforms with variable output ports.

The transform formulas are from Rudich--Rousseau c91cc3f, FactoringAux/
{Transform,MulNet}.lean (Apache-2.0). Integer multiplication is used only
as an independent comparison, never to produce the transform output.
These are finite checks, not a new factoring runtime proof.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from random import Random


@dataclass(frozen=True)
class Params:
    scale: int
    a: int
    length: int
    half: int
    modulus: int
    wire_bits: int

    @classmethod
    def at(cls, scale: int) -> "Params":
        length = 1 << (2 * scale + 1)
        return cls(scale, 1 << (scale + 1), length, 1 << scale,
                   (1 << (4 * length)) + 1, 4 * length + 1)

    @property
    def input_bits(self) -> int:
        return self.length * self.length

    @property
    def digit_count(self) -> int:
        return self.a * self.a


def source_transform(left: int, right: int, scale: int) -> dict:
    """Literal FA/MID/IB layers, followed by deterministic carry ports."""
    p = Params.at(scale)
    assert 0 <= left < 1 << p.input_bits
    assert 0 <= right < 1 << p.input_bits
    digit_mask = (1 << p.length) - 1
    left_digits = [(left >> (p.length * i)) & digit_mask
                   for i in range(p.length)]
    right_digits = [(right >> (p.length * i)) & digit_mask
                    for i in range(p.length)]
    inverse_root = pow(16, p.digit_count - 1, p.modulus)
    # This is the source's explicit tInv, not an assumed field inverse.
    inverse_count = pow(2, (8 * p.length - 1) * (2 * scale + 2), p.modulus)
    assert inverse_count * p.digit_count % p.modulus == 1

    def first_layer(digits):
        return [[pow(16, r * k, p.modulus) * sum(
            digits[r + p.a * i] * pow(16, p.a * i * k, p.modulus)
            for i in range(p.half)) % p.modulus
            for k in range(p.a)] for r in range(p.a)]

    fa_left, fa_right = first_layer(left_digits), first_layer(right_digits)
    mid = [[0] * p.a for _ in range(p.a)]
    for k1 in range(p.a):
        pointwise = []
        for k2 in range(p.a):
            value_left = sum(fa_left[r][k1] * pow(16, p.a * r * k2, p.modulus)
                             for r in range(p.a)) % p.modulus
            value_right = sum(fa_right[r][k1] * pow(16, p.a * r * k2, p.modulus)
                              for r in range(p.a)) % p.modulus
            pointwise.append(value_left * value_right % p.modulus)
        for j1 in range(p.a):
            mid[k1][j1] = pow(inverse_root, k1 * j1, p.modulus) * sum(
                pointwise[k2] * pow(inverse_root, p.a * k2 * j1, p.modulus)
                for k2 in range(p.a)) % p.modulus

    residues = [0] * p.digit_count
    for j1 in range(p.a):
        for j2 in range(p.a):
            residues[j1 + p.a * j2] = inverse_count * sum(
                mid[k1][j1] * pow(inverse_root, p.a * k1 * j2, p.modulus)
                for k1 in range(p.a)) % p.modulus

    output, carries = [], [0]
    for residue in residues:
        total = residue + carries[-1]
        output.append(total & digit_mask)
        carries.append(total >> p.length)
    assert all(0 <= x < 1 << p.wire_bits for x in residues + carries)
    assert carries[-1] == 0
    value = sum(x << (p.length * j) for j, x in enumerate(output))
    return {"params": p, "left_digits": left_digits, "right_digits": right_digits,
            "residues": residues, "output_digits": output,
            "carries": carries, "value": value}


def check_transform(left: int, right: int, scale: int) -> dict:
    record = source_transform(left, right, scale)
    p = record["params"]
    expected = [0] * p.digit_count
    for i, a in enumerate(record["left_digits"]):
        for j, b in enumerate(record["right_digits"]):
            expected[i + j] += a * b
    assert expected == record["residues"]
    assert all(x < p.modulus for x in expected)
    assert record["value"] == left * right
    return record


def copied_bit_wires(source: dict, target_scale: int, copied_bits: int) -> dict:
    """Pack equalities by the actual carry/input-gate pair, preserving order."""
    sp, tp = source["params"], Params.at(target_scale)
    assert 0 <= copied_bits <= min(2 * sp.input_bits, tp.input_bits)
    value = source["value"]
    assert value < 1 << copied_bits if copied_bits else value == 0
    target_digits = [(value >> (tp.length * j)) & ((1 << tp.length) - 1)
                     for j in range(tp.length)]
    by_pair = defaultdict(list)
    for bit in range(copied_bits):
        carry_gate, source_offset = divmod(bit, sp.length)
        digit, target_offset = divmod(bit, tp.length)
        input_gate = digit % tp.a
        input_slot = digit // tp.a
        by_pair[(carry_gate, input_gate)].append(
            (source_offset, input_slot, target_offset))
    total_width = 0
    for (carry_gate, input_gate), coordinates in by_pair.items():
        source_message = sum(((source["output_digits"][carry_gate] >> so) & 1) << j
                             for j, (so, _, _) in enumerate(coordinates))
        target_message = sum(((target_digits[input_gate + tp.a * slot] >> to) & 1) << j
                             for j, (_, slot, to) in enumerate(coordinates))
        assert source_message == target_message
        assert source_message < 1 << len(coordinates)
        total_width += len(coordinates)
    assert total_width == copied_bits
    return {"wires": len(by_pair), "width": total_width}


def verify() -> dict:
    rng = Random(0xF710EB9C)
    counts = {"transform_inputs": 0, "convolution_coefficients": 0,
              "copied_port_bits": 0, "joined_products": 0,
              "reject_controls": 0}
    cases = [(a, b, 0) for a in range(16) for b in range(16)]
    for scale, random_count in ((1, 120), (2, 20), (3, 4)):
        bits = Params.at(scale).input_bits
        anchors = sorted({0, 1, 2, (1 << bits) - 1,
                          1 << (bits - 1), (1 << (bits // 2)) - 1})
        cases.extend((a, b, scale) for a in anchors for b in anchors)
        cases.extend((rng.getrandbits(bits), rng.getrandbits(bits), scale)
                     for _ in range(random_count))
    by_scale = defaultdict(int)
    for left, right, scale in cases:
        result = check_transform(left, right, scale)
        counts["transform_inputs"] += 1
        counts["convolution_coefficients"] += len(result["residues"])
        by_scale[scale] += 1
        # An altered fixed output digit must fail the source carry weight.
        j = rng.randrange(result["params"].digit_count)
        wrong_digit = result["output_digits"][j] ^ 1
        total = result["residues"][j] + result["carries"][j]
        assert total % (1 << result["params"].length) != wrong_digit
        counts["reject_controls"] += 1

    # s=1 uses a composite auxiliary modulus; no field assumption is allowed.
    assert Params.at(1).modulus == 641 * 6700417

    # Join small and large digit formats, including carry-to-strided-input maps.
    for a in range(16):
        square = check_transform(a, a, 0)
        wiring = copied_bit_wires(square, 1, 8)
        counts["copied_port_bits"] += wiring["width"]
        for b in (0, 1, 3, 17, 255):
            final = check_transform(square["value"], b, 1)
            assert final["value"] == a * a * b
            counts["joined_products"] += 1
    for _ in range(32):
        # Both directions of digit refinement are tested, with pinned high zeros.
        a, b = rng.randrange(4), rng.randrange(4)
        result = check_transform(a, b, 1)
        wiring = copied_bit_wires(result, 0, 4)
        counts["copied_port_bits"] += wiring["width"]

    # Three repeated squarings with independent final-product comparisons.
    for a in range(1, 8):
        current = a
        for power in (2, 4, 8):
            result = check_transform(current, current, 1)
            current = result["value"]
            assert current == a ** power
            copied_bit_wires(result, 1, max(1, current.bit_length()))
            counts["joined_products"] += 1
        for b in (1, 3, 9):
            result = check_transform(current, b, 1)
            assert result["value"] == a ** 8 * b
            counts["joined_products"] += 1

    # Independently recompute the compiler's exact arithmetic envelope at large n.
    parameter_cases = 0
    for p in range(2, 41):
        n = p ** 40
        # D<=C*n^(3/4) exactly here because n=P^40.
        for c in (1, 2, 7):
            w, d, g = c * n, c * p ** 30, 3 * p ** 8
            q = w // g + d
            label_width = max(1, (p * w + g * q) // p ** 3)
            h = (2 * g // p) * q
            assert p * p * label_width <= (1 + 3 * c) * n
            assert p * h <= 8 * c * n
            assert 3 * p ** 10 <= 3 * n
            parameter_cases += 1
    return {"status": "pass", "scope": "finite source-port and compiler controls; no new factoring exponent",
            "source_commit": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1",
            "counts": counts, "transform_cases_by_scale": dict(sorted(by_scale.items())),
            "exact_parameter_cases": parameter_cases,
            "composite_modulus_control": Params.at(1).modulus,
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
