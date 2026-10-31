#!/usr/bin/env python3
"""Exact finite controls for Fourier-coset fan-out ranks and optimal codecs.

The source fan-out is faV/faOut at factoring revision c91cc3f.
This checker establishes finite identities only; the uniform proof is in
research/factoring-network-structure.md. No factoring score increase.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations, product
import json
from pathlib import Path
import random


def power_of_two(n: int) -> bool:
    return n > 0 and n & (n - 1) == 0


@dataclass(frozen=True)
class Fanout:
    h: int
    digit_bits: int
    gate_index: int = 0

    def __post_init__(self) -> None:
        if not power_of_two(self.h) or self.digit_bits <= 0:
            raise ValueError("positive power-of-two h and positive digit width required")
        if 4 * self.digit_bits % self.h:
            raise ValueError("h must divide four times the digit width")
        if not 0 <= self.gate_index < self.a:
            raise ValueError("gate index outside source layer")
        if self.h * (self.base - 1) >= self.modulus:
            raise ValueError("zero-coset sum can wrap")

    @property
    def a(self) -> int:
        return 2 * self.h

    @property
    def base(self) -> int:
        return 1 << self.digit_bits

    @property
    def modulus(self) -> int:
        return (1 << (4 * self.digit_bits)) + 1

    @property
    def root(self) -> int:
        return pow(2, 4 * self.digit_bits // self.h, self.modulus)

    @property
    def source_parameters(self) -> bool:
        return self.digit_bits == 2 * self.h * self.h

    def phase(self, j: int) -> int:
        return pow(16, self.gate_index * j, self.modulus)

    def messages(self, digits: tuple[int, ...]) -> tuple[int, ...]:
        if len(digits) != self.h or any(not 0 <= x < self.base for x in digits):
            raise ValueError("invalid source digits")
        # Direct original evaluation, separate from the folded codec.
        return tuple(self.message_at(digits, j) for j in range(self.a))

    def message_at(self, digits: tuple[int, ...], j: int) -> int:
        return (self.phase(j)
                * sum(x * pow(self.root, i * j, self.modulus)
                      for i, x in enumerate(digits)) % self.modulus)


@dataclass(frozen=True)
class Bundle:
    fan: Fanout
    p: int
    c: int

    def __post_init__(self) -> None:
        if not power_of_two(self.p) or not 2 <= self.p <= self.fan.a:
            raise ValueError("P must be a power of two between two and A")
        if not 0 <= self.c < self.p:
            raise ValueError("invalid coset")

    @property
    def k(self) -> int:
        return self.fan.a // self.p

    @property
    def t(self) -> int:
        return self.p // 2

    @property
    def q(self) -> int:
        return self.t if self.c == 0 else self.c & -self.c

    @property
    def r(self) -> int:
        return 1 if self.c == 0 else self.t // self.q

    @property
    def small_base(self) -> int:
        return 1 << (4 * self.fan.digit_bits // self.r)

    @property
    def height(self) -> int:
        return self.q * (self.fan.base - 1) + 1

    @property
    def image_per_coordinate(self) -> int:
        if self.height <= self.small_base:
            return self.height ** self.r
        return self.fan.modulus

    @property
    def rank(self) -> int:
        return self.image_per_coordinate ** self.k

    @property
    def indices(self) -> tuple[int, ...]:
        return tuple(self.c + self.p * j for j in range(self.k))

    def signed_position(self, a: int) -> tuple[int, int]:
        if self.c == 0:
            return 0, 1
        exponent = (self.c // self.q) * a
        return exponent % self.r, -1 if exponent // self.r % 2 else 1

    def translation(self) -> int:
        return sum((self.fan.base - 1) * self.small_base ** j
                   for a in range(self.t)
                   for j, sign in (self.signed_position(a),) if sign < 0)

    def row_code(self, digits: tuple[int, ...]) -> int:
        if len(digits) != self.fan.h or any(not 0 <= x < self.fan.base for x in digits):
            raise ValueError("invalid source digits")
        answer = 0
        for ell in reversed(range(self.k)):
            reduced = [0] * self.r
            for a in range(self.t):
                j, sign = self.signed_position(a)
                value = digits[ell + self.k * a]
                reduced[j] += value if sign > 0 else self.fan.base - 1 - value
            assert all(0 <= value < self.height for value in reduced)
            radix = self.height if self.height <= self.small_base else self.small_base
            code = sum(value * radix ** j for j, value in enumerate(reduced))
            if self.height > self.small_base:
                code %= self.fan.modulus
            answer = answer * self.image_per_coordinate + code
        return answer

    def column_code(self, messages: tuple[int, ...]) -> int | None:
        if len(messages) != self.k or any(not 0 <= x < self.fan.modulus for x in messages):
            return None
        modulus = self.fan.modulus
        stripped = [x * pow(self.fan.phase(j), -1, modulus) % modulus
                    for j, x in zip(self.indices, messages)]
        xi_inverse = pow(pow(self.fan.root, self.p, modulus), -1, modulus)
        gamma_inverse = pow(pow(self.fan.root, self.c, modulus), -1, modulus)
        k_inverse = pow(self.k, -1, modulus)
        translation = self.translation()
        codes = []
        for ell in range(self.k):
            value = (k_inverse * sum(x * pow(xi_inverse, ell * j, modulus)
                                      for j, x in enumerate(stripped))
                     * pow(gamma_inverse, ell, modulus) + translation) % modulus
            if self.height <= self.small_base:
                code, power = 0, 1
                for _ in range(self.r):
                    value, digit = divmod(value, self.small_base)
                    if digit >= self.height:
                        return None
                    code += digit * power
                    power *= self.height
                if value:
                    return None
            else:
                code = value
            codes.append(code)
        return sum(code * self.image_per_coordinate ** ell
                   for ell, code in enumerate(codes))

    def messages_from_code(self, code: int) -> tuple[int, ...]:
        if not 0 <= code < self.rank:
            raise ValueError("channel outside its exact range")
        folded = []
        for _ in range(self.k):
            code, item = divmod(code, self.image_per_coordinate)
            if self.height <= self.small_base:
                value = 0
                for j in range(self.r):
                    item, digit = divmod(item, self.height)
                    value += digit * self.small_base ** j
            else:
                value = item
            folded.append((value - self.translation()) % self.fan.modulus)
        gamma = pow(self.fan.root, self.c, self.fan.modulus)
        xi = pow(self.fan.root, self.p, self.fan.modulus)
        return tuple(
            self.fan.phase(index)
            * sum(value * pow(gamma, ell, self.fan.modulus)
                  * pow(xi, ell * j, self.fan.modulus)
                  for ell, value in enumerate(folded)) % self.fan.modulus
            for j, index in enumerate(self.indices)
        )

    def source_from_code(self, code: int) -> tuple[int, ...]:
        """Construct a source digit vector realizing every channel value."""
        if not 0 <= code < self.rank:
            raise ValueError("channel outside its exact range")
        digits = [0] * self.fan.h
        for ell in range(self.k):
            code, item = divmod(code, self.image_per_coordinate)
            reduced = []
            if self.height <= self.small_base:
                for _ in range(self.r):
                    item, digit = divmod(item, self.height)
                    reduced.append(digit)
            elif item == self.fan.modulus - 1:
                # C^r is represented by C in the last retained position.
                reduced = [0] * (self.r - 1) + [self.small_base]
            else:
                for _ in range(self.r):
                    item, digit = divmod(item, self.small_base)
                    reduced.append(digit)
            for a in range(self.t):
                j, sign = self.signed_position(a)
                value = min(self.fan.base - 1, reduced[j])
                reduced[j] -= value
                digits[ell + self.k * a] = value if sign > 0 else self.fan.base - 1 - value
            assert not any(reduced)
        return tuple(digits)


def matrix_rank(rows: list[list[int]], prime: int) -> int:
    rows = [[x % prime for x in row] for row in rows]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [x * inverse % prime for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            multiplier = rows[i][column]
            if multiplier:
                rows[i] = [(x - multiplier * y) % prime
                           for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def interpolate_coefficients(fan: Fanout, indices: tuple[int, ...],
                             messages: tuple[int, ...]) -> tuple[int, ...]:
    """Independent Lagrange interpolation over the composite residue ring."""
    modulus = fan.modulus
    roots = [pow(fan.root, j, modulus) for j in indices]
    values = [value * pow(fan.phase(j), -1, modulus) % modulus
              for j, value in zip(indices, messages)]
    result = [0] * len(indices)
    for i, (root, value) in enumerate(zip(roots, values)):
        coefficients = [1]
        denominator = 1
        for j, other in enumerate(roots):
            if j == i:
                continue
            next_coefficients = [0] * (len(coefficients) + 1)
            for degree, coefficient in enumerate(coefficients):
                next_coefficients[degree] -= other * coefficient
                next_coefficients[degree + 1] += coefficient
            coefficients = [x % modulus for x in next_coefficients]
            denominator = denominator * (root - other) % modulus
        multiplier = value * pow(denominator, -1, modulus) % modulus
        result = [(a + multiplier * b) % modulus
                  for a, b in zip(result, coefficients)]
    return tuple(result)


def middle_output(fan: Fanout, position: int, left: tuple[int, ...],
                  right: tuple[int, ...]) -> tuple[int, ...]:
    """Original source hatOut, pointwise product, and faOut inverse stage."""
    if not fan.source_parameters:
        raise ValueError("this implementation follows the actual source parameters")
    modulus, a = fan.modulus, fan.a
    omega_inverse = pow(16, -1, modulus)
    transformed_left = [sum(x * pow(16, a * r * j, modulus)
                            for r, x in enumerate(left)) % modulus for j in range(a)]
    transformed_right = [sum(x * pow(16, a * r * j, modulus)
                             for r, x in enumerate(right)) % modulus for j in range(a)]
    return tuple(
        pow(omega_inverse, position * j, modulus)
        * sum(x * y * pow(omega_inverse, a * r * j, modulus)
              for r, (x, y) in enumerate(zip(transformed_left, transformed_right))) % modulus
        for j in range(a)
    )


def inverse_output(fan: Fanout, values: tuple[int, ...]) -> tuple[int, ...]:
    if not fan.source_parameters:
        raise ValueError("this implementation follows the actual source parameters")
    modulus, a = fan.modulus, fan.a
    omega_inverse = pow(16, -1, modulus)
    t_inverse = pow(a * a, -1, modulus)
    return tuple(t_inverse * sum(x * pow(omega_inverse, a * r * j, modulus)
                                 for r, x in enumerate(values)) % modulus
                 for j in range(a))


def run_checks() -> dict:
    rng = random.Random(20261031)
    report: dict = {
        "claim": "exact local fan-out ranks and codecs; no factoring score increase",
        "source_revision": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1",
        "exhaustive_image_cases": [],
    }
    controls = {"unit_differences": 0, "arbitrary_set_interpolations": 0,
                "codec_roundtrips": 0, "channel_preimages": 0,
                "complete_column_domains": 0, "equality_matrix_ranks": 0,
                "larger_source_cosets": 0, "radix_boundary_cases": 0,
                "middle_gate_vectors": 0, "inverse_gate_vectors": 0,
                "consistent_forward_surjections": 0}
    # Independent finite proof controls for the digit-image lemma, including
    # H=C and H=C+1 even when a boundary is absent from a source schedule.
    for radix in (2, 4, 8):
        for length in (1, 2, 3):
            modulus = radix ** length + 1
            for height in sorted({1, radix - 1, radix, radix + 1, radix + 2}):
                image = {sum(x * radix ** j for j, x in enumerate(digits)) % modulus
                         for digits in product(range(height), repeat=length)}
                expected = height ** length if height <= radix else modulus
                assert len(image) == expected
                controls["radix_boundary_cases"] += 1
    models = [Fanout(1, 2, 1), Fanout(2, 1, 1), Fanout(4, 1, 3),
              Fanout(4, 2, 7), Fanout(8, 2, 11), Fanout(2, 8, 3)]
    for fan in models:
        assert pow(fan.root, fan.h, fan.modulus) == fan.modulus - 1
        for left, right in combinations(range(fan.a), 2):
            difference = (pow(fan.root, left, fan.modulus)
                          - pow(fan.root, right, fan.modulus)) % fan.modulus
            assert difference * pow(difference, -1, fan.modulus) % fan.modulus == 1
            controls["unit_differences"] += 1
        for _ in range(20):
            indices = tuple(sorted(rng.sample(range(fan.a), fan.h)))
            digits = tuple(rng.randrange(fan.base) for _ in range(fan.h))
            all_messages = fan.messages(digits)
            assert interpolate_coefficients(fan, indices, tuple(all_messages[j] for j in indices)) == digits
            controls["arbitrary_set_interpolations"] += 1

        # One direct enumeration is reused only as raw input data. The rank
        # formula and the codec use separate folding and inverse transforms.
        all_digits = list(product(range(fan.base), repeat=fan.h))
        all_messages = [fan.messages(digits) for digits in all_digits]
        assert len(set(all_messages)) == fan.base ** fan.h
        for p in (1 << j for j in range(1, fan.a.bit_length())):
            exact = []
            for c in range(p):
                bundle = Bundle(fan, p, c)
                image = {tuple(values[j] for j in bundle.indices) for values in all_messages}
                assert len(image) == bundle.rank, (fan, p, c, len(image), bundle.rank)
                exact.append(len(image))
                sample_indices = range(len(all_digits)) if len(all_digits) <= 256 else rng.sample(range(len(all_digits)), 40)
                for index in sample_indices:
                    digits = all_digits[index]
                    values = tuple(all_messages[index][j] for j in bundle.indices)
                    code = bundle.row_code(digits)
                    assert bundle.column_code(values) == code
                    assert bundle.messages_from_code(code) == values
                    controls["codec_roundtrips"] += 1
                code_samples = list(range(bundle.rank)) if bundle.rank <= 256 else [0, bundle.rank - 1] + [rng.randrange(bundle.rank) for _ in range(30)]
                for code in code_samples:
                    digits = bundle.source_from_code(code)
                    assert bundle.row_code(digits) == code
                    values = tuple(fan.messages(digits)[j] for j in bundle.indices)
                    assert bundle.column_code(values) == code
                    assert bundle.messages_from_code(code) == values
                    assert values in image
                    controls["channel_preimages"] += 1
                if fan.modulus ** bundle.k <= 1000:
                    for values in product(range(fan.modulus), repeat=bundle.k):
                        assert (bundle.column_code(values) is not None) == (values in image)
                        controls["complete_column_domains"] += 1
                # Literal equality matrices, independently eliminated over
                # two unrelated coefficient fields in bounded cases.
                if len(all_digits) <= 16 and fan.modulus ** bundle.k <= 300:
                    columns = list(product(range(fan.modulus), repeat=bundle.k))
                    rows = [[int(tuple(values[j] for j in bundle.indices) == column)
                             for column in columns] for values in all_messages]
                    for prime in (2, 65521):
                        assert matrix_rank(rows, prime) == bundle.rank
                        controls["equality_matrix_ranks"] += 1
            report["exhaustive_image_cases"].append({
                "h": fan.h, "digit_bits": fan.digit_bits,
                "source_parameters": fan.source_parameters, "gate_index": fan.gate_index,
                "P": p, "source_vectors": len(all_digits), "coset_ranks": exact,
            })

    # Exact all-but-four classification and product lower bound, using large
    # integer products rather than floating logarithms.
    larger = []
    for h in (8, 16, 32, 64):
        fan = Fanout(h, 2 * h * h)
        for p in (1 << j for j in range(4, fan.a.bit_length())):
            rank_product = 1
            full = 0
            for c in range(p):
                bundle = Bundle(fan, p, c)
                if c in (0, p // 2):
                    expected = (p // 2) * (fan.base - 1) + 1
                elif c in (p // 4, 3 * p // 4):
                    expected = ((p // 4) * (fan.base - 1) + 1) ** 2
                else:
                    expected = fan.modulus
                    full += 1
                assert bundle.image_per_coordinate == expected
                rank_product *= bundle.rank
                controls["larger_source_cosets"] += 1
            assert full == p - 4
            floor_exponent = (8 * p - 20) * h * fan.digit_bits // p
            assert rank_product >= 1 << floor_exponent
            larger.append({"h": h, "digit_bits": fan.digit_bits, "P": p,
                           "full_message_cosets": full,
                           "sum_log_rank_floor": rank_product.bit_length() - 1,
                           "proved_lower_exponent": floor_exponent})
    report["larger_source_checks"] = larger
    assert (1 << 32) + 1 == 641 * 6700417
    for h in (1, 2, 8, 16):
        fan = Fanout(h, 2 * h * h)
        one = (1,) + (0,) * (fan.a - 1)
        for _ in range(8):
            position = rng.randrange(fan.a)
            left = tuple(rng.randrange(fan.modulus) for _ in range(fan.a))
            actual = middle_output(fan, position, left, one)
            expected = tuple(fan.a * pow(16, -position * j, fan.modulus) * x % fan.modulus
                             for j, x in enumerate(left))
            assert actual == expected
            controls["middle_gate_vectors"] += 1
            transformed = inverse_output(fan, left)
            recovered = tuple(
                fan.a * sum(x * pow(16, fan.a * r * j, fan.modulus)
                            for r, x in enumerate(transformed)) % fan.modulus
                for j in range(fan.a)
            )
            assert recovered == left
            controls["inverse_gate_vectors"] += 1
        if h < 8:
            continue
        positions = [c for c in range(fan.a)
                     if c not in (0, fan.a // 4, fan.a // 2, 3 * fan.a // 4)]
        for position in positions:
            desired = tuple(rng.randrange(fan.modulus) for _ in range(fan.a))
            incoming = tuple(pow(fan.a, -1, fan.modulus)
                             * pow(16, position * j, fan.modulus) * value % fan.modulus
                             for j, value in enumerate(desired))
            all_input_gates = []
            for gate_index, value in enumerate(incoming):
                source_fan = Fanout(h, 2 * h * h, gate_index)
                bundle = Bundle(source_fan, fan.a, position)
                assert bundle.rank == fan.modulus
                code = bundle.column_code((value,))
                assert code is not None
                digits = bundle.source_from_code(code)
                assert source_fan.message_at(digits, position) == value
                all_input_gates.append(digits)
            # The disjoint digit vectors constitute one consistent complete
            # factor input. Other MID recipients receive their own forced
            # evaluations; none is assigned an inconsistent independent value.
            realized = tuple(Fanout(h, 2 * h * h, r).message_at(digits, position)
                             for r, digits in enumerate(all_input_gates))
            assert middle_output(fan, position, realized, one) == desired
            controls["consistent_forward_surjections"] += 1
    report["controls"] = controls
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run_checks()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"exhaustive_image_cases": len(report["exhaustive_image_cases"]),
                      **report["controls"]}, sort_keys=True))


if __name__ == "__main__":
    main()
