#!/usr/bin/env python3
"""Exact finite checks for the parity/power-divisor factoring reduction.

The counting oracles below deliberately enumerate small divisor sets. They
are independent finite reference oracles, not an implementation of the
asymptotic clique evaluator. Only Python's standard library is used.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path
import random


def factor_reference(value: int) -> list[int]:
    """Independent ordinary trial division, including repeated factors."""
    result: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result.append(divisor)
            value //= divisor
        divisor += 1
    if value > 1:
        result.append(value)
    return result


@lru_cache(maxsize=None)
def divisors_reference(value: int) -> tuple[int, ...]:
    """Enumerate by divisibility, without using the exponent formula."""
    low: list[int] = []
    high: list[int] = []
    for divisor in range(1, isqrt(value) + 1):
        if value % divisor == 0:
            low.append(divisor)
            if divisor * divisor != value:
                high.append(value // divisor)
    return tuple(low + high[::-1])


def power_divisors_reference(value: int, exponent_bit: int) -> tuple[int, ...]:
    power = 1 << exponent_bit
    result: list[int] = []
    for divisor in divisors_reference(value):
        if divisor <= 1:
            continue
        raised = divisor**power
        if raised <= value and value % raised == 0:
            result.append(divisor)
    return tuple(result)


def prefix_blocks(limit: int, width: int) -> list[tuple[int, int]]:
    """Disjoint (free-low-bits, high-prefix) blocks covering [0,limit)."""
    assert 0 <= limit < 1 << width
    fixed = 0
    result: list[tuple[int, int]] = []
    for position in range(width - 1, -1, -1):
        bit = (limit >> position) & 1
        if bit:
            result.append((position, fixed * 2))
        fixed = fixed * 2 + bit
    return result


def query_parity(values: tuple[int, ...], free: int, prefix: int) -> int:
    return sum((value >> free) == prefix for value in values) & 1


def recover_odd_prefix(
    values: tuple[int, ...], free: int, prefix: int
) -> tuple[int, int]:
    assert query_parity(values, free, prefix) == 1
    calls = 0
    while free:
        free -= 1
        candidate = prefix * 2
        left = query_parity(values, free, candidate)
        calls += 1
        prefix = candidate if left else candidate + 1
        assert query_parity(values, free, prefix) == 1
    assert prefix in values
    return prefix, calls


def split_or_prime(value: int) -> tuple[int | None, int, str]:
    """The proved reduction, using finite enumerated parity oracles."""
    width = value.bit_length()
    calls = 0
    last_bit = (width - 1).bit_length() - 1
    for exponent_bit in range(1, last_bit + 1):
        values = power_divisors_reference(value, exponent_bit)
        base_width = (width + (1 << exponent_bit) - 1) >> exponent_bit
        calls += 1
        if query_parity(values, base_width, 0):
            divisor, used = recover_odd_prefix(values, base_width, 0)
            assert 1 < divisor < value and value % divisor == 0
            return divisor, calls + used, "power"

    values = tuple(divisor for divisor in divisors_reference(value) if divisor > 1)
    for free, prefix in prefix_blocks(isqrt(value) + 1, width):
        calls += 1
        if query_parity(values, free, prefix):
            divisor, used = recover_odd_prefix(values, free, prefix)
            assert 1 < divisor <= isqrt(value) and value % divisor == 0
            return divisor, calls + used, "interval"
    return None, calls, "prime"


def factor_by_parities(value: int) -> tuple[list[int], int, int]:
    pending = [value]
    primes: list[int] = []
    calls = 0
    nodes = 0
    while pending:
        current = pending.pop()
        nodes += 1
        divisor, used, _ = split_or_prime(current)
        calls += used
        if divisor is None:
            primes.append(current)
        else:
            pending.extend((divisor, current // divisor))
    return sorted(primes), calls, nodes


def power_chain_widths(width: int, exponent_bit: int) -> list[int]:
    return [
        (width + (1 << (exponent_bit - index)) - 1)
        >> (exponent_bit - index)
        for index in range(exponent_bit + 1)
    ]


def least_source_scale(width: int) -> tuple[int, int, int, int]:
    scale = 0
    while 1 << (4 * scale + 2) < width:
        scale += 1
    root = 1 << scale
    digit_width = 2 * root * root
    return scale, digit_width, 2 * root, root


def digit_location(bit: int, digit_width: int, layer_width: int) -> tuple[int, int, int]:
    digit, within = divmod(bit, digit_width)
    gate, slot = digit % layer_width, digit // layer_width
    return gate, slot, within


def product_output_digits(left: int, right: int, digit_width: int) -> list[int]:
    """Literal integer digit convolution and the source's carry rule."""
    base = 1 << digit_width
    digits_left = [(left >> (digit_width * i)) % base for i in range(digit_width)]
    digits_right = [(right >> (digit_width * i)) % base for i in range(digit_width)]
    carry = 0
    output: list[int] = []
    for column in range(2 * digit_width):
        coefficient = sum(
            digits_left[i] * digits_right[column - i]
            for i in range(digit_width)
            if 0 <= column - i < digit_width
        )
        assert coefficient < (1 << (4 * digit_width)) + 1
        total = carry + coefficient
        output.append(total % base)
        carry = total // base
    assert carry == 0
    return output


def output_bit(output: list[int], digit_width: int, bit: int) -> int:
    digit, within = divmod(bit, digit_width)
    return (output[digit] >> within) & 1


def check_arithmetic(limit: int) -> dict[str, int]:
    power_cases = 0
    squarefree_cases = 0
    prime_cases = 0
    route_counts: Counter[str] = Counter()
    maximum_calls = 0
    for value in range(2, limit + 1):
        width = value.bit_length()
        exponents = Counter(factor_reference(value))
        last_bit = (width - 1).bit_length() - 1
        parities: list[int] = []
        for exponent_bit in range(1, last_bit + 1):
            values = power_divisors_reference(value, exponent_bit)
            expected = 1
            for exponent in exponents.values():
                expected *= (exponent >> exponent_bit) + 1
            expected -= 1
            assert len(values) == expected
            parity = len(values) & 1
            assert parity == int(any((a >> exponent_bit) & 1 for a in exponents.values()))
            parities.append(parity)
            power_cases += 1
        squarefree = all(a == 1 for a in exponents.values())
        assert (not any(parities)) == squarefree
        if squarefree:
            count = sum(1 < d <= isqrt(value) for d in divisors_reference(value))
            assert count == (1 << (len(exponents) - 1)) - 1
            squarefree_cases += 1
        divisor, _, route = split_or_prime(value)
        route_counts[route] += 1
        if divisor is None:
            assert factor_reference(value) == [value]
            prime_cases += 1
        actual, calls, nodes = factor_by_parities(value)
        assert actual == factor_reference(value)
        assert nodes <= 2 * width - 1
        assert calls <= (2 * width - 1) * (last_bit + 2 * width)
        maximum_calls = max(maximum_calls, calls)
    return {
        "integers": limit - 1,
        "power_divisor_cases": power_cases,
        "squarefree_cases": squarefree_cases,
        "prime_declarations": prime_cases,
        "power_route_inputs": route_counts["power"],
        "interval_route_inputs": route_counts["interval"],
        "maximum_complete_factorization_queries": maximum_calls,
    }


def check_prefixes(limit: int) -> dict[str, int]:
    blocks = 0
    for endpoint in range(1, limit + 1):
        width = endpoint.bit_length()
        covered: list[int] = []
        for free, prefix in prefix_blocks(endpoint, width):
            covered.extend(range(prefix << free, (prefix + 1) << free))
            blocks += 1
        assert covered == list(range(endpoint))
        assert len(prefix_blocks(endpoint, width)) <= width
    selected, calls = recover_odd_prefix((2, 3, 5), 3, 0)
    assert selected == 5 and selected != 2
    return {"endpoints": limit, "blocks": blocks, "nonminimum_witness": selected, "witness_queries": calls}


def check_chain_budgets() -> dict[str, int]:
    widths = list(range(2, 1025)) + [10**k for k in range(4, 13)]
    cases = 0
    for width in widths:
        last_bit = (width - 1).bit_length() - 1
        for exponent_bit in range(1, last_bit + 1):
            sizes = power_chain_widths(width, exponent_bit)
            assert sizes[-1] == width
            requested = width + sum(sizes[1:])
            wires = sizes[0] + 2 * sum(sizes[1:-1]) + sizes[-1]
            assert requested < 3 * width + exponent_bit <= 4 * width
            assert wires <= 3 * width + 2 * exponent_bit
            capacities = [least_source_scale(s)[1] ** 2 for s in [width, *sizes[1:]]]
            assert all(cap <= 16 * size for cap, size in zip(capacities, [width, *sizes[1:]]))
            assert sum(capacities) <= 64 * width
            cases += 1
    return {"bit_lengths": len(widths), "chain_cases": cases}


def check_ports() -> dict[str, int]:
    generator = random.Random(20261030)
    products = 0
    port_bits = 0
    unequal_layouts = 0
    for left_width in (1, 2, 3, 4, 5, 8, 13, 17, 31, 49, 65):
        for right_width in (1, 2, 4, 7, 16, 33):
            result_width = left_width + right_width
            _, source_digit, _, _ = least_source_scale(result_width)
            for _ in range(5):
                left = generator.randrange(1 << left_width)
                right = generator.randrange(1 << right_width)
                output = product_output_digits(left, right, source_digit)
                recovered = sum(digit << (source_digit * index) for index, digit in enumerate(output))
                assert recovered == left * right
                # Next gadget can have a different digit width and layout.
                _, target_digit, target_layer, target_slots = least_source_scale(2 * result_width)
                if target_digit != source_digit:
                    unequal_layouts += 1
                target = [[0] * target_slots for _ in range(target_layer)]
                for bit in range(result_width):
                    gate, slot, within = digit_location(bit, target_digit, target_layer)
                    target[gate][slot] |= output_bit(output, source_digit, bit) << within
                    port_bits += 1
                flat = sum(
                    target[gate][slot] << (target_digit * (gate + target_layer * slot))
                    for gate in range(target_layer)
                    for slot in range(target_slots)
                )
                assert flat == left * right
                assert all(output_bit(output, source_digit, bit) == 0 for bit in range(result_width, 2 * source_digit**2))
                products += 1
    # Literal arithmetic-chain witness counts on small positive inputs.
    chain_cases = 0
    for value in range(2, 129):
        width = value.bit_length()
        for exponent_bit in range(1, (width - 1).bit_length()):
            sizes = power_chain_widths(width, exponent_bit)
            successful: list[int] = []
            for divisor in range(2, 1 << sizes[0]):
                current = divisor
                valid = True
                for target_width in sizes[1:]:
                    _, digit, _, _ = least_source_scale(target_width)
                    output = product_output_digits(current, current, digit)
                    squared = sum(entry << (digit * i) for i, entry in enumerate(output))
                    assert squared == current * current
                    current = squared
                    if current >= 1 << target_width:
                        valid = False
                        break
                if valid and current <= value and value % current == 0:
                    successful.append(divisor)
            assert tuple(successful) == power_divisors_reference(value, exponent_bit)
            chain_cases += 1
    return {"product_cases": products, "wired_bits": port_bits, "unequal_digit_layouts": unequal_layouts, "arithmetic_chain_cases": chain_cases}


def check_negative_controls() -> dict[str, object]:
    # Keeping the divisor one destroys the prime predicate: it is an odd
    # power-divisor count even on a prime input.
    assert sum(11 % (d * d) == 0 for d in divisors_reference(11)) & 1
    assert power_divisors_reference(11, 1) == ()
    # Missing the equality between the two inputs of the first square
    # allows d*a*e=N. Its nonunit-d count is odd at N=6 although no square
    # greater than one divides six.
    omitted_copy = [(d, a, e) for d, a, e in product(range(1, 7), repeat=3) if d > 1 and d * a * e == 6]
    assert len(omitted_copy) == 5
    assert power_divisors_reference(6, 1) == ()
    # Missing the link between d*d and the final first operand allows
    # d=2, square output 4, independent final first operand 3, cofactor2.
    assert 2 * 2 == 4 and 3 * 2 == 6 and 4 != 3
    return {
        "unremoved_one_false_positive": 11,
        "missing_square_copy_input": 6,
        "missing_square_copy_witnesses": len(omitted_copy),
        "missing_output_link": {"d": 2, "square_output": 4, "final_first_operand": 3, "cofactor": 2, "N": 6},
        "parity_search_is_not_minimum": {"set": [2, 3, 5], "returned": 5},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=8192)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert args.max_n >= 2
    report = {
        "status": "PASS",
        "scope": "Finite exact recovery and port checks; no asymptotic evaluator is executed.",
        "arithmetic": check_arithmetic(args.max_n),
        "prefixes": check_prefixes(1024),
        "chain_budgets": check_chain_budgets(),
        "ports": check_ports(),
        "negative_controls": check_negative_controls(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
