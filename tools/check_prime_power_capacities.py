#!/usr/bin/env python3
"""Exact finite checks for prime-power capacity composition.

The all-length bound is proved in research/integral-schur-crossing.md.
This checker tests integer capacities, the greedy construction, binary
label splitting, and complete mixed-radix piece-index encoding. It does
not infer a crossing rank or an asymptotic factoring rate from these cases.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path


class PartitionTable:
    """Euler's pentagonal recurrence, extended only as needed."""

    def __init__(self) -> None:
        self.values = [1]

    def extend(self, end: int) -> None:
        for n in range(len(self.values), end + 1):
            total = 0
            k = 1
            while True:
                left = k * (3 * k - 1) // 2
                if left > n:
                    break
                sign = 1 if k % 2 else -1
                total += sign * self.values[n - left]
                right = k * (3 * k + 1) // 2
                if right <= n:
                    total += sign * self.values[n - right]
                k += 1
            assert total > 0
            self.values.append(total)

    def __getitem__(self, n: int) -> int:
        self.extend(n)
        return self.values[n]


@dataclass(frozen=True)
class Block:
    capacity: int
    energy: int
    restricted_to_one_bit: bool = False

    @property
    def pieces(self) -> int:
        return 4 * self.energy + 1


def catalog(prime: int, maximum: int, table: PartitionTable) -> list[Block]:
    """Include one capacity above maximum, plus a restricted one-bit block."""
    by_capacity = {1: Block(1, prime, True)}
    energy = prime
    while True:
        dimension = table[energy]
        capacity = dimension.bit_length() - 1
        assert (1 << capacity) <= dimension < (1 << (capacity + 1))
        by_capacity[capacity] = Block(capacity, energy)
        if capacity > maximum:
            break
        energy *= prime
    return [by_capacity[key] for key in sorted(by_capacity)]


def decompose(length: int, blocks: list[Block]) -> list[Block]:
    capacities = [block.capacity for block in blocks]
    remaining = length
    result = []
    while remaining:
        position = bisect_right(capacities, remaining) - 1
        if position < 0:
            raise ValueError("no block fits the remaining positive capacity")
        block = blocks[position]
        result.append(block)
        remaining -= block.capacity
    return result


def validate_decomposition(length: int, chosen: list[Block], blocks: list[Block]) -> None:
    assert sum(block.capacity for block in chosen) == length
    assert all(block in blocks for block in chosen)
    counts = Counter(block.capacity for block in chosen)
    # Independently check the exact consequence of a remainder smaller
    # than the next capacity, including the fixed one-bit block.
    for position, block in enumerate(blocks[:-1]):
        next_capacity = blocks[position + 1].capacity
        assert counts[block.capacity] <= (next_capacity - 1) // block.capacity


def split_label(value: int, chosen: list[Block]) -> list[int]:
    parts = []
    for block in chosen:
        parts.append(value & ((1 << block.capacity) - 1))
        value >>= block.capacity
    assert value == 0
    return parts


def join_label(parts: list[int], chosen: list[Block]) -> int:
    value = 0
    for part, block in reversed(list(zip(parts, chosen))):
        assert 0 <= part < (1 << block.capacity)
        value = (value << block.capacity) | part
    return value


def decode_piece(value: int, chosen: list[Block]) -> list[int]:
    shifts = []
    for block in chosen:
        value, digit = divmod(value, block.pieces)
        shifts.append(digit - 2 * block.energy)
    assert value == 0
    return shifts


def encode_piece(shifts: list[int], chosen: list[Block]) -> int:
    value = 0
    for shift, block in reversed(list(zip(shifts, chosen))):
        assert -2 * block.energy <= shift <= 2 * block.energy
        value = value * block.pieces + shift + 2 * block.energy
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-label-length", type=int, default=256)
    parser.add_argument("--output", type=Path, default=Path("evidence/prime-power-capacities.json"))
    args = parser.parse_args()
    assert args.max_label_length >= 37
    table = PartitionTable()

    # A different partition algorithm checks the pentagonal recurrence.
    check_end = 128
    coin_counts = [1] + [0] * check_end
    for part in range(1, check_end + 1):
        for weight in range(part, check_end + 1):
            coin_counts[weight] += coin_counts[weight - part]
    table.extend(check_end)
    assert table.values[: check_end + 1] == coin_counts

    decomposition_cases = label_cases = piece_cases = 0
    negative_controls = []
    records = []
    for prime in (2, 3, 5, 7):
        blocks = catalog(prime, args.max_label_length, table)
        max_blocks = max_entropy_upper = max_index_bits = 0
        for length in range(1, args.max_label_length + 1):
            chosen = decompose(length, blocks)
            validate_decomposition(length, chosen, blocks)
            decomposition_cases += 1
            max_blocks = max(max_blocks, len(chosen))

            number_of_pieces = 1
            entropy_upper = 0
            for block in chosen:
                number_of_pieces *= block.pieces
                entropy_upper += (block.pieces - 1).bit_length()
            assert number_of_pieces <= (1 << entropy_upper)
            max_entropy_upper = max(max_entropy_upper, entropy_upper)
            max_index_bits = max(max_index_bits, (number_of_pieces - 1).bit_length())

            labels = {0, (1 << length) - 1}
            labels.update(1 << position for position in range(length))
            # Alternating patterns cross every variable block boundary.
            labels.add(sum(1 << position for position in range(0, length, 2)))
            labels.add(sum(1 << position for position in range(1, length, 2)))
            for label in labels:
                assert join_label(split_label(label, chosen), chosen) == label
                label_cases += 1

            indices = {0, number_of_pieces - 1, number_of_pieces // 2}
            if length <= 8:
                # Complete enumeration is checked on these small products.
                indices.update(range(number_of_pieces))
            for index in indices:
                shifts = decode_piece(index, chosen)
                assert encode_piece(shifts, chosen) == index
                piece_cases += 1

        # A mistaken capacity total must not pass the independent validator.
        broken = decompose(37, blocks)[1:]
        try:
            validate_decomposition(37, broken, blocks)
        except AssertionError:
            negative_controls.append({"prime": prime, "control": "removed_block_rejected"})
        else:
            raise AssertionError("negative control accepted a wrong total capacity")

        natural = [block for block in blocks if not block.restricted_to_one_bit]
        if natural[0].capacity > 1:
            try:
                decompose(1, natural)
            except ValueError:
                negative_controls.append({"prime": prime, "control": "missing_one_bit_block_rejected"})
            else:
                raise AssertionError("negative control unexpectedly filled capacity one")

        records.append({
            "prime": prime,
            "max_label_length": args.max_label_length,
            "max_blocks_used": max_blocks,
            "max_sum_of_local_index_widths": max_entropy_upper,
            "max_mixed_radix_index_width": max_index_bits,
            "catalog": [{
                "energy": block.energy,
                "dimension": str(table[block.energy]),
                "capacity": block.capacity,
                "restricted_to_one_bit": block.restricted_to_one_bit,
                "piece_count": block.pieces,
            } for block in blocks],
        })

    result = {
        "claim": "Finite checks of the exact prime-power block-capacity construction; no rank inference.",
        "status": "passed",
        "partition_cross_checks": check_end + 1,
        "decomposition_cases": decomposition_cases,
        "label_codec_cases": label_cases,
        "piece_index_codec_cases": piece_cases,
        "negative_controls": negative_controls,
        "maximum_partition_energy_computed": len(table.values) - 1,
        "prime_records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "prime_records"}, indent=2))


if __name__ == "__main__":
    main()
