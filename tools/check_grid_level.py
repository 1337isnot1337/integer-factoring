"""Exact finite integration certificate for one 2-by-2 grid step.

Compares the original four-node clique, the expanded signed pieces, and
independent row/column network contractions. This is not a fast factoring
implementation. Full grid checks are capped at lambda <= 2 by default.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
from itertools import product
import json
from math import comb
from pathlib import Path
import random
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.check_crossing import expanded_coefficient, psi, residual_left, residual_right


def _natural(name: str, value: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def scaled_piece(width: int, index: int, x: int, y: int, u: int, v: int) -> int:
    """4^t A_index, evaluated by the new rank-56 expansion."""

    t = (width + 2) // 3
    blocks = tuple(tuple((value >> (3 * s)) & 7 for value in (x, y, u, v)) for s in range(t))
    return expanded_coefficient(blocks, index)


def _digits(number: int, radix: int, length: int) -> tuple[int, ...]:
    result = []
    for _ in range(length):
        result.append(number % radix)
        number //= radix
    return tuple(result)


def _expanded_channel(width: int, index: int, label: int, first: int, second: int, left: bool) -> int:
    """Decode the actual 16^t 7^index channel, zero-padding binary labels."""

    t = (width + 2) // 3
    if label >= 16 ** t * 7 ** index:
        return 0
    mask = label % (2 ** t)
    rest = label // (2 ** t)
    classes = _digits(rest % (8 ** t), 8, t)
    labels = _digits(rest // (8 ** t), 7, index)
    chosen = mask.bit_count()
    if chosen > index or any(labels[chosen:]):
        return 0
    value = comb(t - chosen, index - chosen) * 4 ** (t - index) * (-3) ** (index - chosen) if left else 1
    selected = 0
    for s in range(t):
        p, q = (first >> (3 * s)) & 7, (second >> (3 * s)) & 7
        if mask >> s & 1:
            a = labels[selected] + 1
            selected += 1
            value *= residual_left(classes[s], a, p, q) if left else residual_right(classes[s], a, p, q)
        else:
            value *= psi(classes[s], p, q) if left else psi(classes[s], q, p)
    return value


@dataclass(frozen=True)
class PieceChannels:
    width: int
    index: int
    orientation: str
    channel_bits: int
    fallback: bool
    # Indexed by channel label and then by the pair first*B+second.
    left: tuple[tuple[int, ...], ...]
    right: tuple[tuple[int, ...], ...]


@lru_cache(maxsize=32)
def make_channels(width: int, index: int, orientation: str) -> PieceChannels:
    """Build tiny explicit integer channel tables, using source block factors."""

    _natural("width", width)
    _natural("index", index)
    if width > 3:
        raise ValueError("finite channel table construction is limited to width <= 3")
    if orientation not in ("row", "column"):
        raise ValueError("orientation must be row or column")
    t, b = (width + 2) // 3, 1 << width
    if index > t:
        raise ValueError("piece index exceeds the number of blocks")
    j = index if orientation == "row" else t - index
    expanded_bits = 4 * t + (73 * j + 25) // 26
    bits = min(2 * width, expanded_bits)
    fallback = 2 * width <= expanded_bits
    left, right = [], []
    for label in range(1 << bits):
        if fallback:
            p0, q0 = divmod(label, b)
            left.append(tuple(int(first == p0 and second == q0) for first, second in product(range(b), repeat=2)))
            right.append(tuple(scaled_piece(width, j, p0, q0, first, second) for first, second in product(range(b), repeat=2)))
        else:
            left.append(tuple(_expanded_channel(width, j, label, first, second, True) for first, second in product(range(b), repeat=2)))
            right.append(tuple(_expanded_channel(width, j, label, first, second, False) for first, second in product(range(b), repeat=2)))
    return PieceChannels(width, index, orientation, bits, fallback, tuple(left), tuple(right))


def channel_value(channels: PieceChannels, x: int, y: int, u: int, v: int) -> int:
    b = 1 << channels.width
    a, c = (x * b + y, u * b + v) if channels.orientation == "row" else (x * b + u, y * b + v)
    return sum(left[a] * right[c] for left, right in zip(channels.left, channels.right))


def compressed_channel_value(width: int, index: int, orientation: str, x: int, y: int, u: int, v: int,
                             *, max_channels: int = 4096) -> int:
    """Check a compressed channel without materializing every pair of strings."""

    _natural("width", width)
    _natural("index", index)
    _natural("max_channels", max_channels)
    if orientation not in ("row", "column"):
        raise ValueError("orientation must be row or column")
    t = (width + 2) // 3
    if index > t:
        raise ValueError("piece index exceeds the number of blocks")
    j = index if orientation == "row" else t - index
    bits = 4 * t + (73 * j + 25) // 26
    if not max_channels or bits >= max_channels.bit_length():
        raise ValueError("compressed channel count exceeds the explicit cap")
    p, q, r, s = (x, y, u, v) if orientation == "row" else (x, u, y, v)
    return sum(_expanded_channel(width, j, label, p, q, True)
               * _expanded_channel(width, j, label, r, s, False) for label in range(1 << bits))


def _entry(table: tuple[int, ...], b: int, first: int, second: int, third: int) -> int:
    return table[(first * b + second) * b + third]


def original_contraction(tables: tuple[tuple[int, ...], ...], width: int) -> int:
    """Sum over undirected edge labels (01,02,03,12,13,23)."""

    b = 1 << width
    t0, t1, t2, t3 = tables
    total = 0
    for a, c, x, y, e, f in product(range(b), repeat=6):
        total += (_entry(t0, b, a, c, x) * _entry(t1, b, a, y, e)
                  * _entry(t2, b, c, y, f) * _entry(t3, b, x, e, f))
    return total


def expanded_contractions(tables: tuple[tuple[int, ...], ...], width: int) -> tuple[int, ...]:
    """Independent direct sum after freeing the two diagonal edge copies."""

    b, t = 1 << width, (width + 2) // 3
    t0, t1, t2, t3 = tables
    pieces = {
        (x, y, u, v): tuple(scaled_piece(width, i, x, y, u, v) for i in range(t + 1))
        for x, y, u, v in product(range(b), repeat=4)
    }
    totals = [0] * (t + 1)
    for a, c, e, f, x, y, u, v in product(range(b), repeat=8):
        weight = (_entry(t0, b, a, c, x) * _entry(t1, b, a, y, e)
                  * _entry(t2, b, c, u, f) * _entry(t3, b, v, e, f))
        if weight:
            for i, value in enumerate(pieces[x, y, u, v]):
                totals[i] += weight * value
    return tuple(totals)


def _mod(value: int, modulus: int | None) -> int:
    return value if modulus is None else value % modulus


def network_contraction(tables: tuple[tuple[int, ...], ...], channels: PieceChannels, modulus: int | None = None) -> int:
    """Contract actual row or column nodes through their factor channels."""

    if modulus is not None:
        _natural("modulus", modulus)
        if modulus < 2:
            raise ValueError("modulus must be at least 2")
    b = 1 << channels.width
    t0, t1, t2, t3 = tables
    total = 0
    for first_shared, second_shared in product(range(b), repeat=2):
        for left_channel, right_channel in zip(channels.left, channels.right):
            left_value, right_value = 0, 0
            if channels.orientation == "row":
                c, e = first_shared, second_shared
                for a, x, y in product(range(b), repeat=3):
                    left_value += _entry(t0, b, a, c, x) * _entry(t1, b, a, y, e) * left_channel[x * b + y]
                for f, u, v in product(range(b), repeat=3):
                    right_value += _entry(t2, b, c, u, f) * _entry(t3, b, v, e, f) * right_channel[u * b + v]
            else:
                a, f = first_shared, second_shared
                for c, x, u in product(range(b), repeat=3):
                    left_value += _entry(t0, b, a, c, x) * _entry(t2, b, c, u, f) * left_channel[x * b + u]
                for e, y, v in product(range(b), repeat=3):
                    right_value += _entry(t1, b, a, y, e) * _entry(t3, b, v, e, f) * right_channel[y * b + v]
            total = _mod(total + _mod(left_value, modulus) * _mod(right_value, modulus), modulus)
    return total


@dataclass(frozen=True)
class GridCertificate:
    width: int
    original_states: int
    expanded_states: int
    original_value: int
    multiplier: int
    expected_scaled_value: int
    direct_pieces: tuple[int, ...]
    row_pieces: tuple[int, ...]
    column_pieces: tuple[int, ...]
    channel_bits_row: tuple[int, ...]
    channel_bits_column: tuple[int, ...]
    fallback_row: tuple[bool, ...]
    fallback_column: tuple[bool, ...]
    modular_values: tuple[tuple[int, tuple[int, ...], tuple[int, ...]], ...]
    valid: bool


def verify_grid_level(tables: tuple[tuple[int, ...], ...], width: int, *, max_expanded_states: int = 100_000,
                      moduli: tuple[int, ...] = (4, 6, 8, 9, 15)) -> GridCertificate:
    """Run an exact bounded one-level certificate; reject oversized inputs."""

    _natural("width", width)
    _natural("max_expanded_states", max_expanded_states)
    if not max_expanded_states:
        raise ValueError("max_expanded_states must be positive")
    # Compare the exponent to the cap's bit length before constructing 2^(8w).
    if 8 * width >= max_expanded_states.bit_length() or width > 3:
        raise ValueError("expanded state count exceeds the explicit finite-check cap")
    expanded_states, b, t = 1 << (8 * width), 1 << width, (width + 2) // 3
    if len(tables) != 4 or any(len(table) != b ** 3 for table in tables):
        raise ValueError("provide four integer tables, each with (2^width)^3 entries")
    if any(not isinstance(value, int) or isinstance(value, bool) for table in tables for value in table):
        raise ValueError("table entries must be integers")
    for modulus in moduli:
        if not isinstance(modulus, int) or isinstance(modulus, bool) or modulus < 2:
            raise ValueError("moduli must be integers at least 2")
    original = original_contraction(tables, width)
    direct = expanded_contractions(tables, width)
    rows = tuple(make_channels(width, i, "row") for i in range(t + 1))
    columns = tuple(make_channels(width, i, "column") for i in range(t + 1))
    row_values = tuple(network_contraction(tables, channel) for channel in rows)
    column_values = tuple(network_contraction(tables, channel) for channel in columns)
    modular_values = tuple((m, tuple(network_contraction(tables, c, m) for c in rows),
                            tuple(network_contraction(tables, c, m) for c in columns)) for m in moduli)
    multiplier = 8 ** t
    valid = direct == row_values == column_values and sum(direct) == multiplier * original
    valid = valid and all(r == c == tuple(value % m for value in direct) for m, r, c in modular_values)
    return GridCertificate(
        width, b ** 6, expanded_states, original, multiplier, multiplier * original, direct,
        row_values, column_values, tuple(c.channel_bits for c in rows), tuple(c.channel_bits for c in columns),
        tuple(c.fallback for c in rows), tuple(c.fallback for c in columns), modular_values, valid,
    )


def signed_tables(width: int, seed: int) -> tuple[tuple[int, ...], ...]:
    rng = random.Random(seed)
    return tuple(tuple(rng.randrange(-3, 4) for _ in range(1 << (3 * width))) for _ in range(4))


def run_suite() -> dict:
    cases = []
    for width, seeds in ((1, (100, 101, 102, 103)), (2, (200, 201))):
        b = 1 << width
        fixtures = [("all_ones", tuple((1,) * (b ** 3) for _ in range(4)))]
        fixtures += [(f"signed_seed_{seed}", signed_tables(width, seed)) for seed in seeds]
        for name, tables in fixtures:
            certificate = verify_grid_level(tables, width)
            if not certificate.valid:
                raise AssertionError((name, certificate))
            cases.append({"name": name, "tables": tables, "certificate": asdict(certificate)})
    return {
        "status": "finite one-level integration checks only; no full factoring machine is certified",
        "grid_side": 2,
        "vertices": 4,
        "case_count": len(cases),
        "composite_moduli": (4, 6, 8, 9, 15),
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="include input tables and exact certificate values")
    args = parser.parse_args()
    report = run_suite()
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"{report['case_count']} exact one-level grid certificates passed")
        for case in report["cases"]:
            certificate = case["certificate"]
            print(f"lambda={certificate['width']} {case['name']}: Z={certificate['original_value']}, "
                  f"pieces={certificate['direct_pieces']}, multiplier={certificate['multiplier']}")


if __name__ == "__main__":
    main()
