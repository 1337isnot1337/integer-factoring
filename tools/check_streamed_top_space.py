#!/usr/bin/env python3
"""Exact checks for packed, constant-depth original top-table queries."""

import argparse
from itertools import combinations
import json
from pathlib import Path
from random import Random

from check_polynomial_top_queries import Network, interpolate
from streamed_top_queries import BitTape, StreamedTop, edges, lagrange_inplace, pack_boundary


def f4_mul(left, right):
    value = 0
    for bit in range(2):
        if (right >> bit) & 1:
            value ^= left << bit
    if value & 4:
        value ^= 7
    return value


def ring_mul(left, right):
    a, b = left & 3, left >> 2
    c, d = right & 3, right >> 2
    return f4_mul(a, c) | ((f4_mul(a, d) ^ f4_mul(b, c)) << 2)


def bit_weight(target_kind, target_index):
    def query(kind, index, bit):
        if not bit:
            return 1
        phase = (1, 2, 3)[(index + (kind == "q")) % 3]
        return phase | ((phase if (kind, index) == (target_kind, target_index) else 0) << 2)
    return query


def reference_weight(network, cluster, setting, weight):
    if setting is None:
        return None
    result = 1
    for kind, row in network.gates:
        if (kind, row) not in cluster or kind not in ("p", "q"):
            continue
        for position, digit in enumerate(setting[(kind, row)]):
            for bit in range(network.digit_bits):
                index = network.digit_bits * (row + network.a * position) + bit
                result = ring_mul(result, weight(kind, index, (digit >> bit) & 1))
    return result


def verify():
    rng = Random(2026103008)
    interpolation_cases = boundary_addresses = valid_queries = malformed_queries = 0
    prefix_queries = padding_queries = weighted_second_factor = maximum_depth = 0
    recomputed_interpolations = 0
    for s in range(4):
        network = Network(s)
        for trial in range(8):
            selected = tuple(sorted(rng.sample(range(network.a), network.h)))
            points = [pow(network.xi, index, network.modulus) for index in selected]
            values = [rng.randrange(network.modulus) for _ in selected]
            actual = lagrange_inplace(points, values, network.modulus)
            assert tuple(actual) == interpolate(points, values, network.modulus)
            for point, value in zip(points, values):
                assert sum(coefficient * pow(point, degree, network.modulus)
                           for degree, coefficient in enumerate(actual)) % network.modulus == value
            interpolation_cases += 1

    for s in range(3):
        network = Network(s, toy_base=2)
        for trial in range(2):
            first = [rng.randrange(network.base) for _ in range(network.a * network.h)]
            second = [rng.randrange(network.base) for _ in range(network.a * network.h)]
            first[0] = second[0] = 1
            settings, outgoing, number, first_number = network.global_setting(first, second)
            number_tape, zero_tape = BitTape.from_int(number), BitTape.from_int(0)
            for groups in (1 << power for power in range(1, s + 2)):
                grouping = network.grouping(groups)
                for group in range(groups):
                    cluster = {gate for gate in network.gates if grouping[gate] == group}
                    boundary = network.boundary(cluster, settings, outgoing)
                    tape = pack_boundary(network.a, groups, group, network.f, boundary)
                    streamed = StreamedTop(s, groups, group, tape, toy_base=2)
                    assert not hasattr(streamed, "gates") and not hasattr(streamed, "edges")
                    for pair, value in boundary.items():
                        assert streamed.read_boundary(*pair) == value
                        boundary_addresses += 1
                    assert streamed.boundary_bits == (groups - 1) * (
                        8 * (network.a // groups) ** 2 + 1) * network.f
                    kind = "q" if (trial + group) % 2 else "p"
                    target = rng.randrange(network.input_bits)
                    weight = bit_weight(kind, target)
                    expected = reference_weight(network, cluster, settings, weight)
                    assert streamed.query(number_tape, network.input_bits, zero_tape, weight, ring_mul) == expected
                    assert streamed.max_arithmetic_depth <= 3
                    maximum_depth = max(maximum_depth, streamed.max_arithmetic_depth)
                    recomputed_interpolations += streamed.interpolations
                    valid_queries += 1
                    weighted_second_factor += kind == "q"

                    # A valid prefix, then a changed prefix. A group can accept
                    # the latter if its particular input digits are unaffected.
                    cutoff = network.input_bits // 2
                    prefix = first_number >> cutoff
                    for candidate in (prefix, prefix ^ 1):
                        expected_setting = network.query(cluster, boundary, number, cutoff, candidate)
                        expected = reference_weight(network, cluster, expected_setting, weight)
                        assert streamed.query(number_tape, cutoff, BitTape.from_int(candidate), weight, ring_mul) == expected
                        prefix_queries += 1

                    # Neither mutation is assumed to be invalid by fiat. The
                    # independent local-setting query decides its true value.
                    keys = list(boundary)
                    for pair, replacement in (
                            (keys[rng.randrange(len(keys))], None),
                            (keys[rng.randrange(len(keys))], network.modulus)):
                        changed = dict(boundary)
                        changed[pair] = changed[pair] ^ 1 if replacement is None else replacement
                        changed_tape = pack_boundary(network.a, groups, group, network.f, changed)
                        altered = StreamedTop(s, groups, group, changed_tape, toy_base=2)
                        expected_setting = network.query(cluster, changed, number, network.input_bits, 0)
                        expected = reference_weight(network, cluster, expected_setting, weight)
                        assert altered.query(number_tape, network.input_bits, zero_tape, weight, ring_mul) == expected
                        malformed_queries += 1

                    for other in range(groups):
                        if other == group or abs(other - group) == 1:
                            continue
                        offset = streamed.block_position(other) * streamed.edge_bits + (
                            8 * streamed.length ** 2) * streamed.f
                        altered = StreamedTop(s, groups, group, tape | (1 << offset), toy_base=2)
                        assert altered.query(number_tape, network.input_bits, zero_tape) is None
                        padding_queries += 1

    # The edge generator has exactly the original un-oriented edge set.
    edge_sets = 0
    for s in range(5):
        network = Network(s)
        assert set(edges(network.a)) == set(network.edges)
        edge_sets += 1

    # Exact integer estimates, without floating-point asymptotic tests.
    parameter_cases = 0
    for n in range(2, 10000):
        source_s = 0
        while n > 4 * (1 << source_s) ** 4:
            source_s += 1
        a = 2 << source_s
        f = 8 * (1 << source_s) ** 2 + 1
        assert a * a * f <= 144 * n
        for power in range(1, source_s + 2):
            groups = 1 << power
            length = a // groups
            boundary_size = (groups - 1) * (8 * length * length + 1) * f
            assert boundary_size * groups <= 9 * a * a * f
            parameter_cases += 1

    return {
        "status": "pass",
        "interpolation_cases": interpolation_cases,
        "boundary_addresses": boundary_addresses,
        "valid_queries": valid_queries,
        "malformed_queries": malformed_queries,
        "prefix_queries": prefix_queries,
        "nonzero_padding_rejections": padding_queries,
        "second_factor_tag_queries": weighted_second_factor,
        "maximum_arithmetic_dependency_depth": maximum_depth,
        "recomputed_interpolations_in_valid_queries": recomputed_interpolations,
        "source_edge_set_comparisons": edge_sets,
        "integer_parameter_bounds": parameter_cases,
        "scope": "Exact finite transcriptions; asymptotic space follows from the manuscript.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = verify()
    rendered = json.dumps(evidence, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
