#!/usr/bin/env python3
"""Exact checks for polynomial-time queries of grouped multiplication gates.

The production query uses composite-ring interpolation.  The independent
small-instance reference enumerates every local input digit, then propagates
the remaining gates.  Neither method materializes a full clique table.
"""

from itertools import combinations, product
import json
from math import gcd
from random import Random

from check_subexponential_piece_transfer import schedule_logs


def interpolate(points, values, modulus):
    """Lagrange interpolation, with an explicit unit check on every divisor."""
    degree = len(points)
    result = [0] * degree
    for selected in range(degree):
        polynomial = [1]
        denominator = 1
        for other in range(degree):
            if other == selected:
                continue
            following = [0] * (len(polynomial) + 1)
            for power, coefficient in enumerate(polynomial):
                following[power] -= points[other] * coefficient
                following[power + 1] += coefficient
            polynomial = [coefficient % modulus for coefficient in following]
            denominator = denominator * (points[selected] - points[other]) % modulus
        assert gcd(denominator, modulus) == 1
        scale = values[selected] * pow(denominator, -1, modulus) % modulus
        for power, coefficient in enumerate(polynomial):
            result[power] = (result[power] + scale * coefficient) % modulus
    return tuple(result)


class Network:
    def __init__(self, s, toy_base=None):
        self.h = 1 << s
        self.a = 2 * self.h
        self.d = 2 * self.h * self.h
        self.f = 4 * self.d + 1
        self.modulus = (1 << (4 * self.d)) + 1
        self.base = (1 << self.d) if toy_base is None else toy_base
        assert self.base > 1 and self.base & (self.base - 1) == 0
        self.digit_bits = self.base.bit_length() - 1
        self.input_bits = self.a * self.h * self.digit_bits
        self.carry_limit = 1 << self.f
        self.omega = 16 % self.modulus
        self.inverse = pow(self.omega, -1, self.modulus)
        self.xi = pow(self.omega, self.a, self.modulus)
        assert pow(self.xi, self.h, self.modulus) == self.modulus - 1
        self.gates = [(kind, index) for kind in ("p", "q", "mid", "ib")
                      for index in range(self.a)]
        self.gates += [("carry", index) for index in range(self.a ** 2)]
        self.order = {gate: index for index, gate in enumerate(self.gates)}
        self.edges = []
        for kind in ("p", "q"):
            self.edges += [((kind, row), ("mid", column))
                           for row in range(self.a) for column in range(self.a)]
        self.edges += [(("mid", row), ("ib", column))
                       for row in range(self.a) for column in range(self.a)]
        self.edges += [(("ib", index % self.a), ("carry", index))
                       for index in range(self.a ** 2)]
        self.edges += [(("carry", index), ("carry", index + 1))
                       for index in range(self.a ** 2 - 1)]
        self.edges.sort(key=lambda edge: (self.order[edge[0]], self.order[edge[1]]))

    def fa_outputs(self, row, digits):
        return tuple(pow(self.omega, row * column, self.modulus) * sum(
            value * pow(self.xi, position * column, self.modulus)
            for position, value in enumerate(digits)) % self.modulus
            for column in range(self.a))

    def mid_outputs(self, column, first, second):
        transformed = []
        for frequency in range(self.a):
            left = sum(value * pow(self.xi, row * frequency, self.modulus)
                       for row, value in enumerate(first)) % self.modulus
            right = sum(value * pow(self.xi, row * frequency, self.modulus)
                        for row, value in enumerate(second)) % self.modulus
            transformed.append(left * right % self.modulus)
        return tuple(pow(self.inverse, column * row, self.modulus) * sum(
            value * pow(self.inverse, self.a * frequency * row, self.modulus)
            for frequency, value in enumerate(transformed)) % self.modulus
            for row in range(self.a))

    def ib_outputs(self, incoming):
        normalization = pow(self.a ** 2, -1, self.modulus)
        return tuple(normalization * sum(
            value * pow(self.inverse, self.a * row * column, self.modulus)
            for row, value in enumerate(incoming)) % self.modulus
            for column in range(self.a))

    def outputs(self, gate, setting):
        kind, index = gate
        if kind in ("p", "q"):
            return self.fa_outputs(index, setting)
        if kind == "mid":
            return self.mid_outputs(index, *setting)
        if kind == "ib":
            return self.ib_outputs(setting)
        return ((setting[0] + setting[1]) // self.base,)

    def message(self, source, destination, settings, outgoing):
        kind, index = source
        other_kind, other_index = destination
        if kind in ("p", "q"):
            assert other_kind == "mid"
            return outgoing[source][other_index]
        if kind == "mid":
            if other_kind in ("p", "q"):
                return settings[source][other_kind == "q"][other_index]
            assert other_kind == "ib"
            return outgoing[source][other_index]
        if kind == "ib":
            if other_kind == "mid":
                return settings[source][other_index]
            assert other_kind == "carry" and other_index % self.a == index
            return outgoing[source][other_index // self.a]
        if other_kind == "ib":
            return settings[source][0]
        assert other_kind == "carry" and abs(other_index - index) == 1
        return outgoing[source][0] if other_index == index + 1 else settings[source][1]

    def global_setting(self, first, second):
        settings, outgoing = {}, {}
        for kind, digits in (("p", first), ("q", second)):
            for row in range(self.a):
                gate = (kind, row)
                settings[gate] = tuple(digits[row + self.a * i] for i in range(self.h))
                outgoing[gate] = self.outputs(gate, settings[gate])
        for column in range(self.a):
            gate = ("mid", column)
            settings[gate] = tuple(tuple(outgoing[(kind, row)][column]
                                        for row in range(self.a)) for kind in ("p", "q"))
            outgoing[gate] = self.outputs(gate, settings[gate])
        for row in range(self.a):
            gate = ("ib", row)
            settings[gate] = tuple(outgoing[("mid", column)][row]
                                   for column in range(self.a))
            outgoing[gate] = self.outputs(gate, settings[gate])
        carry = 0
        for index in range(self.a ** 2):
            gate = ("carry", index)
            coefficient = outgoing[("ib", index % self.a)][index // self.a]
            expected = sum(first[i] * second[index - i]
                           for i in range(len(first)) if 0 <= index - i < len(second))
            assert coefficient == expected
            settings[gate] = (coefficient, carry)
            outgoing[gate] = self.outputs(gate, settings[gate])
            carry = outgoing[gate][0]
        assert carry == 0
        first_value = sum(value * self.base ** i for i, value in enumerate(first))
        second_value = sum(value * self.base ** i for i, value in enumerate(second))
        return settings, outgoing, first_value * second_value, first_value

    def boundary(self, cluster, settings, outgoing):
        result = {}
        for first, second in self.edges:
            if (first in cluster) != (second in cluster):
                inside, outside = (first, second) if first in cluster else (second, first)
                result[(inside, outside)] = self.message(inside, outside, settings, outgoing)
        return result

    def grouping(self, count):
        assert 2 <= count <= self.a and self.a % count == 0
        block = self.a // count
        result = {}
        for gate in self.gates:
            kind, index = gate
            result[gate] = (index // self.a if kind == "carry" else index) // block
        return result

    def prefix_ok(self, row, digits, cutoff, prefix):
        if row == 0 and prefix >= 1 << max(0, self.input_bits - cutoff):
            return False
        for i, value in enumerate(digits):
            for bit in range(self.digit_bits):
                position = self.digit_bits * (row + self.a * i) + bit
                if position >= cutoff and ((value >> bit) ^ (prefix >> (position - cutoff))) & 1:
                    return False
        return True

    def propagate(self, cluster, boundary, first_settings, number, cutoff, prefix):
        """Forward propagation from fixed FA digits; no interpolation is used."""
        expected_keys = {(first, second) if first in cluster else (second, first)
                         for first, second in self.edges
                         if (first in cluster) != (second in cluster)}
        if boundary.keys() != expected_keys:
            return None
        if any(not 0 <= value < 1 << self.f for value in boundary.values()):
            return None
        settings, outgoing = {}, {}

        def read(destination, source):
            if source in cluster:
                return self.message(source, destination, settings, outgoing)
            return boundary[(destination, source)]

        for gate in self.gates:
            if gate not in cluster:
                continue
            kind, index = gate
            if kind in ("p", "q"):
                setting = first_settings[gate]
                if len(setting) != self.h or any(not 0 <= value < self.base for value in setting):
                    return None
                if kind == "p" and not self.prefix_ok(index, setting, cutoff, prefix):
                    return None
            elif kind == "mid":
                setting = tuple(tuple(read(gate, (input_kind, row))
                                      for row in range(self.a)) for input_kind in ("p", "q"))
                if any(not 0 <= value < self.modulus for vector in setting for value in vector):
                    return None
            elif kind == "ib":
                setting = tuple(read(gate, ("mid", row)) for row in range(self.a))
                if any(not 0 <= value < self.modulus for value in setting):
                    return None
            else:
                coefficient = read(gate, ("ib", index % self.a))
                carry = 0 if index == 0 else read(gate, ("carry", index - 1))
                if not 0 <= coefficient < self.modulus or not 0 <= carry < self.carry_limit:
                    return None
                total = coefficient + carry
                if total % self.base != number // (self.base ** index) % self.base:
                    return None
                if index == self.a ** 2 - 1 and total // self.base:
                    return None
                setting = (coefficient, carry)
            settings[gate] = setting
            outgoing[gate] = self.outputs(gate, setting)
            # This is only early rejection of explicitly enumerated digits.
            # It does not solve for or omit any candidate input assignment.
            for (inside, outside), value in boundary.items():
                if inside == gate and self.message(inside, outside, settings, outgoing) != value:
                    return None
        for (inside, outside), value in boundary.items():
            if self.message(inside, outside, settings, outgoing) != value:
                return None
        for first, second in self.edges:
            if first in cluster and second in cluster:
                assert self.message(first, second, settings, outgoing) == (
                    self.message(second, first, settings, outgoing))
        return settings

    def query(self, cluster, boundary, number, cutoff, prefix):
        assert sum(kind == "mid" for kind, _ in cluster) <= self.h
        first_settings = {}
        for gate in self.gates:
            kind, row = gate
            if gate not in cluster or kind not in ("p", "q"):
                continue
            selected = [column for column in range(self.a) if ("mid", column) not in cluster][:self.h]
            if any((gate, ("mid", column)) not in boundary for column in selected):
                return None
            values = [boundary[(gate, ("mid", column))] for column in selected]
            if any(not 0 <= value < self.modulus for value in values):
                return None
            points = [pow(self.xi, column, self.modulus) for column in selected]
            values = [value * pow(self.inverse, row * column, self.modulus) % self.modulus
                      for column, value in zip(selected, values)]
            first_settings[gate] = interpolate(points, values, self.modulus)
        return self.propagate(cluster, boundary, first_settings, number, cutoff, prefix)

    def enumerate_query(self, cluster, boundary, number, cutoff, prefix):
        inputs = [gate for gate in self.gates if gate in cluster and gate[0] in ("p", "q")]
        solutions = []
        for digits in product(range(self.base), repeat=len(inputs) * self.h):
            first_settings = {gate: tuple(digits[index * self.h:(index + 1) * self.h])
                              for index, gate in enumerate(inputs)}
            setting = self.propagate(cluster, boundary, first_settings, number, cutoff, prefix)
            if setting is not None:
                solutions.append(setting)
        return solutions


def verify():
    rng = Random(20261030)
    unit_cases = interpolation_cases = 0
    composite_interpolations = 0
    for s in range(5):
        network = Network(s)
        points = [pow(network.xi, column, network.modulus) for column in range(network.a)]
        for first, second in combinations(points, 2):
            assert gcd(first - second, network.modulus) == 1
            unit_cases += 1
        if network.a <= 8:
            selections = list(combinations(range(network.a), network.h))
        else:
            selections = [tuple(sorted(rng.sample(range(network.a), network.h))) for _ in range(8)]
        for trial in range(4):
            digits = tuple(rng.randrange(network.base) for _ in range(network.h))
            rows = range(network.a) if network.a <= 8 else (0, 1, network.a // 2, network.a - 1)
            for row in rows:
                outputs = network.fa_outputs(row, digits)
                for selected in selections:
                    recovered = interpolate([points[column] for column in selected],
                        [outputs[column] * pow(network.inverse, row * column, network.modulus)
                         % network.modulus for column in selected], network.modulus)
                    assert recovered == digits
                    interpolation_cases += 1
                    composite_interpolations += s == 1
        if s == 1:
            assert network.modulus % 641 == 0

    edge_pairs = volume_cases = 0
    for s in range(6):
        network = Network(s)
        for power in range(1, s + 2):
            count = 1 << power
            groups = network.grouping(count)
            block = network.a // count
            counts = {(a, b): 0 for a in range(count) for b in range(a + 1, count)}
            for first, second in network.edges:
                a, b = sorted((groups[first], groups[second]))
                if a != b:
                    counts[a, b] += 1
            for (a, b), value in counts.items():
                assert value == 8 * block ** 2 + (b == a + 1)
                edge_pairs += 1
            label = (8 * block ** 2 + 1) * network.f
            assert count ** 2 * label <= 9 * network.a ** 2 * network.f
            # Smallest bit length selecting this s (except the finite s=0 case).
            minimum_n = 2 if s == 0 else 4 * (1 << (s - 1)) ** 4 + 1
            assert count ** 2 * label <= 1296 * minimum_n
            volume_cases += 1

    cluster_cases = malformed_cases = prefix_cases = weighted_cases = 0
    arbitrary_cluster_cases = toy_cluster_cases = 0

    def compare(network, cluster, boundary, number, cutoff, prefix, toy=False):
        nonlocal cluster_cases, weighted_cases, toy_cluster_cases
        fast = network.query(cluster, boundary, number, cutoff, prefix)
        direct = network.enumerate_query(cluster, boundary, number, cutoff, prefix)
        assert len(direct) <= 1
        assert direct == ([] if fast is None else [fast])
        if fast is not None:
            fast_weight = sum(value * network.base ** (row + network.a * position)
                              for (kind, row), digits in fast.items() if kind == "p"
                              for position, value in enumerate(digits))
            direct_weight = sum(value * network.base ** (row + network.a * position)
                                for (kind, row), digits in direct[0].items() if kind == "p"
                                for position, value in enumerate(digits))
            assert fast_weight == direct_weight
            weighted_cases += 1
        cluster_cases += 1
        toy_cluster_cases += toy

    network = Network(0)
    groups = network.grouping(2)
    clusters = [{gate for gate in network.gates if groups[gate] == group} for group in range(2)]
    for first_value, second_value in product(range(network.base ** 2), repeat=2):
        first = tuple(first_value // network.base ** i % network.base for i in range(2))
        second = tuple(second_value // network.base ** i % network.base for i in range(2))
        settings, outgoing, number, value = network.global_setting(first, second)
        for cluster in clusters:
            boundary = network.boundary(cluster, settings, outgoing)
            for cutoff in (0, network.input_bits // 2, network.input_bits):
                compare(network, cluster, boundary, number, cutoff, value >> cutoff)
                prefix_cases += 1
            if (first_value + 3 * second_value) % 8 == 0:
                for key in boundary:
                    changed = dict(boundary)
                    changed[key] = (changed[key] + 1) % network.modulus
                    compare(network, cluster, changed, number, network.input_bits, 0)
                    malformed_cases += 1
                compare(network, cluster, boundary, number, 0, value ^ 1)
                prefix_cases += 1

    # Disconnected and empty clusters test the stronger arbitrary-cluster lemma.
    for _ in range(48):
        cluster = {gate for gate in network.gates if rng.randrange(2)}
        middles = [gate for gate in cluster if gate[0] == "mid"]
        if len(middles) > network.h:
            cluster.remove(sorted(middles)[0])
        first = tuple(rng.randrange(network.base) for _ in range(2))
        second = tuple(rng.randrange(network.base) for _ in range(2))
        settings, outgoing, number, value = network.global_setting(first, second)
        boundary = network.boundary(cluster, settings, outgoing)
        compare(network, cluster, boundary, number, network.input_bits, 0)
        arbitrary_cluster_cases += 1
        if boundary:
            changed = dict(boundary)
            changed[next(iter(changed))] = network.carry_limit
            compare(network, cluster, changed, number, network.input_bits, 0)
            malformed_cases += 1

    # Source s=1 modulus/root, but a deliberately smaller radix permits exhaustive
    # enumeration of all four FA digits in each cluster. It is labelled as a toy
    # radix test, not a source-parameter instance.
    toy = Network(1, toy_base=4)
    groups = toy.grouping(4)
    for _ in range(12):
        first = tuple(rng.randrange(toy.base) for _ in range(toy.a * toy.h))
        second = tuple(rng.randrange(toy.base) for _ in range(toy.a * toy.h))
        settings, outgoing, number, value = toy.global_setting(first, second)
        for group in range(4):
            cluster = {gate for gate in toy.gates if groups[gate] == group}
            boundary = toy.boundary(cluster, settings, outgoing)
            compare(toy, cluster, boundary, number, toy.input_bits // 2,
                    value >> (toy.input_bits // 2), toy=True)
            changed = dict(boundary)
            key = list(changed)[rng.randrange(len(changed))]
            changed[key] = (changed[key] + 1) % toy.modulus
            compare(toy, cluster, changed, number, toy.input_bits, 0, toy=True)
            malformed_cases += 1

    # A legal residue polynomial can have an illegal digit coefficient.
    one_gate = {("p", 0)}
    invalid_digits = (network.base,)
    invalid_outputs = network.fa_outputs(0, invalid_digits)
    invalid_boundary = {(("p", 0), ("mid", column)): invalid_outputs[column]
                        for column in range(network.a)}
    compare(network, one_gate, invalid_boundary, 0, network.input_bits, 0)
    assert network.query(one_gate, invalid_boundary, 0, network.input_bits, 0) is None

    # Negative control: an FA together with all middle gates has four extensions
    # when the second-factor boundary is zero, so the midpoint cap matters.
    settings, outgoing, _, _ = network.global_setting((0, 0), (0, 0))
    oversized = {("p", 0)} | {("mid", index) for index in range(network.a)}
    boundary = network.boundary(oversized, settings, outgoing)
    extensions = network.enumerate_query(oversized, boundary, 0, network.input_bits, 0)
    assert len(extensions) == network.base == 4

    # The larger-top cutoff inequalities use integers only.
    schedule_cases = 0
    for logarithm in schedule_logs(6)[1:]:
        p = 1 << logarithm
        n = (p ** 4) // 4
        for volume in (p ** 2 * (2 * p ** 2 + 1) * 9,
                       p ** 2 * ((n.bit_length() + 1) ** 2)):
            if volume ** 4 <= n ** 3:
                continue
            label = volume // (p ** 2)
            assert 4096 * label ** 4 > n
            for child_log in schedule_logs(6):
                if 2 * child_log > logarithm:
                    continue
                child_squared = 1 << (2 * child_log)
                assert child_squared <= p
                assert child_squared ** 2 * n <= 2 * volume ** 2
                schedule_cases += 1

    return {
        "scope": "exact transcription checks for composite interpolation, unique cluster extension, and explicit grouping; no improved logarithmic exponent",
        "unit_vandermonde_differences": unit_cases,
        "interpolation_recoveries": interpolation_cases,
        "certified_composite_modulus_recoveries": composite_interpolations,
        "exact_intergroup_edge_counts": edge_pairs,
        "information_volume_cases": volume_cases,
        "independent_cluster_enumerations": cluster_cases,
        "malformed_boundary_cases": malformed_cases,
        "prefix_cases": prefix_cases,
        "weighted_unique_extensions": weighted_cases,
        "arbitrary_cluster_cases": arbitrary_cluster_cases,
        "toy_radix_cluster_cases": toy_cluster_cases,
        "larger_top_integer_cutoff_checks": schedule_cases,
        "negative_control": {
            "cluster": "one first-factor input gate and every middle gate",
            "same_boundary_extension_count": len(extensions),
            "illegal_digit_residue_rejected": True,
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
