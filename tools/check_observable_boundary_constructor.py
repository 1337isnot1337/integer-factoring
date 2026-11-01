#!/usr/bin/env python3
"""Exact relay grouping, unit interpolation, and local table decoding checks.

The full written table still enumerates its boundary labels.  These controls
compare the new per-entry decoder with an independent enumeration of every
local bounded input setting, including input length/prefix restrictions,
invalid residue labels, and all forced carry states.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
import random

from check_adaptive_transform_constructor import TransformModel, ceil_div
from check_operand_length_transform import LengthTransform, support_count


def balanced_blocks(a, p):
    if a < 2 * p:
        raise ValueError("every contiguous block must contain at least two indices")
    sizes = [a // p + int(group < a % p) for group in range(p)]
    blocks, start = [], 0
    for size in sizes:
        blocks.append(tuple(range(start, start + size)))
        start += size
    assert start == a and max(sizes) == ceil_div(a, p)
    return blocks


def relay_colors(a, p):
    blocks = balanced_blocks(a, p)
    colors = [0] * a
    starts = set()
    for group, block in enumerate(blocks):
        starts.add(block[0])
        for index in block:
            colors[index] = group
    carry = [row % p if column in starts else colors[column]
             for row in range(a) for column in range(a)]
    return blocks, colors, carry


def lagrange_coefficients(nodes, values, modulus):
    answer = [0] * len(nodes)
    for index, (node, value) in enumerate(zip(nodes, values)):
        polynomial, denominator = [1], 1
        for j, other in enumerate(nodes):
            if j == index:
                continue
            denominator = denominator * (node - other) % modulus
            following = [0] * (len(polynomial) + 1)
            for k, coefficient in enumerate(polynomial):
                following[k] = (following[k] - other * coefficient) % modulus
                following[k + 1] = (following[k + 1] + coefficient) % modulus
            polynomial = following
        assert gcd(denominator, modulus) == 1
        scalar = value * pow(denominator, -1, modulus) % modulus
        for k, coefficient in enumerate(polynomial):
            answer[k] = (answer[k] + scalar * coefficient) % modulus
    return tuple(answer)


def fa_message(model, residue, digits, recipient):
    return (model.powers[(residue * recipient) % model.length]
            * sum(value * model.powers[(model.a * i * recipient) % model.length]
                  for i, value in enumerate(digits))) % model.modulus


def interpolate_input(model, residue, height, recipients, prescribed):
    """Enumerate only the last h-min(h,external) bounded coefficients."""
    if any(not 0 <= value < model.modulus for value in prescribed):
        return []
    count = min(height, len(recipients))
    selected = recipients[:count]
    nodes = [model.powers[(model.a * recipient) % model.length] for recipient in selected]
    candidates = []
    for free in product(range(model.base), repeat=height - count):
        values = []
        for position, recipient in enumerate(selected):
            value = prescribed[position] * pow(model.powers[(residue * recipient) % model.length],
                                                -1, model.modulus) % model.modulus
            value -= sum(digit * model.powers[(model.a * (count + i) * recipient) % model.length]
                         for i, digit in enumerate(free))
            values.append(value % model.modulus)
        initial = lagrange_coefficients(nodes, values, model.modulus)
        candidate = initial + tuple(free)
        if all(value < model.base for value in initial) and all(
                fa_message(model, residue, candidate, recipient) == value
                for recipient, value in zip(recipients, prescribed)):
            candidates.append(candidate)
    assert len(candidates) == len(set(candidates))
    return candidates


class BoundaryNetwork:
    def __init__(self, model, groups, target, prefix_bits, prefix):
        self.model, self.groups, self.target = model, groups, target
        self.n = target.bit_length()
        self.prefix_bits, self.prefix = prefix_bits, prefix
        self.blocks, self.colors, self.carry_colors = relay_colors(model.a, groups)
        if isinstance(model, LengthTransform):
            assert self.n == model.n
            self.length_conditioned = True
            self.supports = (model.left_digits, model.right_digits)
            self.bit_caps = (model.ell, model.right_bits)
            self.prefix_capacity = model.n
        else:
            self.length_conditioned = False
            self.supports = (model.digits, model.digits)
            self.bit_caps = (model.capacity, model.capacity)
            self.prefix_capacity = model.capacity
        self.heights = {(family, r): support_count(self.supports[family], r, model.a)
                        for family in range(2) for r in range(model.a)}
        assert all(self.heights[0, r] + self.heights[1, r] <= model.a for r in range(model.a))
        self.edges = []
        for family, r, k in product(range(2), range(model.a), range(model.a)):
            self.edges.append((("fa", family, r, k), self.colors[r], self.colors[k], "residue"))
        for k, j in product(range(model.a), repeat=2):
            self.edges.append((("mid", k, j), self.colors[k], self.colors[j], "residue"))
        for j in range(model.length):
            self.edges.append((("ib", j), self.colors[j % model.a], self.carry_colors[j], "residue"))
        for j in range(model.length - 1):
            self.edges.append((("carry", j), self.carry_colors[j], self.carry_colors[j + 1], "carry"))

    def input_valid(self, family, residue, values):
        model = self.model
        if len(values) != self.heights[family, residue] or any(not 0 <= a < model.base for a in values):
            return False
        for i, value in enumerate(values):
            index = residue + model.a * i
            width = max(0, min(model.digit_bits, self.bit_caps[family] - index * model.digit_bits))
            if value >= 1 << width:
                return False
        if family == 0:
            if self.length_conditioned:
                high = self.bit_caps[0] - 1
                digit, bit = divmod(high, model.digit_bits)
                if digit % model.a == residue and not (values[digit // model.a] >> bit & 1):
                    return False
            if residue == 0 and not 0 <= self.prefix < 1 << max(0, self.n - self.prefix_bits):
                return False
            for position in range(self.prefix_bits, self.prefix_capacity):
                digit, bit = divmod(position, model.digit_bits)
                if digit % model.a != residue:
                    continue
                value = values[digit // model.a] if digit // model.a < len(values) else 0
                if value >> bit & 1 != self.prefix >> (position - self.prefix_bits) & 1:
                    return False
        return True

    def middle(self, k, left, right):
        model, a, modulus = self.model, self.model.a, self.model.modulus
        products = []
        for frequency in range(a):
            x = sum(left[r] * model.powers[(a * r * frequency) % model.length] for r in range(a))
            y = sum(right[r] * model.powers[(a * r * frequency) % model.length] for r in range(a))
            products.append(x * y % modulus)
        return tuple(model.inverse_powers[(k * j) % model.length] * sum(
            products[frequency] * model.inverse_powers[(a * frequency * j) % model.length]
            for frequency in range(a)) % modulus for j in range(a))

    def inverse_block(self, values):
        model = self.model
        return tuple(model.inverse_length * sum(
            values[k] * model.inverse_powers[(model.a * k * j2) % model.length]
            for k in range(model.a)) % model.modulus for j2 in range(model.a))

    def trace(self, left, right):
        model = self.model
        digits = [[(value >> (model.digit_bits * i)) & (model.base - 1)
                   for i in range(self.supports[family])]
                  for family, value in enumerate((left, right))]
        fa = {(family, r): tuple(fa_message(model, r, tuple(digits[family][r::model.a]), k)
                                 for k in range(model.a))
              for family in range(2) for r in range(model.a)}
        middle = {k: self.middle(k, [fa[0, r][k] for r in range(model.a)],
                                  [fa[1, r][k] for r in range(model.a)]) for k in range(model.a)}
        inverse = {j: self.inverse_block([middle[k][j] for k in range(model.a)]) for j in range(model.a)}
        coefficients = [inverse[j % model.a][j // model.a] for j in range(model.length)]
        value, carries, valid = model.carry(coefficients)
        assert valid and value == left * right
        labels = {}
        for edge, _, _, kind in self.edges:
            if edge[0] == "fa":
                _, family, r, k = edge
                labels[edge] = fa[family, r][k]
            elif edge[0] == "mid":
                labels[edge] = middle[edge[1]][edge[2]]
            elif edge[0] == "ib":
                labels[edge] = coefficients[edge[1]]
            else:
                labels[edge] = carries[edge[1] + 1]
        return labels

    def boundary(self, group, labels):
        return {edge: labels[edge] for edge, a, b, _ in self.edges if a != b and group in (a, b)}

    def boundary_valid(self, group, boundary, padding):
        if padding:
            return False
        for edge, a, b, kind in self.edges:
            if a != b and group in (a, b):
                value = boundary[edge]
                bound = self.model.modulus if kind == "residue" else self.model.carry_bound + 1
                if not 0 <= value < bound:
                    return False
        return True

    def accepts(self, group, boundary, input_values):
        model, inside = self.model, set(self.blocks[group])
        if any(not self.input_valid(family, r, values) for (family, r), values in input_values.items()):
            return False
        fa = {(family, r): tuple(fa_message(model, r, values, k) for k in range(model.a))
              for (family, r), values in input_values.items()}
        for (family, r), row in fa.items():
            if any(row[k] != boundary["fa", family, r, k] for k in range(model.a) if k not in inside):
                return False
        middle = {}
        for k in inside:
            operands = [[fa[family, r][k] if r in inside else boundary["fa", family, r, k]
                         for r in range(model.a)] for family in range(2)]
            middle[k] = self.middle(k, *operands)
            if any(middle[k][j] != boundary["mid", k, j] for j in range(model.a) if j not in inside):
                return False
        inverse = {}
        for j in inside:
            incoming = [middle[k][j] if k in inside else boundary["mid", k, j] for k in range(model.a)]
            inverse[j] = self.inverse_block(incoming)
            for j2 in range(model.a):
                index = j + model.a * j2
                if self.carry_colors[index] != group and inverse[j][j2] != boundary["ib", index]:
                    return False
        outgoing = {}
        for j in range(model.length):
            if self.carry_colors[j] != group:
                continue
            coefficient = inverse[j % model.a][j // model.a] if j % model.a in inside else boundary["ib", j]
            incoming = (0 if j == 0 else outgoing[j - 1] if self.carry_colors[j - 1] == group
                        else boundary["carry", j - 1])
            if not 0 <= incoming <= model.carry_bound:
                return False
            next_carry, digit = divmod(coefficient + incoming, model.base)
            if digit != self.target >> (model.digit_bits * j) & (model.base - 1):
                return False
            if next_carry > model.carry_bound:
                return False
            outgoing[j] = next_carry
            if j == model.length - 1:
                if next_carry != 0:
                    return False
            elif self.carry_colors[j + 1] != group and next_carry != boundary["carry", j]:
                return False
        return True

    def entry(self, group, boundary, decode, padding=False):
        if not self.boundary_valid(group, boundary, padding):
            return 0, 0
        inside = set(self.blocks[group])
        recipients = tuple(k for k in range(self.model.a) if k not in inside)
        keys, choices = [], []
        for family, r in product(range(2), self.blocks[group]):
            height = self.heights[family, r]
            if decode:
                prescribed = tuple(boundary["fa", family, r, k] for k in recipients)
                candidates = interpolate_input(self.model, r, height, recipients, prescribed)
            else:
                candidates = product(range(self.model.base), repeat=height)
            candidates = tuple(values for values in candidates if self.input_valid(family, r, values))
            keys.append((family, r))
            choices.append(candidates)
        answer, trials = 0, 0
        for values in product(*choices):
            trials += 1
            answer += self.accepts(group, boundary, dict(zip(keys, values)))
        return answer, trials


def check_relay_edges(counts, records):
    cases = [(a, p) for a in (4, 8, 16, 32, 64) for p in range(2, a // 2 + 1)]
    cases += [(a, p) for a in (128, 256) for p in (2, 3, 5, 7, 8, 13, 16, 31)]
    for a, p in cases:
        blocks, colors, carry = relay_colors(a, p)
        transform, ib, chain = Counter(), Counter(), Counter()
        for left, right in product(range(a), repeat=2):
            x, y = colors[left], colors[right]
            if x != y:
                transform[tuple(sorted((x, y)))] += 3
        for j in range(a * a):
            x, y = colors[j % a], carry[j]
            if x != y:
                ib[tuple(sorted((x, y)))] += 1
        for j in range(a * a - 1):
            x, y = carry[j], carry[j + 1]
            if x != y:
                chain[tuple(sorted((x, y)))] += 1
        k = ceil_div(a, p)
        assert max(transform.values(), default=0) <= 6 * k * k
        assert max(ib.values(), default=0) <= 2 * k
        assert max(chain.values(), default=0) <= 4 * k
        assert all(len(block) >= 2 for block in blocks)
        counts["complete_relay_colorings"] += 1
        counts["explicitly_counted_original_edges"] += 4 * a * a + a * a - 1
        if p in (2, 3, 5) and a <= 32:
            records.append({"A": a, "P": p, "maximum_block_size": k,
                            "max_transform_edges": max(transform.values(), default=0),
                            "max_IB_carry_edges": max(ib.values(), default=0),
                            "max_carry_chain_edges": max(chain.values(), default=0)})


def check_vandermonde_decoder(counts):
    rng = random.Random(3110202608)
    for model in (TransformModel(1, 2, 6), TransformModel(2, 1), TransformModel(2, 1, 24),
                  TransformModel(4, 1, 96)):
        roots = [model.powers[(model.a * k) % model.length] for k in range(model.a)]
        for i, j in combinations(range(model.a), 2):
            assert gcd(roots[i] - roots[j], model.modulus) == 1
            counts["composite_ring_unit_root_differences"] += 1
        recipient_sets = [tuple(x) for size in range(model.a + 1)
                          for x in combinations(range(model.a), size)] if model.a <= 4 else [
                              tuple(sorted(rng.sample(range(model.a), size))) for size in range(model.a + 1)]
        for recipients in recipient_sets:
            for height in range(model.a + 1):
                all_vectors = tuple(product(range(model.base), repeat=height))
                selected = all_vectors if len(all_vectors) <= 32 else rng.sample(all_vectors, 10)
                for residue in range(model.a):
                    for vector in selected:
                        prescribed = tuple(fa_message(model, residue, vector, k) for k in recipients)
                        expected = {other for other in all_vectors if all(
                            fa_message(model, residue, other, k) == value
                            for k, value in zip(recipients, prescribed))}
                        decoded = set(interpolate_input(model, residue, height, recipients, prescribed))
                        assert decoded == expected
                        assert len(decoded) <= model.base ** max(0, height - len(recipients))
                        counts["complete_bounded_input_fibers"] += 1


def check_table_entries(counts, records):
    rng = random.Random(3110202609)
    examples = [
        (TransformModel(2, 1), 3, 5),
        (TransformModel(2, 1), 7, 9),
        (TransformModel(2, 2), 17, 23),
        (LengthTransform(2, 9, 1), 1, 269),
        (LengthTransform(2, 10, 10), 617, 1),
        (LengthTransform(2, 19, 1), 1, 267013),
        (LengthTransform(2, 19, 19), 267013, 1),
        (LengthTransform(2, 20, 10), 719, 1451),
    ]
    for model, left, right in examples:
        target = left * right
        assert target.bit_length() == model.n if isinstance(model, LengthTransform) else True
        for cut in sorted({0, target.bit_length() // 2, target.bit_length()}):
            network = BoundaryNetwork(model, 2, target, cut, left >> cut)
            trace = network.trace(left, right)
            for group in range(2):
                boundary = network.boundary(group, trace)
                actual, decoded_trials = network.entry(group, boundary, True)
                expected, full_trials = network.entry(group, boundary, False)
                assert actual == expected == 1
                free = sum(max(0, network.heights[family, r] - (model.a - len(network.blocks[group])))
                           for family, r in product(range(2), network.blocks[group]))
                assert free <= len(network.blocks[group]) ** 2
                assert decoded_trials <= model.base ** free
                if not network.length_conditioned:
                    assert free == 0 and actual == 1 and decoded_trials == 1
                assert decoded_trials <= 1
                counts["valid_boundary_entry_comparisons"] += 1
                counts["reference_local_input_settings"] += full_trials
                counts["decoder_local_input_settings"] += decoded_trials
                counts["strictly_positive_free_digit_cases"] += free > 0
                for _ in range(2):
                    edge = rng.choice(tuple(boundary))
                    changed = dict(boundary)
                    bound = model.carry_bound + 1 if edge[0] == "carry" else model.modulus
                    changed[edge] = (changed[edge] + rng.randrange(1, bound)) % bound
                    answer, _ = network.entry(group, changed, True)
                    reference, trials = network.entry(group, changed, False)
                    assert answer == reference
                    counts["mutated_boundary_entry_comparisons"] += 1
                    counts["reference_local_input_settings"] += trials
                for kind in ("residue", "carry"):
                    choices = [edge for edge in boundary if (edge[0] == "carry") == (kind == "carry")]
                    if choices:
                        changed = dict(boundary)
                        changed[choices[0]] = model.carry_bound + 1 if kind == "carry" else model.modulus
                        assert network.entry(group, changed, True)[0] == 0
                        assert network.entry(group, changed, False)[0] == 0
                        counts["out_of_range_boundary_rejections"] += 1
                assert network.entry(group, boundary, True, padding=True)[0] == 0
                counts["nonzero_padding_rejections"] += 1
                if cut == target.bit_length():
                    records.append({"kind": "operand length" if network.length_conditioned else "first constructor",
                                    "N": target, "p": left, "q": right, "group": group,
                                    "digit_width": model.digit_bits, "free_digits": free,
                                    "entry_value": actual, "reference_settings": full_trials,
                                    "decoded_settings": decoded_trials})


def check_parameter_bounds(counts, records):
    for p in list(range(2, 130)) + [255, 256, 257, 511, 512, 1023, 1024, 4095, 4096]:
        h = 1 << ((p ** 9 - 1).bit_length())
        a, d = 2 * h, 2 * h * h
        t, k, n0 = d.bit_length() - 1, ceil_div(a, p), p ** 40
        denominator = a * a - a - 1
        preliminary0 = ceil_div(n0 + 1, denominator)
        candidates = {n0, n0 + 1, n0 + d - 1, n0 + d, 2 * n0, 8 * n0}
        candidates.update(denominator * (preliminary0 + shift) + delta
                          for shift in range(4) for delta in (-2, -1, 0, 1))
        for n in sorted(x for x in candidates if x >= n0):
            for even_round in (False, True):
                w = ceil_div(n + 1, denominator)
                if even_round:
                    w += w % 2
                m = d * ceil_div(2 * w + t, d)
                f, c = m + 1, w + t
                width = 6 * k * k * f + 2 * k * f + 4 * k * c
                free = k * k * w
                search_bits = free + 1
                bits = 1
                assert p * p * width < 17 * n
                assert 100 * p * p * search_bits < 101 * n
                assert p * ((p - 1) * width + free) < 18 * n
                assert 100000 * d * w < 50101 * n
                assert k <= a // 2
                counts["exact_large_parameter_cases"] += 1
                if n == n0 and p in (2, 3, 4, 8, 16, 128) and even_round:
                    records.append({"P": p, "n": str(n), "A": a, "D": d,
                                    "k": k, "w": str(w), "m": str(m),
                                    "lambda": str(width), "F": str(free), "B_tab": str(bits),
                                    "enumeration_bit_bound": str(search_bits),
                                    "top_volume_ratio": [str(p * p * width), str(n)],
                                    "scaled_entry_bits_ratio": [str(p * p * bits), str(n)]})


def check_contiguous_injectivity(counts, records):
    for model in (TransformModel(2, 1), TransformModel(2, 2),
                  TransformModel(2, 3), TransformModel(4, 1)):
        a, b, modulus = model.a, model.base, model.modulus
        assert modulus > model.digits * b * b
        for hidden_size in range(1, a // 2 + 1):
            for start in range(a - hidden_size + 1):
                hidden = set(range(start, start + hidden_size))
                recipients = tuple(i for i in range(a) if i not in hidden)
                for height in range(a):
                    for residue in range(a):
                        images = set()
                        for vector in product(range(b), repeat=height):
                            image = tuple(fa_message(model, residue, vector, k) for k in recipients)
                            assert image not in images
                            images.add(image)
                            counts["contiguous_boundary_distinct_images"] += 1
                        assert len(images) == b ** height
                        counts["complete_contiguous_injective_maps"] += 1
        zero = (0,) * a
        ones = (1,) * a
        first_hidden = set(range(a // 2))
        assert all(fa_message(model, 0, zero, k) == fa_message(model, 0, ones, k)
                   for k in range(a) if k not in first_hidden)
        even = tuple(int(i % 2 == 0) for i in range(a))
        odd = tuple(int(i % 2 == 1) for i in range(a))
        second_hidden = set(range(a // 2, a))
        assert all(fa_message(model, 0, even, k) == fa_message(model, 0, odd, k)
                   for k in range(a) if k not in second_hidden)
        counts["sharp_excluded_height_A_collisions"] += 2
        # An interleaved hidden set violates the interval hypothesis.
        left = (0,) * (a // 2) + (1,)
        right = (1,) + (0,) * (a // 2)
        assert len(left) <= a - 1
        assert all(fa_message(model, 0, left, k) == fa_message(model, 0, right, k)
                   for k in range(0, a, 2))
        counts["noncontiguous_hidden_set_collisions"] += 1
        records.append({"A": a, "digit_width": model.digit_bits,
                        "M": str(modulus), "largest_proved_height": a - 1})
    # The principal roots alone do not give bounded injectivity when the
    # modulus is below the coefficient-size threshold.
    invalid = TransformModel(2, 5, 8)
    assert not invalid.modulus > invalid.digits * invalid.base ** 2
    collision, zero = (16, 17, 1), (0, 0, 0)
    assert all(value < invalid.base for value in collision)
    assert all(fa_message(invalid, 0, collision, k) == fa_message(invalid, 0, zero, k)
               for k in (2, 3))
    counts["undersized_modulus_observability_collisions"] += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts, coloring, entries, parameters, injectivity = Counter(), [], [], [], []
    check_relay_edges(counts, coloring)
    check_vandermonde_decoder(counts)
    check_table_entries(counts, entries)
    check_contiguous_injectivity(counts, injectivity)
    check_parameter_bounds(counts, parameters)
    evidence = {"status": "all exact checks passed",
                "scope": "explicit group widths, Boolean grouped tables, and per-entry input decoding; all boundary labels remain written; factoring exponent unchanged",
                "counts": dict(sorted(counts.items())), "relay_colorings": coloring,
                "table_entry_examples": entries, "large_parameter_examples": parameters,
                "bounded_injectivity_models": injectivity}
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": evidence["status"], "counts": evidence["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
