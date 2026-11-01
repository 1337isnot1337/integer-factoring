#!/usr/bin/env python3
"""Exact controls for quadratic-pole jets and inverse residue-block sums."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random
from time import perf_counter

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at, smith_orders
from interval_unit_jets import QuadraticJets
from pole_cauchy_quadratic_selectors import affine_quadratic_strata, affine_quadratic_jet_product
from pole_cauchy_quadratic_selectors import inverse_quadratic_pole_contraction, inverse_residue_block_parity
from pole_cauchy_rational_poles import inverse_two_pole_contraction

ROOT = Path(__file__).resolve().parents[1]


def quad(pair):
    return Quad(*(Fraction(value) for value in pair))


def power(value, exponent):
    result = Quad(Fraction(1))
    for _ in range(exponent):
        result *= value
    return result


def controls():
    rng = random.Random(2026110111)
    strata_checks = affine_jet_checks = zero_product_rejections = 0
    for _ in range(700):
        first = tuple(rng.randrange(-40, 41) for _ in range(2))
        step = tuple(rng.randrange(-12, 13) for _ in range(2))
        length = rng.randrange(1, 15)
        factors = [quad(tuple(a + i * b for a, b in zip(first, step))) for i in range(length)]
        if not all(factors):
            try:
                affine_quadratic_strata(first, step, length)
            except ValueError:
                zero_product_rejections += 1
            else:
                raise AssertionError("A zero affine factor was accepted")
            continue
        strata = affine_quadratic_strata(first, step, length)
        actual_indices = []
        for offset, stride, count, valuation in strata:
            indices = [offset + stride * j for j in range(count)]
            actual_indices += indices
            assert all(factors[i].valuation() == valuation for i in indices)
        assert sorted(actual_indices) == list(range(length))
        assert sum(count * valuation for _, _, count, valuation in strata) == sum(x.valuation() for x in factors)
        strata_checks += 1
        shift = max(x.valuation() for x in factors) + 1
        first_derivative = tuple(rng.randrange(-4, 5) for _ in range(2))
        step_derivative = tuple(rng.randrange(-3, 4) for _ in range(2))
        ring = QuadraticJets(rng.randrange(3, 14), rng.randrange(4))
        expected = ring.constant(1)
        for i, factor in enumerate(factors):
            valuation = factor.valuation()
            term = list(ring.constant(int(factor.a) >> valuation, int(factor.b) >> valuation))
            if ring.size > 1:
                term[1] = tuple((a + i * b) * (1 << (shift - valuation)) % ring.modulus
                                for a, b in zip(first_derivative, step_derivative))
            expected = ring.multiply(expected, tuple(term))
        actual = affine_quadratic_jet_product(ring, first, step, first_derivative,
                                              step_derivative, strata, shift)
        assert actual == expected
        affine_jet_checks += 1
    selector_values = 0
    for bits in range(8):
        for residue in range(1 << bits):
            for offset in range(-4, 5):
                a = residue + offset
                denominator = quad((a - residue, 1 << bits))
                selector = power((1 << bits) / denominator, 3)
                assert selector.valuation() is None or selector.valuation() >= 0
                assert mod_pair(selector, 1) == (int((a - residue) % (1 << bits) == 0), 0)
                selector_values += 1
    contractions = block_checks = nonzero_blocks = integer_specializations = 0
    examples = []
    for index in range(60):
        size = rng.randrange(1, 5)
        left = tuple((a, 1) for a in sorted(rng.sample(list(range(1, 22, 2)), size)))
        right = ((rng.randrange(1, 8, 2), size),)
        aa, bb = expand(left), expand(right)
        n = rng.randrange(31, 200, 2)
        k = max(abs(n - a * b) for a in aa for b in bb).bit_length() + 1
        dense = inverse(matrix_at(n, aa, bb, k))
        c = (rng.randrange(-4, 6), -(1 << rng.randrange(5)))
        d = (rng.randrange(-4, 6), -(1 << rng.randrange(5)))
        r, s = 1 + index % 3, 1 + (index // 3) % 3
        expected = sum((dense[j][i] / (power(a - quad(c), r) * power(b - quad(d), s))
                        for i, a in enumerate(aa) for j, b in enumerate(bb)), Quad())
        actual = inverse_quadratic_pole_contraction(n, left, right, k, c, d, r, s, 5)
        assert actual["scale_exponent"] + actual["working_precision"] == 5
        assert tuple(actual["residue"]) == mod_pair(expected * (1 << (-actual["scale_exponent"])),
                                                    actual["working_precision"])
        contractions += 1
        row_bits, column_bits = index % 5, (index // 5) % 5
        row_residue, column_residue = rng.randrange(16), rng.randrange(16)
        block = inverse_residue_block_parity(n, left, right, k, row_residue, row_bits,
                                             column_residue, column_bits)
        expected_block = sum((dense[j][i] for i, a in enumerate(aa) for j, b in enumerate(bb)
                              if (a - row_residue) % (1 << row_bits) == 0
                              and (b - column_residue) % (1 << column_bits) == 0), Quad())
        expected_block *= 1 << block["inverse_normalization_exponent"]
        assert tuple(block["parity"]) == mod_pair(expected_block, 1)
        block_checks += 1
        nonzero_blocks += any(block["parity"])
        if index < 12:
            old = inverse_two_pole_contraction(n, left, right, k, -1, -3, 5)
            new = inverse_quadratic_pole_contraction(n, left, right, k, (-1, 0), (-3, 0), 1, 1, 5)
            assert old["residue"] == new["residue"]
            assert old["scale_exponent"] == new["scale_exponent"]
            integer_specializations += 1
            examples.append({"n": n, "size": size, "row_pole": c, "column_pole": d,
                             "row_order": r, "column_order": s,
                             "jet_shift_exponent": actual["jet_shift_exponent"],
                             "working_precision": actual["working_precision"],
                             "scale_exponent": actual["scale_exponent"],
                             "block_parity": block["parity"],
                             "inverse_normalization_exponent": block["inverse_normalization_exponent"]})
    rejected = 0
    invalid = (((3, 0), (-1, -1), 1, 1),
               ((1, 0), (101, 1 << 10), 1, 1),
               ((-1, -1), (-3, -1), 0, 1))
    for c, d, r, s in invalid:
        try:
            inverse_quadratic_pole_contraction(101, ((3, 2),), ((3, 2),), 10, c, d, r, s, 5)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("An invalid quadratic-pole contract was accepted")
    try:
        inverse_quadratic_pole_contraction(101, ((1, 2),), ((3, 2),), 10,
                                           (-1, 0), (101, 1 << 10), 1, 1, 5)
    except ValueError:
        rejected += 1
    else:
        raise AssertionError("A zero cross-product factor was accepted")
    normalizing_domain = ((1, 3),)
    normalizing_orders = smith_orders(matrix_at(15, (1, 3, 5), (1, 3, 5), 5))
    assert normalizing_orders == [0, 0, 4]
    normalizing_blocks = []
    for a, b, expected_parity in ((1, 1, [1, 1]), (3, 5, [0, 0]), (5, 3, [0, 0])):
        block = inverse_residue_block_parity(15, normalizing_domain, normalizing_domain, 5, a, 3, b, 3)
        assert block["inverse_normalization_exponent"] == 4
        assert block["parity"] == expected_parity
        normalizing_blocks.append({"a": a, "b": b, "is_factor_pair": a * b == 15,
                                   "inverse_block_parity": block["parity"]})
    return {"affine_valuation_partitions": strata_checks, "literal_affine_jet_products": affine_jet_checks,
            "zero_affine_product_rejections": zero_product_rejections,
            "quadratic_residue_selector_values": selector_values,
            "dense_repeated_quadratic_pole_contractions": contractions,
            "normalized_inverse_residue_blocks": block_checks, "nonzero_block_controls": nonzero_blocks,
            "integer_simple_pole_specializations": integer_specializations,
            "invalid_contract_rejections": rejected, "examples": examples,
            "original_unit_directions_counterexample": {"n": 15, "nodes": [1, 3, 5], "k": 5,
                "smith_orders": normalizing_orders, "inverse_normalization_exponent": 4,
                "singleton_blocks": normalizing_blocks}}


def large_control():
    n, size = (1 << 401) + 45, (1 << 200) + 13
    left, right = ((1, size),), ((3, size),)
    k = max(abs(n - 3), abs(n - (2 * size - 1) * (2 * size + 1))).bit_length() + 12
    started = perf_counter()
    block = inverse_residue_block_parity(n, left, right, k, 5, 7, 11, 9)
    elapsed = perf_counter() - started
    result = block["quadratic_pole_contraction"]
    higher = inverse_quadratic_pole_contraction(n, left, right, k,
                                                result["row_pole"], result["column_pole"], 3, 3, 5)
    assert result["scale_exponent"] == higher["scale_exponent"]
    assert result["residue"] == [coordinate % (1 << result["working_precision"])
                                  for coordinate in higher["residue"]]
    return {"n": n, "input_bits": n.bit_length(), "domain_size": size, "left": left,
            "right": right, "k": k, "residue_block": block,
            "higher_precision_result": higher, "lower_precision_seconds": elapsed,
            "enumerated_matrix_nodes": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--small-only", action="store_true")
    args = parser.parse_args()
    checked = controls()
    print(json.dumps(checked, indent=2), flush=True)
    large = None if args.small_only else large_control()
    report = {"claim": "Exact quadratic-pole jets give normalized inverse residue-block parity in polynomial bit time",
              "controls": checked, "large_control": large,
              "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "pole_cauchy_quadratic_selectors.py",
                                             "pole_cauchy_domain_precision.py", "modular_interval_jets.py")},
              "factoring_exponent_improved": False,
              "scope": "Block sums concern the normalized inverse, not original unit Smith directions or factor-pair counts."}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    if large is not None:
        result = large["residue_block"]["quadratic_pole_contraction"]
        print(json.dumps({"input_bits": large["input_bits"], "domain_size": large["domain_size"],
                          "block_parity": large["residue_block"]["parity"],
                          "working_precision": result["working_precision"],
                          "jet_shift_exponent": result["jet_shift_exponent"],
                          "affine_valuation_strata": result["affine_valuation_strata"],
                          "lower_precision_seconds": large["lower_precision_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
