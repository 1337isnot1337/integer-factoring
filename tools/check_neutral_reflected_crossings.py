#!/usr/bin/env python3
"""Exact controls for reflected modular duals and neutral Pfaffian seeds."""

from itertools import combinations_with_replacement
from math import comb
from pathlib import Path
import json

from check_modular_crossing import (
    column_bits, determinant_operator_two, homogeneous_piece_two,
    modular_operator, rank_two, regular,
)
from check_partition_channels import partition_denominator
from neutral_reflected_crossing import (
    binomial_polynomial, crossing_entry, dense_polynomial_seed,
    kernel_entry, sector_labels, tensor_coefficients,
)


def reflected_duality():
    entries = 0
    records = []
    for prime, period, length in ((2, 8, 2), (2, 16, 2), (2, 16, 4),
                                  (2, 32, 4), (3, 3, 2), (3, 9, 2),
                                  (3, 9, 4), (5, 5, 2), (5, 5, 4),
                                  (7, 7, 2)):
        weight = period * length // 2
        allowed = tuple(k for k in range(1, period) if k % prime)
        labels = tuple(mu for mu in combinations_with_replacement(allowed, length)
                       if sum(mu) == weight and regular(mu, prime))
        for xi, mu in enumerate(labels):
            parameter = tuple(sorted(period - k for k in mu))
            scale = pow(partition_denominator(mu) % prime, -1, prime)
            for yi, nu in enumerate(labels):
                assert scale * modular_operator(prime, parameter, (), nu) % prime == (xi == yi)
                entries += 1
        records.append({"prime": prime, "period": period, "length": length,
                        "weight": weight, "labels": len(labels)})
    return entries, records


def contraction_checks():
    entries = 0
    for period in (4, 8, 16, 32):
        for k in range(1, period, 2):
            for v in range(1, period, 2):
                # Expansion at t=z^-1 of t^v(1+t)^(k-v)/(1+t^Q).
                numerator = (binomial_polynomial(k - v) << v
                             if v & ~k == 0 else 0)
                for j in range(1, 4 * period + 1):
                    direct = (j * comb(j + period - k - 1, period - k - 1)
                              * (comb(j - 1, v - 1) if j >= v else 0)) % 2
                    via_series = sum((numerator >> (j - multiple * period)) & 1
                                     for multiple in range(j // period + 1)) % 2
                    assert direct == via_series
                    entries += 1
    return entries


def verify():
    dual_entries, dual_records = reflected_duality()
    contraction_entries = contraction_checks()
    profile_records = []
    source_entries = endpoints = 0
    for period, length, weight in ((8, 2, 8), (16, 2, 16), (32, 2, 32),
                                   (16, 4, 32), (16, 4, 24)):
        labels = sector_labels(period, length, weight)
        dimension = len(labels)
        matrices = tensor_coefficients(period, labels)
        total = [0] * dimension ** 2
        for rows in matrices.values():
            total = [a ^ b for a, b in zip(total, rows)]
        assert total == [1 << i for i in range(dimension ** 2)]
        endpoints += dimension ** 4
        assert min(matrices) >= -period * length
        assert max(matrices) <= 2 * weight
        profile = [[shift, rank_two(rows), rank_two(column_bits(rows, dimension))]
                   for shift, rows in matrices.items()]
        profile_records.append({"period": period, "length": length,
                                "weight": weight, "labels": dimension,
                                "row_column_profile": profile})
        if 2 * weight == period * length:
            duals = tuple((tuple(sorted(period - k for k in mu)),) for mu in labels)
            # Full source-piece comparison at small weights. At length four,
            # test nonzero and vanishing coefficients without the large tails.
            shifts = (range(-2 * weight, 2 * weight + 1)
                      if period <= 16 and length == 2 else
                      (-2 * weight, -3 * weight // 2, -weight, -weight // 2, 0))
            for shift in shifts:
                actual = matrices.get(shift, [0] * dimension ** 2)
                expected = homogeneous_piece_two(weight, shift, labels, duals)
                assert actual == expected
                source_entries += dimension ** 4
        if (period, length, weight) == (16, 4, 32):
            assert profile == [[-64, 1, 64], [-32, 8, 30], [-16, 6, 6],
                               [0, 29, 29], [16, 6, 6], [32, 30, 8], [64, 64, 1]]
            values, normalization = dense_polynomial_seed(period, labels)
            saved = json.loads((Path(__file__).resolve().parents[1]
                                / "evidence/neutral-reflected-q16-l4-h32-seed.json").read_text())
            assert saved["values"] == [[value, 0] for value in values]
            assert saved["normalization"] == normalization

    # Equal length without fixed weight does not force a Pascal minor diagonal.
    # At z=1, singletons x=(3),v=(1) have bitwise inclusion and hence endpoint 1.
    endpoint = lambda entry: entry[1].bit_count() % 2
    assert endpoint(kernel_entry(8, 0, 3, 1, 1)) == 1
    assert endpoint(kernel_entry(8, 2, 1, 3, 1)) == 1
    assert (3,) != (1,)
    # The source reflected dual is not obtained by merely keeping the same part.
    mu = (1, 7, 11, 13)
    assert determinant_operator_two(mu, (), mu) == 0
    assert determinant_operator_two((3, 5, 9, 15), (), mu) == 1
    return {"status": "exact controls passed; no factoring improvement",
            "reflected_duality_entries": dual_entries,
            "reflected_duality_ranges": dual_records,
            "periodic_contraction_coefficients": contraction_entries,
            "literal_source_piece_entries": source_entries,
            "complete_endpoint_entries": endpoints,
            "complete_profiles": profile_records,
            "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5"}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
