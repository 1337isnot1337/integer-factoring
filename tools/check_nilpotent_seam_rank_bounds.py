#!/usr/bin/env python3
"""Exact counter, rational window certificates and complete finite cut checks."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from check_quantum_braid_path_profiles import cut_ranks
from nilpotent_seam_rank_bounds import central_weights, rank_envelope, vacuum_dimension
from quantum_braid_nilpotent_seam import coefficient_tensor

ROOT = Path(__file__).resolve().parents[1]


def admissible(word):
    # Independent literal definition: erase zeros, then inspect two-runs.
    compact = "".join(str(d) for d in word if d)
    return all(len(run) <= 3 for run in compact.split("1"))


def counter_controls():
    vacuum = [len(paths(2 * f)) for f in range(9)]
    assert [vacuum_dimension(f) for f in range(9)] == vacuum
    assert not admissible((2, 0, 2, 0, 2, 0, 2))
    assert admissible((2, 2, 2, 1, 2, 2, 2))
    words, accepted, histograms = 0, 0, []
    for length in range(9):
        literal = [0] * (2 * length + 1)
        for word in product(range(3), repeat=length):
            words += 1
            if admissible(word):
                literal[sum(word)] += vacuum[length - word.count(2)]
                accepted += 1
        assert central_weights(length) == literal
        histograms.append({"length": length, "weighted_degree_histogram": literal})
    for length in range(65):
        weights = central_weights(length)
        assert all(v >= 0 for v in weights)
        assert max(a for a, v in enumerate(weights) if v) == 2 * length - length // 4
        top = 2 * length - length // 4
        for gap in range(min(4, top) + 1):
            bound = (
                (gap + 1)
                * (length + 1) ** (5 * gap + 3)
                * vacuum_dimension(length // 4 + gap)
            )
            assert weights[top - gap] <= bound
    return {
        "literal_ternary_words": words,
        "admissible_words": accepted,
        "independent_path_dimensions": vacuum,
        "degree_bound_lengths": 65,
        "histograms": histograms,
    }


def rational_certificates():
    c, z, bound = Fraction(21, 8), 256, 24000
    w = [Fraction(11, 4) ** h for h in range(4)]
    image = [c * w[0] + c * z * sum(w)] + [
        c * w[h] + z * z * w[h - 1] for h in range(1, 4)
    ]
    slacks = [bound * x - y for x, y in zip(w, image)]
    assert all(v > 0 for v in slacks)
    assert sum(w) == Fraction(2055, 64)
    inequalities = {
        "middle_row": (24000**28 * 5**35, 256**48 * 8**35),
        "middle_column": (147 * 5**9, 64 * 7**8),
        "shoulder_column": (105 * 6**6 * 5**2, 8**9),
    }
    assert all(lhs < rhs for lhs, rhs in inequalities.values())
    # Rational bounds on phi, proved by squaring positive numbers.
    assert Fraction(11, 5) ** 2 < 5 < Fraction(9, 4) ** 2
    # Check the transfer majorant against exact counter evaluations too.
    for length in (0, 1, 2, 4, 8, 16, 32, 64):
        weights = central_weights(length)
        assert (
            sum(value * z**a for a, value in enumerate(weights))
            <= sum(w) * bound**length
        )
    return {
        "transfer_free_weight": str(c),
        "tilt": z,
        "transfer_upper_bound": bound,
        "positive_vector": [str(v) for v in w],
        "componentwise_slacks": [str(v) for v in slacks],
        "rate_inequalities": {
            name: {"lhs": lhs, "rhs": rhs} for name, (lhs, rhs) in inequalities.items()
        },
        "middle_window": "a>=ceil(12*n/7), b<=floor(n/8): ranks <= (L^(5/4+o(1)),L^(9/8+o(1)))",
        "shoulder_window": "a=2*n-floor(n/4), b<=floor(n/6): ranks <= (L^(1/2+o(1)),L^(4/3+o(1)))",
    }


def reflection_controls():
    random, field = Random(781204), Field(0b111)
    for _ in range(12):
        d = random.choice((2, 3, 4))
        tensor = np.array(
            [random.randrange(4) for _ in range(d**4)], dtype=np.uint8
        ).reshape(d, d, d, d)
        ranks = cut_ranks(tensor)
        assert ranks == [
            field.rank(tensor.reshape(d * d, d * d).tolist()),
            field.rank(tensor.transpose(0, 2, 1, 3).reshape(d * d, d * d).tolist()),
        ]
        assert cut_ranks(tensor.transpose(0, 2, 1, 3)) == ranks[::-1]
    return 12


def check():
    counter = counter_controls()
    certificates = rational_certificates()
    reflections = reflection_controls()
    profiles, coefficient_checks = [], 0
    for length in (4, 6, 8):
        tensor = coefficient_tensor(length)
        a_count, b_count, d = tensor.shape[:3]
        weights = central_weights(length)
        swap = np.zeros(tensor.shape[2:], dtype=np.uint8)
        for p, q in product(range(d), repeat=2):
            swap[p, q, q, p] = 1
        assert np.array_equal(
            np.bitwise_xor.reduce(
                tensor.reshape(a_count * b_count, d, d, d, d), axis=0
            ),
            swap,
        )
        assert np.array_equal(swap.transpose(0, 2, 1, 3), swap)
        row_ranks, column_ranks = [], []
        for a in range(a_count):
            rr, cc = [], []
            for b in range(b_count):
                row, column = cut_ranks(tensor[a, b])
                row_bound, column_bound = rank_envelope(length, a, b, weights)
                assert row <= min(d * d, row_bound)
                assert column <= min(d * d, column_bound)
                rr.append(row)
                cc.append(column)
                coefficient_checks += 1
            row_ranks.append(rr)
            column_ranks.append(cc)
        profiles.append(
            {
                "length": length,
                "labels": d,
                "central_degree_bound": a_count - 1,
                "boundary_degree_bound": b_count - 1,
                "central_weight_histogram": weights,
                "row_ranks_by_a_b": row_ranks,
                "column_ranks_by_a_b": column_ranks,
                "all_joint_bounds_pass": True,
            }
        )
    integer_examples = []
    for length in (32, 64, 128, 256):
        a, b = (12 * length + 6) // 7, length // 8
        row_bound, column_bound = rank_envelope(length, a, b)
        integer_examples.append(
            {
                "length": length,
                "central_degree": a,
                "boundary_degree": b,
                "row_upper_bound": row_bound,
                "column_upper_bound": column_bound,
                "scope": "Exact uncapped bounds with polynomial prefactors; not measured tensor ranks.",
            }
        )
    sources = (
        "research/nilpotent-seam-rank-envelope.md",
        "research/nilpotent-seam-crossings.md",
        "research/fusion-projector-seams.md",
        "tools/nilpotent_seam_rank_bounds.py",
        "tools/check_nilpotent_seam_rank_bounds.py",
        "tools/quantum_braid_nilpotent_seam.py",
        "tools/quantum_braid_path.py",
        "tools/check_quantum_braid_path_profiles.py",
        "tools/check_fixed_field_descent.py",
        "tools/search_auxiliary_extension_span.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_block_braid.py",
    )
    return {
        "status": "pass",
        "counter_controls": counter,
        "exact_rational_certificates": certificates,
        "independent_reflection_and_rank_controls": reflections,
        "complete_coefficient_rank_pairs_checked": coefficient_checks,
        "profiles": profiles,
        "large_exact_bound_examples": integer_examples,
        "scope": "All-degree uniform rank envelopes and certified candidate windows; no favorable complete spanning identity, rank-sum optimum or stronger factoring score.",
        "source_sha256": {
            f: sha256((ROOT / f).read_bytes()).hexdigest() for f in sources
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
