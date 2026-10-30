"""Exact counterexample to the literal low-mode product-span bound.

The original Appendix B operators are retained. Only the input label set
is restricted to partitions with parts <= L. This checks a finite span
obstruction, not an asymptotic tensor-rank or factoring lower bound.
"""
from hashlib import sha256
from itertools import permutations, product
from math import comb
from pathlib import Path
import json

try:
    from .check_half_piece_ranks import afam, opcoef, partitions
except ImportError:
    from check_half_piece_ranks import afam, opcoef, partitions


def determinant(matrix):
    """Leibniz determinant over integers; only a 4 by 4 matrix is used."""
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    answer = 0
    for order in permutations(range(size)):
        inversions = sum(order[i] > order[j]
                         for i in range(size) for j in range(i + 1, size))
        term = (-1) ** inversions
        for row, column in enumerate(order):
            term *= matrix[row][column]
        answer += term
    return answer


def check():
    weight = maximum_part = 2
    angle_index = 2 * weight
    shift = weight
    labels = tuple(sorted(part for part in partitions(weight)
                          if max(part, default=0) <= maximum_part))
    output_basis = tuple(sorted(partitions(2 * weight)))
    low_mode_basis = tuple(part for part in output_basis
                          if max(part, default=0) <= maximum_part)
    rows = []
    for operator_label, state_label in product(labels, repeat=2):
        # The source's a^<i> multiplies each weight-j component by C(h+j,i).
        coefficients = {
            target: sum(coefficient * comb(sum(monomial) + weight, angle_index)
                        * opcoef(monomial, target, state_label)
                        for monomial, coefficient in afam(operator_label)
                        if sum(monomial) + weight >= angle_index)
            for target in output_basis
        }
        rows.append((operator_label, state_label, coefficients))
    minor_columns = ((1, 1, 1, 1), (4,), (2, 1, 1), (2, 2))
    minor = [[coefficients[column] for column in minor_columns]
             for _, _, coefficients in rows]
    value = determinant(minor)
    assert len(labels) == 2 and len(low_mode_basis) == 3
    assert value != 0
    assert len(rows) > len(low_mode_basis)
    return {
        "source_repository": "anthropics/formal-math",
        "source_commit": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "claim": "The literal p_{<=L}(2h) intermediate-span bound fails.",
        "scope": "Exact finite counterexample; no asymptotic rank lower bound.",
        "operators": "Original unprojected Appendix B operators",
        "weight": weight,
        "maximum_label_part": maximum_part,
        "central_piece_index": 0,
        "angle_index": angle_index,
        "operator_shift": shift,
        "labels": [list(part) for part in labels],
        "rows": [
            {"operator_label": list(operator_label),
             "state_label": list(state_label),
             "output_terms": [
                 {"partition": list(part), "coefficient": coefficient}
                 for part, coefficient in coefficients.items() if coefficient]}
            for operator_label, state_label, coefficients in rows
        ],
        "low_mode_output_basis": [list(part) for part in low_mode_basis],
        "low_mode_output_dimension": len(low_mode_basis),
        "minor_columns": [list(part) for part in minor_columns],
        "minor": minor,
        "minor_determinant": value,
        "intermediate_span_dimension": len(rows),
        "program_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "shared_operator_checker_sha256": sha256(
            Path(__file__).with_name("check_half_piece_ranks.py").read_bytes()
        ).hexdigest(),
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2, sort_keys=True))
