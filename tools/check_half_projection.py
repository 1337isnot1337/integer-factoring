"""Exact counterexample to projecting Appendix B operators to low modes."""
from math import comb
import json

try:
    from .check_half_piece_ranks import opcoef, partitions
except ImportError:
    from check_half_piece_ranks import opcoef, partitions


def check():
    part = (1, 1)
    full, projected = [], []
    for shift in range(-6, 7):
        value = projected_value = 0
        for intermediate in partitions(2 - shift):
            term = (opcoef(part, part, intermediate) *
                    opcoef(part, intermediate, part))
            value += term
            if max(intermediate, default=0) <= 2:
                projected_value += term
        full.append(value)
        projected.append(projected_value)
    full_piece = [sum((-1) ** i * comb(4, i) * full[s-2+i+6]
                      for i in range(5)) for s in range(-4, 5)]
    projected_piece = [sum((-1) ** i * comb(4, i) * projected[s-2+i+6]
                           for i in range(5)) for s in range(-4, 5)]
    assert full_piece == [2, -8, 28, -56, 70, -56, 28, -8, 2]
    assert projected_piece == [-36, 80, -64, 8, 46, -56, 28, -8, 2]
    assert sum(full_piece) == 2
    assert sum(projected_piece) == 0
    return {
        "source_commit": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "status": "exact finite counterexample; not a general restriction theorem",
        "weight": 2,
        "maximum_mode": 2,
        "label_domain": [[1, 1], [2]],
        "x_y_u_v": [[1, 1]] * 4,
        "s": list(range(-4, 5)),
        "original_piece_values": full_piece,
        "projected_piece_values": projected_piece,
        "expected_crossing": 2,
        "original_sum": sum(full_piece),
        "projected_sum": sum(projected_piece),
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2, sort_keys=True))
