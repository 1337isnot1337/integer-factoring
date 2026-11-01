#!/usr/bin/env python3
"""A finite varying-image rank-one space whose affordable span contains E."""

import argparse
from collections import Counter
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

from explore_cubic_complement_outputs import ComplementOutputs


def rank_binary(rows, low=False):
    pivots = {}
    for value in rows:
        while value:
            pivot = ((value & -value).bit_length() - 1) if low else value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def xor_matrices(*matrices):
    result = [0] * 49
    for matrix in matrices:
        result = [a ^ b for a, b in zip(result, matrix)]
    return result


def rank_pair(pair):
    result = [rank_binary(matrix) for matrix in pair]
    assert result == [rank_binary(matrix, low=True) for matrix in pair]
    return result


def add_basis(pivots, matrix):
    value = sum(row << (49 * i) for i, row in enumerate(matrix))
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            pivots[pivot] = value
            return True
        value ^= pivots[pivot]
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(7)
    map_pairs = sorted({(u, v) for i, j in combinations(range(7), 2)
                        for u, v in product((1 << i, 1 << j, (1 << i) | (1 << j)), repeat=2)})
    assert len(map_pairs) == 154
    map_index = {pair: i for i, pair in enumerate(map_pairs)}
    q3, all_raw, records = [], [], []
    q3_basis, full_basis, selected_full = {}, {}, []
    for index, (u, v) in enumerate(map_pairs):
        complement = [v if (u >> i) & 1 else 0 for i in range(7)]
        assert rank_binary(complement) == 1
        rows, columns = family.matrices(complement)
        assert rows[0] == columns[0] == [1 << i for i in range(49)]
        u_square = sum(((u >> a) & 1) * ((u >> b) & 1) << (7 * a + b)
                       for a in range(7) for b in range(7))
        v_square = sum(((v >> a) & 1) * ((v >> b) & 1) << (7 * a + b)
                       for a in range(7) for b in range(7))
        c_vector = sum(((u >> a) & 1) * ((v >> b) & 1) << (7 * a + b)
                       for a in range(7) for b in range(7))
        transpose_vector = sum(((v >> a) & 1) * ((u >> b) & 1) << (7 * a + b)
                               for a in range(7) for b in range(7))
        direct_rows = [transpose_vector if (c_vector >> i) & 1 else 0 for i in range(49)]
        direct_columns = [v_square if (u_square >> i) & 1 else 0 for i in range(49)]
        assert rows[4] == direct_rows and columns[4] == direct_columns
        pair = (rows[4], columns[4])
        assert rank_pair(pair) == [1, 1]
        q3.append(pair)
        all_raw.append((rows, columns))
        add_basis(q3_basis, rows[4])
        for coefficient, matrix in enumerate(rows):
            if add_basis(full_basis, matrix):
                selected_full.append([index, coefficient])
        records.append({"image_mask": u, "covector_mask": v,
                        "complement_rows": complement, "q3_rank_pair": [1, 1]})

    identity = [1 << i for i in range(49)]
    diagonal = sum(1 << (8 * i) for i in range(7))
    rank_one_diagonal = [diagonal if i // 7 == i % 7 else 0 for i in range(49)]
    swap = [1 << (7 * (i % 7) + i // 7) for i in range(49)]
    endpoint = (identity, identity)
    h_pairing = (rank_one_diagonal, swap)
    n_pairing = (swap, rank_one_diagonal)
    assert rank_pair(h_pairing) == [1, 49]
    assert rank_pair(n_pairing) == [49, 1]

    edge_pieces, parity = [], Counter()
    edge_total = [[0] * 49, [0] * 49]
    for i, j in combinations(range(7), 2):
        vectors = (1 << i, 1 << j, (1 << i) | (1 << j))
        indices = [map_index[(u, v)] for u, v in product(vectors, repeat=2)]
        parity.update(indices)
        piece = tuple(xor_matrices(*(q3[index][cut] for index in indices)) for cut in range(2))
        alternating = (1 << (7 * i + j)) | (1 << (7 * j + i))
        assert piece[1] == [alternating if (alternating >> k) & 1 else 0 for k in range(49)]
        assert rank_pair(piece) == [4, 1]
        edge_total = [xor_matrices(old, new) for old, new in zip(edge_total, piece)]
        edge_pieces.append({"labels": [i, j], "q3_map_indices": indices,
                            "rank_pair": [4, 1], "cost": 2})
    assert edge_total == [xor_matrices(e, h) for e, h in zip(endpoint, h_pairing)]
    active_q3 = sorted(index for index, count in parity.items() if count % 2)
    assert len(active_q3) == 147
    h_from_space = tuple(xor_matrices(endpoint[cut], *(q3[index][cut] for index in active_q3))
                         for cut in range(2))
    assert h_from_space == h_pairing
    affordable_sum = tuple(xor_matrices(h_pairing[cut], edge_total[cut]) for cut in range(2))
    assert affordable_sum == endpoint

    # C=e_0*(e_0+e_1)^T has C*C^T=0, so A(C)+D(C)+H_C=N.
    n_index = map_index[(1, 3)]
    n_from_space = tuple(xor_matrices(cut[1], cut[3], cut[4]) for cut in all_raw[n_index])
    assert n_from_space == n_pairing
    balanced = tuple(xor_matrices(e, h, n) for e, h, n in zip(endpoint, h_pairing, n_pairing))
    assert rank_pair(balanced) == [22, 22]
    assert tuple(xor_matrices(h, n, b) for h, n, b in zip(h_pairing, n_pairing, balanced)) == endpoint
    assert 7 + 21 * 2 == 49 and 7 + 7 + 22 == 36
    assert 16 ** 70 > 7 ** 99 and 99 ** 2 > 2 * 70 ** 2
    assert 36 > 16 and 49 > 16
    result = {
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "explore_cubic_complement_outputs.py",
                                       "explore_cubic_multirow_outputs.py")},
        "field": "F2", "labels": 7, "rank_one_output_maps": records,
        "map_count": 154, "physical_raw_coefficient_matrices_regenerated": 770,
        "q3_formula_entries_checked_in_both_cuts": 154 * 2 * 2401,
        "q3_span_dimension": len(q3_basis), "full_raw_coefficient_span_dimension": len(full_basis),
        "full_basis_source_indices": selected_full,
        "affordable_edge_pieces": edge_pieces,
        "h_membership": {"crossing_coefficient": 1, "q3_map_indices": active_q3,
                         "rank_pair": [1, 49], "cost": 7},
        "affordable_decomposition": {"pieces": "H plus the 21 displayed edge pieces",
                                     "piece_count": 22, "cost": 49,
                                     "every_piece_cost_strictly_below_seed_threshold": True},
        "n_membership": {"map_index": n_index, "raw_coefficient_indices": [1, 3, 4],
                         "rank_pair": [49, 1], "cost": 7},
        "three_piece_decomposition": {"pieces": ["H", "N", "E+H+N"],
                                      "rank_pairs": [[1, 49], [49, 1], [22, 22]], "cost": 36},
        "affordable_piece_span_contains_crossing": True,
        "recorded_decompositions_pass_strict_total_cost_gate": False,
        "optimality_claim": False,
        "scope": "The explicit full mixed coefficient space of 154 binary rank-one complementary output maps; both image and covector vary. The affine-symmetric span escape is proved, but no strict seed or factoring exponent improvement is claimed.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("rank_one_output_maps", "full_basis_source_indices", "affordable_edge_pieces",
                                     "h_membership")}), flush=True)


if __name__ == "__main__":
    main()
