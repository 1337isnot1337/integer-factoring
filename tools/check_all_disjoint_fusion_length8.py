#!/usr/bin/env python3
"""Exact thirteen-label all-placement span and separating functional."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import tempfile

import numpy as np

from check_all_disjoint_fusion_placements import add_piece, crossing, digest, identity
from check_quantum_braid_fusion_seam import rank as field_rank
from quantum_braid_fusion_seam import FIELD
from quantum_braid_graded_seam import (
    _layout,
    all_disjoint_projector_sites,
    coefficient_tensor_sites,
    free_counts,
    projector_sites,
)

ROOT = Path(__file__).resolve().parents[1]
LENGTH = 8
EXPECTED = {
    "labels": 13,
    "placements": 277,
    "nonzero_pieces": 4086,
    "span_dimension": 4014,
    "dual_support": 3663,
}


def build_matrix(cache: Path | None):
    labels, _, _, _ = _layout(LENGTH)
    if cache is not None:
        data = cache.read_bytes()
        rows, columns = struct.unpack("<QQ", data[:16])
        matrix = np.frombuffer(data[16:], dtype=np.uint8).reshape(rows, columns)
        return matrix, labels
    pieces = []
    for placement in all_disjoint_projector_sites(LENGTH):
        tensor = coefficient_tensor_sites(LENGTH, placement)
        for degree in range(LENGTH + 1):
            add_piece(pieces, labels, tensor[degree], placement, degree, "row")
    add_piece(pieces, labels, identity(labels), None, None, "identity")
    matrix = np.array([piece[0] for piece in pieces], dtype=np.uint8)
    return matrix, labels


def prefix_matrix(labels):
    pieces = []
    for left in free_counts(LENGTH):
        for right in free_counts(LENGTH):
            placement = projector_sites(LENGTH, left, right)
            tensor = coefficient_tensor_sites(LENGTH, placement)
            for degree in range(LENGTH + 1):
                add_piece(pieces, labels, tensor[degree], placement, degree, "row")
    add_piece(pieces, labels, identity(labels), None, None, "identity")
    return np.array([piece[0] for piece in pieces], dtype=np.uint8)


def write_matrix(path: Path, matrix: np.ndarray) -> None:
    with path.open("wb") as handle:
        handle.write(struct.pack("<QQ", matrix.shape[0], matrix.shape[1]))
        handle.write(matrix.tobytes())


def compile_tool(destination: Path) -> None:
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    subprocess.run(
        compiler
        + [
            "-O3",
            "-std=c++17",
            "-Wall",
            "-Wextra",
            "-Werror",
            str(ROOT / "tools/rank_f16_matrix_dual.cpp"),
            "-o",
            str(destination),
        ],
        check=True,
    )


def verify_dual(matrix: np.ndarray, target: np.ndarray, dual: np.ndarray):
    nonzero = np.flatnonzero(dual)
    values = dual[nonzero]
    products = FIELD[values[None, :], matrix[:, nonzero]]
    dots = np.bitwise_xor.reduce(products, axis=1)
    if np.any(dots):
        raise RuntimeError("separating functional does not vanish on the span")
    endpoint_dot = 0
    for value, coordinate in zip(values.tolist(), nonzero.tolist()):
        endpoint_dot ^= int(FIELD[value, target[coordinate]])
    if endpoint_dot == 0:
        raise RuntimeError("separating functional vanishes on the crossing")
    return endpoint_dot


def native_dual(matrix: np.ndarray, target: np.ndarray, work: Path):
    matrix_path = work / "matrix.bin"
    target_path = work / "target.bin"
    dual_path = work / "dual.bin"
    tool_path = work / "rank_dual"
    write_matrix(matrix_path, matrix)
    target_path.write_bytes(target.tobytes())
    compile_tool(tool_path)
    process = subprocess.run(
        [str(tool_path), str(matrix_path), str(target_path), str(dual_path)],
        check=True,
        text=True,
        capture_output=True,
    )
    match = re.fullmatch(
        r"rank=(\d+) free_column=(\d+) target_dot=(\d+)\n?", process.stdout
    )
    if match is None:
        raise RuntimeError(f"unexpected native output: {process.stdout!r}")
    rank, free_column, target_dot = (int(value) for value in match.groups())
    dual = np.frombuffer(dual_path.read_bytes(), dtype=np.uint8)
    computed_dot = verify_dual(matrix, target, dual)
    if computed_dot != target_dot:
        raise RuntimeError("native and Python target dots disagree")
    return rank, free_column, target_dot, dual
    nonzero = np.flatnonzero(dual)
    values = dual[nonzero]
    products = FIELD[values[None, :], matrix[:, nonzero]]
    dots = np.bitwise_xor.reduce(products, axis=1)
    if np.any(dots):
        raise RuntimeError("separating functional does not vanish on the span")
    endpoint_dot = 0
    for value, coordinate in zip(values.tolist(), nonzero.tolist()):
        endpoint_dot ^= int(FIELD[value, target[coordinate]])


def check(output: Path, cache: Path | None, dual_cache: Path | None) -> dict:
    matrix, labels = build_matrix(cache)
    if matrix.shape[1] != labels**4:
        raise AssertionError("matrix width does not match labels")
    if matrix.shape[0] != EXPECTED["nonzero_pieces"]:
        raise AssertionError("unexpected piece count")
    endpoint = crossing(labels).reshape(-1)
    prefix = prefix_matrix(labels)
    prefix_dimension = field_rank(prefix)
    prefix_contains = field_rank(np.concatenate([prefix, endpoint.reshape(1, -1)])) == prefix_dimension
    if dual_cache is None:
        with tempfile.TemporaryDirectory(prefix="all-disjoint-n8-") as directory:
            work = Path(directory)
            rank, free_column, target_dot, dual = native_dual(
                matrix, endpoint, work
            )
        rank_source = "native elimination"
    else:
        dual = np.frombuffer(dual_cache.read_bytes(), dtype=np.uint8)
        target_dot = verify_dual(matrix, endpoint, dual)
        rank = EXPECTED["span_dimension"]
        free_column = 3910
        rank_source = "cached native elimination"
    nonzero = np.flatnonzero(dual)
    if rank != EXPECTED["span_dimension"]:
        raise AssertionError(f"unexpected span dimension {rank}")
    if len(nonzero) != EXPECTED["dual_support"]:
        raise AssertionError(f"unexpected dual support {len(nonzero)}")
    if target_dot == 0:
        raise AssertionError("crossing is not separated")
    axis0 = np.count_nonzero(dual.reshape(labels, -1), axis=1)
    support_axis0 = np.flatnonzero(axis0).tolist()
    record = {
        "length": LENGTH,
        "labels": labels,
        "placements": EXPECTED["placements"],
        "prefix_pieces": int(prefix.shape[0]),
        "prefix_span_dimension": prefix_dimension,
        "prefix_contains_crossing": bool(prefix_contains),
        "nonzero_pieces": int(matrix.shape[0]),
        "span_dimension": rank,
        "crossing_in_span": False,
        "separating_functional": {
            "free_column": free_column,
            "target_dot": target_dot,
            "support": len(nonzero),
            "rank_source": rank_source,
            "first_coordinate_support": support_axis0,
            "nonzero": [
                [int(coordinate), int(value)]
                for coordinate, value in zip(nonzero.tolist(), dual[nonzero].tolist())
            ],
        },
    }
    receipt = {
        "field": "F16 with x^4+x+1",
        "status": "completed exact finite check",
        "record": record,
        "matrix_sha256": sha256(matrix.tobytes()).hexdigest(),
        "target_sha256": sha256(endpoint.tobytes()).hexdigest(),
        "dual_sha256": sha256(dual.tobytes()).hexdigest(),
        "source_sha256": {
            "tools/check_all_disjoint_fusion_length8.py": digest(
                ROOT / "tools/check_all_disjoint_fusion_length8.py"
            ),
            "tools/rank_f16_matrix_dual.cpp": digest(
                ROOT / "tools/rank_f16_matrix_dual.cpp"
            ),
            "tools/quantum_braid_graded_seam.py": digest(
                ROOT / "tools/quantum_braid_graded_seam.py"
            ),
        },
        "scope": (
            "All 277 disjoint three-strand f3 placements at thirteen labels, "
            "both cut orientations and the two identity pairings, over F16. "
            "The 4086-piece span has exact dimension 4014 and does not contain "
            "the crossing, while the 146-piece prefix span has dimension 146 "
            "and also misses it. A stored separating functional vanishes on "
            "every piece and evaluates to 4 on the crossing. No rank-sum optimum, "
            "uniform family, or stronger factoring score is claimed."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "evidence/all-disjoint-fusion-placements-length8.json",
    )
    parser.add_argument("--matrix-cache", type=Path)
    parser.add_argument("--dual-cache", type=Path)
    args = parser.parse_args()
    receipt = check(args.output, args.matrix_cache, args.dual_cache)
    record = receipt["record"]
    print(
        json.dumps(
            {
                "labels": record["labels"],
                "placements": record["placements"],
                "nonzero_pieces": record["nonzero_pieces"],
                "span_dimension": record["span_dimension"],
                "crossing_in_span": record["crossing_in_span"],
                "dual_support": record["separating_functional"]["support"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
