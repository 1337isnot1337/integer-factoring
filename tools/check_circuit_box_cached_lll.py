#!/usr/bin/env python3
"""Independent cache audits and complete saved-receipt equivalence checks."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from random import Random
from time import monotonic, process_time
from types import FunctionType

import circuit_box_cached_lll as cached
import circuit_box_centered_windows as windows
import circuit_box_ratio_small_roots as roots
import check_circuit_box_centered_windows as window_checker
import check_circuit_box_ratio_small_roots as root_checker


ROOT = Path(__file__).resolve().parents[1]
SOURCES = [
    "tools/circuit_box_cached_lll.py",
    "tools/check_circuit_box_cached_lll.py",
    "tools/circuit_box_centered_windows.py",
    "tools/check_circuit_box_centered_windows.py",
    "tools/circuit_box_ratio_small_roots.py",
    "tools/check_circuit_box_ratio_small_roots.py",
]
RECEIPTS = [
    "evidence/circuit-box-centered-windows.json",
    "evidence/circuit-box-ratio-small-roots.json",
]


def _digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _determinant(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next((i for i in range(column, len(work)) if work[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            result = -result
        value = work[column][column]
        result *= value
        for i in range(column + 1, len(work)):
            multiplier = work[i][column] / value
            for j in range(column, len(work)):
                work[i][j] -= multiplier * work[column][j]
    return result


def _clone(function, replacements):
    namespace = dict(function.__globals__)
    namespace.update(replacements)
    result = FunctionType(function.__code__, namespace, function.__name__,
                          function.__defaults__, function.__closure__)
    result.__kwdefaults__ = function.__kwdefaults__
    return result


def run_checks():
    hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    receipt_hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in RECEIPTS}
    frozen_functions = (roots._short_lattice_vector, roots.residue_short_polynomial,
                        windows.centered_window_factors,
                        window_checker.centered_window_factors,
                        root_checker.residue_short_polynomial)
    rng = Random(202611010850)
    counts = {"full_lll_certificates": 0, "exact_cache_state_recomputations": 0,
              "reference_early_stop_equalities": 0, "half_integer_ties": 0,
              "rectangular_lattices": 0, "invalid_inputs_rejected": 0,
              "complete_saved_suite_equalities": 0}
    traces = []
    cases = []
    for dimension in range(2, 9):
        for trial in range(8):
            width = dimension + int(trial % 4 == 0)
            rows = [[rng.randrange(-15, 16) if j < i or j >= dimension
                     else (rng.randrange(1, 20) if j == i else 0)
                     for j in range(width)] for i in range(dimension)]
            for _ in range(3 * dimension):
                i, j = rng.sample(range(dimension), 2)
                multiple = rng.choice((-2, -1, 1, 2))
                rows[i] = [a + multiple * b for a, b in zip(rows[i], rows[j])]
            cases.append(rows)
    for numerator in (1, -1, 5, -5):
        cases.append([[2, 0], [numerator, 2]])
        counts["half_integer_ties"] += 1
    for bits in (64, 257, 1024, 4096):
        integer = (1 << bits) + 3
        cases.append([[integer, 0, 0], [integer - 1, 3, 0], [integer + 5, -7, 11]])

    for rows in cases:
        untouched = [row[:] for row in rows]
        result = cached.reduce_lattice_cached(rows, audit_cache=True)
        assert rows == untouched
        assert abs(_determinant(result["transform"])) == 1
        dimension = len(rows)
        coefficients, norms = roots._gram_schmidt(result["basis"])
        assert all(abs(coefficients[i][j]) <= Fraction(1, 2)
                   for i in range(dimension) for j in range(i))
        assert all(norms[i] >= (Fraction(3, 4) - coefficients[i][i - 1] ** 2) * norms[i - 1]
                   for i in range(1, dimension))
        assert result["cache_audits"] == result["lattice_swaps"] + result["lattice_reductions"]
        threshold = sum(map(abs, result["basis"][0])) + 1
        accelerated = cached.short_lattice_vector_cached(rows, threshold)
        reference = roots._short_lattice_vector(rows, threshold)
        assert accelerated == reference
        counts["full_lll_certificates"] += 1
        counts["exact_cache_state_recomputations"] += result["cache_audits"]
        counts["reference_early_stop_equalities"] += 1
        counts["rectangular_lattices"] += len(rows[0]) != dimension
        traces.append({"dimension": dimension, "width": len(rows[0]),
                       "largest_entry_bits": max(abs(value).bit_length() for row in rows for value in row),
                       "swaps": result["lattice_swaps"],
                       "reductions": result["lattice_reductions"],
                       "cache_states_recomputed": result["cache_audits"],
                       "result_sha256": _digest(result)})

    invalid = [([], None), ([[]], None), ([[0]], None),
               ([[1, 2], [3]], None), ([[1.0]], None), ([[True]], None),
               ([[1, 2], [2, 4]], None), ([[1]], 0), ([[1]], 1)]
    for rows, threshold in invalid:
        try:
            cached.reduce_lattice_cached(rows, threshold)
        except (ValueError, AssertionError):
            counts["invalid_inputs_rejected"] += 1
        else:
            raise AssertionError("an invalid lattice contract was accepted")

    comparisons, timings = [], []
    for label, checker, replacements, receipt in (
        ("centered_windows", window_checker.run_checks,
         {"centered_window_factors": cached.centered_window_factors,
          "centered_factor_rectangle": cached.centered_factor_rectangle,
          "smallest_factor_by_centered_windows": cached.smallest_factor_by_centered_windows},
         RECEIPTS[0]),
        ("quarter_bit_endpoint", root_checker.run_checks,
         {"residue_short_polynomial": cached.residue_short_polynomial,
          "dyadic_residue_factors": cached.dyadic_residue_factors,
          "smallest_factor_by_small_roots": cached.smallest_factor_by_small_roots},
         RECEIPTS[1]),
    ):
        run = _clone(checker, replacements)
        started_wall, started_cpu = monotonic(), process_time()
        actual = run()
        timings.append({"suite": label, "wall_seconds": monotonic() - started_wall,
                        "cpu_seconds": process_time() - started_cpu})
        expected = json.loads((ROOT / receipt).read_text())
        assert actual == expected, label
        comparisons.append({"suite": label, "complete_dictionary_equality": True,
                            "result_sha256": _digest(actual), "counts": actual["counts"]})
        counts["complete_saved_suite_equalities"] += 1

    reference = json.loads((ROOT / RECEIPTS[0]).read_text())["local_controls"][0]["result"]
    gram_control = {"input_bits": 64, "dimension": reference["dimension"],
                    "lattice_swaps": reference["lattice_swaps"],
                    "lattice_reductions": reference["lattice_reductions"],
                    "reference_full_gram_evaluations": 1 + reference["lattice_swaps"] + reference["lattice_reductions"],
                    "cached_full_gram_evaluations": 1,
                    "complete_result_sha256": _digest(reference)}
    assert frozen_functions == (roots._short_lattice_vector, roots.residue_short_polynomial,
                                windows.centered_window_factors,
                                window_checker.centered_window_factors,
                                root_checker.residue_short_polynomial)
    assert hashes == {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    assert receipt_hashes == {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in RECEIPTS}
    return ({"status": "passed", "seed": 202611010850, "counts": counts,
             "full_reduction_traces": traces, "saved_suite_comparisons": comparisons,
             "large_gram_control": gram_control, "source_sha256": hashes,
             "reference_receipt_sha256": receipt_hashes,
             "scope": "Exact cache and full-output equivalence. Same classical small-root bound; no factoring-score improvement. Native timings are printed separately and excluded from this deterministic receipt."},
            timings)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result, timings = run_checks()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], **result["counts"],
                      "large_gram_control": result["large_gram_control"],
                      "native_suite_timings": timings}))


if __name__ == "__main__":
    main()
