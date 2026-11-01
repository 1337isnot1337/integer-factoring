#!/usr/bin/env python3
"""Independent literal controls and retained F4 boundary catalogue readback."""

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from coherent_f4_boundary import (
    SOURCES, ROOT, boundary_maps, catalogue, compare, counting_record, dot,
    expand_counts, frobenius, multiply, orbit, profiles, transpose,
)


def independent_multiply(a, b):
    result = 0
    for i in range(2):
        if b >> i & 1:
            result ^= a << i
    if result & 4:
        result ^= 7
    return result


def rank4(matrix):
    work, rank = [list(row) for row in matrix], 0
    for column in range(4):
        pivot = next((i for i in range(rank, 4) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = (0, 1, 3, 2)[work[rank][column]]
        work[rank] = [independent_multiply(value, inverse) for value in work[rank]]
        for i in range(rank + 1, 4):
            coefficient = work[i][column]
            work[i] = [a ^ independent_multiply(coefficient, b)
                       for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


@lru_cache(None)
def spread(value):
    answer = 0
    while value:
        bit = value & -value
        answer |= 1 << (2 * (bit.bit_length() - 1))
        value ^= bit
    return answer


def literal_controls():
    from eight_letter_reflection_retraction import ambient_normalized_tensor
    from weighted_smith_oracle_f4 import WeightedSmithOracle
    ambient = ambient_normalized_tensor()
    def reverse(word):
        return sum(((word >> i) & 1) << (3 - i) for i in range(4))
    results = []
    for coefficient in (2, 3):
        embedding = [{word: 1} for word in (0, 1, 2, 4, 8)] + [{6: 1, 9: coefficient}]
        retraction = [{word: 1} for word in (0, 1, 2, 4, 8, 6)]
        native = profiles(embedding, retraction, 4)
        q = len(embedding)
        matrices = [[[0] * (q * q) for _ in range(q * q)] for _ in range(2)]
        for x, y, u, v in product(range(q), repeat=4):
            values = [0, 0]
            for (a, ca), (b, cb), (c, cc), (d, cd) in product(
                    retraction[x].items(), retraction[y].items(), embedding[u].items(), embedding[v].items()):
                scalar = independent_multiply(independent_multiply(ca, cb), independent_multiply(cc, cd))
                a, b, c, d = map(reverse, (a, b, c, d))
                for cut, mask in enumerate((a | b << 4 | c << 8 | d << 12,
                                            a | c << 4 | b << 8 | d << 12)):
                    # The even-bit spread has disjoint supports from its
                    # one-bit shift, so multiplying by 1,2,3 packs F4 scalars.
                    values[cut] ^= spread(ambient[mask]) * scalar
            matrices[0][q * x + y][q * u + v] = values[0]
            matrices[1][q * x + u][q * y + v] = values[1]
            for value in values:
                endpoint = 0
                while value:
                    endpoint ^= value & 3
                    value >>= 2
                assert endpoint == int(x == v and y == u)
        for cut, matrix in enumerate(matrices):
            actual, _ = WeightedSmithOracle(matrix).profile([0] * (q * q), [0] * (q * q))
            expected = native["infinity_orders_t" if cut else "zero_orders_t"]
            assert sorted(actual) == expected
        results.append({"coherent_coefficient": coefficient, "complete_literal_orders": 2 * q * q,
                        "literal_endpoint_entries": 2 * q ** 4, "profiles": counting_record(native)})
    return results


def catalogue_controls():
    representatives = catalogue()
    keys = set().union(*(orbit(key) for key in representatives))
    assert len(representatives) == 1416 and len(keys) == 5440
    fixed = [0, 0, 0]
    for key in keys:
        kernel, covector = key[:4], key[4:]
        assert next(x for x in kernel if x) == 1
        assert dot(kernel, covector) == 1
        assert frobenius(frobenius(key)) == transpose(transpose(key)) == key
        assert frobenius(transpose(key)) == transpose(frobenius(key))
        fixed[0] += frobenius(key) == key
        fixed[1] += transpose(key) == key
        fixed[2] += frobenius(transpose(key)) == key
        matrix = [[int(i == j) ^ independent_multiply(kernel[i], covector[j])
                   for j in range(4)] for i in range(4)]
        assert rank4(matrix) == 3
        for i, j in product(range(4), repeat=2):
            entry = 0
            for k in range(4):
                entry ^= independent_multiply(matrix[i][k], matrix[k][j])
            assert entry == matrix[i][j]
    assert fixed == [120, 64, 40]
    assert (5440 + sum(fixed)) // 4 == 1416
    return {"rank_three_idempotents": 5440, "orbits": 1416,
            "fixed_by_frobenius_transpose_composite": fixed}


def native_controls():
    from coherent_hot_compounds import profiles as binary_profiles
    from check_coherent_q32_pressure import displayed_maps, ratio_interval
    from check_large_interval_pressure import BinaryIntervals
    result, baseline = [], None
    for coherent in (False, True):
        binary_embedding, binary_retraction = displayed_maps(coherent)
        embedding = [{word: 1 for word in range(256) if value >> word & 1}
                     for value in binary_embedding]
        retraction = [{word: 1 for word in range(256) if value >> word & 1}
                      for value in binary_retraction]
        expected = binary_profiles(binary_embedding, binary_retraction, 8)
        actual = profiles(embedding, retraction, 8)
        for key in ("heights", "zero_orders_t", "infinity_orders_t"):
            assert actual[key] == expected[key]
        if coherent:
            baseline = actual
        result.append({"coherent": coherent, "complete_orders": 2048,
                       "profiles": counting_record(actual)})
    lower, upper = BinaryIntervals.fraction_bounds(ratio_interval(baseline, 1781, 2501, 5405, 5406))
    assert 0 < lower <= upper < 1
    key = (1, 1, 2, 2, 3, 1, 1, 3)
    record = None
    controls = []
    for equivalent in sorted(orbit(key)):
        actual = profiles(*boundary_maps(equivalent[:4], equivalent[4:]))
        if record is None:
            record = actual
        assert actual["zero_orders_t"] == record["zero_orders_t"]
        assert actual["infinity_orders_t"] == record["infinity_orders_t"]
        controls.append(equivalent)
    return {"binary_controls": result, "equivariance_keys": controls,
            "equivariance_complete_orders": 2048 * len(controls),
            "reference_exact_rate": "rho < 1781/2501", "reference_ratio_strictly_below_one": True}


def read_catalogue(path):
    saved = json.loads(Path(path).read_text())
    representatives = catalogue()
    expected_hash = sha256(json.dumps(representatives, separators=(",", ":")).encode()).hexdigest()
    assert saved["catalogue_sha256"] == expected_hash
    assert saved["source_sha256"] == {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                       for name in SOURCES}
    keys = [tuple(record["key"]) for record in saved["records"]]
    assert keys == representatives[:len(keys)]
    reference = {"heights": saved["reference"]["heights"],
                 "zero_orders_t": expand_counts(saved["reference"]["zero_counts"]),
                 "infinity_orders_t": expand_counts(saved["reference"]["infinity_counts"])}
    comparisons = {}
    for digest, record in saved["profiles"].items():
        encoded = {key: record[key] for key in ("heights", "zero_counts", "infinity_counts")}
        assert sha256(json.dumps(encoded, sort_keys=True, separators=(",", ":")).encode()).hexdigest() == digest
        actual = {"heights": record["heights"], "zero_orders_t": expand_counts(record["zero_counts"]),
                  "infinity_orders_t": expand_counts(record["infinity_counts"])}
        assert actual["heights"] == reference["heights"]
        assert len(actual["zero_orders_t"]) == len(actual["infinity_orders_t"]) == 1024
        certificate = compare(reference, actual)
        assert certificate == record["comparison_to_reference"]
        comparisons[certificate["comparison"]] = comparisons.get(certificate["comparison"], 0) + 1
    covered = set()
    for key, record in zip(keys, saved["records"]):
        assert record["profile"] in saved["profiles"]
        actual_orbit = orbit(key)
        assert not covered.intersection(actual_orbit)
        covered.update(actual_orbit)
        assert record["orbit_size"] == len(actual_orbit)
        assert record["deletion_checks"] == 512
        assert all(sum(block["size"] for block in side) == 1024 for side in record["blocks"])
        assert all(block["precision_attempts"][-1][1] == 0
                   for side in record["blocks"] for block in side)
    if saved["complete"]:
        assert len(keys) == 1416 and len(covered) == 5440
    return {"complete": saved["complete"], "checked_orbits": len(keys), "covered_projectors": len(covered),
            "distinct_profile_pairs": len(saved["profiles"]), "exact_comparisons": comparisons}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--catalogue-only", action="store_true")
    args = parser.parse_args()
    if args.catalogue_only and not args.catalogue:
        raise ValueError("catalogue-only readback needs a catalogue")
    result = {"status": "passed", "factoring_exponent_improved": False}
    if not args.catalogue_only:
        assert all(multiply(a, b) == independent_multiply(a, b) for a, b in product(range(4), repeat=2))
        result["projector_controls"] = catalogue_controls()
        result["literal_controls"] = literal_controls()
        result["native_controls"] = native_controls()
    if args.catalogue:
        result["catalogue_readback"] = read_catalogue(args.catalogue)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    summary = {key: value for key, value in result.items() if key not in ("literal_controls", "native_controls")}
    if "literal_controls" in result:
        summary["literal_complete_orders"] = sum(row["complete_literal_orders"] for row in result["literal_controls"])
        summary["literal_endpoint_entries"] = sum(row["literal_endpoint_entries"] for row in result["literal_controls"])
        summary["native_binary_complete_orders"] = 4096
        summary["native_equivariance_complete_orders"] = result["native_controls"]["equivariance_complete_orders"]
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
