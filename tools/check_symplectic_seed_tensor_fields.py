#!/usr/bin/env python3
"""Independent weighted-block F2 replay and complete F4 tensor-power screen."""

import argparse
from collections import Counter
import hashlib
from itertools import combinations_with_replacement
import json
from math import factorial, prod
from pathlib import Path
import random
import subprocess

from check_pfaffian_span_f4 import INVERSE, remainder, scale_vector, span_basis


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def kronecker(first, second):
    width = len(second)
    return [sum(value << (width * j) for j in bits(row))
            for row in first for value in second]


def polynomial_tensor(first, second):
    size = len(first[0]) * len(second[0])
    result = [[0] * size for _ in range(len(first) + len(second) - 1)]
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            product = kronecker(left, right)
            result[i + j] = [a ^ b for a, b in zip(result[i + j], product)]
    return result


def tensor_blocks(types, power, maximum_block=0):
    result = []
    omitted = 0
    for indices in combinations_with_replacement(range(len(types)), power):
        counts = Counter(indices)
        multiplicity = factorial(power)
        for index, count in counts.items():
            multiplicity //= factorial(count)
        for index in indices:
            multiplicity *= types[index]["multiplicity"]
        size = prod(types[index]["size"] for index in indices)
        if maximum_block and size > maximum_block:
            omitted += size * multiplicity
            continue
        matrices = [[1]]
        for index in indices:
            matrices = polynomial_tensor(matrices, types[index]["generator_rows"])
        result.append({"size": len(matrices[0]), "multiplicity": multiplicity,
                       "generator_rows": matrices})
    result.sort(key=lambda block: (block["size"], -block["multiplicity"]))
    if omitted:
        # A zero replacement preserves the input layout's dimension while
        # explicitly discarding these positive direct-summand ranks.
        result.append({"size": 1, "multiplicity": omitted,
                       "generator_rows": [[0] for _ in range(3 * power + 1)]})
    return result, omitted


def integer_gate(labels):
    assert 99 ** 2 > 2 * 70 ** 2
    lower, upper = 1, labels * labels
    while upper - lower > 1:
        middle = (lower + upper) // 2
        if middle ** 70 > labels ** 99:
            upper = middle
        else:
            lower = middle
    return upper


def export_native(cuts, labels, path):
    dimension = len(cuts["row"][0]["generator_rows"])
    gate = integer_gate(labels)
    with path.open("w") as stream:
        stream.write(f"{labels} {dimension} {gate}\n")
        for cut in ("row", "column"):
            blocks = cuts[cut]
            assert sum(block["size"] * block["multiplicity"] for block in blocks) == labels * labels
            stream.write(f"{len(blocks)}\n")
            for block in blocks:
                size = block["size"]
                stream.write(f"{size} {block['multiplicity']}\n")
                for matrix in block["generator_rows"]:
                    for row in matrix:
                        stream.write(" ".join(f"{(row >> shift) & ((1 << 64) - 1):x}"
                                              for shift in range(0, size, 64)) + "\n")
    return dimension, gate


def independent_field_rank(low, high):
    """Scalar F4 elimination in bit planes, using low pivots."""
    pivots = {}
    for lo, hi in zip(low, high):
        while lo or hi:
            pivot = ((lo | hi) & -(lo | hi)).bit_length() - 1
            scalar = ((lo >> pivot) & 1) + 2 * ((hi >> pivot) & 1)
            if pivot not in pivots:
                if scalar == 2:
                    lo, hi = lo ^ hi, lo
                elif scalar == 3:
                    lo, hi = hi, lo ^ hi
                pivots[pivot] = lo, hi
                break
            a, b = pivots[pivot]
            if scalar == 1:
                lo, hi = lo ^ a, hi ^ b
            elif scalar == 2:
                lo, hi = lo ^ b, hi ^ a ^ b
            else:
                lo, hi = lo ^ a ^ b, hi ^ a
    return len(pivots)


def exact_ranks(cuts, code):
    result = []
    for name in ("row", "column"):
        total = 0
        for block in cuts[name]:
            lo, hi = [0] * block["size"], [0] * block["size"]
            for i, matrix in enumerate(block["generator_rows"]):
                scalar = (code >> (2 * i)) & 3
                if scalar & 1:
                    lo = [a ^ b for a, b in zip(lo, matrix)]
                if scalar & 2:
                    hi = [a ^ b for a, b in zip(hi, matrix)]
            total += block["multiplicity"] * independent_field_rank(lo, hi)
        result.append(total)
    return result


def run_native(native, input_path, output_path, arguments):
    with input_path.open() as source, output_path.open("w") as target:
        subprocess.run([str(native), *arguments], stdin=source, stdout=target, check=True)
    return [json.loads(line) for line in output_path.read_text().splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--types", type=Path, required=True)
    parser.add_argument("--power", type=int, default=2)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--compare-binary", type=Path)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--maximum-block", type=int, default=0,
                        help="omit larger direct summands to obtain rank lower bounds; zero keeps all")
    parser.add_argument("--field2-only", action="store_true")
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.power < 1 or args.samples < 0 or args.maximum_block < 0:
        raise ValueError("invalid power or sample count")
    types = json.loads(args.types.read_text())
    certificate_check = None
    if "types" in types:
        from check_symplectic_seed_tensor_modules import verify_certificate
        certificate_check = verify_certificate(types)
        types = types["types"]
    built = {cut: tensor_blocks(types[cut], args.power, args.maximum_block)
             for cut in ("row", "column")}
    cuts = {cut: value[0] for cut, value in built.items()}
    omitted = {cut: value[1] for cut, value in built.items()}
    labels = 26 ** args.power
    args.export_dir.mkdir(parents=True, exist_ok=True)
    input_path = args.export_dir / f"omega-z0-power-{args.power}.txt"
    dimension, gate = export_native(cuts, labels, input_path)
    binary_path = args.export_dir / f"omega-z0-power-{args.power}-F2.jsonl"
    binary = run_native(args.native, input_path, binary_path, ("--field2", "--all"))
    assert binary[-1]["projective_tested"] == (1 << dimension) - 1
    binary_records = [record for record in binary if record["event"] == "rank"]
    assert len(binary_records) == (1 << dimension) - 1
    if args.compare_binary:
        assert not any(omitted.values()), "a full-rank comparison cannot use omitted blocks"
        independent = json.loads(args.compare_binary.read_text())
        if "powers" in independent:
            independent = next(record for record in independent["powers"] if record["power"] == args.power)
        independent = independent["rank_tables"]
        for record in binary_records:
            mask = record["index"]
            assert [record["row_rank"], record["column_rank"]] == [
                independent["row"][mask], independent["column"][mask]]
    binary_eligible = [record for record in binary_records if record["eligible"]]
    binary_span = span_basis((record["code"] for record in binary_eligible), dimension)
    binary_residue = remainder(1, binary_span, dimension)
    rng = random.Random(529451)
    binary_selected = {record["index"] for record in binary_eligible}
    while len(binary_selected) < min(len(binary_records), len(binary_eligible) + args.samples):
        binary_selected.add(rng.randrange(1, 1 << dimension))
    for index in sorted(binary_selected):
        record = binary_records[index - 1]
        assert exact_ranks(cuts, record["code"]) == [record["row_rank"], record["column_rank"]]
    print(json.dumps({"event": "binary_replay_complete", "vectors": len(binary_records),
                      "eligible": len(binary_eligible), "span": len(binary_span),
                      "endpoint_remainder": binary_residue}), flush=True)
    result = {
        "power": args.power, "labels": labels, "coefficient_dimension": dimension,
        "rank_lower_bounds_only": bool(any(omitted.values())),
        "maximum_block_requested": args.maximum_block,
        "omitted_direct_summand_dimensions": omitted,
        "direct_sum_certificate_check": certificate_check,
        "weighted_block_counts": {cut: len(blocks) for cut, blocks in cuts.items()},
        "largest_blocks": {cut: max(block["size"] for block in blocks) for cut, blocks in cuts.items()},
        "complete_binary_rank_replay": len(binary_records),
        "binary_independent_table_compared": bool(args.compare_binary),
        "binary_comparison_sha256": (hashlib.sha256(args.compare_binary.read_bytes()).hexdigest()
                                     if args.compare_binary else None),
        "independently_checked_binary_vectors": len(binary_selected),
        "integer_cost_gate": gate, "F2_eligible_vectors": binary_eligible,
        "F2_eligible_span_dimension": len(binary_span),
        "F2_endpoint_remainder": binary_residue, "F2_excluded_by_span": bool(binary_residue),
    }
    if args.field2_only:
        events = None
    else:
        if dimension > 13:
            raise ValueError("use --field2-only for coefficient dimensions above 13")
        events = screen_field4(args, cuts, dimension, gate, input_path)
        result.update(events)
    assert exact_ranks(cuts, 1) == [labels * labels - omitted[cut] for cut in ("row", "column")]
    result.update({
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "screen_weighted_pfaffian_span.cpp",
                                       "screen_pfaffian_span_f4.cpp", "check_pfaffian_span_f4.py",
                                       "check_symplectic_seed_tensor_modules.py")},
        "types_sha256": hashlib.sha256(args.types.read_bytes()).hexdigest(),
        "native_input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "F2_screen_sha256": hashlib.sha256(binary_path.read_bytes()).hexdigest(),
    })
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("F2_eligible_vectors", "F4_eligible_vectors")}, indent=2))


def screen_field4(args, cuts, dimension, gate, input_path):
    screen_path = args.export_dir / f"omega-z0-power-{args.power}-F4.jsonl"
    events = run_native(args.native, input_path, screen_path, ())
    assert events[-1]["projective_tested"] == (4 ** dimension - 1) // 3
    eligible = [record for record in events if record["event"] == "eligible"]
    assert len(eligible) == events[-1]["eligible"]
    records = {record["code"]: record for record in eligible}
    assert len(records) == len(eligible)
    selected = set(records)
    rng = random.Random(529452)
    target_count = min((4 ** dimension - 1) // 3, len(selected) + args.samples)
    while len(selected) < target_count:
        code = rng.randrange(1, 4 ** dimension)
        pivot = (code.bit_length() - 1) // 2
        selected.add(scale_vector(code, INVERSE[(code >> (2 * pivot)) & 3], dimension))
    for code in sorted(selected):
        ranks = exact_ranks(cuts, code)
        assert (code in records) == (ranks[0] * ranks[1] < gate * gate)
        if code in records:
            assert ranks == [records[code]["row_rank"], records[code]["column_rank"]]
    basis = span_basis(records, dimension)
    residue = remainder(1, basis, dimension)
    return {
        "F4_projective_vectors_tested": events[-1]["projective_tested"],
        "integer_cost_gate": gate, "F4_eligible_vectors": eligible,
        "independently_checked_F4_vectors": len(selected),
        "F4_eligible_span_dimension": len(basis), "endpoint_remainder": residue,
        "F4_excluded_by_span": bool(residue),
        "F4_screen_sha256": hashlib.sha256(screen_path.read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    main()
