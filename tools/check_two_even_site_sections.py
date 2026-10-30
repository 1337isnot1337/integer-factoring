#!/usr/bin/env python3
"""Replay the graph-section and five-site two-even-space certificates."""

import argparse
import hashlib
import json
from pathlib import Path

from analyze_pfaffian_binary_screen import analyze, unpack_binary
from check_double_endpoint_retracts import bit_indices
from check_modular_crossing import rank_two
from check_pfaffian_span_f4 import verify
from explore_symplectic_boundary_codes import coefficient_basis
from explore_two_even_site_sections import parent_coefficients, restrict, save_basis, two_even_code


def literal_controls(parent):
    checks = 0
    for choices in ((2, 2, 2, 2, 2, 2), (1, 2, 1, 2, 1, 2), (1, 0, 2, 0, 1, 0)):
        data, rows, columns = restrict(parent, choices)
        selected = data["selected_original_labels"]
        original = parent["code"]
        code = dict(original)
        code["labels"] = len(selected)
        for key in ("input_vectors", "syndromes"):
            code[key] = [original[key][i] for i in selected]
        outputs = []
        for i in selected:
            support = bit_indices(original["input_vectors"][i])
            support.sort(key=lambda index: next(j for j, x in enumerate(
                original["ambient_occupations"][index]) if x == 3))
            outputs.append(support[choices[parent["edge_of"][i]] - 1])
        code["output_indices"] = outputs
        for key in ("input_Gram_rows", "induced_mode_complement_rows"):
            code[key] = [sum(((original[key][i] >> j) & 1) << k
                            for k, j in enumerate(selected)) for i in selected]
        literal, literal_rows, literal_columns = coefficient_basis(code)
        assert data["block_sizes"] == literal["block_sizes"]
        coordinates = {tuple(profile[:2]): mask for profile, mask in zip(
            literal["universal_coefficient_profile"], literal["universal_coefficient_coordinates"])}
        for i, parent_index in enumerate(data["basis_parent_indices"]):
            monomial = tuple(parent["ambient_coefficient_monomials"][parent_index])
            mask = coordinates[monomial]
            for actual, source in ((rows, literal_rows), (columns, literal_columns)):
                for block, size in enumerate(data["block_sizes"]):
                    expected = [0] * size
                    for j in bit_indices(mask):
                        expected = [a ^ b for a, b in zip(expected, source[j][block])]
                    assert expected == actual[i][block]
                    checks += size
    return checks


def audit_graphs(manifest_path, screens_path, work_dir):
    manifest = json.loads(manifest_path.read_text())
    archive = json.loads(screens_path.read_text())
    assert manifest["complete"] and manifest["all_excluded_over_F2"]
    assert len(manifest["cases"]) == len(archive) == 728
    assert len({tuple(x["section_choices"]) for x in manifest["cases"]}) == 728
    for name, digest in manifest["source_sha256"].items():
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest
    parent = parent_coefficients()
    literal_entries = literal_controls(parent)
    tested = eligible = 0
    for record in manifest["cases"]:
        data, rows, columns = restrict(parent, record["section_choices"])
        basis_path = work_dir / "graph-basis.json"
        screen_path = work_dir / "graph-screen.jsonl"
        save_basis(basis_path, data, rows, columns)
        screen_path.write_text(archive[record["name"]])
        replay = analyze(basis_path, screen_path)
        assert all(record[key] == value for key, value in replay.items())
        events = [json.loads(line) for line in archive[record["name"]].splitlines()]
        affordable = {unpack_binary(event["code"], len(rows)): event
                      for event in events if event["event"] == "eligible"}
        gate_squared = record["integer_cost_gate"] ** 2
        work = [[[0] * size for size in data["block_sizes"]] for _ in range(2)]
        previous = 0
        for step in range(1, 1 << len(rows)):
            mask = step ^ (step >> 1)
            changed = (mask ^ previous).bit_length() - 1
            for current, cut in zip(work, (rows, columns)):
                for target, source in zip(current, cut[changed]):
                    for i, value in enumerate(source):
                        target[i] ^= value
            ranks = [sum(rank_two(block) for block in current) for current in work]
            assert (mask in affordable) == (ranks[0] * ranks[1] < gate_squared)
            if mask in affordable:
                assert ranks == [affordable[mask]["row_rank"], affordable[mask]["column_rank"]]
                eligible += 1
            tested += 1
            previous = mask
    assert tested == 42536 and eligible == 3942
    return {"cases": 728, "nonzero_binary_vectors_independently_checked": tested,
            "affordable_rank_pairs_checked": eligible, "literal_control_matrix_rows_checked": literal_entries,
            "every_affordable_span_misses_endpoint": True,
            "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            "screens_sha256": hashlib.sha256(screens_path.read_bytes()).hexdigest()}


def audit_five_site(screen_dir, work_dir):
    code = two_even_code(5)
    data, rows, columns = coefficient_basis(code)
    # Preserve the original export order and encoding for the archived hashes.
    exported = dict(data, code=code, name="omega_two_even_P5",
                    basis_row_blocks=[[[hex(v) for v in b] for b in m] for m in rows],
                    basis_column_blocks=[[[hex(v) for v in b] for b in m] for m in columns])
    path = work_dir / "five-site-basis.json"
    path.write_text(json.dumps(exported, indent=2) + "\n")
    binary = analyze(path, screen_dir / "two-even-site-P5-F2.jsonl")
    assert binary["tested_vectors"] == 255 and binary["eligible_vectors"] == 6
    assert binary["eligible_span_dimension"] == 5 and binary["endpoint_remainder"]
    field4 = verify(path, screen_dir / "two-even-site-P5-F4.jsonl", 16)
    assert field4["projective_vectors_tested"] == 21845
    assert field4["eligible_count"] == 8 and field4["eligible_span_dimension"] == 5
    assert field4["endpoint_remainder"]
    return {"F2": binary, "F4": field4}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=("graphs", "five-site", "all"), default="all")
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    evidence = Path(__file__).resolve().parent.parent / "evidence"
    result = {"scope": args.scope, "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.scope in ("graphs", "all"):
        result["graphs"] = audit_graphs(evidence / "two-even-site-sections.json",
                                       evidence / "two-even-site-section-screens.json", args.work_dir)
        print(json.dumps(result["graphs"]), flush=True)
    if args.scope in ("five-site", "all"):
        result["five_site"] = audit_five_site(evidence, args.work_dir)
        print(json.dumps({key: value for key, value in result["five_site"]["F4"].items()
                          if key != "eligible_vectors"}), flush=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
