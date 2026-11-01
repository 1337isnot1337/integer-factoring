#!/usr/bin/env python3
"""Rebuild five structural eight-feature screens and audit all exact ranks."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
from audit_eight_feature_screen import SCREENS, candidates, audit
from make_eight_feature_candidates import matrices


def complete(path, parameters, boundary, shift):
    if not path.exists():
        return False
    try:
        with path.open() as stream:
            for expected in parameters:
                value = json.loads(next(stream))
                if value["matrices"] != matrices(expected) or value["t"] != 4 or value["boundary"] != boundary or value.get("off_shift", 0) != shift:
                    return False
            return next(stream, None) is None
    except (ValueError, KeyError, StopIteration):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory / "eight_feature_crossings"
    subprocess.run(["g++", "-O3", "-std=c++17", str(root / "eight_feature_crossings.cpp"), "-o", str(binary)], check=True)
    pools = {}
    for pool in ("broad", "sparse", "efficient"):
        parameters = candidates(pool)
        path = directory / f"{pool}-candidates.txt"
        path.write_text("".join(" ".join(map(str, values)) + "\n" for values in parameters))
        pools[pool] = parameters, path
    for filename, pool, boundary, shift in SCREENS:
        parameters, candidate_file = pools[pool]
        output = directory / filename
        if not (args.resume and complete(output, parameters, boundary, shift)):
            command = [str(binary), "screen", "4", str(output), str(boundary), str(candidate_file), str(shift)]
            with (directory / (filename + ".log")).open("w") as log:
                subprocess.run(command, check=True, stdout=log, stderr=log)
        assert complete(output, parameters, boundary, shift), filename
        print(json.dumps({"completed_screen": filename, "verified_records": len(parameters)}), flush=True)
    result = audit(directory)
    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
