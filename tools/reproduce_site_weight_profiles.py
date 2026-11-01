#!/usr/bin/env python3
"""Rebuild every reported site-weight profile and its compact audit."""
from __future__ import annotations

import argparse
from math import isqrt
from pathlib import Path
import subprocess
import sys


def cases():
    result = [(f"site-square-t{t}.json", [i * i for i in range(1, t + 1)])
              for t in (10, 12, 14, 16)]
    result += [(f"site-beatty-t{t}.json", [i + isqrt(2 * i * i) for i in range(1, t + 1)])
               for t in (12, 16)]
    result.append(("site-asymmetric-t12.json", [0, 1, 2, 4, 5, 7, 8, 12, 14, 18, 21, 23]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory / "site-weight-crossings"
    subprocess.run(["g++", "-O3", "-std=c++17", str(root / "site_weight_crossings.cpp"),
                    "-o", str(binary)], check=True)
    for name, weights in cases():
        print(f"Rebuilding {name}", flush=True)
        subprocess.run([str(binary), "profile", str(len(weights)), str(directory / name),
                        *map(str, weights)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    command = [sys.executable, str(root / "audit_site_weight_profiles.py"), "--directory", str(directory)]
    if args.output:
        command += ["--output", str(args.output.resolve())]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
