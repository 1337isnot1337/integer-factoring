#!/usr/bin/env python3
"""Screen exact merged Gram channel bounds across low-mode shell energies."""

import argparse
import hashlib
import json
from math import isqrt, log2
from pathlib import Path
import time

from folded_gram_aggregate import complete_gram_channels, shell_dimension


def one_case(modulus, length, energy):
    started = time.monotonic()
    labels = shell_dimension(modulus, length, energy)
    if labels < 2:
        return None
    channels = complete_gram_channels(modulus, length, energy)
    scale = 10 ** 9
    upper = 0
    profile = []
    for shift, row, column in channels:
        row, column = min(row, labels * labels), min(column, labels * labels)
        value = row * column * scale * scale
        root = isqrt(value)
        upper += root + (root * root != value)
        profile.append([shift, row, column])
    assert upper >= labels * scale
    rho = (log2(upper) - log2(scale)) / (2 * log2(labels))
    return {"modulus": modulus, "length": length, "energy": energy,
            "labels": labels, "row_energy_offset": 2 * energy - length * modulus,
            "rank_sum_upper_numerator": upper, "rank_sum_denominator": scale,
            "bound_rho": rho, "score_from_this_upper_bound": -log2(rho),
            "strict_improved_seed_by_integer_comparison": upper ** 985 < labels ** 1393 * scale ** 985,
            "channels_after_dimension_cap": profile, "seconds": time.monotonic() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="64:8,64:12,64:16,128:16,128:24,128:32")
    parser.add_argument("--energy-step", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in (Path(__file__).name, "folded_gram_aggregate.py", "folded_gram_ranks.py")}
    results = []
    for spec in args.cases.split(","):
        modulus, length = map(int, spec.split(":"))
        step = args.energy_step or max(2, modulus // 16)
        if step < 2 or step % 2:
            raise ValueError("energy step must be positive and even")
        minimum = length * length
        maximum = length * modulus // 2
        maximum -= (maximum - length) % 2
        energies = sorted(set(range(minimum, maximum + 1, step)) | {maximum})
        for energy in energies:
            record = one_case(modulus, length, energy)
            if record is None:
                continue
            results.append(record)
            report = {"status": "exact constructive channel upper bounds; a failed upper bound is not a lower bound on actual tensor ranks",
                      "cases": results, "source_sha256": hashes}
            temporary = args.output.with_suffix(args.output.suffix + ".tmp")
            temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            temporary.replace(args.output)
            print(json.dumps({key: value for key, value in record.items()
                              if key != "channels_after_dimension_cap"}), flush=True)


if __name__ == "__main__":
    main()
