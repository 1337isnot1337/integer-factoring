#!/usr/bin/env python3
"""Verify the completed k=4 diagonal determinant-code profile and entry controls."""
import argparse
import json
import math
from hashlib import sha256
from pathlib import Path

from coherent_hot_compounds_f8 import partial_profiles
from determinant_inverse_codes import code
from probe_determinant_inverse_codes import entry_controls, lse, serializable

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence/determinant-inverse-k4-diagonal-profile.json"
NAMES = [
    "determinant_inverse_codes.py",
    "probe_determinant_inverse_codes.py",
    "coherent_hot_compounds_f8.py",
    "coherent_hot_compounds_f8.cpp",
    "coherent_f8_maps.py",
    "coherent_hot_compounds.py",
    "native_rectangular_smith_f4.cpp",
    "selected_gaussian_coefficients.cpp",
    "binary_extension_smith.py",
]


def without_seconds(value):
    if isinstance(value, dict):
        return {k: without_seconds(v) for k, v in value.items() if k != "seconds"}
    if isinstance(value, list):
        return [without_seconds(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="Reconstruct all 9,800 orders at cutoff 511.")
    args = parser.parse_args()
    data = json.loads(EVIDENCE.read_text())
    assert data["status"] == "complete" and data["k"] == 4 and data["kind"] == "diagonal"
    assert data["cutoff"] == 511 and data["q"] == 70 and data["modes"] == 16
    assert not data["factoring_exponent_improved"]
    assert set(data["source_sha256"]) == set(NAMES)
    for name, digest in data["source_sha256"].items():
        assert sha256((ROOT / "tools" / name).read_bytes()).hexdigest() == digest, name
    acquired = code(4, "diagonal")
    assert acquired["q"] == 70 and acquired["modes"] == 16
    profile = data["profile"]
    assert profile["heights"] == [x["height"] for x in acquired["labels"]]
    assert len(profile["heights"]) == 70
    zero = profile["zero_orders_t"]
    infinity = profile["infinity_orders_t"]
    assert len(zero) == len(infinity) == 70 * 70 == 4900
    assert zero == sorted(zero) and infinity == sorted(infinity)
    assert min(zero + infinity) >= 0 and max(zero + infinity) <= 511
    construction = data["construction"]
    assert construction["all_orders_resolved"]
    assert [c["tail_count"] for c in construction["cuts"]] == [0, 0]
    assert construction["cuts"][0]["retained_orders_t"] == zero
    assert construction["cuts"][1]["retained_orders_t"] == infinity
    fit = data["numerical_fit"]
    rho, scale, heights = fit["rho"], fit["scale"], profile["heights"]
    raw = sum(lse([-scale * v for v in orders]) for orders in (zero, infinity))
    z = lse([-scale * h / rho for h in heights])
    assert abs(raw - 4 * rho * z) < 1e-11
    assert abs(rho - 0.7394802516174483) < 1e-12
    assert rho > 1 / math.sqrt(2)
    assert fit["strictly_below_critical_numerically"] is False
    assert fit["global_minimum_certified"] is False
    assert all(abs(p - math.exp(-scale * h / rho - z)) < 1e-13
               for p, h in zip(fit["physical_probabilities"], heights))
    controls = entry_controls(acquired)
    assert serializable(controls) == data["independent_controls"]
    if args.full:
        recomputed = partial_profiles(
            acquired["embedding"], acquired["retraction"], acquired["modes"], 511
        )
        assert recomputed["all_orders_resolved"]
        assert [c["tail_count"] for c in recomputed["cuts"]] == [0, 0]
        assert serializable(without_seconds(recomputed)) == without_seconds(construction)
    print(json.dumps({
        "status": "PASS",
        "complete_orders": 9800,
        "cutoff": 511,
        "sampled_natural_rho": rho,
        "independent_polynomial_entries": controls["entries"],
        "nonbinary_entries": controls["nonbinary_entries"],
        "full_profile_recomputed_this_run": args.full,
        "factoring_exponent_improved": False,
    }))
    print("This is a complete finite diagnostic; no global scale optimum, all-scale floor or score improvement is claimed.")


if __name__ == "__main__":
    main()
