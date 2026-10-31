#!/usr/bin/env python3
"""Explore endpoint-vanishing exterior perturbations of a polynomial seed.

The first pass uses an explicit Gaussian gauge and several profile scales.
Only the best candidates receive a full-letter floating optimization.
Every reported number is exploratory, with no exact rate claim.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from rectangular_frozen_core_rows import build_frozen_row
from ramified_polynomial_crossings import endpoint_perturbation, occupation_charge
from native_rectangular_smith_f4 import RectangularOracle
from screen_full_alphabet_monomial_pressure import optimize


HEIGHTS = [4, 4, 8, -4, -4, -8, -1, 1, 0, -1, 1, 0,
           -5, -3, -8, 3, 5, 8, 0, 0, 0, 0, 0, 0]
BASE_SCALE = 0.9380048567639383


def candidates(family, coefficients, degrees):
    if family == "single":
        for degree in degrees:
            for positions in itertools.combinations(range(8), degree):
                mask = sum(1 << bit for bit in positions)
                if occupation_charge(mask, 2):
                    continue
                for coefficient in coefficients:
                    yield [(mask, coefficient, 0, 1)]
    elif family == "opposite-cubics":
        minus = [mask for mask in range(256) if mask.bit_count() == 3
                 and occupation_charge(mask, 2) == -1]
        plus = [mask for mask in range(256) if mask.bit_count() == 3
                and occupation_charge(mask, 2) == 1]
        for a, b, coefficient in itertools.product(minus, plus, coefficients):
            yield [(a, 1, 0, 1), (b, coefficient, 0, 1)]
    else:
        raise ValueError("Unknown candidate family.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=("single", "opposite-cubics"), default="single")
    parser.add_argument("--degrees", type=int, nargs="+", default=[4])
    parser.add_argument("--ramifications", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    parser.add_argument("--coefficients", type=int, nargs="+", choices=(1, 2, 3), default=[1, 2])
    parser.add_argument("--scale-factors", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    parser.add_argument("--cuts", type=int, nargs="+", choices=(0, 1), default=[1])
    parser.add_argument("--top", type=int, default=40)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--maxfun", type=int, default=200)
    parser.add_argument("--candidate-seconds", type=float, default=45)
    parser.add_argument("--seconds", type=float, default=3600)
    parser.add_argument("--seed", type=int, default=310733)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if any(d < 1 or d > 16 for d in args.ramifications):
        raise ValueError("This screen uses the 512-coefficient oracle; choose ramification 1 through 16.")
    if any(s <= 0 for s in args.scale_factors):
        raise ValueError("Profile scale factors must be positive.")
    source = build_frozen_row([1]+[0]*11)
    settings = list(itertools.product(args.ramifications,
                                     list(candidates(args.family, args.coefficients, args.degrees)),
                                     args.cuts))
    started = time.monotonic()
    deadline = started+args.seconds
    rng = np.random.default_rng(args.seed)
    result = {"status": "exploratory floating search; no exact rate certificate",
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "ramified_polynomial_crossings.py",
                                             "rectangular_frozen_core_rows.py")},
              "family": args.family, "random_seed": args.seed,
              "ramifications": args.ramifications, "scale_factors": args.scale_factors,
              "cuts": args.cuts, "gaussian_normalized_heights": HEIGHTS,
              "gaussian_optimal_scale": BASE_SCALE, "planned_candidates": len(settings),
              "screen": [], "optimized": [], "screen_completed": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        result["elapsed_seconds"] = time.monotonic()-started
        values = [record["rho"] for record in result["screen"] if "rho" in record]
        values += [record["optimization"]["best"]["rho"] for record in result["optimized"]
                   if record["optimization"]["best"]]
        result["best_rho"] = min(values, default=None)
        temporary = args.output.with_suffix(args.output.suffix+".tmp")
        temporary.write_text(json.dumps(result, indent=2)+"\n")
        temporary.replace(args.output)

    for d, monomials, cut in settings:
        if time.monotonic() >= deadline:
            break
        record = {"ramification": d, "monomials": monomials, "cut": cut}
        seed = endpoint_perturbation(source, monomials, ramification=d)
        record.update({key: seed[key] for key in ("coefficient_sha256", "degree")})
        oracle = RectangularOracle(seed["values"], 4, 4, endpoint_check=False)
        best = None
        errors = []
        for scale in args.scale_factors:
            ell = BASE_SCALE*scale
            theta = [height*ell for height in HEIGHTS]+[ell/d]
            try:
                rho, *_ = oracle.types(theta, cut)
            except ArithmeticError as error:
                errors.append(str(error))
                continue
            if best is None or rho < best[0]:
                best = rho, theta
        if best:
            record["rho"], record["parameters"] = best
        if errors:
            record["unresolved_scales"] = errors
        result["screen"].append(record)
        if len(result["screen"]) % 50 == 0:
            save()
            print(json.dumps({"phase": "screen", "candidates": len(result["screen"]),
                              "planned": len(settings), "best_rho": result["best_rho"]}), flush=True)
    result["screen_completed"] = len(result["screen"]) == len(settings)
    save()
    ordered = sorted((record for record in result["screen"] if "rho" in record), key=lambda r:r["rho"])
    for record in ordered[:args.top]:
        if time.monotonic() >= deadline:
            break
        seed = endpoint_perturbation(source, record["monomials"], ramification=record["ramification"])
        outcome = optimize(seed, cut=record["cut"], rng=rng, restarts=args.restarts,
                           maxfun=args.maxfun, height_bound=24*record["ramification"],
                           deadline=min(deadline, time.monotonic()+args.candidate_seconds),
                           starts=[record["parameters"]])
        result["optimized"].append({"ramification": record["ramification"],
                                    "monomials": record["monomials"], "cut": record["cut"],
                                    "optimization": outcome})
        save()
        print(json.dumps({"phase": "optimize", "candidate": len(result["optimized"]),
                          "ramification": record["ramification"],
                          "monomials": record["monomials"],
                          "rho": outcome["best"]["rho"] if outcome["best"] else None,
                          "best_rho": result["best_rho"]}), flush=True)
    save()


if __name__ == "__main__":
    main()
