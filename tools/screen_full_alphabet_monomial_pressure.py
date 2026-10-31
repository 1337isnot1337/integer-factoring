#!/usr/bin/env python3
"""Explore the exact critical-cubic seeds with floating weighted profiles.

The constructor proves the endpoint and generates every coefficient exactly.
The optimization is exploratory; a score claim needs an exact pressure replay.
Each completed candidate retains its seed hash, weights, types, and both cuts.
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

from check_full_alphabet_monomial_residual import cubic_pair, full_residual_monomial_seed
from native_rectangular_smith_f4 import RectangularOracle

PAIRS = {"A": (0, 2), "B": (0, 7), "C": (1, 3), "D": (1, 6)}


def entropy(probabilities):
    return -sum(p * math.log(p) for p in probabilities if p)


def optimize(seed, *, cut, rng, restarts, maxfun, height_bound, deadline, starts=()):
    # SciPy is optional for the exact checkers that import this module.
    from scipy.optimize import minimize
    oracle = RectangularOracle(seed["values"], seed["q_x"], seed["q_u"],
                               seed["denominator_degree"], endpoint_check=False)
    dimension = oracle.parameters
    best = None
    queries = failures = 0
    records = []

    def objective(x):
        nonlocal best, queries, failures
        if time.monotonic() >= deadline:
            raise TimeoutError
        theta = np.r_[x[:-1] * x[-1], x[-1]]
        queries += 1
        try:
            rho, px, py, orders, gx, gy = oracle.types(theta, cut)
            total_entropy = entropy(px) + entropy(py)
            if total_entropy <= 1e-12:
                raise ArithmeticError("Both paired types have negligible entropy.")
            cost, derivative = oracle.evaluate(theta, px, py, cut)
            if abs(cost / (2 * total_entropy) - rho) > 1e-7 * max(1, abs(rho)):
                raise ArithmeticError("Paired Gibbs and pressure evaluations disagree.")
            derivative = np.asarray(derivative) / (2 * total_entropy)
            gradient = np.r_[derivative[:-1] * x[-1],
                             derivative[-1] + derivative[:-1] @ x[:-1]]
            if not np.all(np.isfinite(gradient)):
                raise ArithmeticError("Nonfinite search gradient.")
        except ArithmeticError:
            failures += 1
            return 100.0, np.zeros(dimension)
        if best is None or rho < best["rho"]:
            best = {"rho": rho, "parameters": theta.tolist(),
                    "normalized_heights": x[:-1].tolist(), "cut": cut,
                    "type_x": px, "type_y": py, "entropy": total_entropy,
                    "orders_times_log_s": orders, "G_X_times_log_s": gx,
                    "G_Y_times_log_s": gy}
        return rho, gradient

    points = [np.r_[np.zeros(dimension - 1), scale] for scale in (0.125, 0.5, 2.0)]
    for theta in starts:
        if len(theta) != dimension or theta[-1] <= 0:
            raise ValueError("A saved start has incompatible dimension or scale.")
        points.append(np.r_[np.asarray(theta[:-1]) / theta[-1], theta[-1]])
    try:
        for point in points:
            objective(point)
        for restart in range(restarts):
            if best is None:
                break
            point = np.r_[best["normalized_heights"], best["parameters"][-1]]
            if restart:
                point[:-1] += rng.normal(0, 0.5 if restart % 3 else 2.0, dimension - 1)
                point[-1] *= math.exp(rng.normal(0, 0.4))
            point[:-1] = np.clip(point[:-1], -height_bound, height_bound)
            point[-1] = np.clip(point[-1], 0.002, 8.0)
            run = minimize(objective, point, jac=True, method="L-BFGS-B",
                           bounds=[(-height_bound, height_bound)] * (dimension - 1)
                                  + [(0.002, 8.0)],
                           options={"maxfun": maxfun, "maxiter": maxfun,
                                    "maxls": 20, "ftol": 1e-13, "gtol": 1e-7})
            records.append({"restart": restart, "success": bool(run.success),
                            "message": str(run.message), "evaluations": int(run.nfev),
                            "rho": float(run.fun)})
    except TimeoutError:
        records.append({"status": "time budget reached"})
    return {"cut": cut, "best": best, "queries": queries,
            "failures": failures, "restarts": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, default=2)
    parser.add_argument("--indices", type=int, nargs="+")
    parser.add_argument("--pairs", choices=tuple(PAIRS), nargs="+", default=list(PAIRS))
    parser.add_argument("--coefficients", type=int, choices=(1, 2, 3), nargs="+", default=[1, 2, 3])
    parser.add_argument("--gammas", type=int, choices=(0, 1, 2, 3), nargs="+", default=[1])
    parser.add_argument("--repairs", choices=("partial", "inverse"), nargs="+", default=["partial", "inverse"])
    parser.add_argument("--minus-powers", type=int, nargs="+", default=[0])
    parser.add_argument("--plus-powers", type=int, nargs="+", default=[0])
    parser.add_argument("--cuts", type=int, choices=(0, 1), nargs="+", default=[0, 1])
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--maxfun", type=int, default=160)
    parser.add_argument("--seconds", type=float, default=1800)
    parser.add_argument("--candidate-seconds", type=float, default=120)
    parser.add_argument("--height-bound", type=float, default=24)
    parser.add_argument("--seed", type=int, default=310731)
    parser.add_argument("--starts", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    indices = tuple(args.indices) if args.indices else tuple(range(args.modes))
    if not 2 <= len(indices) <= 4 or indices[0] != 0 or any(a >= b for a, b in zip(indices, indices[1:])):
        raise ValueError("Use two through four strictly increasing modes, beginning with zero.")
    if args.seconds <= 0 or args.candidate_seconds <= 0 or args.height_bound <= 0:
        raise ValueError("Search time and height bounds must be positive.")
    ambient = 1 << max(1, indices[-1].bit_length())
    rng = np.random.default_rng(args.seed)
    started = time.monotonic()
    deadline = started + args.seconds
    settings = list(itertools.product(args.pairs, args.coefficients, args.gammas,
                                      args.repairs, args.minus_powers, args.plus_powers))
    starts = []
    if args.starts:
        supplied = json.loads(args.starts.read_text())
        starts = supplied if isinstance(supplied, list) else [supplied["best"]["parameters"]]
    result = {"status": "exploratory floating critical-cubic search; no exact rate certificate",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "indices": list(indices), "ambient": ambient, "random_seed": args.seed,
              "settings": [list(setting) for setting in settings], "cuts": args.cuts,
              "restarts": args.restarts, "maxfun": args.maxfun,
              "normalized_height_bound": args.height_bound, "initial_starts": starts,
              "planned_candidates": len(settings), "candidates": [], "completed": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        result["elapsed_seconds"] = time.monotonic() - started
        successful = [cut["best"]["rho"] for record in result["candidates"]
                      for cut in record.get("cuts", []) if cut["best"]]
        result["best_rho"] = min(successful, default=None)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2) + "\n")
        temporary.replace(args.output)

    for pair, coefficient, gamma, repair, minus_power, plus_power in settings:
        if time.monotonic() >= deadline:
            break
        minus, plus = PAIRS[pair]
        specifications = cubic_pair(indices[-1], minus, plus, coefficient, minus_power, plus_power)
        record = {"pair": pair, "coefficient_product": coefficient, "gamma": gamma,
                  "repair": repair, "minus_power": minus_power, "plus_power": plus_power,
                  "residual_monomials": specifications, "cuts": []}
        result["candidates"].append(record)
        try:
            seed = full_residual_monomial_seed(specifications, indices, ambient, gamma, repair,
                                              check_types=False)
        except ValueError as error:
            record["rejected"] = str(error)
            save()
            continue
        record.update({key: seed[key] for key in
                       ("coefficient_sha256", "denominator_degree", "global_w_power",
                        "residual_rank", "residual_determinant", "degree", "nonzero_entries")})
        for cut in args.cuts:
            if time.monotonic() >= deadline:
                break
            outcome = optimize(seed, cut=cut, rng=rng, restarts=args.restarts,
                               maxfun=args.maxfun, height_bound=args.height_bound,
                               deadline=min(deadline, time.monotonic() + args.candidate_seconds),
                               starts=starts)
            record["cuts"].append(outcome)
            save()
            print(json.dumps({"candidate": len(result["candidates"]),
                              "pair": pair, "coefficient": coefficient, "gamma": gamma,
                              "repair": repair, "powers": [minus_power, plus_power],
                              "cut": cut, "rho": outcome["best"]["rho"] if outcome["best"] else None,
                              "best_rho": result["best_rho"], "queries": outcome["queries"]}), flush=True)
    result["completed"] = (len(result["candidates"]) == len(settings)
                           and all("rejected" in record or len(record["cuts"]) == len(args.cuts)
                                   for record in result["candidates"]))
    save()


if __name__ == "__main__":
    main()
