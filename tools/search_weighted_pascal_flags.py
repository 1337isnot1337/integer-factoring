#!/usr/bin/env python3
"""Explore unrestricted mode gauges of weighted Pascal kernels and binary flags.

The polynomial kernel and Clifford identities are exact. Native polynomial
arithmetic is exact, but weighted pivots and pressures here use floating point.
Every proposed improvement therefore needs a separate exact replay.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import expit

from native_matrix_smith_f4 import MatrixOracle
from permuted_pascal_kernels import permuted_pascal_kernel
from weighted_pascal_flags import change_flag, clifford_data, natural_heights, weighted_kernel


class WeightedPascalPressure:
    def __init__(self, seed, precision=None):
        self.seed = seed
        self.n, self.outer = seed["n"], seed["outer"]
        self.size, self.queries, self.last = 4 * self.n, 0, None
        self.cuts = []
        if precision is None:
            precision = max(512, 1 << (16 * self.outer - 1).bit_length())
        for place in (0, 1):
            matrix, shifts, rows, columns = clifford_data(seed, place)
            oracle = MatrixOracle([[(value, 0) for value in row] for row in matrix],
                                  4 * self.n, shifts, precision)
            self.cuts.append((oracle, np.asarray(rows), np.asarray(columns)))

    def natural(self, scale):
        return np.asarray([scale * height for place in natural_heights(self.seed)
                           for height in place] + [scale])

    def evaluate(self, parameters, gradient=False):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (2 * self.size + 1,) or not np.all(np.isfinite(parameters)):
            raise ValueError("Supply one finite height per physical mode and place, plus the scale.")
        scale = parameters[-1]
        if scale <= 0:
            raise ValueError("The profile scale must be positive.")
        modes = 2 * self.n
        log_b = -self.outer * self.n * scale - parameters[:-1].sum() / 2
        derivative = np.r_[np.full(2 * self.size, -.5), -float(self.outer * self.n)]
        records = []
        for place, (oracle, rows, columns) in enumerate(self.cuts):
            heights = parameters[place * self.size:(place + 1) * self.size]
            row, column = heights[rows], heights[columns]
            result = oracle.profile(np.r_[row, -row], np.r_[column, -column],
                                    scale, infinity=bool(place))
            orders = np.asarray(result["orders"])
            if np.max(np.abs(orders + orders[::-1])) > 2e-7:
                raise ArithmeticError("The complete Clifford pairing was not resolved.")
            log_b += .5 * np.logaddexp(orders / 2, -orders / 2).sum()
            coefficients = .25 * np.tanh(orders / 2)
            base_orders = np.diff([0] + result["base_prefix"])
            for coefficient, r, c, base in zip(coefficients, result["rows"], result["cols"], base_orders):
                derivative[place * self.size + rows[r % modes]] += coefficient * (1 if r < modes else -1)
                derivative[place * self.size + columns[c % modes]] += coefficient * (1 if c < modes else -1)
                derivative[-1] += coefficient * base
            records.append(result)
        heights = parameters[:-1].reshape(2, 4, self.n)
        physical = np.asarray([heights[:, 0, :].sum(axis=0) + heights[:, 3, :].sum(axis=0),
                               heights[:, 1, :].sum(axis=0) + heights[:, 2, :].sum(axis=0)])

        def equation(rate):
            return float(log_b - rate * np.logaddexp(0, -physical / rate).sum())

        if equation(1e-8) <= 0:
            raise ArithmeticError("The positive Gibbs root is unresolved near zero entropy.")
        high = 4.
        while equation(high) > 0:
            high *= 2
            if high > 1e8:
                raise ArithmeticError("The paired pressure root was not bracketed.")
        rate = brentq(equation, 1e-8, high, xtol=1e-13)
        probabilities = expit(-physical / rate)
        entropy = float((np.logaddexp(0, -physical / rate)
                         + physical / rate * probabilities).sum())
        if entropy < 1e-8:
            raise ArithmeticError("The physical entropy is too small for a stable rate.")
        for place in range(2):
            for leg in range(4):
                start = place * self.size + leg * self.n
                derivative[start:start + self.n] += probabilities[0 if leg in (0, 3) else 1]
        self.queries += 1
        self.last = {"rho": float(rate / 2), "entropy": entropy,
                     "critical_pressure": equation(math.sqrt(2)),
                     "parameters": parameters.tolist(), "clifford_profiles": records}
        return (rate / 2, derivative / (2 * entropy)) if gradient else rate / 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    family = parser.add_mutually_exclusive_group(required=True)
    family.add_argument("--weights", type=int, nargs="+")
    family.add_argument("--images", type=int, nargs="+",
                        help="Images of an arbitrary power-of-two exponent permutation.")
    parser.add_argument("--outer", type=int,
                        help="Required with --weights; equals the permutation length with --images.")
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=600)
    parser.add_argument("--seconds", type=float, default=1800)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--flag-shears", type=int, default=0)
    parser.add_argument("--precision", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if min(args.restarts, args.iterations, args.seconds) <= 0 or args.flag_shears < 0:
        raise ValueError("Use positive search budgets and a nonnegative flag-shear count.")
    if args.images is None:
        if args.outer is None:
            parser.error("--weights requires --outer")
        seed = weighted_kernel(args.weights, args.outer)
    else:
        if args.outer is not None and args.outer != len(args.images):
            parser.error("--outer must equal the permutation length when supplied with --images")
        seed = permuted_pascal_kernel(args.images)
        args.outer = seed["outer"]
    rng = np.random.default_rng(args.seed)
    if args.flag_shears:
        constant = np.eye(seed["n"], dtype=np.uint8)
        for _ in range(args.flag_shears):
            row, source = rng.choice(seed["n"], size=2, replace=False)
            constant[row] ^= constant[source]
        seed = change_flag(seed, constant.tolist())
    pressure = WeightedPascalPressure(seed, args.precision)
    fitted_scale = minimize_scalar(lambda x: pressure.evaluate(pressure.natural(math.exp(x))),
                                  bounds=(-8., math.log(2.)), method="bounded")
    scale = math.exp(fitted_scale.x)
    natural = pressure.natural(scale)
    pressure.evaluate(natural)
    baseline, best = pressure.last.copy(), pressure.last.copy()
    control = natural + rng.normal(size=len(natural)) * .00003
    control[-1] = scale
    _, analytic = pressure.evaluate(control, True)
    gradient_error = 0.
    for coordinate in list(range(min(10, len(control) - 1))) + [len(control) - 1]:
        step = 1e-7 if coordinate + 1 < len(control) else 1e-7 / (args.outer * seed["n"])
        delta = np.zeros_like(control)
        delta[coordinate] = step
        finite = (pressure.evaluate(control + delta) - pressure.evaluate(control - delta)) / (2 * step)
        gradient_error = max(gradient_error, abs(float(finite - analytic[coordinate])))
    if gradient_error > 1e-4:
        raise ArithmeticError("The local gradient control was unresolved: " + str(gradient_error))
    start, runs = time.monotonic(), []
    deadline = start + args.seconds
    root = Path(__file__).resolve().parent.parent
    sources = ["tools/weighted_pascal_flags.py", "tools/search_weighted_pascal_flags.py",
               "research/weighted-pascal-flag-crossings.md"]
    if args.images is not None:
        sources += ["tools/permuted_pascal_kernels.py", "research/permuted-pascal-zero-spectra.md"]
    report = {"status": "exploratory floating search", "weights": args.weights,
              "exponent_images": args.images,
              "outer": args.outer, "constant_flag": seed.get("constant_flag"),
              "random_seed": args.seed, "baseline": baseline, "best": best, "runs": runs,
              "gradient_max_error": gradient_error, "rejected_queries": 0,
              "budget_exhausted": False,
              "source_sha256": {name: sha256((root / name).read_bytes()).hexdigest() for name in sources},
              "scope": "Complete Clifford profiles are used; finite floating optimization is neither a strict improvement certificate nor a universal exclusion."}

    def save():
        report.update({"best": best, "queries": pressure.queries,
                       "seconds": time.monotonic() - start})
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(args.output)

    save()
    print(json.dumps({"weights": args.weights, "exponent_images": args.images, "outer": args.outer,
                      "baseline_rho": baseline["rho"], "gradient_error": gradient_error}), flush=True)
    for trial in range(args.restarts):
        if time.monotonic() >= deadline:
            report["budget_exhausted"] = True
            break
        initial = np.r_[natural[:-1] / (scale * args.outer), math.log(scale)]
        if trial:
            initial[:-1] += rng.normal(size=len(initial) - 1) * (.05 if trial < args.restarts // 2 else .5)

        def objective(parameters):
            nonlocal best
            if time.monotonic() >= deadline:
                raise TimeoutError("The requested search budget is exhausted.")
            current_scale = math.exp(parameters[-1])
            current = np.r_[current_scale * args.outer * parameters[:-1], current_scale]
            try:
                value, derivative = pressure.evaluate(current, True)
            except ArithmeticError:
                report["rejected_queries"] += 1
                return 2., np.zeros_like(parameters)
            if value < best["rho"] - 1e-12:
                best = pressure.last.copy()
                save()
            transformed = np.r_[derivative[:-1] * current_scale * args.outer, derivative @ current]
            return value, transformed

        try:
            fitted = minimize(objective, initial, jac=True, method="SLSQP",
                              bounds=[(-8., 8.)] * (len(initial) - 1)
                              + [(math.log(.0001), math.log(2.))],
                              options={"maxiter": args.iterations, "ftol": 1e-11})
            runs.append({"trial": trial, "rho": float(fitted.fun), "iterations": int(fitted.nit),
                         "success": bool(fitted.success), "message": str(fitted.message)})
        except TimeoutError:
            report["budget_exhausted"] = True
            break
        save()
        print(json.dumps({"trial": trial, "best_rho": best["rho"],
                          "queries": pressure.queries, "seconds": report["seconds"]}), flush=True)
    save()


if __name__ == "__main__":
    main()
