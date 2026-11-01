#!/usr/bin/env python3
"""Exact minimum-row separation and verified common lower potentials.

The row oracle is exact and accepts nonnegative potentials, including
zero coordinates. The optional cutting-plane driver uses floating LP
solutions only to propose rational potentials; every returned lower
certificate is verified with exact fractions and the full row oracle.
"""

from __future__ import annotations

from fractions import Fraction


def minimum_transfer(graph, potential, point):
    """Exact minimum row costs for the subspace-channel action family."""
    potential = tuple(Fraction(value) for value in potential)
    point = Fraction(point)
    assert len(potential) == len(graph.states) and any(potential)
    assert all(value >= 0 for value in potential) and point > 1
    endpoint_degree = graph.degree if graph.infinity else 0
    rows, image = [], []
    for source, actions in enumerate(graph.actions):
        rank = len(graph.local_factors[source])
        weighted = sorted(actions, key=lambda action: (
            point ** (endpoint_degree-action["order"]) * potential[action["target"]],
            action["order"], action["target"],
        ))
        pivots, row = {}, {}
        for action in weighted:
            increment = 0
            for residue, _ in action["basis"]:
                value = residue
                for pivot in sorted(pivots, reverse=True):
                    if value >> pivot & 1:
                        value ^= pivots[pivot]
                if value:
                    pivots[value.bit_length()-1] = value
                    increment += 1
            if increment:
                target = action["target"]
                value = increment * point ** (endpoint_degree-action["order"])
                row[target] = row.get(target, Fraction(0)) + value
            if len(pivots) == rank:
                break
        assert len(pivots) == rank
        rows.append(row)
        image.append(sum(value * potential[target] for target, value in row.items()))
    return tuple(rows), tuple(image)


def verify_lower_potential(row_oracle, potential, rate):
    """Certify rho(T)>=rate for every policy in the row-product family."""
    potential = tuple(Fraction(value) for value in potential)
    rate = Fraction(rate)
    if rate < 0 or not any(potential) or any(value < 0 for value in potential):
        raise ValueError("expected a nonzero nonnegative potential and rate")
    rows, image = row_oracle(potential)
    margins = tuple(value-rate*height for value, height in zip(image, potential))
    if any(value < 0 for value in margins):
        return None
    positive = [i for i, value in enumerate(potential) if value]
    certificate = {
        "rate": str(rate),
        "potential": tuple(map(str, potential)),
        "minimum_image": tuple(map(str, image)),
        "margins": tuple(map(str, margins)),
        "support": positive,
        "minimum_positive_coordinate_ratio": str(min(image[i]/potential[i] for i in positive)),
        "minimizing_transfer": tuple({target: str(value) for target, value in row.items()} for row in rows),
        "scope": "All-state spectral lower bound for every admissible stationary policy; not a boundary-word or coefficient-rank lower bound",
    }
    if len(positive) == len(potential):
        certificate["selected_policy_upper_rate"] = str(max(value/height for value, height in zip(image, potential)))
    return certificate


def rational_proposals(values):
    """Several rationalizations; acceptance always requires exact replay."""
    base = [max(0.0, float(value)) for value in values]
    candidates = []
    for scale in (10**8, 10**12, 10**16):
        candidates.append(tuple(Fraction(round(value*scale), scale) for value in base))
    candidates.append(tuple(Fraction(str(value)).limit_denominator(10**9) for value in base))
    seen = set()
    for candidate in candidates:
        total = sum(candidate)
        if not total:
            continue
        normalized = tuple(value/total for value in candidate)
        if normalized not in seen:
            seen.add(normalized)
            yield normalized


def find_lower_potential(row_oracle, dimension, rate, maximum_iterations=100):
    """Search by cutting-plane LP, returning only exact lower certificates.

    This practical driver has an iteration cap, not a polynomial iteration
    guarantee. A no-certificate return is inconclusive. SciPy is used only
    here; verify_lower_potential and minimum_transfer use exact arithmetic.
    """
    from scipy.optimize import linprog

    rate = Fraction(rate)
    if dimension < 1 or rate < 0 or maximum_iterations < 1:
        raise ValueError("invalid separation-search parameters")
    cuts = {}

    def add_rows(rows):
        additions = 0
        for source, row in enumerate(rows):
            values = tuple(Fraction(row.get(target, 0)) for target in range(dimension))
            key = source, values
            if key not in cuts:
                cuts[key] = values
                additions += 1
        return additions

    rows, _ = row_oracle((Fraction(1),) * dimension)
    add_rows(rows)
    history = []
    for iteration in range(1, maximum_iterations+1):
        constraints = []
        for (source, _), row in cuts.items():
            constraints.append([float((rate if target == source else 0)-value)
                                for target, value in enumerate(row)] + [1.0])
        # The common margin is also bounded by every coordinate. At margin
        # zero this leaves the entire nonnegative simplex feasible.
        constraints.extend([[-float(i == j) for j in range(dimension)] + [1.0]
                            for i in range(dimension)])
        result = linprog(
            [0.0] * dimension + [-1.0],
            A_ub=constraints, b_ub=[0.0] * len(constraints),
            A_eq=[[1.0] * dimension + [0.0]], b_eq=[1.0],
            bounds=[(0.0, 1.0)] * dimension + [(None, None)],
            method="highs",
        )
        if not result.success:
            return {"status": "no_exact_certificate", "reason": "floating LP did not produce a candidate",
                    "iterations": iteration, "cuts": len(cuts), "history": history}
        proposals = tuple(rational_proposals(result.x[:dimension]))
        for potential in proposals:
            certificate = verify_lower_potential(row_oracle, potential, rate)
            if certificate is not None:
                return {"status": "certified", "iterations": iteration,
                        "cuts": len(cuts), "certificate": certificate, "history": history}
        potential = proposals[-1]
        rows, image = row_oracle(potential)
        additions = add_rows(rows)
        history.append({"iteration": iteration, "cuts": len(cuts),
                        "new_cuts": additions, "floating_margin": float(result.x[-1]),
                        "worst_exact_margin": str(min(value-rate*height
                                                        for value, height in zip(image, potential)))})
        if additions == 0:
            return {"status": "no_exact_certificate", "reason": "no new row and no verified rational candidate",
                    "iterations": iteration, "cuts": len(cuts), "history": history}
    return {"status": "no_exact_certificate", "reason": "cutting-plane iteration budget exhausted",
            "iterations": maximum_iterations, "cuts": len(cuts), "history": history}
