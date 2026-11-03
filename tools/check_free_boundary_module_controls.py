#!/usr/bin/env python3
"""Free-boundary q=2 module controls for the ordered auxiliary certificate."""
import json
from fractions import Fraction
from itertools import product

import check_cancellation_module_automata as module
from check_cancellation_module_crossing import local_crossing

FAMILIES = {
    "elementary": module.elementary(2)[0],
    "pairing": module.pairing()[0],
    "mixed_loop": module.mixed_loop()[0],
    "module_crossing": local_crossing(),
}
SCALES = tuple(map(Fraction, ("1", "3/2", "2", "3", "5")))
PRECISIONS = (1, 2, 3, 4, 5, 6)
HEIGHTS = ((0, 0), (0, 1), (0, -1))


def best_free_boundary(local):
    candidates = {}
    maximum_states = 0
    for cut, side, height, infinity in product(range(2), range(2), HEIGHTS, (False, True)):
        factors = module.local_factors(local, cut, side, height, infinity)
        for precision in PRECISIONS:
            automaton = module.ModuleAutomaton(
                factors["channels"], precision, free_boundaries=True
            )
            maximum_states = max(maximum_states, len(automaton.states))
            for s in SCALES:
                matrix, _ = automaton.weighted_matrix(factors["exponents"], s)
                beta, _ = module.potential(matrix)
                key = (s, cut, infinity)
                if key not in candidates or beta < candidates[key][0]:
                    candidates[key] = (beta, precision, len(automaton.states))
    best = None
    for s in SCALES:
        for cut in (0, 1):
            high = candidates.get((s, cut, True))
            low = candidates.get((s, 1 - cut, False))
            if high is None or low is None:
                continue
            value = high[0] * low[0]
            if best is None or value < best[0]:
                best = value, s, cut, high, low
    return {
        "best_B": str(best[0]),
        "best_scale": str(best[1]),
        "infinity_cut": best[2],
        "infinity_beta": str(best[3][0]),
        "zero_beta": str(best[4][0]),
        "precisions": sorted({best[3][1], best[4][1]}),
        "maximum_reachable_states": maximum_states,
    }


def main():
    report = {
        "alphabet": 2,
        "auxiliary_dimension": 2,
        "free_boundaries": True,
        "threshold": "2^(2*sqrt2)",
        "families": {name: best_free_boundary(local) for name, local in FAMILIES.items()},
        "scope": (
            "Finite q=2 D=2 module controls and sampled precisions/scales. "
            "Not a theorem and not a score certificate."
        ),
    }
    assert report["families"]["elementary"]["best_B"] == "256/27"
    assert report["families"]["pairing"]["best_B"] == "16"
    assert report["families"]["mixed_loop"]["best_B"] == "16"
    assert Fraction(report["families"]["module_crossing"]["best_B"]) > Fraction(710, 100)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
