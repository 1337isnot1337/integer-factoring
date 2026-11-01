#!/usr/bin/env python3
"""Remove an exact common parameter ramification from an explicit crossing."""
from __future__ import annotations

from functools import reduce
from math import gcd, prod

from weyl_lattice_filter_synthesis import polynomial_degrees


def deflate(seed, ramification=None):
    values = seed['values']
    powers = [power for pair in values for plane in pair for power in polynomial_degrees(plane) if power]
    available = reduce(gcd, powers) if powers else 1
    ramification = available if ramification is None else ramification
    if not isinstance(ramification, int) or ramification < 1 or available % ramification:
        raise ValueError('The positive ramification must divide every nonzero polynomial exponent.')
    actual = max((plane.bit_length()-1 for pair in values for plane in pair), default=0)
    if actual != seed['degree']:
        raise ValueError('The input degree must be its actual coefficient degree.')
    result = []
    for pair in values:
        output = tuple(sum(1 << (power//ramification) for power in polynomial_degrees(plane)) for plane in pair)
        recovered = tuple(sum(1 << (power*ramification) for power in polynomial_degrees(plane)) for plane in output)
        assert recovered == tuple(pair)
        result.append(output)
    dimensions = seed.get('dimensions', [seed['q_x'], seed['q_u'], seed['q_u'], seed['q_x']])
    assert len(dimensions) == 4 and prod(dimensions) == len(result)
    for entry, pair in enumerate(result):
        remainder, labels = entry, []
        for size in dimensions:
            labels.append(remainder % size)
            remainder //= size
        x, y, u, v = labels
        assert ((pair[0].bit_count() & 1) | ((pair[1].bit_count() & 1) << 1)) == int(x == v and y == u)
    return seed | {'values': result, 'degree': actual//ramification,
                   'construction': f'Exact parameter deflation by {ramification}',
                   'parameter_ramification_removed': ramification}
