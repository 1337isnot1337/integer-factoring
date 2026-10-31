#!/usr/bin/env python3
"""Replay the stored large relation matrix, determinant and torsion certificate."""
import json
from pathlib import Path

from arithmetic_norm_sieve import prime_base, roots_x, verify
from binary_quadratic_forms import reduced_forms


def check_roots():
    controls = 0
    for discriminant in (-23, -84, -348, -66548):
        base = prime_base(discriminant, 23)
        for a, b, c in reduced_forms(discriminant):
            for y in range(4):
                for prime, root, _ in base:
                    expected = {x for x in range(prime) if (a*x*x+b*x*y+c*y*y) % prime == 0}
                    assert set(roots_x((a, b, c), y, prime, root)) == expected
                    controls += 1
    return controls


if __name__ == '__main__':
    path = Path(__file__).resolve().parents[1]/'evidence/arithmetic-norm-sieve-60bit.json'
    certificate = json.loads(path.read_text())
    assert certificate['status'] == 'factor_found'
    controls = check_roots()
    result = verify(certificate)
    result['small_residue_controls'] = controls
    print(json.dumps(result, indent=2))
