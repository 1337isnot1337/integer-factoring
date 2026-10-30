#!/usr/bin/env python3
"""Construct exact projectors in the odd-characteristic Weyl coefficient algebra."""

import argparse
import json
from math import comb
from pathlib import Path


def make_algebra(p, h):
    assert h % p == 0
    m = h//p
    weights = [1]
    w = 1
    while w <= m:
        weights.append(w)
        w *= p
    H = 2*m+sum(weights)
    D = m+1
    table = [[None]*D for _ in range(D)]
    for i in range(D):
        for j in range(D):
            value = [0]*D
            for k in range(min(i, j)+1):
                n = i+j-k
                if n <= m:
                    value[n] = (comb(n, i)*comb(n, j)*comb(H-n-1, k)) % p
            table[i][j] = [(k, c) for k, c in enumerate(value) if c]
    for i in range(D):
        for j in range(D):
            assert table[i][j] == table[j][i]
            for k in range(D):
                left, right = [0]*D, [0]*D
                for n, c in table[i][j]:
                    for r, value in table[n][k]:
                        left[r] = (left[r]+c*value) % p
                for n, c in table[j][k]:
                    for r, value in table[i][n]:
                        right[r] = (right[r]+c*value) % p
                assert left == right

    def multiply(a, b):
        result = [0]*D
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if y:
                    for k, c in table[i][j]:
                        result[k] = (result[k]+x*y*c) % p
        return result

    one = [1]+[0]*m
    def power(a, n):
        result = one[:]
        while n:
            if n & 1:
                result = multiply(result, a)
            a = multiply(a, a)
            n >>= 1
        return result

    basis = [[int(i == j) for i in range(D)] for j in range(D)]
    exponent = 1
    while exponent < D:
        exponent *= p
    projectors = [one]
    separating = []
    for i, x in enumerate(basis[1:], 1):
        semisimple = power(x, exponent)
        assert power(semisimple, p) == semisimple
        splitters = []
        for scalar in range(p):
            shifted = semisimple[:]
            shifted[0] = (shifted[0]-scalar) % p
            value = power(shifted, p-1)
            value = [(-a) % p for a in value]
            value[0] = (value[0]+1) % p
            splitters.append(value)
        refined = [value for e in projectors for s in splitters
                   if any(value := multiply(e, s))]
        if len(refined) != len(projectors):
            separating.append(i)
        projectors = refined
    assert [sum(e[i] for e in projectors) % p for i in range(D)] == one
    for i, e in enumerate(projectors):
        for j, f in enumerate(projectors):
            assert multiply(e, f) == (e if i == j else [0]*D)
    raw = [[((-1)**j*comb(i, j)) % p if i >= j else 0
            for i in range(D)] for j in range(D)]
    assert [sum(a[i] for a in raw) % p for i in range(D)] == one
    zero_projector = projectors[0]
    assert zero_projector[0] == 1
    assert all(not any(multiply(zero_projector, power(x, exponent))) for x in basis[1:])
    assert power(raw[0], exponent) == zero_projector
    assert multiply(raw[0], zero_projector) == raw[0]
    masks = [(f"raw_{j}", a) for j, a in enumerate(raw)]
    masks += [(f"projector_{j}", a) for j, a in enumerate(projectors)]
    masks += [(f"refined_{e}_{j}", multiply(pr, a))
              for e, pr in enumerate(projectors) for j, a in enumerate(raw)
              if any(multiply(pr, a))]
    return {"characteristic": p, "original_energy": h, "reduced_energy": m,
            "weights": weights, "cartan_weight": H,
            "abstract_dimension": D, "projector_count": len(projectors),
            "basis_associativity_checks": D**3,
            "frobenius_exponent": exponent,
            "zero_character_local_unit_verified": True,
            "separating_indices": separating, "projectors": projectors,
            "masks": masks}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("--energy", type=int, default=27)
    parser.add_argument("--mask-output", type=Path)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()
    result = make_algebra(args.prime, args.energy)
    if args.mask_output:
        args.mask_output.write_text(str(len(result["masks"]))+"\n"+"\n".join(
            name+" "+" ".join(map(str, value)) for name, value in result["masks"])+"\n")
    if args.json_output:
        args.json_output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key: value for key, value in result.items() if key != "masks"}))


if __name__ == "__main__":
    main()
