#!/usr/bin/env python3
"""Independent finite checks of weighted Weyl and coefficient interfaces."""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import comb, factorial
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
if not (HERE / "check_small_space_berkowitz.py").exists():
    sys.path.insert(0, str(HERE.parent / "integer-factoring" / "tools"))
from check_small_space_berkowitz import gaussian  # noqa: E402


@lru_cache(None)
def labels(weights, energy):
    if energy < 0:
        return ()
    return tuple(alpha for alpha in product(*(range(energy//w+1) for w in weights))
                 if sum(a*w for a, w in zip(alpha, weights)) == energy)


def binomial(top, bottom):
    result = 1
    for n, k in zip(top, bottom):
        if k < 0 or k > n:
            return 0
        result *= comb(n, k)
    return result


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x-y for x, y in zip(a, b))


def nonnegative(a):
    return all(x >= 0 for x in a)


def coefficient_vector(alpha, y, u, beta, weights, h, p):
    result = [0] * (h+1)
    if add(alpha, u) != add(beta, y):
        return result
    for gamma in product(*(range(min(a,b)+1) for a,b in zip(alpha,y))):
        scalar = binomial(beta, subtract(alpha,gamma)) * binomial(y,gamma)
        exponent = sum(w*g for w,g in zip(weights,gamma))
        for i in range(exponent+1):
            result[i] = (result[i] + scalar * (-1)**i * comb(exponent,i)) % p
    return result


def b_entry(alpha, y, u, beta, weights, e, p):
    if add(alpha,u) != add(beta,y):
        return 0
    return sum(binomial(beta,k) * binomial(y,subtract(alpha,k))
               for k in labels(weights,e)) % p


def row_factor_entry(alpha, y, u, beta, weights, e, p):
    result = 0
    for k, ell in product(labels(weights,e), repeat=2):
        v, other_v = subtract(alpha,k), subtract(y,ell)
        w, other_w = subtract(beta,k), subtract(u,ell)
        if v == other_v and w == other_w and nonnegative(v) and nonnegative(w):
            result += binomial(y,ell) * binomial(beta,k)
    return result % p


def column_factor_entry(alpha, y, u, beta, weights, h, i, p):
    left_total, right_total = add(alpha,u), add(beta,y)
    result = 0
    for gamma in labels(weights,i):
        for eta in labels(weights,2*h-i):
            if add(gamma,eta) == left_total == right_total:
                result += binomial(eta,subtract(alpha,gamma)) * (-1)**sum(gamma) * binomial(y,gamma)
    return result % p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    row_checks = column_checks = coefficient_checks = rank_checks = crossings = 0
    profiles = []
    for p, weights, h in [(2,(1,2),4),(2,(1,2,4),4),(3,(1,1),4),
                           (3,(1,3),5),(5,(1,5),6)]:
        domain = labels(weights,h)
        d = len(domain)
        rows = [[[0]*(d*d) for _ in range(d*d)] for _ in range(h+1)]
        cols = [[[0]*(d*d) for _ in range(d*d)] for _ in range(h+1)]
        for a,b,c,v in product(range(d),repeat=4):
            alpha,y,u,beta=domain[a],domain[b],domain[c],domain[v]
            coefficients=coefficient_vector(alpha,y,u,beta,weights,h,p)
            assert sum(coefficients)%p == int(a==v and b==c)
            crossings += 1
            b_values=[]
            for e in range(h+1):
                value=b_entry(alpha,y,u,beta,weights,e,p)
                assert value == row_factor_entry(alpha,y,u,beta,weights,e,p)
                b_values.append(value)
                row_checks += 1
            for i,coefficient in enumerate(coefficients):
                expanded=sum((-1)**i*comb(h-e,i)*b_values[e] for e in range(h-i+1))%p
                assert coefficient==expanded
                coefficient_checks += 1
                assert coefficient==column_factor_entry(alpha,y,u,beta,weights,h,i,p)
                column_checks += 1
                rows[i][a*d+b][c*d+v]=coefficient
                cols[i][a*d+c][b*d+v]=coefficient
        row_ranks=[gaussian(matrix,p)[0] for matrix in rows]
        col_ranks=[gaussian(matrix,p)[0] for matrix in cols]
        for i in range(h+1):
            row_bound=sum(len(labels(weights,e))**2 for e in range(h-i+1) if comb(h-e,i)%p)
            col_bound=len(labels(weights,i))*len(labels(weights,2*h-i))
            assert row_ranks[i]<=row_bound and col_ranks[i]<=col_bound
            rank_checks += 2
        profiles.append({"prime":p,"weights":weights,"energy":h,"labels":d,
                         "row_ranks":row_ranks,"col_ranks":col_ranks})

    # The one-channel identity over a multivariate truncated box is exact,
    # but its common top scalar squares to zero.
    top_identities = 0
    for weights,h in [((1,2),2),((1,2),4),((1,2,4),4)]:
        domain=labels(weights,h)
        bounds=tuple(max(a[j] for a in domain) for j in range(len(weights)))
        for alpha,y,u,beta in product(domain,repeat=4):
            t=tuple(a+k-b for a,k,b in zip(alpha,bounds,beta))
            s=tuple(a+k-b for a,k,b in zip(y,bounds,u))
            survives=all(a<=k for a,k in zip(t,bounds)) and all(a<=k for a,k in zip(s,bounds))
            assert survives == (alpha==beta and y==u)
            if survives:
                assert t==s==bounds
            top_identities += 1
        assert any(k>0 for k in bounds)
        assert any(2*k>k for k in bounds)
        box_size=1
        for k in bounds:box_size*=k+1
        assert box_size**2>=len(domain)**2

    divided_power_checks=0
    for p in (2,3,5):
        for alpha in product(range(5),repeat=2):
            if alpha==(0,0):continue
            coefficient=1
            for a in alpha:coefficient*=factorial(p*a)//factorial(a)**p
            assert coefficient%p==0
            divided_power_checks+=1

    result={"status":"passed",
            "scope":"Finite coefficient/factor identities and exact ranks; no improved seed or asymptotic rank claim.",
            "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "crossing_entries":crossings,"row_factor_entries":row_checks,
            "column_factor_entries":column_checks,"basis_coefficient_entries":coefficient_checks,
            "rank_bound_comparisons":rank_checks,"profiles":profiles,
            "nonunit_top_coefficient_entries":top_identities,
            "divided_power_nilpotence_checks":divided_power_checks}
    rendered=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:args.output.write_text(rendered)
    print(rendered,end="")


if __name__=="__main__":
    main()
