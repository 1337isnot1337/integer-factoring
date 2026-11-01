#!/usr/bin/env python3
"""Specialize the published rectangular cofactor-jet cap to our q6 point."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THEOREM_COMMIT = '365f319fbb62c734bb409c93d39323d566b87c49'
THEOREM_PATH = 'research/rectangular-adjugate-first-jet.md'


def certificate():
    path = ROOT/'evidence/hereditary-q6-integer-target.json'
    point = json.loads(path.read_text())
    d, r, energies = point['deficit'], point['target'], point['energies']
    p, b = 2, 6
    x, y = energies[:p], energies[:b]
    n, h = p*b, min(p, b)**2
    c = n-h
    weight = b*sum(x)+p*sum(y)
    # Equation (5) at k=n: L_n=c(n-1), TopCell_n=W.
    coefficient = c*(n-1)
    rhs = n*((n-1)*d+weight)-sum(a+v for a in x for v in y)
    lhs = n*sum(r[:n-1])+coefficient*r[1]
    assert (n, c, coefficient, weight) == (12, 8, 88, 140)
    assert (sum(r[:n-1]), r[1], lhs, rhs) == (162, 6, 2472, 2332)
    assert lhs-rhs == 140
    return dict(scope='Published all-field necessary cap rejects the retained q6 point and its dominating average-prefix sector.',
                theorem_commit=THEOREM_COMMIT, theorem_path=THEOREM_PATH,
                theorem_equation=5, physical_alphabet=[p, b], cofactor_nullity=c,
                first_prefix=n-1, first_coefficient=n, second_prefix=2,
                second_coefficient=coefficient, lhs=lhs, rhs=rhs, violation=lhs-rhs,
                target_sha256=sha256(path.read_bytes()).hexdigest(),
                checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/hereditary-q6-rectangular-rejection.json').read_text())
    print(json.dumps(result))


if __name__ == '__main__':
    main()
