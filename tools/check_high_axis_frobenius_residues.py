#!/usr/bin/env python3
"""Exact higher-form identities and independent prefix obstruction for one kernel."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

from check_high_axis_plucker import toeplitz_ranks
from high_axis_frobenius_residues import crossblock,residue
from high_axis_plucker import MUL,power
from native_matrix_smith_f4 import MatrixOracle
from regularized_transpose_words import field_rank
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import pack

ROOT = Path(__file__).resolve().parents[1]


def symbolic_form(size,exponent):
    result = {}
    for indices in combinations(range(8),size):
        result[sum(1 << i for i in indices)] = {tuple(exponent if i in indices else 0 for i in range(8))}
    return result


def symbolic_product(first,second):
    result = {}
    for a,terms in first.items():
        for b,others in second.items():
            if a & b:
                continue
            target = result.setdefault(a | b,set())
            for x in terms:
                for y in others:
                    value = tuple(i+j for i,j in zip(x,y))
                    if value in target:
                        target.remove(value)
                    else:
                        target.add(value)
    return {mask:terms for mask,terms in result.items() if terms}


def leading_form(size,base,raised,extra):
    result = {}
    for indices in combinations(range(8),size):
        terms = set()
        for high in combinations(indices,raised):
            terms.add(tuple((base+extra if i in high else base) if i in indices else 0 for i in range(8)))
        result[sum(1 << i for i in indices)] = terms
    return result


def evaluate(form,values):
    result = {}
    for mask,terms in form.items():
        coefficient = 0
        for exponents in terms:
            value = 1
            for x,e in zip(values,exponents):
                value = MUL[value][power(x,e)]
            coefficient ^= value
        result[mask] = coefficient
    return result


def multiplication(form):
    return [[form.get(i ^ j,0) if j & ~i == 0 else 0 for j in range(256)] for i in range(256)]


def literal_entry(subset,values):
    selected = [v for i,v in enumerate(values) if subset >> i & 1]
    e,j,f = residue(len(selected)+1)
    factor = 1
    for v in selected:
        factor = MUL[factor][power(v,e)]
    coefficients = {}
    # Position zero is the free axis. Other positions are field constants.
    for indices in combinations(range(len(selected)+1),j):
        value = factor
        for i in indices:
            if i:
                value = MUL[value][power(selected[i-1],f)]
        exponent = e+(f if 0 in indices else 0)
        coefficients[exponent] = coefficients.get(exponent,0) ^ value
    low,high = 0,0
    for exponent,value in coefficients.items():
        low |= (value & 1) << exponent
        high |= ((value >> 1) & 1) << exponent
    return low,high


def run():
    expected = [(64,0,0),(32,0,0),(21,1,1),(16,0,0),(12,1,4),(10,2,2),(9,1,1),(8,0,0),(7,1,1)]
    assert [residue(k) for k in range(1,10)] == expected
    quartic = leading_form(4,12,1,4)
    qright = symbolic_product(symbolic_product(symbolic_form(1,16),symbolic_form(1,12)),symbolic_form(2,12))
    quintic = leading_form(5,10,2,2)
    fright = symbolic_product(symbolic_product(symbolic_form(2,12),symbolic_form(1,10)),symbolic_form(2,10))
    assert quartic == qright and quintic == fright
    old = multiplication({mask:1 for mask in range(256) if mask.bit_count() == 4})
    generic_fixture = multiplication(evaluate(quartic,[2]+[1]*7))
    assert field_rank(old) == 46 and field_rank(generic_fixture) == 28
    entries = 0
    for values in ([1]*8,[2]+[1]*7):
        matrix = crossblock(9,values)
        for i in range(256):
            for j in range(256):
                if j & ~i:
                    assert matrix[i][j] == (0,0)
                else:
                    assert matrix[i][j] == literal_entry(i ^ j,values)
                    entries += 1
    controls = []
    for dimension in (5,6):
        for nonbinary in (False,True):
            values = [1]*(dimension-1)
            if nonbinary:
                values[0] = 2
            matrix = crossblock(dimension,values)
            exact,_ = WeightedSmithOracle([[pack(v) for v in row] for row in matrix]).profile([0]*len(matrix),[0]*len(matrix))
            native = MatrixOracle(matrix,len(matrix),precision=256).profile()['orders']
            assert list(map(int,native)) == list(exact)
            controls.append(dict(dimension=dimension,values=values,orders=exact))
    path = ROOT/'evidence/high-axis-frobenius-residues.json'
    diagnostic = json.loads(path.read_text())
    assert diagnostic['status'] == 'stopped_by_independent_all_scale_prefix' and len(diagnostic['cases']) == 1
    record = diagnostic['cases'][0]
    assert record['label'] == 'all-one'
    matrix = crossblock(9,[1]*8)
    assert MatrixOracle(matrix,256,precision=256).profile()['orders'] == record['orders']
    ranks = toeplitz_ranks(matrix,29)
    assert ranks == record['independent_toeplitz_ranks'] and ranks[-2:] == [840,960]
    m,total = ranks[-1]-ranks[-2],29*(ranks[-1]-ranks[-2])-ranks[-1]
    assert (m,total) == (120,2520) and 3*m > 256 and 3*total < 64*m
    names = ['tools/check_high_axis_frobenius_residues.py','tools/high_axis_frobenius_residues.py',
             'tools/high_axis_plucker.py','tools/high_axis_cayley.py','tools/check_high_axis_plucker.py',
             'tools/weighted_smith_oracle_f4.py','tools/native_matrix_smith_f4.py']
    return dict(scope='Generic quartic factor identity and rank change, literal polynomial controls, and an independent all-scale stopping prefix; no stronger factoring bound.',
                quartic_symbolic_monomials=sum(map(len,quartic.values())),
                quintic_symbolic_monomials=sum(map(len,quintic.values())),
                old_quartic_rank=46,new_quartic_generic_rank=28,
                independently_expanded_nonzero_pattern_entries=entries,
                small_complete_profiles=controls,independent_untruncated_orders=sum(len(c['orders']) for c in controls),
                toeplitz_ranks_precisions_1_through_29=ranks,stopping_channels=m,stopping_prefix_sum=total,
                input_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    begun = time.monotonic()
    result = json.loads(json.dumps(run()))
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/high-axis-frobenius-residue-controls.json').read_text())
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','small_complete_profiles','toeplitz_ranks_precisions_1_through_29')}))
    print(json.dumps({'seconds':time.monotonic()-begun}))


if __name__ == '__main__':
    main()
