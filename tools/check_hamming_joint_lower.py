#!/usr/bin/env python3
"""Check the finite joint-span lower bound using a retained RUP refutation.

Default replay needs only Python's standard library and a C11 compiler.
--regenerate-proof additionally needs python-sat; it is not used for checking.
"""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

from check_modular_recombination import (
    CROSSING_SYMMETRIES, candidates, crossing, decode, encode, insert, permute, represent,
)
from hamming_sectors import binary_sector_labels, sector_matrices

ROOT = Path(__file__).resolve().parents[1]
ELL = 0x686610c0
SEPARATOR_ENTRIES = (20, 24, 66)


class CNF:
    def __init__(self):
        self.variables, self.clauses, self.definitions, self.gates = 0, [], [], {}
        self.false = self.new()
        self.definitions.append((self.false, 'false'))
        self.add(-self.false)

    def new(self):
        self.variables += 1
        return self.variables

    def add(self, *literals):
        self.clauses.append(list(literals))

    def xor(self, a, b):
        if a == self.false:
            return b
        if b == self.false:
            return a
        if a == b:
            return self.false
        key = tuple(sorted((a, b)))
        if key not in self.gates:
            t = self.new()
            self.gates[key] = t
            self.definitions.append((t, 'xor', a, b))
            self.add(-a, -b, -t)
            self.add(a, b, -t)
            self.add(a, -b, t)
            self.add(-a, b, t)
        return self.gates[key]

    def conj(self, a, b):
        t = self.new()
        self.definitions.append((t, 'and', a, b))
        self.add(-a, -b, t)
        self.add(a, -t)
        self.add(b, -t)
        return t

    def equal(self, a, b):
        self.add(a, -b)
        self.add(-a, b)

    def dimacs(self):
        header = f'p cnf {self.variables} {len(self.clauses)}\n'
        return (header+''.join(' '.join(map(str, c))+' 0\n' for c in self.clauses)).encode()


def source_basis():
    basis = {}
    for count in (1, 2):
        free, domain = binary_sector_labels(2, count)
        for vector in candidates(sector_matrices(free, count, 4), 4, False):
            insert(basis, vector, 0)
    assert len(basis) == 33
    return [v for v, mask in basis.values()]


def rank_formula(generators, rank=4, separator_value=1):
    f = CNF()
    coeff = [f.new() for _ in generators]
    left = [[f.new() for _ in range(rank)] for _ in range(16)]
    right = [[f.new() for _ in range(16)] for _ in range(rank)]
    pivot = [[f.new() for _ in range(17)] for _ in range(rank)]
    cache = {0: f.false, **{1 << j: v for j, v in enumerate(coeff)}}
    def linear(form):
        if form not in cache:
            bit = form & -form
            cache[form] = f.xor(linear(bit), linear(form ^ bit))
        return cache[form]
    f.add(linear(ELL) if separator_value else -linear(ELL))
    for i in range(16):
        for j in range(16):
            form = sum(((g >> (16*i+j)) & 1) << k for k, g in enumerate(generators))
            value = f.false
            for k in range(rank):
                value = f.xor(value, f.conj(left[i][k], right[k][j]))
            f.equal(linear(form), value)
    for k in range(rank):
        f.add(*pivot[k])
        for i in range(17):
            for j in range(i):
                f.add(-pivot[k][i], -pivot[k][j])
            if i == 16:
                for j in range(16):
                    f.add(-pivot[k][i], -right[k][j])
            else:
                f.add(-pivot[k][i], right[k][i])
                for j in range(i):
                    f.add(-pivot[k][i], -right[k][j])
                for h in range(rank):
                    if h != k:
                        f.add(-pivot[k][i], -right[h][i])
        if k+1 < rank:
            for i in range(17):
                for j in range(i+1):
                    if i == j == 16:
                        continue
                    f.add(-pivot[k][i], -pivot[k+1][j])
    return f, coeff, left, right, pivot


def concrete_assignment(generators, coefficient, data):
    f, coeff, left, right, pivot = data
    vector = 0
    for j, g in enumerate(generators):
        if (coefficient >> j) & 1:
            vector ^= g
    rows = decode(vector, 4)
    reduced, pivots, rank = rows[:], [], 0
    for column in range(16):
        selected = next((i for i in range(rank, 16) if (reduced[i] >> column) & 1), None)
        if selected is None:
            continue
        reduced[rank], reduced[selected] = reduced[selected], reduced[rank]
        for i in range(16):
            if i != rank and (reduced[i] >> column) & 1:
                reduced[i] ^= reduced[rank]
        pivots.append(column)
        rank += 1
    assert rank <= len(right)
    values = [None]*(f.variables+1)
    for j, var in enumerate(coeff):
        values[var] = (coefficient >> j) & 1
    for k in range(len(right)):
        p = pivots[k] if k < rank else 16
        for j, var in enumerate(pivot[k]):
            values[var] = int(j == p)
        for j, var in enumerate(right[k]):
            values[var] = (reduced[k] >> j) & 1 if k < rank else 0
        for i in range(16):
            values[left[i][k]] = (rows[i] >> p) & 1 if k < rank else 0
    for var, operation, *inputs in f.definitions:
        if operation == 'false':
            values[var] = 0
        elif operation == 'xor':
            values[var] = values[inputs[0]] ^ values[inputs[1]]
        else:
            values[var] = values[inputs[0]] & values[inputs[1]]
    assert all(v is not None for v in values[1:])
    assert all(any(values[abs(lit)] == int(lit > 0) for lit in clause)
               for clause in f.clauses)
    return rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--regenerate-proof', action='store_true')
    args = parser.parse_args()
    generators = source_basis()
    for j, vector in enumerate(generators):
        assert sum((vector >> i) & 1 for i in SEPARATOR_ENTRIES)%2 == ((ELL >> j) & 1)
    coordinates = {}
    for j, vector in enumerate(generators):
        insert(coordinates, vector, 1 << j)
    target = represent(coordinates, encode(crossing(4), 4))
    assert target is not None and (target & ELL).bit_count()%2 == 1
    for permutation in CROSSING_SYMMETRIES:
        for j, vector in enumerate(generators):
            changed = encode(permute(decode(vector, 4), 4, permutation), 4)
            mask = represent(coordinates, changed)
            assert mask is not None and (mask & ELL).bit_count()%2 == ((ELL >> j) & 1)
    data = rank_formula(generators)
    formula = data[0]
    assert (formula.variables, len(formula.clauses)) == (2357, 9802)
    zero_separator = rank_formula(generators, separator_value=0)
    for coefficient in (0, 4, 0x80000, 0x80200, 0x904):
        concrete_assignment(generators, coefficient, zero_separator)
    assert concrete_assignment(generators, 0x1200, rank_formula(generators, rank=5)) == 5
    proof_path = ROOT/'evidence/hamming-joint-rank4.rup.gz'
    if args.regenerate_proof:
        from pysat.solvers import Glucose4
        with Glucose4(bootstrap_with=formula.clauses, with_proof=True) as solver:
            assert not solver.solve()
            proof = '\n'.join(line for line in solver.get_proof() if not line.startswith('d'))+'\n'
        with proof_path.open('wb') as out:
            with gzip.GzipFile(filename='', fileobj=out, mode='wb', mtime=0) as packed:
                packed.write(proof.encode())
    proof = gzip.decompress(proof_path.read_bytes())
    source = ROOT/'tools/check_rup.c'
    with tempfile.TemporaryDirectory(prefix='factoring-joint-lower-') as name:
        temporary = Path(name)
        binary, cnf = temporary/'check-rup', temporary/'rank.cnf'
        subprocess.run(['cc', '-D_FORTIFY_SOURCE=3', '-O3', '-std=c11', '-Wall', '-Wextra',
                        '-Werror', str(source), '-o', str(binary)], check=True)
        cnf.write_bytes(formula.dimacs())
        # Incomplete, out-of-range, and unproved-empty certificates must fail.
        for invalid in (b'0\n', b'1', f'{formula.variables+1} 0\n'.encode()):
            check = subprocess.run([str(binary), str(cnf)], input=invalid, capture_output=True)
            assert check.returncode != 0
        small = temporary/'small.cnf'
        small.write_text('p cnf 2 4\n1 2 0\n-1 2 0\n1 -2 0\n-1 -2 0\n')
        check = subprocess.run([str(binary), str(small)], input=b'1 0\n0\n',
                               capture_output=True, check=True)
        small.write_text('p cnf 1 1\n1 0\n')
        assert subprocess.run([str(binary), str(small)], input=b'0\n',
                              capture_output=True).returncode != 0
        checked = subprocess.run([str(binary), str(cnf)], input=proof,
                                 capture_output=True, check=True).stdout.decode().strip()
    additions = proof.count(b'\n')
    assert checked == (f'VERIFIED UNSAT additions {additions} variables '
                       f'{formula.variables} initial_clauses {len(formula.clauses)}')
    paths = [Path(__file__), source, proof_path]
    result = {'status': 'exact finite optimum 25 over F2 in the declared 33-dimensional joint span',
              'separator': hex(ELL), 'separator_tensor_entries': SEPARATOR_ENTRIES,
              'symmetry_checks': 8*33,
              'encoding_satisfying_controls': 6, 'proof_negative_controls': 4,
              'cnf_variables': formula.variables, 'cnf_clauses': len(formula.clauses),
              'cnf_sha256': sha256(formula.dimacs()).hexdigest(),
              'rup_additions': additions, 'rup_payload_sha256': sha256(proof).hexdigest(),
              'native_result': checked,
              'sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
