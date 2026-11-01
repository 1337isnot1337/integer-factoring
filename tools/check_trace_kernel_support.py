#!/usr/bin/env python3
"""Exact coefficient-support certificate for the first eligible trace-kernel size."""
import argparse
from array import array
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess
import sys
import tempfile
from time import monotonic

from check_trace_bilinear_kernel import Field, irreducibility_by_division, trace_polynomial


def read_coefficients(path, size):
    values = array('I')
    assert values.itemsize == 4
    values.frombytes(path.read_bytes())
    if sys.byteorder != 'little':
        values.byteswap()
    assert len(values) == size
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = monotonic()
    here = Path(__file__).resolve().parent
    source = here/'trace_kernel_support.cpp'
    with tempfile.TemporaryDirectory(prefix='trace-kernel-support-') as temp:
        temp = Path(temp)
        executable = temp/'support'
        subprocess.run(['g++', '-std=c++17', '-O2', str(source), '-o', str(executable)], check=True)

        def run(t, polynomial, k, name):
            path = temp/name
            completed = subprocess.run([str(executable), str(t), hex(polynomial), str(k), str(path)],
                                       check=True, capture_output=True, text=True)
            profiles = [json.loads(line) for line in completed.stdout.splitlines()]
            values = read_coefficients(path, (1 << t)-1)
            return profiles, values, hashlib.sha256(path.read_bytes()).hexdigest()

        small_profiles = 0
        small_coefficients = 0
        for t, polynomial in ((2, 7), (3, 11), (4, 19), (5, 37), (6, 67), (7, 131)):
            irreducibility_by_division(polynomial)
            field = Field(polynomial)
            k = min(t-1, 3)
            profiles, values, _ = run(t, polynomial, k, f'small-{t}.bin')
            for profile in profiles:
                expected = trace_polynomial(field, profile['k'])
                assert profile['rank'] == len(expected)
                small_profiles += 1
            expected = trace_polynomial(field, k)
            assert list(values) == [expected.get(e, 0) for e in range(field.order-1)]
            small_coefficients += len(values)

        t, k, polynomial = 19, 7, 0x80027
        divisors = irreducibility_by_division(polynomial)
        assert divisors == 511
        field = Field(polynomial)
        profiles, values, digest = run(t, polynomial, k, 'large.bin')
        size = (1 << t)-1
        assert all(value < field.order for value in values)
        for exponent, value in enumerate(values):
            assert bool(value) == (1 <= exponent.bit_count() <= k)
            if value:
                assert field.mul(value, value) == values[(2*exponent) % size]
        rank = sum(bool(value) for value in values)
        assert rank == sum(comb(t, j) for j in range(1, k+1)) == 94183
        for profile in profiles:
            kk = profile['k']
            assert profile['weights'] == [comb(t, j) if 1 <= j <= kk else 0 for j in range(t+1)]

        # Solve the independent trace equations to include two positive values,
        # rather than relying only on the much more frequent zero values.
        pivots = {}
        for j in range(k):
            a = field.power(field.alpha, j)
            mask = sum(field.trace(field.mul(a, 1 << i)) << i for i in range(t))
            rhs = 1
            while mask:
                pivot = mask.bit_length()-1
                if pivot in pivots:
                    other, value = pivots[pivot]
                    mask ^= other
                    rhs ^= value
                else:
                    pivots[pivot] = (mask, rhs)
                    break
            else:
                assert not rhs
        assert len(pivots) == k
        positive = 0
        free = next(i for i in range(t) if i not in pivots)
        kernel = 1 << free
        for pivot in sorted(pivots):
            mask, rhs = pivots[pivot]
            positive ^= (rhs ^ ((mask & positive).bit_count() & 1)) << pivot
            kernel ^= ((mask & kernel).bit_count() & 1) << pivot
        assert kernel
        # Independent field arithmetic evaluates the entire supplied polynomial
        # against the original product of binary traces at five field points.
        evaluations = []
        for z in (1, 2, 3, positive, positive ^ kernel):
            actual = 0
            for coefficient in reversed(values):
                actual = field.mul(actual, z) ^ coefficient
            expected = 1
            for j in range(k):
                expected &= field.trace(field.mul(field.power(field.alpha, j), z))
            assert actual == expected
            evaluations.append([z, actual])
        assert evaluations[-2][1] == evaluations[-1][1] == 1

        # The cyclic-carry support theorem gives this upper bound for every
        # member of the stated family, independently of its defining field.
        excluded_parameters = []
        for tt in range(2, 19):
            dd = 1 << tt
            for kk in range(1, tt):
                upper = sum(comb(tt, j) for j in range(1, kk+1))
                rr = (dd-1)*(1 << kk)
                assert rr*rr >= dd*upper*upper
                excluded_parameters.append([tt, kk])
        d = 1 << t
        cut_rank = (d-1)*(1 << k)
        assert cut_rank*cut_rank < d*rank*rank
        report = {
            'status': 'passed', 't': t, 'k': k, 'labels': d,
            'field_polynomial': polynomial, 'polynomial_hex': hex(polynomial),
            'irreducibility_divisors_checked': divisors,
            'ordinary_kernel_rank': rank, 'exact_row_cut_rank': cut_rank,
            'exact_column_cut_rank': cut_rank,
            'cut_rank_product': cut_rank*cut_rank, 'rejected_bound': d*rank*rank,
            'coefficient_slots_checked': len(values), 'coefficient_sha256': digest,
            'nonzero_coefficient_frobenius_checks': rank,
            'independent_small_rank_profiles': small_profiles,
            'independent_small_coefficients': small_coefficients,
            'independent_full_polynomial_evaluations': evaluations,
            'native_profiles': profiles,
            'smaller_family_parameter_pairs_excluded': len(excluded_parameters),
            'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (Path(__file__), source, here/'check_trace_bilinear_kernel.py')},
            'scope': 'Exact rank-inequality counterexample at t=19; minimal t only in the stated trace-product family. No crossing synthesis or improved factoring score.',
            'seconds': monotonic()-start,
        }
    rendered = json.dumps(report, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
