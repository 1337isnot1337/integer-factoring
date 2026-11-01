#!/usr/bin/env python3
"""Independent low-jet identities, cofactor ranks and graded word controls."""
from itertools import product,combinations
from random import Random
import json

def detpoly(A,P,prime):
 n=len(A);dp=[[0]*P for _ in range(1<<n)];dp[0][0]=1
 for mask in range((1<<n)-1):
  i=mask.bit_count()
  for j in range(n):
   if mask>>j&1:continue
   sign=-1 if (mask>>(j+1)).bit_count()%2 else 1;out=dp[mask|1<<j]
   for a,v in enumerate(dp[mask]):
    if v:
     for b,w in enumerate(A[i][j][:P-a]):out[a+b]=(out[a+b]+sign*v*w)%prime
 return dp[-1]
def determinant(A,p):return detpoly([[[v] for v in row] for row in A],1,p)[0]
def adjpoly(A,P,p):
 n=len(A);out=[]
 for i in range(n):
  line=[]
  for j in range(n):
   poly=detpoly([[A[a][b] for b in range(n) if b!=i] for a in range(n) if a!=j],P,p)
   line.append([(-v if(i+j)%2 else v)%p for v in poly])
  out.append(line)
 return out
def adj(A,p):return [[v[0] for v in row] for row in adjpoly([[[a] for a in row] for row in A],1,p)]
def tr(A):return list(map(list,zip(*A)))
def mm(A,B,p):return [[sum(a*b for a,b in zip(row,col))%p for col in zip(*B)] for row in A]
def channel(A,B):
 q=len(A);return [[A[x][u]*B[v][y] for v,u in product(range(q),repeat=2)] for x,y in product(range(q),repeat=2)]
def basis(rng,q,p):
 A=[[int(i==j) for j in range(q)] for i in range(q)]
 for _ in range(5*q):
  i,j=rng.sample(range(q),2);c=rng.randrange(1,p);A[i]=[(a+c*b)%p for a,b in zip(A[i],A[j])]
 return A
import argparse
from fractions import Fraction as F
from hashlib import sha256
import lzma
from pathlib import Path
from check_cell_minor_caps import PRIME, literal_case, polynomial_minor, random_crossing, smith_orders
from higher_jet_word_caps import inequalities, word_cost_profile
from check_q12_second_jet_target import certificate as old_target_certificate
ROOT=Path(__file__).resolve().parents[1]


def homogeneous_controls():
    rng, records, count = Random(0xF0F7), [], 0
    for prime, q, policy in product((2, 10007), (2, 3), (0, 1, 2)):
        n = q*q
        A, B = basis(rng, q, prime), basis(rng, q, prime)
        if policy:
            A[-1] = [0]*q
        if policy == 2:
            B[-1] = [0]*q
        R = channel(A, B)
        S = channel(tr(adj(B, prime)), tr(adj(A, prime)))
        D = determinant(A, prime)*determinant(B, prime) % prime
        E = [[rng.randrange(prime) for _ in range(n)] for _ in range(n)]
        actual = adjpoly([[[R[i][j] % prime, E[i][j]] for j in range(n)] for i in range(n)], q, prime)
        W = mm(S, E, prime)
        powers = [[[int(i == j) for j in range(n)] for i in range(n)]]
        for _ in range(1, q):
            powers.append(mm(powers[-1], W, prime))
        elementary = [1]+[sum(determinant([[W[i][j] for j in ids] for i in ids], prime)
                             for ids in combinations(range(n), a)) % prime for a in range(1, q)]
        epsilon = (-1)**(q*(q-1)//2)
        for ell in range(q):
            expected = [[0]*n for _ in range(n)]
            for b in range(ell+1):
                coefficient = epsilon*pow(D, q-1-ell, prime)*elementary[ell-b]*(-1)**b
                word = mm(powers[b], S, prime)
                for i, j in product(range(n), repeat=2):
                    expected[i][j] = (expected[i][j]+coefficient*word[i][j]) % prime
            assert expected == [[actual[i][j][ell] for j in range(n)] for i in range(n)]
            count += n*n
        records.append(dict(prime=prime, alphabet=q, singularity_policy=policy, A=A, B=B,
                            E=E, low_adjugate=actual))
    return dict(coefficient_identities=count, cases=records)


def rank_mod(A):
    rows, rank = [row[:] for row in A], 0
    for j in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][j] % PRIME), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][j], PRIME-2, PRIME)
        rows[rank] = [v*inverse % PRIME for v in rows[rank]]
        for i in range(rank+1, len(rows)):
            scalar = rows[i][j]
            if scalar:
                rows[i] = [(a-scalar*b) % PRIME for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def toeplitz_profile(A, cap):
    n, ranks = len(A), [0]
    for P in range(1, cap+1):
        block = [[(A[i][j][a-b] if 0 <= a-b < len(A[i][j]) else 0)
                  for b in range(P) for j in range(n)] for a in range(P) for i in range(n)]
        ranks.append(rank_mod(block))
    below = [0]+[b-a for a, b in zip(ranks, ranks[1:])]
    orders = []
    for k in range(cap):
        orders.extend([k]*(below[k+1]-below[k]))
    return orders+[cap]*(n-len(orders)), ranks


def extend(vector):
    n, anchor = len(vector), next(i for i, v in enumerate(vector) if v)
    return [[vector[i]]+[int(i == j) for j in range(n) if j != anchor] for i in range(n)]


def cofactor_controls():
    rng, records, count = Random(0xBA6E7), [], 0
    cases = [(2, [0, 1, 2, 4], False), (3, [0, 1]+[3]*7, False),
             (3, [0, 1, 2]+[4]*6, False), (3, [0, 1]+[3]*7, True)]
    for q, orders, singular in cases:
        n = q*q
        U, V = basis(rng, n, PRIME), basis(rng, n, PRIME)
        if singular:
            a = [int(i == j and i < q-1) for i, j in product(range(q), repeat=2)]
            b = [int(i == j) for i, j in product(range(q), repeat=2)]
            U, V = extend(a), tr(extend(b))
        source = [[[sum(U[i][k]*V[k][j] for k in range(n) if orders[k] == e) % PRIME
                    for e in range(max(orders)+1)] for j in range(n)] for i in range(n)]
        assert smith_orders(source, n*max(orders)+1) == orders
        row = [[source[x*q+z][y*q+v] for v, z in product(range(q), repeat=2)]
               for x, y in product(range(q), repeat=2)]
        entries = {i: dict(enumerate(line)) for i, line in enumerate(row)}
        cofactors = [[polynomial_minor(entries, [a for a in range(n) if a != j],
                                      [a for a in range(n) if a != i]) for j in range(n)] for i in range(n)]
        for i, j in product(range(n), repeat=2):
            if (i+j) % 2:
                cofactors[i][j] = [-v % PRIME for v in cofactors[i][j]]
        crossed = [[cofactors[x*q+y][v*q+z] for y, v in product(range(q), repeat=2)]
                   for x, z in product(range(q), repeat=2)]
        actual, ranks = toeplitz_profile(crossed, q*orders[1])
        lower = word_cost_profile(q, orders)
        assert all(a >= b for a, b in zip(actual, lower))
        if not singular:
            assert tuple(actual) == lower
        records.append(dict(alphabet=q, singular_leading_factor=singular, crossed_source=source,
                            source_orders=orders, adjugate=cofactors, toeplitz_ranks=ranks,
                            truncated_crossed_orders=actual, word_lower=list(lower)))
        count += n*n
    return dict(interpolated_cofactors=count, cases=records)


def combinatorial_controls():
    rng, heap_checks, shifts = Random(0xA0D17), 0, 0
    for q in (2, 3, 4):
        n = q*q
        for trial in range(12):
            source = [F(0)]+sorted(F(rng.randrange(1, 15), 4) for _ in range(n-1))
            bound = word_cost_profile(q, source)
            if trial < 2:
                costs = [F(0)]+[sum(word) for length in range(1, q) for word in product(source[1:], repeat=length)]
                expected = tuple(sorted(min(v, q*source[1]) for v in costs)[:n])
                assert bound == expected
                heap_checks += 1
            b = F(rng.randrange(10), 3)
            raw, local = word_cost_profile(q, [v+b for v in source]), bound
            for k in range(1, n+1):
                delta = sum(raw[:k])-sum(local[:k])
                assert 0 <= delta <= (k-1)*q*b
                assert k*(n-1)*b-delta >= (k*(n-1)-(k-1)*q)*b >= 0
                shifts += 1
    return dict(independent_word_enumerations=heap_checks, local_shift_prefixes=shifts)


def graded_controls():
    rng, count, records = Random(0xA0DE2), 0, []
    for q in (2, 3, 4):
        for _ in range(3):
            tensor = random_crossing(rng, q, 2, .65)
            heights = [[[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2)]
            data = literal_case(rng, q, tensor, heights, cell_samples=0)
            for pair in data['profiles'].values():
                caps = inequalities(data['deficit'], pair['r'], pair['s'], data['paired_energies_X'], data['paired_energies_Y'])
                assert min(cap['slack'] for cap in caps) >= 0
                count += len(caps)
            records.append(data)
    return dict(tensors=9, complete_profiles=36, inequalities=count,
                records_sha256=sha256(json.dumps(records, sort_keys=True, default=str).encode()).hexdigest())


def q12_rejection():
    previous = old_target_certificate()
    path = ROOT/'evidence/q12-second-jet-target.json'
    point = json.loads(path.read_text())
    r, d, energies = point['orders'], point['deficit'], point['paired_energies']
    rejected = [cap for cap in inequalities(d, r, r, energies, energies) if cap['slack'] < 0]
    assert len(rejected) == 258
    worst = min(rejected, key=lambda cap: cap['slack'])
    assert (worst['physical_size'], worst['size'], worst['slack']) == (12, 144, -244667)
    return dict(previous_relaxation_checked=previous['checks'], previous_pressure_upper_bound=previous['pressure_upper_bound'],
                oriented_word_violations=len(rejected), worst=json.loads(json.dumps(worst, default=str)),
                target_sha256=sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    names = ['tools/higher_jet_word_caps.py', 'tools/check_higher_jet_word_caps.py',
             'tools/hereditary_adjugate_caps.py', 'tools/check_q12_second_jet_target.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py',
             'tools/inherited_second_jet_caps.py', 'tools/multichannel_cofactor_caps.py',
             'tools/rectangular_adjugate_jet_caps.py']
    result = dict(scope='All-field higher-jet word lower profile and graded square restrictions; necessary constraints, no stronger factoring algorithm.',
                  homogeneous=homogeneous_controls(), cofactor_profiles=cofactor_controls(),
                  combinatorics=combinatorial_controls(), graded=graded_controls(), q12=q12_rejection(),
                  source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    if args.output:
        args.output.write_bytes(lzma.compress((json.dumps(result, sort_keys=True, default=str)+'\n').encode()))
    else:
        saved = json.loads(lzma.decompress((ROOT/'evidence/higher-jet-word-controls.json.xz').read_bytes()))
        assert json.loads(json.dumps(result, default=str)) == saved
    print(json.dumps(dict(homogeneous_identities=result['homogeneous']['coefficient_identities'],
                         interpolated_cofactors=result['cofactor_profiles']['interpolated_cofactors'],
                         combinatorics=result['combinatorics'], graded=result['graded'], q12=result['q12'])))


if __name__ == '__main__':
    main()
