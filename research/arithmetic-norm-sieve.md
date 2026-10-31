# A deterministic norm sieve with an actual factoring certificate

The arithmetic framework now has a concrete search experiment beyond
class enumeration. For the 60-bit integer

```text
1000000016000000063 = 1000000007 * 1000000009
```

the deterministic search constructs 88 strictly diagonally dominant relation
rows on the prime-form factor base through 1000. Its determinant has 10,560
bits and two-adic valuation 723. Determinant projection and the generic
two-group algorithm give one cyclic basis element of order 8; its top is
the ambiguous form (1000000007,0,1000000009), giving a proper factor.

This run does not enumerate the class group or find generator orders by
cycling. It constructs norm witnesses, verifies relations, computes the
determinant and extracts torsion. It also does not prove a new factoring
exponent: success within the search box is unproved uniformly, and the
factor-base subgroup is not assumed to be the full class group. The
strongest complete score remains 1/2. This is a finite certificate for the
new mechanism, not a practical factoring benchmark or a novelty claim over
existing class-group factoring methods.

## The deterministic search

Use Delta=-4N and the prime forms g_i at unramified splitting primes <=1000.
Choose K=2^(2*bits(N)). For each i, reduce the target class g_i^K. First
test its leading coefficient as a possible smooth represented norm. Then
search Q_i(x,y) for 0<=x<=4096 and 1<=y<=64, one row in y at a time.

For each splitting prime p, the discriminant root gives all x residue
classes where Q_i(x,y) is divisible by p. Divide all powers of p from the
corresponding candidate residuals. A residual of one gives a smooth norm.
The [bounded-norm converter](bounded-smooth-norm.md) checks balanced
divisibility and returns two smooth forms in the target quotient class.

Their prime-form words v_a,v_b give the relation

```math
r_i=K e_i-v_a+v_b.
```

The correction norm is much smaller than K, and every row is checked for
strict diagonal dominance. Thus the square matrix has nonzero determinant
regardless of relations among the supplied generators. No random sample
or parity-distribution assumption is used in this successful finite run.

All 88 targets found witnesses within the stated box; the maximum y used
was 5 and the maximum x was 3713. Three targets used the leading-coefficient
test. Those observations describe this input, not a uniform search bound.
The sieved row is computed in full before candidate inspection stops, so
an early witness does not mean only that many values were evaluated.

## Certificate and reproduction

The [certificate](../evidence/arithmetic-norm-sieve-60bit.json) stores only
the input parameters, the 88 integer witness vectors, a relation-matrix
hash, the exact determinant in hexadecimal, and the extracted basis and
factor. The matrix and all prime forms are reconstructed. No timing data
or precomputed factor is used by the search.

Run the complete algebraic replay:

```sh
make verify-arithmetic-sieve
```

The checker reconstructs every target and norm pair, independently checks
that each row evaluates to the identity, checks diagonal dominance,
recomputes the exact determinant, and repeats projection and full two-group
basis extraction on the factor-base subgroup. It checks the recovered
proper divisor. It does not simply trust the stored determinant or factor.

To repeat the deterministic witness search as well:

```sh
python3 tools/arithmetic_norm_sieve.py 1000000016000000063 > /tmp/norm-sieve.json
```

The [search implementation](../tools/arithmetic_norm_sieve.py) reports
progress on stderr and writes JSON on stdout. A search failure reports the
first exhausted target instead of claiming a relation certificate. A full
matrix can also fail to yield a split if its subgroup misses the needed
ambiguous directions. Both limits are part of the research contract.

The experiment uses the simple exact composition routine from the finite
controls, not an optimized class-group package. Default `make verify`
continues to check the small controls; this larger determinant/torsion
replay has a separate scoped target. Both passed locally for this change.
Substantial OpenAI Codex assistance was used. No Lean source changed.
