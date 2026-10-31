# A deterministic norm sieve with actual factoring certificates

The arithmetic framework now has concrete search experiments beyond
class enumeration. The two certified inputs are

```text
1000000016000000063 = 1000000007 * 1000000009
906694351050364540708897 = 824633708473 * 1099511627689
```

The deterministic search constructs strictly diagonally dominant relation
rows on the prime-form factor base through 1000. The second input has
factors separated by 274,877,919,216.

| Input bits | Relation rows | Determinant bits | Two-adic valuation | Cyclic two-part orders | Recovered factor |
| --- | --- | --- | --- | --- | --- |
| 60 | 88 | 10,560 | 723 | 8 | 1000000007 |
| 80 | 79 | 12,640 | 165 | 4, 2 | 824633708473 |

Determinant projection and the generic two-group algorithm produce these
bases and their ambiguous tops. For the first input the top is
(1000000007,0,1000000009). Both runs yield proper factors.

These runs do not enumerate the class group or find generator orders by
cycling. They construct norm witnesses, verify relations, compute the
determinant and extract torsion. They do not prove a new factoring
exponent: success within the search box is unproved uniformly, and the
factor-base subgroup is not assumed to be the full class group. The
strongest complete score remains 1/2. These are finite certificates for the
new mechanism, not a practical factoring benchmark or a novelty claim over
existing class-group factoring methods.

## The deterministic search

Use Delta=-4N and the prime forms g_i at unramified splitting primes <=1000.
Choose K=2^(2*bits(N)). For each i, reduce the target class g_i^K. First
test its leading coefficient as a possible smooth represented norm. Then
search Q_i(x,y) for 0<=x<=4096, one row in y at a time. The height bounds
are 64 for the stored 60-bit run and 1024 for the separated-factor 80-bit run.

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
or parity-distribution assumption is used in these successful finite runs.

All 88 targets in the 60-bit run found witnesses; the maximum y used
was 5 and the maximum x was 3713. Three targets used the leading-coefficient
test. Those observations describe this input, not a uniform search bound.
The 80-bit run found all 79 witnesses, with maximum y=359 and x=3855.
The sieved row is computed in full before candidate inspection stops, so
an early witness does not mean only that many values were evaluated.

## Certificate and reproduction

The [60-bit certificate](../evidence/arithmetic-norm-sieve-60bit.json) and
[80-bit certificate](../evidence/arithmetic-norm-sieve-80bit.json) store
the input parameters, integer witness vectors, a relation-matrix
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
The checker also replays two [fixed-box exhaustion controls](../evidence/arithmetic-norm-search-exhaustion.json).
The nearby-factor 80-bit input 1208925818180866012192151 exhausts height
64 at generator prime 197. The separated-factor input exhausts height 256
at generator prime 127. Thus the smaller boxes do not suffice even for
these finite controls. The larger-box certificates establish no uniform
height bound.

To repeat the deterministic witness search as well:

```sh
python3 tools/arithmetic_norm_sieve.py 1000000016000000063 > /tmp/norm-sieve.json
python3 tools/arithmetic_norm_sieve.py 906694351050364540708897 --y-bound 1024 > /tmp/norm-sieve-80.json
```

The [search implementation](../tools/arithmetic_norm_sieve.py) reports
progress on stderr and writes JSON on stdout. A search failure reports the
first exhausted target instead of claiming a relation certificate. A full
matrix can also fail to yield a split if its subgroup misses the needed
ambiguous directions. Both limits are part of the research contract.

## Related methods and the uniformity gap

Sieving quadratic-form values is standard in class-group methods and
quadratic-sieve variants. Bradford--Wagstaff's
[SQUFOF2](https://homes.cerias.purdue.edu/~ssw/squfof2.pdf), Section 8.2,
explicitly bases its complexity analysis on heuristic hypotheses, including
smoothness in its sieve rectangle. We do not import those hypotheses as
theorems. The certificates here check dominant relation matrices and full
two-part extraction on their actual factor-base subgroups.

[Balog--Blomer--Dartyge--Tenenbaum](https://tenenb.perso.math.cnrs.fr/PPP/friableforms.pdf),
Theorem 1, proves friable-value results with constants depending on the
fixed binary form. That is insufficient by itself for a short-box bound
uniform in this growing discriminant. The ideal-in-class sampling theorem
of [de Boer--Pellet-Mary--Wesolowski](https://koendeboer.com/publication/rigorousmethods/RigorousMethods.pdf),
Theorem 1.1, is randomized and assumes ERH. Neither statement supplies the
deterministic uniform witness-search guarantee left open here.

The experiment uses the simple exact composition routine from the finite
controls, not an optimized class-group package. Default `make verify`
continues to check the small controls; this larger determinant/torsion
replay has a separate scoped target and its own CI job. Both targets passed
locally for this change.
Substantial OpenAI Codex assistance was used. No Lean source changed.
