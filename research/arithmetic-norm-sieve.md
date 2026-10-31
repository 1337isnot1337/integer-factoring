# A deterministic norm sieve with actual factoring certificates

The arithmetic framework now has concrete search experiments beyond
class enumeration. The two certified inputs are

```text
1000000016000000063 = 1000000007 * 1000000009
906694351050364540708897 = 824633708473 * 1099511627689
```

The deterministic search constructs strictly diagonally dominant relation
rows on the prime-form factor base through 1000. The current word route
uses supplied smooth norms directly, polynomial ideal-lattice composition,
and a diagonal chosen from the proved word-length bound. The original
pair-form certificates remain reproducible. The second input has factors
separated by 274,877,919,216.

| Input bits | Relation correction | Rows | Diagonal K | Determinant bits | Two-adic valuation | Cyclic two-part orders |
| --- | --- | --- | --- | --- | --- | --- |
| 60 | Pair forms | 88 | 2^120 | 10,560 | 723 | 8 |
| 60 | Smooth word | 88 | 60 | 520 | 19 | 8 |
| 80 | Pair forms | 79 | 2^160 | 12,640 | 165 | 4, 2 |
| 80 | Smooth word | 79 | 80 | 500 | 16 | 4, 2 |

Determinant projection and the generic two-group algorithm produce these
bases and their ambiguous tops. For the first input the top is
(1000000007,0,1000000009). Both inputs yield proper factors.

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
For the original pair format, K=2^(2*bits(N)). For the direct word format,
let R be the strict radius in the norm interface, L=floor(log2(R^2)) and
K=L+1. For each i, reduce the target class g_i^K. First
test its leading coefficient as a possible smooth represented norm. Then
search Q_i(x,y), one row in y at a time. The old pair format uses
0<=x<=4096; the word format uses -4096<=x<=4096. Positive y and the
primitive leading-coefficient test cover the signed box up to simultaneous
negation, since Q(-x,-y)=Q(x,y). The height bounds are 64 for the 60-bit
certificates and 1024 for the separated-factor 80-bit certificates.

For each splitting prime p, the discriminant root gives all x residue
classes where Q_i(x,y) is divisible by p. Divide all powers of p from the
corresponding candidate residuals. A residual of one gives a smooth norm.
The [bounded-norm converter](bounded-smooth-norm.md) checks balanced
divisibility and returns two smooth forms in the target quotient class.
The word mode instead uses [direct prime-form word extraction](smooth-norm-relations.md),
accepting every eligible smooth norm without balanced divisibility.

Their prime-form words v_a,v_b give the relation

```math
r_i=K e_i-v_a+v_b.
```

For the word format the row is simply r_i=K e_i-v_n. Since the supplied
norm is at most R^2, ||v_n||_1<=L<K. The diagonal-dominance proof is
identical. The integration checker includes the existing form (24,1,50)
of discriminant -4799: its norm 24 gives a short word, while no pair at
radius 5 exists. Thus the word search does not silently retain balancing.

The correction norm is much smaller than K, and every row is checked for
strict diagonal dominance. Thus the square matrix has nonzero determinant
regardless of relations among the supplied generators. No random sample
or parity-distribution assumption is used in these successful finite runs.

All 88 targets in the 60-bit run found witnesses; the maximum y used
was 5 and the maximum x was 3713. Three targets used the leading-coefficient
test. Those observations describe this input, not a uniform search bound.
The 80-bit run found all 79 witnesses, with maximum y=359 and x=3855.
The word certificates have maximum y=5 and 247, respectively. Their target
powers and signed-coordinate ordering differ from the original certificates;
the failed boxes retain their own original target definitions.
The sieved row is computed in full before candidate inspection stops, so
an early witness does not mean only that many values were evaluated.

## Certificate and reproduction

The [60-bit certificate](../evidence/arithmetic-norm-sieve-60bit.json) and
[80-bit certificate](../evidence/arithmetic-norm-sieve-80bit.json) store
the input parameters, integer witness vectors, a relation-matrix
hash, the exact determinant in hexadecimal, and the extracted basis and
factor. The matrix and all prime forms are reconstructed. No timing data
or precomputed factor is used by the search.
The direct word certificates are stored separately for
[60 bits](../evidence/arithmetic-norm-word-sieve-60bit.json) and
[80 bits](../evidence/arithmetic-norm-word-sieve-80bit.json). The format tag
selects the exact row reconstruction used by the verifier.

Run the complete algebraic replay:

```sh
make verify-arithmetic-sieve
```

The checker reconstructs every target, norm witness and relation, independently checks
that each row evaluates to the identity, checks diagonal dominance,
recomputes the exact determinant, and repeats projection and full two-group
basis extraction on the factor-base subgroup. It checks the recovered
proper divisor. It does not simply trust the stored determinant or factor.
It replays all four certificates, checks the word-only control and compares
signed sieve residuals with direct division on independent small instances.
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
python3 tools/arithmetic_norm_sieve.py 1000000016000000063 --word-relations > /tmp/norm-word-sieve.json
python3 tools/arithmetic_norm_sieve.py 906694351050364540708897 --word-relations --y-bound 1024 > /tmp/norm-word-sieve-80.json
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

The main group operations now use the explicit polynomial
[ideal-lattice composer](factored-norm-search.md), retaining its source
attribution and independent controls. The legacy pair converter remains
available for exact replay. Default `make verify`
continues to check the small controls; this larger determinant/torsion
replay has a separate scoped target and its own CI job. Both targets passed
locally for this change.
Substantial OpenAI Codex assistance was used. No Lean source changed.
