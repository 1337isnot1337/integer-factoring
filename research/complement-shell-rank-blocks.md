# Exact rank blocks and protected minors for complementary shells

## Result

The complementary-shell construction has a conserved binary label that splits
**both complete tensor flattenings** into exact blocks. This substantially
reduces the cost of measuring a full crossing decomposition. For example, the
496-label shell below has 60,523,872,256 potential tensor entries, but only
5,127,880 lie in permitted blocks.

The same local tensor gives explicit identity minors in its central
coefficient. These prove an all-size limitation: for every fixed length
\(L\geq3\), the full central shells fail the finite seed criterion once the
mode cutoff is sufficiently large. This statement does not cover lengths
growing with the cutoff, selected code families, or modified local tensors.

The new exact profiles reach dimensions 86, 120, 190, and 496. They give no
improved factoring exponent.

## 1. Setup and exact finite tensor

Use the notation and proof of
[complementary fixed-length shells](complement-shell-crossings.md).
Let \(M\geq4\) be a power of two. A label is an \(L\)-element subset of the
odd modes below \(M\), with common weight \(j\leq LM/2\).
Parameter and state label sets can be different subsets of the same shell.
They must have the same cardinality \(d\).

The partner of parameter label \(v\) is the mode complement

\[
 \bar v=\{M-r:r\in v\}.
\]

Its weight is \(h=LM-j\). The raw vacuum pairing is the identity. The state
energy is \(j\).

Set \(t=z^M\). The exact complete crossing polynomial, multiplied by \(t^L\),
is a product over the \(M/4\) complementary pairs \(\{r,M-r\}\).
At one such pair, a local occupation is one of the four masks
\(0,1,2,3\). For local labels \(x,y,u,v\), take the Pfaffian on the occupied
vertices of four blocks \(U,A,V,Y\). Its only edges are:

| Blocks | Mode relation | Edge polynomial |
|---|---|---|
| \(U,A\) | same | \(t+t^2\) |
| \(A,Y\) | complementary | \(1+t\) |
| \(A,V\) | same | \(t\) |
| \(U,V\) | complementary | \(1\) |
| \(V,Y\) | same | \(1\) |
| \(U,Y\) | same | \(1\) |

All arithmetic is over \(\mathbb F_2\). The complete product has degree at
most \(2L\). Its degree-\(k\) coefficient is the original crossing piece
at shift \((k-L)M\). No high intermediate modes are removed: their exact
geometric sums are already incorporated in this local tensor.

The two rank cuts, in the repository's indexing convention, are

\[
 (x,y)\mid(v,u),\qquad (x,u)\mid(v,y).
\]

## 2. Conserved parity of single occupation

For a label \(a\), define its syndrome

\[
 \sigma(a)_r=
 \bigl|a\cap\{r,M-r\}\bigr|\pmod2
 \quad\in\mathbb F_2^{M/4}.
\]

Thus a complementary pair contributes one precisely when it has one occupied
mode. A local Pfaffian with an odd number of vertices is zero. Every nonzero
local tensor entry therefore satisfies

\[
 \sigma(x)_r+\sigma(y)_r+\sigma(u)_r+\sigma(v)_r=0.
\]

This is an identity of complete polynomials, and hence holds coefficient by
coefficient.

**Block theorem.** For every coefficient, the first cut is a direct sum over
syndromes \(s\), with row and column labels satisfying

\[
 \sigma(x)+\sigma(y)=s=\sigma(v)+\sigma(u).
\]

The second cut has the same block sizes, using

\[
 \sigma(x)+\sigma(u)=s=\sigma(v)+\sigma(y).
\]

Both statements remain valid when the parameter and state domains differ.
The rank of each complete cut is the sum of its block ranks. There is no
inequality or relaxation in this reduction.

If \(n_s\) is the number of parameter-state label pairs with syndrome \(s\),
then the full tensor can be assembled with

\[
 \sum_s n_s^2
\]

entry evaluations, instead of \(d^4\). A single evaluation updates both cuts.
The checker below uses exact binary elimination on every block.

The syndrome is coarser than the XOR of the individual mode indicators.
The latter is not conserved by the shell tensor; complementary-mode edges
give explicit counterexamples in the earlier shell evidence.

## 3. Protected central identity minors

Call a label *simple* if it contains at most one mode from each complementary
pair. For two simple labels \(a,b\), call the ordered pair *disjoint* when
their complementary-pair supports are disjoint.

Let \(\mathcal D\) be all such ordered pairs, with \(a\) in the parameter
domain and \(b\) in the state domain. Define two sets of ordinary mode unions:

\[
 \mathcal U_R=\{a\cup b:(a,b)\in\mathcal D\},\qquad
 \mathcal U_C=\{a\cup\bar b:(a,b)\in\mathcal D\}.
\]

Every union has \(2L\) distinct modes and at most one mode from any
complementary pair.

**Central-minor theorem.** If \(r_0,c_0\) are the two ranks of the original
shift-zero coefficient, then

\[
 r_0\geq|\mathcal U_R|,\qquad
 c_0\geq|\mathcal U_C|.
\]

In fact, each bound is realized by an explicit identity submatrix.

### First cut

Restrict rows \((x,y)\) and columns \((v,u)\) to \(\mathcal D\).
Each occupied complementary pair has one vertex on each side of the cut.
A nonzero local edge requires the same actual mode on both sides. Consequently
the entry is zero unless

\[
 x\cup y=v\cup u.
\]

If the unions agree, let \(k\) be the number of modes sent from \(A\) to \(U\).
All \(L\) vertices of \(A\) supply a factor \(t\); exactly those \(k\)
vertices supply an additional factor \(1+t\). Vertices of \(Y\) supply \(1\).
The complete entry is therefore

\[
 t^L(1+t)^k.
\]

Its degree-\(L\) coefficient is one. The restricted central matrix is thus a
direct sum of all-ones matrices, one for each union in \(\mathcal U_R\).
Choose one pair representing each union. The resulting square submatrix is
the identity.

### Second cut

Restrict rows \((x,u)\) and columns \((v,y)\) to \(\mathcal D\), and complement
the state modes when recording their union. The local edges are nonzero
precisely when

\[
 x\cup\bar u=v\cup\bar y.
\]

Suppose \(k\) vertices of \(A\) are sent to \(Y\). The complete entry is

\[
 t^{L-k}(1+t)^k.
\]

Its leading coefficient, at degree \(L\), is one. Again the central matrix
on the selected rows and columns is a direct sum of all-ones matrices.
One representative of each union gives the second identity minor.

### Consequence for the complete seed score

The criterion used in this repository is

\[
 S=\sum_s\sqrt{r_sc_s}<d^{\sqrt2}.
\]

The central coefficient alone gives the exact lower bound

\[
 S\geq
 \sqrt{|\mathcal U_R|\,|\mathcal U_C|}.
\]

Moreover, each union has at most \(\binom{2L}{L}\) decompositions into its two
label subsets. If \(D=|\mathcal D|\), then

\[
 |\mathcal U_R|,\ |\mathcal U_C|
 \geq \frac{D}{\binom{2L}{L}},
 \qquad
 S\geq\frac{D}{\binom{2L}{L}}.
\]

These bounds are valid for arbitrary equally sized parameter and state
subsets of the shell. They can reject a proposed code before computing the
whole tensor.

## 4. An all-size fixed-length limitation

Fix \(L\geq3\). Take the full shell at central admissible weight

\[
 j_M=
 \begin{cases}
 LM/2,&L\text{ even},\\
 LM/2-1,&L\text{ odd}.
 \end{cases}
\]

Let its size be \(d_M\). As \(M\) tends through powers of two,

\[
 d_M=\Theta_L(M^{L-1}).
\]

Here is an elementary count. An ordered choice of the first \(L-1\) modes
determines the last, which proves the upper bound. For the lower bound,
choose those modes in a fixed small interval about \(M/2\). Make the interval
short enough, depending only on \(L\), that the determined last mode remains
strictly between zero and \(M\). The parity of \(j_M\) makes it odd.
Repeated modes exclude only \(O_L(M^{L-2})\) ordered choices. Dividing by at
most \(L!\) gives the lower bound.

There are only \(O_L(M^{L-2})\) nonsimple labels. To count them, choose a
complementary pair in \(O(M)\) ways and then count the \(L-2\) remaining modes
of prescribed sum. Since \(L-2\geq1\), this count is
\(O_L(M^{L-3})\).

For a fixed simple label, a second label that meets one of its complementary
pairs has one of at most \(2L\) specified modes. The remaining \(L-1\) modes
have a prescribed sum, giving \(O_L(M^{L-2})\) possibilities. Therefore

\[
 D_M=d_M^2-O_L(M^{2L-3})
     =(1-o(1))d_M^2.
\]

The protected-minor theorem now yields

\[
 S_M\geq
 \frac{1-o(1)}{\binom{2L}{L}}d_M^2.
\]

Since \(2>\sqrt2\), this exceeds \(d_M^{\sqrt2}\) for all sufficiently large
\(M\). Increasing only the cutoff at fixed \(L\) cannot improve the
factoring exponent through these full-shell decompositions.

The fixed-length hypothesis matters. The constant
\(\binom{2L}{L}\) grows with \(L\), and the disjoint-pair count can collapse
for dense or deliberately overlapping label families.

## 5. Complete finite profiles

The companion evidence uses the original local tensor and computes every
nonzero coefficient in both cuts:

| \(M,L,j\) | Labels \(d\) | Central ranks \((r_0,c_0)\) | Complete \(S\) | \(d^{\sqrt2}\) |
|---|---:|---:|---:|---:|
| \(16,4,32\) | 8 | \(29,29\) | 87.983867 | 18.930501 |
| \(32,3,47\) | 28 | \(194,392\) | 835.960654 | 111.324439 |
| \(32,4,64\) | 86 | \(2795,2795\) | 5745.474886 | 544.244371 |
| \(64,3,95\) | 120 | \(4408,9510\) | 12394.889300 | 871.778826 |
| \(128,3,191\) | 496 | \(93584,198336\) | 198624.513187 | 6486.146356 |
| \(32,5,79\) | 190 | \(8550,9424\) | 39465.065913 | 1669.722301 |

For the 120-label shell, the protected identity minors have sizes 2719 and
4566. They already imply \(S\geq3523.486058\), exceeding the threshold.
For the 496-label shell their sizes are 78565 and 145090, giving
\(S\geq106766.080054\).

The block reduction is particularly effective for low fixed length:

| \(M,L\) | Syndrome blocks | Largest block | Evaluated entries | Dense entries |
|---|---:|---:|---:|---:|
| \(32,3\) | 103 | 70 | 13792 | 614656 |
| \(32,4\) | 121 | 900 | 1534992 | 54700816 |
| \(64,3\) | 3216 | 330 | 301272 | 207360000 |
| \(128,3\) | 81583 | 1426 | 5127880 | 60523872256 |
| \(32,5\) | 128 | 1058 | 12569432 | 1303210000 |

The evidence records full profiles and per-degree distributions of block
ranks. It is not a sample of random restrictions.

## 6. Reproduction and reusable interface

Files:

- [Exact blocked checker](../tools/check_complement_shell_ranks.py).
- [Complete evidence](../evidence/complement-shell-ranks.json).
- [Shared local tensor](../tools/complement_shell_tensor.py).

Run:

    python3 tools/check_complement_shell_ranks.py \
      --cases 16:4:32,32:3:47,32:4:64,64:3:95,128:3:191,32:5:79 \
      --output evidence/complement-shell-ranks.json

The checker verifies all 256 local tensor entries for the parity symmetry and
all 32 local identities used by the protected minors. It independently
compares dense and blocked ranks for the eight-label case. The earlier shell
checker compares the local tensor with the full homogeneous Pfaffian, so the
two rank implementations do not supply their own unverified crossing input.

For selected domains, import the function

    rank_labels(modulus, parameters, states=None, compare_dense=False)

Each domain is a sequence of distinct tuples of modes. Both domains must have
the same positive cardinality, length, and energy. Omitting the state domain
uses the parameter domain. The function returns both complete rank profiles,
the correct rank sum \(S=\sum_s\sqrt{r_sc_s}\), the threshold, syndrome-block
statistics, and the protected central-minor sizes.

The cheaper function

    protected_central_minors(modulus, parameters, states=None)

returns the proven central lower bound without constructing the coefficient
matrices.
