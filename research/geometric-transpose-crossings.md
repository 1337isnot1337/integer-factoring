# Geometric transpose channels and their two-channel inverses

## Status

This is a new structured crossing construction and a finite exact checkpoint.
It does **not** improve the factoring exponent or the speedrun score.

For label sizes `q = 2, 4, 8, 16`, the construction below gives a nilpotent
superoperator of index `2q`, with `q/2` Jordan blocks of that size. Its associated
crossing polynomial is unimodular in the row cut. Inverting the corresponding
unipotent superoperator reduces its reshuffled rank from `q` to exactly `2`,
while retaining the same Jordan type and row-unimodularity.

The formulas defining both channels are valid for every size. Nilpotence with
the stated Jordan type, and row-unimodularity, are certified here only for the
four listed sizes. A general proof for all binary sizes remains open in this
checkpoint.

An exhaustive normalized polynomial-change screen adds two precise conclusions:

- For the `q=8` geometric channel, inversion is the unique one of `16,384`
  changes whose leading reshuffled rank is at most `5`.
- For the earlier regular `q=4` two-channel control, the original channel is
  the unique one of `16,384` changes whose leading reshuffled rank is at most
  `3`.

These are exhaustive statements about the stated finite polynomial spaces,
not about all crossing constructions.

A protected subleading-coefficient argument below also excludes the remaining
candidate in each space. Consequently neither complete polynomial space can
improve the target by ordinary tensor-power coefficient decompositions. This
conclusion uses actual rank lower bounds and does not depend on the sufficient
Smith bound.

## 1. Construction

Work over a field of characteristic two. Let `N` be the `q` by `q` lower
Jordan shift, and put

\[
 P=(I+N^T)^{-1},\qquad C=I+N^2.
\]

For left and right multiplication operators write
`L_A(Y)=AY`, `R_A(Y)=YA`, and define

\[
 \tau_P(Y)=PY^TP^{-1},\qquad
 U=(I+L_NR_N)^{-1}R_C\tau_P.
\]

Since `N` is nilpotent, every inverse in this formula is a finite sum. More
explicitly,

\[
 U(Y)=\sum_{i=0}^{q-1}N^iPY^TP^{-1}N^i(I+N^2).                 \tag{1}
\]

This is an invertible linear map on `Mat_q`. It fixes the identity and preserves
the trace. For the first assertion, set `Y=I` in (1) and telescope

\[
 \sum_iN^{2i}(I+N^2)=I.
\]

For the second, use cyclicity of the trace to move both powers of `N` and
`I+N²` together, and use the same identity. These arguments do not assume the
finite nilpotence checks stated later. Also, the two lists of matrices
`N^i P` and `P^{-1} N^i C`, for `0<=i<q`, are linearly independent.
Under reshuffling, (1) is a sum of outer products of these two lists.
Thus `rank Phi(U)=q` for every `q`; this rank identity is an all-size
algebraic fact.

Set `S=I+U`. Whenever `S` is nilpotent of index `m`, apply the general
nilpotent-superoperator construction from `nilpotent-flag-crossings.md`:

\[
 K_S(z)=\Phi\left(\sum_{k=0}^{m-1}
             z^{m-1-k}S^k(I+S)\right),                       \tag{2}
\]

where

\[
 \Phi(M)_{(x,y),(u,v)}=M(E_{vy})_{x,u}.
\]

Then `K_S(1)` is the crossing tensor. If the Jordan blocks of `S` have lengths
`b_j`, the column-zero Smith orders are one `m-b_j` per block and `m` repeated
`b_j-1` times per block.

## 2. Inversion gives an exact rank-two channel

The inverse of (1) has a particularly short formula:

\[
\begin{aligned}
 U^{-1}(Y)
 ={}&P^TC^{-T}Y^TP^{-T}\\
   &+P^TC^{-T}N^TY^TN^TP^{-T}.                                \tag{3}
\end{aligned}
\]

This follows by reversing the three factors defining `U`:

\[
 U^{-1}=\tau_P^{-1}R_{C^{-1}}(I+L_NR_N).
\]

Each map `Y -> AY^TB` becomes an outer product under `Φ`. The two left factors
in (3) are linearly independent, as are the two right factors: after cancelling
invertible factors, this is the independence of `I` and the nonzero nilpotent
matrix `N^T`. Consequently

\[
 \operatorname{rank}\Phi(U^{-1})=2\qquad(q\ge2).              \tag{4}
\]

Equation (4) is a general theorem about this formula, not merely a finite
observation.

If `S` is nilpotent, the new nilpotent associated to `U^{-1}` is

\[
 S'=I+U^{-1}=\frac{S}{1+S}=S+S^2+\cdots+S^{m-1}.              \tag{5}
\]

Multiplication by the unit `(1+S)^{-1}` preserves the kernels of every power,
so `S'` has exactly the same Jordan type as `S`.

### A useful transfer identity

Suppose `m` is a power of two. The operator polynomial inside the
reshuffling has the closed expression

\[
 F_{S'}(z)=(z+(1+z)S)^{m-1}.                  \tag{6}
\]

The literal column flattening is this operator matrix followed by the
fixed transpose permutation. Equation (6) refers to the operator matrix.

Indeed, `U^m=I`, and the finite geometric sum in (2), after using (5), becomes
the characteristic-two identity

\[
 \sum_{j=0}^{m-1}a^jb^{m-1-j}=(a+b)^{m-1}.
\]

The row constant coefficient is the same for both curves, namely
`Φ(S^(m-1))`. In the algebra spanned by `I,S,...,S^(m-1)`, the positive
coefficients of **both** curves span

\[
 \{p(S):\deg p<m,\ p(1)=0\}.                                \tag{7}
\]

For the original curve this follows from its coefficients
`(1+S)S^j`, `0 <= j <= m-2`. For (6), substitution `S=1` gives the constant
polynomial `1`; hence every positive coefficient lies in the subspace (7).
The coefficient change is invertible, so these `m-1` coefficients span it.

Thus a common nilpotent-algebra certificate for the positive row coefficients
of one curve transfers to the other. This explains the identical dimension
chains found below; it is not an accidental numerical equality.

## 3. Exact finite certificates

For each `q=2,4,8,16`, the checker verifies

\[
 \operatorname{rank}S^k=\frac q2(2q-k)
 \quad(0\le k\le2q).                                        \tag{8}
\]

This proves nilpotence, its exact index, and the stated Jordan type. The
column-zero orders are therefore

\[
 0^{[q/2]},\quad(2q)^{[q^2-q/2]}.                            \tag{9}
\]

The checker constructs the literal polynomial (2), verifies every entry of
the crossing endpoint, and verifies that its constant row coefficient `K_0`
is invertible. For the normalized positive coefficients

\[
 A_i=K_0^{-1}K_i\quad(i>0),
\]

it computes the descending sequence

\[
 W_0=\mathbb F_2^{q^2},\qquad W_{j+1}=\sum_{i>0}A_iW_j.
\]

Every sequence reaches zero. Hence all sufficiently long products of the
`A_i` vanish. It follows that

\[
 I+\sum_{i>0}z^iA_i
\]

has a polynomial inverse, and its determinant is the constant `1`. This is an
exact unimodularity certificate; it does not infer a determinant from a
truncated Smith computation.

| q | Index m | Jordan blocks | Steps until W becomes zero | rank Φ(U) | rank Φ(U⁻¹) |
|---:|---:|---|---:|---:|---:|
| 2 | 4 | one block of size 4 | 4 | 2 | 2 |
| 4 | 8 | two blocks of size 8 | 6 | 4 | 2 |
| 8 | 16 | four blocks of size 16 | 15 | 8 | 2 |
| 16 | 32 | eight blocks of size 32 | 36 | 16 | 2 |

For example, the complete `q=8` chain is

\[
 64,62,60,58,53,47,41,34,26,20,14,8,5,3,1,0.
\]

The audit JSON records all four complete chains. The same chain holds for each
inverse curve by (7), and the implementation checks it independently.

### Smith profiles

For a polynomial of degree `L=m-1`, the following are the local Smith orders
of `w^L K(1/w)` at `w=0`. Superscripts in square brackets give multiplicities.

| q | Original geometric row-infinity orders |
|---:|---|
| 2 | `0^[2], 1, 11` |
| 4 | `0^[4], 1^[6], 8^[2], 15^[2], 30^[2]` |
| 8 | `0^[8], 1^[28]`, then each of `7,8,13,18,20,28,30,31,33,34,46,56,64,78` twice |

For the inverse channel at `q=8`, the orders are

\[
 0^{[2]},\ 1^{[49]},\ 2^{[2]},\ 5,9,13,
 96^{[4]},\ 120^{[2]},\ 128^{[2]}.                            \tag{10}
\]

Their sum is `960=64*15`, consistent with row-unimodularity. The finite
Smith computations use precision strictly greater than every reported order.

The existing two-cut sufficient bound gives approximate diagnostic exponents
`rho=0.911927, 0.933965, 0.947360` for the original `q=2,4,8` examples. The
inverse `q=8` example gives `rho≈0.978692`. These are numerical diagnostics of
that sufficient bound, not factoring exponents and not exact certificates of
an improvement.

## 4. Complete polynomial-change screens

Fix a nilpotent `S` of index `m`. Every normalized polynomial

\[
 f(S)=S+\sum_{j=2}^{m-1}c_jS^j,\qquad c_j\in\mathbb F_2,
\]

has the same Jordan type. For `m=16` there are exactly `2^14=16,384` choices.
The checker visits every coefficient mask in Gray order and computes the
exact binary rank of

\[
 \Phi(I+f(S)).                                               \tag{11}
\]

It records the entire rank histogram. No random sampling is used.

For the geometric `q=8` channel the minimum rank is `2`, attained uniquely by
the all-one coefficient mask in (5). The next rank is `6`, attained once.
For the independent regular `q=4` channel from the nilpotent-flag checkpoint,
the minimum is `2`, attained uniquely by `f(S)=S`; the next rank is `4`.

These searches have a direct endpoint interpretation. The highest coefficient
of (2) is `Φ(I+f(S))` in the row cut, and an invertible matrix in the column
cut. If its row rank is `r`, its contribution to an ordinary univariate
tensor-power coefficient decomposition at tensor exponent `t` is

\[
 (q\sqrt r)^t.
\]

Thus an improvement from that decomposition requires

\[
 r<q^{2\sqrt2-2}.                                           \tag{12}
\]

At `q=4`, every rank at least `4` violates (12). At `q=8`, every rank at least
`6` violates it; an exact comparison is `6^35 > 8^29`, combined with
`sqrt(2)<99/70`. Therefore the two searches leave exactly one
endpoint-admissible polynomial change in each stated space.

This endpoint argument excludes the **original** geometric curves from
improving the target through ordinary tensor-power coefficient grouping:
their leading row rank is `q`, giving the protected contribution
`(q^(3/2))^t`. The inverse curves require the following stronger argument.

## 5. A protected subleading-coefficient lemma

Let `K(z)=sum_(i=0)^L z^i K_i` be any polynomial crossing. Write `A=K_L`
and `B=K_(L-1)` in the row cut. Define

\[
 r_* = \operatorname{rank}(B|_{\ker A})
     = \operatorname{rank}\begin{pmatrix}A\\B\end{pmatrix}
        -\operatorname{rank}A.
\]

In particular, `r_* >= rank(B)-rank(A)`.

Restrict the inputs of the row cut of `K(z)^tensor(t)` to
`(ker A)^tensor(t)`. Every summand using the top coefficient `A` vanishes.
The coefficient of degree `(L-1)t` then has only one remaining tensor term:
`B^tensor(t)`. Hence its row rank is at least `r_*^t`.

For a nilpotent-superoperator curve, the top two column coefficients are
`U` and `SU`. Right multiplication by `U^{-1}` makes them `I` and `S`.
On each Jordan block of `S`, choose alternating input and output positions.
The selected rectangle of `I` is zero, while the selected rectangle of `S`
has rank the floor of half the block length. Therefore, with

\[
 h=\sum_j\left\lfloor\frac{b_j}{2}\right\rfloor,
\]

the same degree `(L-1)t` coefficient has column rank at least `h^t`.
The two cut restrictions are independent rank witnesses; they do not need to
come from the same four-tensor restriction. Consequently this single
coefficient contributes at least

\[
 (r_*h)^{t/2}                                               \tag{13}
\]

to the rank cost. This is an actual-rank obstruction valid for every tensor
exponent.

The finite controls give the following values:

| Curve | rank A | rank B | r_* | h | r_* h |
|---|---:|---:|---:|---:|---:|
| Geometric inverse, q=4 | 2 | 11 | 9 | 8 | 72 |
| Geometric inverse, q=8 | 2 | 53 | 51 | 32 | 1632 |
| Geometric inverse, q=16 | 2 | 233 | 231 | 128 | 29568 |
| Earlier regular control, q=4 | 2 | 14 | 13 | 8 | 104 |

In every row, `(r_* h)^35 > q^99`. Since `sqrt(2)<99/70`, this proves
`r_* h > q^(2 sqrt(2))`; the lower bound (13) is already too large.

For `q=2` the simpler leading-coefficient obstruction applies. Combining
(13) with the complete screens in Section 4 excludes every one of the
`16,384` normalized polynomial changes of the `q=8` geometric channel, and
every one of the `16,384` changes of the earlier regular `q=4` channel, from
succeeding through ordinary univariate tensor-power coefficient grouping.

This does not exclude arbitrary linear recombinations of coefficient tensors,
different curves, or different Jordan types. It supplies a fast constructive
filter for new candidates: both the leading reshuffled map and a suitable
restriction of the next coefficient must be small.

## 6. Inverting the row polynomial matrix

This operation is distinct from replacing the superoperator `U` by its inverse.
Since the row polynomial `K(z)` is unimodular and `K(1)=E`, its row-matrix
inverse is again a polynomial crossing, with endpoint `E^{-1}=E`.

The checker constructs this inverse from the common nilpotent flag. A basis
adapted to the spaces `W_j` makes `K_0^{-1}K(z)` unit upper triangular. Polynomial
back-substitution then gives the inverse without adjugates or rational-function
reconstruction. Both matrix products with the original are checked entry by
entry.

For `q=4,8`, inversion gives the following exact column-zero spectra:

| q | Source curve | Row-inverse degree | Column-zero orders | Diagnostic rho |
|---:|---|---:|---|---:|
| 4 | geometric | 23 | `0^[4], 8^[12]` | 0.975629 |
| 4 | inverse superoperator | 25 | `0^[4], 8^[12]` | 0.978449 |
| 8 | geometric | 63 | `0^[14], 8^[12], 16^[26], 24^[4], 32^[4], 40^[4]` | 0.982008 |
| 8 | inverse superoperator | 113 | same as preceding row | 0.993267 |

The row Cartan degrees negate under inversion. Equivalently, if the source
has degree `L` and row-infinity orders `nu_i`, the inverse has degree
`D=max(nu_i)-L` and row-infinity orders `max(nu_i)-nu_i`.

The leading and next row coefficients of the four inverse curves have rank
`2`, so the preceding subleading test does not close this branch. Their
constant coefficients do: the constant row matrix is invertible, while its
column ranks are `4` at `q=4` and `14` at `q=8`. Thus the constant tensor-power
coefficient has protected rank products `64` and `896`, respectively. Both
exceed the target `q^(2 sqrt(2))` by exact integer comparisons.

### All finite basepoints at q=4

An affine change of the parameter can move a different finite point to zero.
For the two `q=4` row inverses, the complete column determinant factorizations
are

\[
 \det K_{\mathrm{geom}}^{-1,\mathrm{col}}(z)
 =z^{96}(z^{14}+z^{12}+z^{10}+z^8+z^4+z+1)^{16},
\]

and

\[
 \det K_{\mathrm{inverse}}^{-1,\mathrm{col}}(z)
 =z^{96}(z^{12}+z^{10}+z^4+z+1)^{16}.
\]

The displayed degree-14 and degree-12 factors are irreducible over `F2`.
The checker verifies the determinants exactly, proves irreducibility by
finite-field polynomial arithmetic, and computes ranks in the corresponding
extension fields. At a root of either nonzero factor, the row rank is `16`
and the column rank is `14`.

Consequently the column ranks at all finite points over the algebraic closure
are `4` at zero, `14` at the other singular points, and `16` elsewhere. Every
finite basepoint therefore fails the protected constant-coefficient test for
ordinary tensor-power coefficient grouping. This is a complete finite-point
statement for these two `q=4` curves; no analogous all-root factorization is
claimed here at `q=8`.

### The adjacent affine-pencil operation

In operator/trace-dual column bases, the inverse of the nilpotent polynomial
satisfies

\[
 z^m F_S(z)^{-1}=U^{-1}(zI-S).
\]

After reshuffling this is the affine pencil

\[
 P(z)=E+(z+1)R,\qquad R=\Phi(U^{-1}).
\]

For literal column flattening, the inverse differs from this display by a
row transpose, because that flattening is `F_S` followed by the transpose
permutation. This distinction preserves all ranks and Smith spectra.

For the geometric formula and every `q` divisible by four, the rank-two
matrix `R` satisfies `(ER)^2=0`. To see this, write the two summands of (3)
as `A_i Y^T B_i`. The two-by-two Gram matrix governing `(ER)^2` has entries
`Tr(B_i^T A_j)`. Counting the entries in the upper triangular all-ones matrix
`P` gives, for even `q`,

\[
 G=\left(\frac q2\bmod2\right)
       \begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]

One way to verify the trace calculation is to put `J=N^T`, so
`C^{-T}=P^2`, `PJ=P+I`, and `NP^T=P^T+I`. The first Gram entry is
`Tr(PP^T)=q(q+1)/2`; the other three differ from it only by `Tr(P)=q`
or `Tr(PJ)=0`. All four therefore equal `q/2` in characteristic two
when `q` is even.

Thus `G=0` when four divides `q`. In particular, the row inverse of the pencil
stays affine:

\[
 P(z)^{-1}=E+(z+1)ERE.
\]

The checker includes literal controls at `q=4,8,16`. This short adjacent
inversion does not by itself produce a new favorable spectrum. Longer
birational sequences with parameter changes are not excluded by this check.

### Constant twists and the column-degree budget

An arbitrary constant row twist of the form

\[
 K'(z)=EF^{-1}E K(z)F
\]

preserves the crossing endpoint and the row Smith data, while it can change
the column data. The following gate allows the most favorable possible column
orders, without restricting `F`.

Let `n=q²`, let `L` be the polynomial degree, and let `d_i` be the fixed row
Cartan degrees. At any finite point the column Smith orders are nonnegative;
at least one is zero because the row matrix is invertible there. Their sum is
at most `nL`, the degree bound for the column determinant. Convexity gives

\[
 \sum_i s^{-\nu_i}\ge
 1+(n-1)s^{-nL/(n-1)}\qquad(s\ge1).
\]

Thus every such twist's present two-cut sufficient bound is at least

\[
 \left(\sum_i s^{d_i}\right)
 \left(1+(n-1)s^{-nL/(n-1)}\right).                          \tag{14}
\]

The checker proves that (14) exceeds the target for both original source
curves and all four row-inverse curves in the preceding table. The proof is
exact. Substitute `s=t^b`, where `b` clears the denominator of `nL/(n-1)`,
subtract the integer target `ceil(q^(99/35))`, and multiply by a power of `t`
to remove negative exponents. Expanding at `t=1+u` gives a polynomial with
nonnegative coefficients and positive constant. In two cases a further
positive factor `t^k` is used before the expansion:

| Row-inverse case | Integer target | Extra power k | Certificate degree |
|---|---:|---:|---:|
| q=4, geometric source | 51 | 0 | 818 |
| q=4, inverse-superoperator source | 51 | 7 | 183 |
| q=8, geometric source | 359 | 32 | 174 |
| q=8, inverse-superoperator source | 359 | 0 | 15296 |

Every coefficient is checked with integer arithmetic. Since
`2 sqrt(2)<99/35`, these certificates exclude a favorable sufficient Smith
bound for all constant twists of the listed curves, even if a twist reduces
the column constant rank to one. This is a statement about that sufficient
bound; it does not replace the actual-rank obstructions proved separately.

## 7. Reproduction

Place the checker beside `check_nilpotent_flag_crossings.py` and run:

```sh
python3 tools/check_geometric_transpose_crossings.py \
  --polynomial-screen \
  --output evidence/geometric-transpose-audit.json

python3 tools/check_geometric_row_inverses.py \
  --output evidence/geometric-row-inverse-audit.json
```

The implementation uses only the Python standard library. Small `q=2,4`
controls compare the packed channel construction against independent literal
matrix products. The checks cover the exact power ranks, both inverse
identities, unitality, trace preservation, every crossing endpoint entry,
the nilpotent-algebra certificates, the displayed finite Smith profiles, and
both complete polynomial screens. It also verifies every rank and exact
integer comparison in the subleading-coefficient table.

The files for this checkpoint are:

- `tools/check_geometric_transpose_crossings.py`
- `research/geometric-transpose-crossings.md`
- `evidence/geometric-transpose-audit.json`
- `tools/check_geometric_row_inverses.py`
- `evidence/geometric-row-inverse-audit.json`

The remaining constructive targets are a proof of the binary-size family,
a change that retains long Jordan blocks while improving the first two
reshuffled coefficients, or a coefficient-recombination method that evades
the protected terms identified here.
