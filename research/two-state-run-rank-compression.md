# Exact rank compression for two-state path crossings

This note accelerates the two-state subfamily of
[finite-state family](finite-state-pairing-crossings.md). It proves the compression and rank
identities used by `two_state_run_profile.cpp`. It does not infer an
asymptotic contraction from a finite profile.

## 1. The local nilpotent representation

Over a field of characteristic two, put `S=R+H=E+C`. In the row cut,
`C=I_4` and `S=I_4+P`, where `P` exchanges the two bit coordinates. Thus
`S` is a square-zero map of rank one. Its Jordan decomposition is

\[
 V\simeq J_2\oplus J_1\oplus J_1.                         \tag{1}
\]

In the column cut, `C=Q` is a square-zero rank-one map and `E=P` satisfies
`PQ=QP=Q`. Consequently `S=P+Q` has `S^2=I`. Multiplying every local
matrix on the left by `S` preserves ranks and sends `S` to `I` and `C`
to `Q`. Therefore the same representation (1) applies, with the two
symbols exchanged.

Let

\[
 A_m=K[\epsilon_1,\ldots,\epsilon_m]/(\epsilon_i^2).
\]

The `J_2` summand is the regular representation of `A_1`; the two `J_1`
summands are trivial representations. For any polynomial `f` in the
commuting square-zero local maps, (1) gives the exact rank identity

\[
 \operatorname{rank}_{V^{\otimes t}} f
 =\sum_{J\subseteq[t]}2^{t-|J|}
   \operatorname{rank}_{A_{|J|}}m_{f(\epsilon_J,0)}.        \tag{2}
\]

Here `m_g` denotes multiplication by `g` in the indicated regular algebra.
Formula (2) is a direct-sum identity with explicit multiplicities.

For a path word with vertex costs `a_S,a_C` and transition costs
`b_SS,b_SC,b_CS,b_CC`, form

\[
 f(z,\epsilon)=\sum_w z^{e(w)}\prod_{i:w_i=S}\epsilon_i.
\]

Its coefficient multiplication maps give the row ranks in (2). Replacing
the product by `prod_{i:w_i=C} epsilon_i` gives the column ranks. Inactive
positions in (2) are therefore forced to `C` for row ranks and to `S`
for column ranks. Each word contributes one square-free monomial, so the
coefficient polynomials can be generated directly without forming an
XOR-circulant kernel or a full crossing matrix.

## 2. Grouping subsets by their active runs

Fix the inactive symbol. Decompose an active subset `J` into its maximal
consecutive intervals. Each interval has a length and two flags recording
whether an inactive position lies immediately to its left and to its
right. An interval's polynomial depends only on these three parameters.
The contributions of separate active intervals multiply in disjoint
variables. All edges lying entirely in inactive intervals contribute a
single scalar monomial.

At fixed total length `t`, that scalar depends only on the number of
inactive positions and the number of inactive intervals. Both are
determined by the multiset of active interval lengths and boundary flags.
Permuting active intervals merely permutes their algebra variables.
Therefore all subsets with the same multiset have identical coefficient
ranks in (2).

The program enumerates all `2^t` subsets to count these types, keeps one
representative per type, and multiplies its rank contribution by the type
count and `2^{t-|J|}`. The all-active interval has both flags absent; the
empty active subset is its own type. These cases are included explicitly.

| Positions `t` | Active subsets | Run types | Total assignments over representatives |
|---:|---:|---:|---:|
| 8 | 256 | 121 | 4,707 |
| 10 | 1,024 | 285 | 32,691 |
| 12 | 4,096 | 619 | 215,763 |
| 14 | 16,384 | 1,264 | 1,368,083 |
| 16 | 65,536 | 2,456 | 8,393,363 |

The last column is `sum_types 2^{|J|}`. No claim that the number of types
is bounded by a constant is needed. This is an exact finite computation.

## 3. Exact algebraic rank reductions

### 3.1 Units, linear terms, factors, and absent variables

An element with nonzero constant term is a unit, so its multiplication
rank is `2^m`. If its constant is zero and its linear part is nonzero,
an invertible change of variables sends it to `epsilon_1`; the rank is
`2^{m-1}`. These are the same local-algebra arguments used for the
XOR-circulant shortcuts in the preceding note.

If every monomial contains a common set `B` of variables, factor out
`prod_{i in B} epsilon_i`. Multiplication by this monomial on its own
variables has rank one. Variables absent from all monomials contribute
identity factors of rank two. The remaining polynomial therefore gives
the rank after multiplying by `2` for each absent variable. This also
handles a polynomial consisting of one monomial.

### 3.2 Quadratic pivot lemma

**Lemma.** Suppose `f` has zero constant and linear terms, and the
coefficient of `xy` is one. In the algebra on the remaining `m-2`
variables write uniquely

\[
 f=xy\,a+x\,b+y\,c+d.
\]

Then `a(0)=1`, while `b,c,d` have zero constants. Put

\[
 g=d+abc.
\]

There is an algebra automorphism sending `f` to `x'y'+g`, and

\[
 \boxed{\operatorname{rank}m_f
 =2^{m-2}+2\operatorname{rank}m_g.}                       \tag{3}
\]

**Proof.** Every element of the maximal ideal squares to zero in
characteristic two. Hence `a^2=1`, and the change

\[
 x'=ax+c,\qquad y'=y+ab
\]

preserves the square-zero relations. Its inverse is
`x=a(x'+c), y=y'+ab`, with the other variables fixed. Expanding gives

\[
 x'y'=axy+xb+yc+abc,
\]

so `f=x'y'+g`.

Let `B=A_{m-2}`. The two subspaces `Bx'` and `By'` each contribute a
copy of multiplication by `g`. On `B+Bx'y'`, multiplication by `x'y'+g`
has block matrix

\[
 \begin{pmatrix}m_g&0\\I&m_g\end{pmatrix}.
\]

Because `g` has zero constant, `g^2=0`. Its kernel consists exactly of
`(gv,v)` for `v in B`, and therefore this block has rank `dim B`.
Adding the other two contributions proves (3).

The zero-constant condition is essential to the displayed block rank.
For `g=1`, the element `xy+1` is a unit and has rank four on two
variables; substituting rank one into (3) would incorrectly give three.
The implementation tests the unit case before attempting a quadratic
pivot. A residual produced by the lemma still has zero constant and
linear parts. If a square-zero residual with a nonzero linear part is
considered separately, the block proof remains valid and gives the
expected half-dimensional rank.

For a homogeneous quadratic polynomial, an alternating change of basis
gives a sum of `r` disjoint terms `x_i y_i`. Iterating (3) gives

\[
 \operatorname{rank}m_f
 =\frac{2^m-2^{m-r}}2,                                  \tag{4}
\]

where `2r` is the rank of its alternating coefficient matrix. The program
uses (4) for this special case and (3) for an inhomogeneous quadratic
pivot. It falls back to exact binary Gaussian elimination when the
polynomial has no suitable reduction.

### 3.3 Homogeneous degree blocks

If `f` is homogeneous of degree `k`, multiplication splits into maps

\[
 (A_m)_w\longrightarrow(A_m)_{w+k}.
\]

Their ranks add. The pairing that takes the coefficient of
`epsilon_1...epsilon_m` in the product `uv` is perfect, and multiplication
by `f` is self-adjoint for this pairing. Thus the block with source
degree `w` is the transpose of the block with source degree `m-k-w`
after complementary bases are chosen. The program computes only
`w <= (m-k)/2`, doubling the rank except at the middle block.

This is particularly effective for the count/run bigrading. Each
coefficient fixes the symbol count, so every regular-algebra coefficient
is homogeneous even though its monomial support depends on the order of
the positions.

## 4. Verification and scope

`check_square_zero_rank_reduction.py` independently checks the quadratic
recurrence against full multiplication matrices for every one of the
65,536 elements of `A_4`, and for 60 additional random examples in five
through nine variables. It explicitly checks the constant-term exception
and a linear residual. Its evidence is
`square-zero-rank-reduction.json`.

`check_two_state_run_profile.py` compiles the native sources, replays the
65,536-element test in the native engine, compares twelve three-position
examples with direct expansions of the full four-leg tensors, and
compares twelve six-position examples with the independent
XOR-circulant engine. It also replays stored independent eight- and
twelve-position profiles, including the count/run bigrading. Every
comparison uses both full integer rank lists. Its evidence is
`two-state-run-verification.json`.

The output uses 64-bit total ranks and 128-bit rank products. This avoids
the overflow that a `2^{16}`-dimensional regular block or a `2^{32}`-rank
full cut would cause in narrower integer types. Floating-point logarithms
and square roots appear only in reported descriptive metrics.

For example:

```sh
g++ -O3 -std=c++17 tools/two_state_run_profile.cpp -o /tmp/two_state_run_profile
/tmp/two_state_run_profile profile 16 /tmp/two-state-domainwall-t16.json 0 5 0 1 1 0
/tmp/two_state_run_profile profile 16 /tmp/two-state-run-count-bigrading-t16.json 0 9 0 1 0 0
python3 tools/check_square_zero_rank_reduction.py
python3 tools/check_two_state_run_profile.py
```

The six numerical parameters after the output path are, in order,
`a_S,a_C,b_SS,b_SC,b_CS,b_CC`. These algorithms establish the exact
finite profiles and make longer constructive searches practical. An
asymptotic rank theorem for an improved family is a separate requirement.
