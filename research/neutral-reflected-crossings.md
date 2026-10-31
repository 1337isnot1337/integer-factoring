# Reflected monomial duals and a neutral-fermion crossing compiler

## Results and scope

There is an explicit homogeneous dual family inside the modular partition
construction, obtained by reflecting every part below a prime power.
It requires no Gram-matrix inverse.  In characteristic two the complete
crossing has a finite Pfaffian formula that avoids enumerating intermediate
partitions.

The Pfaffian formula also defines an exact Laurent-polynomial crossing on
every fixed-length, fixed-weight subset of distinct odd indices below a
power of two.  A common monomial makes it polynomial.  Its endpoint is
proved directly, including sectors that are not balanced under reflection.
These sectors restrict the physical labels at each copy, rather than
requiring one occupation from each reflection pair as in
[reflection-pair compression](reflection-dual-rail-compression.md).
Both constructions use Pascal identities; no inequivalence of all their
possible gauge-transformed tensors is asserted.

An eight-letter example has a complete degree-eight polynomial seed after
an exactly checked variable compression.  Its seven nonzero coefficients
and both ranks are given below.  This is a constructive search family, not
an improved factoring bound.  The supported exponent remains `1/2`.

The source operators and modular quotient come from
[modular-crossing.md](modular-crossing.md), based on Claude's Appendix B
at [`anthropics/formal-math`, revision `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5).
The separate [parameter-descent proof](modular-parameter-descent.md)
justifies interpreting both operator parameters and states in the quotient.
All new results here are ordinary proofs with exact finite controls.

## 1. Reflected homogeneous duality in every prime characteristic

Fix a prime `p`, a power `Q=p^r`, and an even positive integer `l`.
Let `h=Ql/2`.  Choose any family `L` of partitions of weight `h` and
length `l`, with every part strictly between zero and `Q`, no part divisible
by `p`, and every multiplicity less than `p`.

For `mu` in `L`, put

\[
 \mu^*=\{Q-k:k\in\mu\},\qquad
 z_\mu=\prod_k k^{m_k(\mu)}m_k(\mu)!,\qquad
 a^\mu=z_\mu^{-1}t^{\mu^*},\quad b^\mu=t^\mu.
 \tag{1}
\]

Both parameters have weight `h`.  The scalar `z_mu` is a unit modulo `p`.

**Theorem.** For all `mu,nu` in `L`,

\[
 [1]a^\mu[-h]b^\nu=[\mu=\nu].                         \tag{2}
\]

**Proof.** Only annihilation terms can contribute to the vacuum.  For
operator part `Q-k` and input part `j`, the source scalar is

\[
 (-1)^{Q-k-1}j\binom{j+Q-k-1}{Q-k-1}.                  \tag{3}
\]

Here `0<j,k<Q`.  Lucas's identity says that a nonzero binomial coefficient
in (3) requires every base-`p` digit of `j` to be at most the corresponding
digit of `k`.  In particular, `j<=k`.  Every contributing matching of the
parts therefore has `sum j<=sum k`.  The two sums are both `h`, so equality
holds in each matched inequality: `j=k` for every pair.  Thus the multisets
must agree.

On the diagonal, `binom(Q-1,k)=(-1)^k mod p`.  The two signs in (3) cancel
because `(-1)^(Q-1)=1` in `F_p`.  The one-part scalar is `k`; repeated-part
matchings give the additional factorial.  The unnormalized diagonal is
exactly `z_mu`, proving (2).

The same proof applies before balancing the reflected weight: if every
label has length `l` and weight `H`, its reflected parameter has weight
`Ql-H` and its mode `-H` still has pairing (2).  Balanced weight is the
case in which both parameters are homogeneous of the same weight, as in
the existing homogeneous crossing interface.

## 2. The characteristic-two quotient is a Clifford module

In `W_2`, the surviving variables are the odd-index `t_j`, with `t_j^2=0`.
Write `X_j` for multiplication and `D_j` for differentiation.  They satisfy

\[
 X_j^2=D_j^2=0,\qquad
 D_jX_k+X_kD_j=[j=k]I.                                 \tag{4}
\]

Operators with distinct indices commute, which is also anticommutation
in characteristic two.  Thus normally ordered source fields can be
evaluated by the usual matching expansion for Clifford generators.
All signs become one; the matching sum is the Pfaffian of the symmetric
zero-diagonal contraction matrix.  This assertion can also be proved
directly by moving the first generator through the others using (4)
and inducting on the number of generators.

For odd `m`, use the field with its weight shift removed:

\[
 \Phi_m(z)=
 \sum_{j\text{ odd}}\binom{j-1}{m-1}X_jz^j
 +\sum_{j\text{ odd}}\binom{j+m-1}{m-1}D_jz^{-j}.
 \tag{5}
\]

Each binomial coefficient is reduced modulo two.  Products belonging to
one source parameter remain normally ordered; there are no contractions
inside that parameter.

Fix `Q` a power of two and odd `k,v<Q`.  The contraction between
`Phi_(Q-k)(z)` and `Phi_v(1)` is

\[
 \begin{aligned}
 C_{k,v}(z)
 &=\sum_{j\text{ odd}>0}
   \binom{j+Q-k-1}{Q-k-1}\binom{j-1}{v-1}z^{-j}\\
 &=[v\subseteq k]\,
   \frac{z^{-v}(1+z^{-1})^{k-v}}{1+z^{-Q}}\\
 &=[v\subseteq k]\,z^{Q-k}(1+z)^{k-v-Q}.
 \end{aligned}                                          \tag{6}
\]

The symbol `v subseteq k` denotes bitwise inclusion.  To prove the second
line, Lucas reduces the coefficient to
`[v subseteq (j mod Q) subseteq k]`.  Sum first over the free bits between
`v` and `k`, then over the nonnegative multiples of `Q`.  This is an
identity of the expansion at infinity of a rational function.  It does
not use an analytic convergence argument.

## 3. An explicit finite Laurent kernel

Let `L` now be any family of distinct-odd-part subsets of `{1,...,Q-1}`,
all with the same length `l` and the same weight `H`.  No balancing
condition `2H=Ql` is imposed in this section.

For labels `x,y,v,u` in `L`, form four groups of nodes:

- `A_k` for each `k` in `x`, representing the reflected parameter;
- `B_v` for each part `v` in the second parameter;
- `Y_j` for the input parts in `y`;
- `U_j` for the output parts in `u`.

Put `c=1+z`.  The symmetric kernel has zero entries within each group,
and these entries between groups:

| Pair | Entry |
| --- | --- |
| `A_k,B_v` | `[v subseteq k]` |
| `A_k,Y_j` | `[j subseteq k] z^(k-Q-j)c^(Q-k)` |
| `A_k,U_j` | `[(Q-k-1) subseteq (j-1)] z^(k-Q+j)c^(Q-k)` |
| `B_v,Y_j` | `[(v-1) & j = 0] c^v` |
| `B_v,U_j` | `[(v-1) subseteq (j-1)] c^v` |
| `Y_j,U_k` | `[j=k]` |

Define `T(z;x,y,v,u)` to be the Pfaffian of this `4l`-by-`4l` matrix.
Its entries are finite Laurent polynomials, so `T` is finite as well.

For the source interpretation, start with (5)--(6), and multiply every
`A_k` node by `z^(k-Q)c^(Q-k)` and every `B_v` node by `c^v`.  Their total
product is

\[
 z^{H-Ql}c^{Ql}.                                         \tag{7}
\]

Equation (6) becomes the constant Pascal entry in the first row of the
table.  All other entries follow directly from the creation, annihilation,
and unchanged-state terms in (5).  The matching expansion proves the
table and its source relation.  In the balanced case `H=h=Ql/2`, (7)
is exactly the generating prefactor for the homogeneous source pieces:

\[
 \sum_s Q_s(a^x,b^v)z^s
 =z^{-h}(z-1)^{2h}\sum_r a^x[r]b^v[-r]z^r.
 \tag{8}
\]

The right side of (8) is read entrywise between the specified states.
There are no factorial denominators at `p=2`, since all parts are odd
and distinct.  The independent checker compares the Pfaffian formula
with the existing source recurrence before using it for new sectors.

## 4. Endpoint and polynomiality

**Endpoint theorem.** Every fixed-length, fixed-weight sector satisfies

\[
 T(1;x,y,v,u)=[x=v][y=u].                                \tag{9}
\]

At `z=1`, every positive power of `c` vanishes.  Only the `A,B` Pascal
block and the `Y,U` identity block remain.  The Pfaffian is therefore
`[y=u] det([v subseteq k])_(k in x,v in v)`.

Every nonzero determinant term matches each part of `v` to a numerically
greater or equal part of `x`.  Their total weights agree, so every match
is equality.  Thus the determinant is zero unless `x=v`, and is one
when they agree.  This proves (9), including arbitrary subfamilies of
the sector.

Every nonzero `A,Y` entry has lowest exponent at least `-Q`, because
`j subseteq k` implies `j<=k`.  All other entries have nonnegative
lowest exponents.  There are only `l` nodes in group `A`, so

\[
 z^{Ql}T(z;x,y,v,u)\in\mathbf F_2[z].                    \tag{10}
\]

The highest degree of an `A,U` entry is `j`; that of a `B,Y` or `B,U`
entry is `v`; that of an `A,Y` entry is `-j`.  A matching therefore
has highest degree at most `2H`.  The polynomial in (10) has degree
at most `Ql+2H`.  Its endpoint remains (9).

The common monomial shifts the two local profiles in opposite directions;
it is retained in the literal polynomial seed.  No local denominator is
dropped.  Every finite fixed choice supplies an effective seed for the
existing [graded-power transfer](graded-even-power-transfer.md), if its
complete pressure ever meets the required strict inequality.

## 5. An eight-letter complete seed

Take `Q=16`, `l=4`, and `H=32`.  There are eight labels.  The compiler
constructs every one of the `8^4=4096` Laurent entries and checks every
nonzero exponent.  For this finite tensor, every exponent is divisible
by 16.  Replacing `z^16` by a new variable and translating the common
lowest exponent produces a degree-eight polynomial crossing.

The complete coefficient profile before translation is:

| Original exponent | Compressed exponent | Row rank | Column rank |
| ---: | ---: | ---: | ---: |
| -64 | -4 | 1 | 64 |
| -32 | -2 | 8 | 30 |
| -16 | -1 | 6 | 6 |
| 0 | 0 | 29 | 29 |
| 16 | 1 | 6 | 6 |
| 32 | 2 | 30 | 8 |
| 64 | 4 | 64 | 1 |

The coefficient at compressed exponent zero already has both ranks 29;
these individual ranks do not establish a better asymptotic contraction.
The [literal seed](../evidence/neutral-reflected-q16-l4-h32-seed.json)
retains all coefficients and the verified translation/compression data.
Exponent compression is tested from the complete tensor for each seed;
no general divisibility theorem is inferred from these examples.

The simpler length-two balanced sectors at `Q=8,16,32` have three nonzero
coefficients with rank pairs `(1,D^2)`, `(1+binom(D,2),1+binom(D,2))`,
and `(D^2,1)`, where `D=2,4,8`, respectively.  These are finite profile
statements here, not an all-size rank theorem.

## Verification and next test

```sh
python3 tools/check_neutral_reflected_crossings.py
```

The checker verifies reflected duality in several prime characteristics,
the periodic contraction formula, exact Pfaffian/source comparisons,
all eight-letter endpoint entries, complete coefficient ranks, and the
literal seed reconstruction.  It tests an unbalanced sector separately.
A negative control removes the fixed-weight condition and exhibits an
off-diagonal Pascal endpoint.

The constructive next test is a search over complete physical label
weights and the two paired label distributions for these fixed-sector
seeds.  Such a search must retain both complete Smith profiles.  Any
floating candidate needs exact profile and pressure replay before it can
be used as a factoring bound.  No such strict witness is asserted here.

This construction, derivation, and checker used substantial OpenAI
assistant support.
