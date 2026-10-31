# Large linear families of polynomial reflection frames

## Constructive result

For every power-of-two reflection interval `M>=8`, put `m=M/2` and
`q=2^m`. Split the physical letters according to their last logical
bit. In Bell coordinates, every matrix

```math
G=\begin{pmatrix}0&F\\0&0\end{pmatrix},
\qquad
F_{a,b}=0\ \text{unless }\
\operatorname{popcount}(a)\ne\operatorname{popcount}(b)\pmod2
                                                               \tag{1}
```

gives a polynomial coupled frame `A=I+w^-k G` for each
`k=1,2,3,4`. The two blocks have size `q/2`; the indices `a,b` in
(1) are the lower `m-1` logical bits.

The entries allowed in (1) are independent parameters over any
characteristic-two field. The family has dimension `q^2/8`, can
contain generators of rank `q/2`, and retains the exact crossing
endpoint and the original degree. It is not restricted to
one-particle or exterior-power frame generators.

More generally, within the entire upper-block space
`G=[[0,F],[0,0]]`, admissibility for every `k<=M/2` reduces to a
system of **linear** equations in `F`. The quadratic terms in the
[uniform physical-jet criterion](uniform-reflection-jet-criterion.md)
vanish identically on this space.

These are constructive polynomial-path results. Complete weighted
profiles and common physical types still have to be charged. No
improved pressure or factoring exponent is claimed.

Independent Laurent coefficients are also allowed. At interval
`M>=8`, the frame
`A=I+w^-1 G_1+w^-2 G_2+w^-3 G_3+w^-4 G_4`
is polynomial after the coupled action whenever all four matrices
have this upper-block support, with `G_1,G_2` arbitrary and only
`G_3,G_4` parity preserving. This gives `3q^2/4` independent field
parameters. Section 5 proves the complete joint linear criterion.

## 1. A triangular one-particle representative

Use the matrices and literal conventions of the uniform jet theorem:

```math
Q(w)=\rho((1+w)H_M(w^2)),\qquad
B(w)=Q(0)^{-1}Q(w),
```

where `rho` takes the selected `M/2` by `M/2` reflection minors.
Modulo `w^(2M)`, selected-minor multiplicativity gives

```math
B(w)=\rho(C_M(w)),\qquad
C_M(w)=G_M(0)(1+w)H_M(w^2).                     \tag{2}
```

Here `G_M(0)` is an exact symplectic matrix and the second factor is
symplectic over the truncated ring. Thus the ring argument from the
uniform theorem applies directly.

The bitwise formulas factor the one-particle product completely:

```math
C_M(w)=(1+w)\bigotimes_{r=0}^{\log_2M-1}
\begin{pmatrix}
1+w^{2^{r+1}}&0\\
w^{2^{r+1}}&1
\end{pmatrix}.                                  \tag{3}
```

Indeed, on one bit the required multiplication is

```math
\begin{pmatrix}1&1\\1&0\end{pmatrix}
\begin{pmatrix}z&1\\1&1\end{pmatrix}
=\begin{pmatrix}1+z&0\\z&1\end{pmatrix},
```

with `z=w^(2^(r+1))`. Tensoring proves (3) at every interval.

## 2. The second physical jet is literal parity

Equation (3) also removes the scalar ambiguity in the second jet:

```math
B_2=P,\qquad P_{ss}=\operatorname{popcount}(s)\pmod2
                                                               \tag{4}
```

for every power-of-two `M>=4`.

To prove (4), remove the scalar `(1+w)` before taking selected minors.
At order `w^2`, only the lowest one-particle bit changes. Its diagonal
entry is `1+w^2` on an even index and one on an odd index; its only
off-diagonal transition is from an even index `j` to `j+1`.

A single such off-diagonal transition cannot take one retained
reflection letter to another. It removes a mode from one reflection
pair and adds a mode from a different pair. The two indices would be
reflected partners only if `j=m-1`, but `m-1` is odd while this
transition requires even `j`. Therefore every off-diagonal entry of
the second selected-minor coefficient vanishes.

For a letter `s`, the diagonal coefficient counts selected even
one-particle indices modulo two. Starting from the all-low letter
there are `m/2` such indices; each flipped logical bit toggles that
parity. This contribution is
`m/2+popcount(s) mod2`. The scalar on an `m` by `m` minor is

```math
(1+w)^m=1+w^m.
```

Its coefficient at `w^2` equals `m/2 mod2`: both are one when `m=2`
and zero for the larger power-of-two values. The constant terms
cancel, leaving (4).

The normalized parity matrix `H=Q(0)^-1 P Q(0)` from the first-jet
proof is also exactly `P`, uniformly. That proof and (4) imply
`H=P+cI`. The entry `Q(0)_(0,q-1)` is one: its selected minor is
`H_m(0)`, which is invertible with determinant one. Both endpoint
letters here have even parity because `m` is even. The relation
`P Q(0)=Q(0)(P+cI)` at this entry forces `c=0`.

Thus the literal parity descriptions previously checked at small
alphabets hold throughout the reflection family.

## 3. A common invariant half of the physical alphabet

Let `L` be the span of all letters whose last logical bit is zero.
The selected one-particle set of such an input letter contains index
`m-1`. A letter outside `L` does not contain that index; it contains
the reflected index `m` instead.

Column `m-1` of (3) is particularly simple. Its only nonzero entries
are

```math
(C_M)_{m-1,m-1}=(1+w)(1+w^M),
\qquad
(C_M)_{M-1,m-1}=(1+w)w^M.                       \tag{5}
```

Any determinant term taking an input letter in `L` to an output
outside `L` therefore uses the second entry in (5), or is zero.
It is divisible by `w^M`. By (2),

```math
B_j L\subseteq L\qquad(0\le j<M).               \tag{6}
```

This establishes the invariant subspace for every low jet needed by
the bounded-pole theorem. It is not inferred from finite matrices.

## 4. Quadratic admissibility becomes linear

Take any `G` with image in `L` and with `G|L=0`; these are exactly
the upper-block matrices in the result. They satisfy `G^2=0`.
Equation (6) gives

```math
G B_j G=0\qquad(0\le j<M).                       \tag{7}
```

For `k<=M/2`, all quadratic terms used by the uniform criterion have
indices below `M`. Consequently its necessary and sufficient
conditions reduce to

```math
[G,B_j]=0\qquad(0\le j<k).                       \tag{8}
```

These equations are linear in the entries of `G`.
Since `B_0=I`, all odd jets vanish, and `B_2=P`, the consequences are:

| Pole range | Conditions within the full upper-block space |
|---|---|
| `k=1,2` | none |
| `k=3,4` | `[G,P]=0` |
| `k=5,6` | `[G,P]=[G,B_4]=0` |
| `k=7,8` | `[G,P]=[G,B_4]=[G,B_6]=0` |

Each row is used only when its pole satisfies `k<=M/2`. No condition
on a higher pole is extrapolated beyond that range.

Commuting with `P` means that a nonzero matrix entry joins physical
letters of equal total parity. The source half has its last bit one,
while the target half has that bit zero. Thus its lower-bit labels
must have opposite parity, exactly the constraint in (1).

There are `q/4` even and `q/4` odd letters in each half. Hence (1)
has `2(q/4)^2=q^2/8` independent entries. Choosing an invertible
map between the appropriate lower-bit parity classes gives rank
`q/2`. Arbitrary nonbinary coefficients and arbitrary smaller ranks
are allowed.

## 5. Independent Laurent coefficients and a complete linear criterion

Fix `K<=M/2` and let

```math
N(w)=\sum_{j=1}^{K}w^{-j}G_j,\qquad A=I+N,
                                                               \tag{9}
```

where every `G_j` has image in the same half `L` and kills `L`.
All products `G_i G_j` vanish, so `N^2=0` and `A^-1=A`.
Also (6) gives

```math
G_i B_r G_j=0\qquad(r<M)                         \tag{10}
```

for every pair of coefficient matrices, whether or not they preserve
parity.

The inner-automorphism approximation and matrix-unit argument in the
uniform theorem apply with maximum pole `K`. Its error has order at
least `2M-4K>=0`, while

```math
ABA^{-1}=B+NB+BN+NBN.                            \tag{11}
```

Every possible negative coefficient of `NBN` uses a jet
`r<i+j<=2K<=M`; it vanishes by (10). Thus the complete coupled
tensor is polynomial **if and only if**

```math
\sum_{\substack{1\le j\le K\\e+j\ge0}}
 [G_j,B_{e+j}]=0
\qquad(-K\le e<0).                              \tag{12}
```

These are joint linear equations in all the independent coefficient
matrices. They include any cancellation between different Laurent
orders; no separate regularity requirement is imposed on each term.

For `K=2`, (12) uses only `B_0=I` and `B_1=0`, so both upper-block
matrices are arbitrary. In the four-letter case this is eight free
`F4` coefficients, giving exactly `4^8=65536` coefficient choices.
Distinct coefficients define distinct physical frames; the resulting
tensors are not asserted to be distinct.

For `K=4`, (12) says precisely

```math
[G_3,P]=[G_4,P]=0.
```

There is no restriction on `G_1,G_2`. The total parameter dimension is

```math
2\frac{q^2}{4}+2\frac{q^2}{8}=\frac{3q^2}{4}.
```

At alphabet sixteen this is 192 independent `F4` coefficients.
The parity restriction at the two highest poles is essential: an
odd-parity entry in `G_3` produces a nonzero negative linear jet.

At one the matrix `I+sum_j G_j` is still its own inverse, and at
infinity it is identity. Thus endpoint and degree preservation hold
for (9), just as for the single-coefficient frames.

## 6. Effective construction and limitations

To choose a pole-three or pole-four frame, fill any allowed entries
of `F`, set the other entries to zero, and use the general coupled
constructor. The polynomiality proof applies to every choice.
The construction does not need a nonlinear search or a separate
full-tensor pole test for its mathematical validity. The implemented
finite searches nevertheless retain that exact check as a control.

For larger permitted poles, construct the physical jets as in the
uniform theorem and solve the linear equations (8) on this known
upper-block space. One can then choose arbitrary coefficients in
the resulting linear subspace. This does not assert that different
chosen frames can be composed while preserving the same jet
hypotheses; each frame is applied to the original reflection seed.

The family can substantially change physical gauge costs. In
particular, an allowed frame need not preserve the natural weighted
infinity profile. Polynomiality and a large parameter space do not
give a rank bound for free.

The independent-coefficient family is obtained by solving (12)
instead. Its input has at most `K q^2/4` field coordinates and the
displayed conditions are explicit finite linear equations. Ordinary
Gaussian elimination supplies a basis effectively. A concrete full
tensor or pressure computation is still charged separately.

## 7. Exact controls

The checker verifies the one-particle factorization (3) entry by entry
through interval 64, including the two-entry column (5). It constructs
complete physical jets at intervals 4, 8, and 16, checks `B_2=P`,
checks that `Q(0)` preserves parity, and verifies every forbidden
lower-left entry below order `M`.

At intervals four and eight, full-rank and dense nonbinary upper-block
examples are also compared against the complete four-leg Laurent
test through the claimed pole range. The interval-16 check uses
256-letter physical matrices and does not construct `256^4`
tensor entries.

Additional independent exponent-dictionary controls construct
four-letter frames with two unrelated upper-block coefficient
matrices, including nonbinary and parity-changing coefficients.
They check polynomiality, endpoint, degree, and the involution.
Sixteen-letter matrix-jet controls test all mixed terms in (11),
including a rejected odd-parity highest-pole coefficient.

```sh
python3 tools/check_triangular_reflection_jet_families.py \
  --output evidence/triangular-reflection-jet-families.json
```

The [checker](../tools/check_triangular_reflection_jet_families.py) and
[evidence](../evidence/triangular-reflection-jet-families.json) retain
the exact matrices, counts, and source digests. These controls test
the implementation and indexing. Equations (2)--(12) prove the
uniform families.
