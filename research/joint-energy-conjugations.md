# Joint energy conjugations and the five-letter reference

Conjugating the crossing operator on a **joint pair of physical labels**
gives a class of strict crossings that preserves one complete weighted
Smith profile. This is broader than changing the basis on each physical
leg separately. The other profile can change, so it is a possible route
from a strong seed without sacrificing both places at once.

For the retained five-letter Gaussian reference, the complete F4 family
that also fixes the leading weighted infinity coefficient has 1,296
distinct tensors. Exact fraction-free profiles and prefix majorization
show that none improves the reference's weighted Smith envelope at its
natural heights, at **any real scale `s>=1`**. Fixing its next coefficient
as well forces the conjugation to leave the entire reference unchanged.

These are a construction and a scoped obstruction, not a better factoring
algorithm or a lower bound on arbitrary crossings. The score stays `1/2`.

## 1. A transformation preserving one place exactly

Let `T(t)` be a scalar polynomial crossing on `q` letters, with `T(1)=S`,
where `S(x,y,u,v)=[x=v][y=u]`. Regard it as a matrix with row `(x,y)`
and column `(u,v)`, and put `A(t)=T(t)S`. Thus the columns of `A` are
in `(v,u)` order and `A(1)=I`.

Choose physical half-energies `h_x`, and let

```math
H_{(x,y)}=h_x+h_y,\qquad D_H(z)=\operatorname{diag}(z^{H_{(x,y)}}).
```

For any constant invertible matrix `U` on the `q^2` joint labels that
preserves the exact **integer** `H` blocks, define

```math
A_U(t)=U A(t)U^{-1},\qquad T_U(t)=A_U(t)S.               \tag{1}
```

Then `T_U(1)=S`, with no change in polynomial degree. At zero, use the
row cut and the gauges `h_x+h_y` on the row side and zero on the other.
Up to the fixed column swap, its matrix is `D_H(z)A(z)`. Since `U`
commutes with `D_H(z)`,

```math
D_H(z)A_U(z)=U\bigl(D_H(z)A(z)\bigr)U^{-1}.             \tag{2}
```

Both factors are constant units. Consequently the **entire zero-place
Smith profile**, including every determinantal order, is unchanged.
The condition uses integer energy equality, not energy equality modulo
the coefficient characteristic.

At infinity, use the other cut with entry gauge `h_u+h_v`. Before
realignment into rows `(x,u)` and columns `(y,v)`, this is

```math
B(z)=A(z^{-1})D_H(z).
```

Equation (1) replaces it by `U B(z)U^-1`. Realignment is not a left/right
matrix equivalence, so its Smith profile can change. This is the useful
asymmetry of the construction. For a fixed finite seed, complete entries
and factors retain the existing `2^O(label width)` generation bound.

## 2. The retained reference and its coefficient algebras

Use the exact degree-18 five-letter reference from
[the aligned Gaussian construction](aligned-gaussian-complement-perturbations.md),
computed by `five_letter_source` in
[`search_f4_recentered_crossings.py`](../tools/search_f4_recentered_crossings.py).
Its physical letters are `(0,2,4,8,12)` in the reflection construction,
and its natural half-energies in reference order are

```text
h = (0,5,3,1,4).
```

The 25 joint-label energies have block sizes
`1,4,2,2,4,1,3,3,2,2,1` in first-occurrence order. Thus the full energy
block algebra has dimension `sum d_i^2=69`. Write

```math
B(z)=\sum_k z^k B_k,\qquad \min k=-8.
```

Exact binary linear algebra gives two nested algebras:

```math
\mathcal L=\{U:[U,H]=0,\ [U,B_{-8}]=0\},\quad\dim\mathcal L=10,
```

```math
\mathcal Z=\{U\in\mathcal L:[U,B_{-7}]=0\},\quad\dim\mathcal Z=4.
                                                               \tag{3}
```

Here `[U,H]=0` denotes the integer energy block condition; it is not a
commutator with an energy matrix reduced modulo two. Each basis matrix
of `Z` commutes with **every** coefficient `B_k`, not just the first two.
Since it also commutes with `D_H`, it commutes with `A(t)`. Thus `Z` is
the full stabilizer algebra of the reference within `L`.

All the equations have binary coefficients. Their kernels extend by base
change to every field of characteristic two with these same dimensions.
It follows in particular that a constant energy-preserving conjugation
fixing both complete matrices `B_-8` and `B_-7` is trivial on `T(t)`.
Preserving their ranks or Smith orders is a weaker condition and is not
excluded by this statement.

## 3. Counting the complete orbit without a million-case assumption

On each energy block, the ten binary basis matrices of `L` admit one
common complete invariant flag over `F2`. In that basis all are upper
triangular. The distinct diagonal entries of a general element are six
independent binary linear forms in its ten coordinates. It is invertible
exactly when all six forms are nonzero.

The four-dimensional stabilizer algebra `Z` has two independent diagonal
forms in the same flags. Therefore, over any finite field `F_Q` of
characteristic two,

```math
|\mathcal L^\times|=Q^4(Q-1)^6,\qquad
|\mathcal Z^\times|=Q^2(Q-1)^2.                          \tag{4}
```

Orbit-stabilizer gives the exact number of distinct conjugated tensors

```math
Q^2(Q-1)^4.                                            \tag{5}
```

In particular, the binary orbit has four tensors. Over `F4` there are
186,624 units, 144 stabilizer units, and **1,296** tensors.

The implementation also constructs the orbit directly without testing
every algebra element. The bimodule spanned by `L B_-7 L` has dimension
16. Twenty-one invertible matrices of the form `I+alpha L_i` act on
that module and generate 1,296 distinct images of `B_-7`. The two-coefficient
stabilizer result in (3) means these images distinguish the complete
tensors. Reaching the cardinality in (5) proves completeness; it is not
necessary to assume that the chosen generators generate the whole unit
group. The corresponding units are retained as orbit representatives.

## 4. Exact profile dominance at every scale

At the natural heights, the reference's sorted infinity orders are

```text
(-8,-7,-5,-4,-3,-2,-1,0,0,1,1,2,2,3,3,4,4,4,5,5,6,7,7,8,12).
```

Its zero orders are this list plus 16. Equation (2) preserves the latter
profile exactly. The 1,296 F4 conjugates have 21 distinct infinity
profiles. For **every** conjugate with sorted orders `mu`, the exact
calculation gives

```math
\sum_{i=1}^k\mu_i\le\sum_{i=1}^k\mu_i^{\rm ref}
\qquad(1\le k\le25).                                  \tag{6}
```

This is sufficient for a uniform scale comparison. For sorted `mu`,

```math
\sum_i(t-\mu_i)^+=\max_{0\le k\le25}
 \left(kt-\sum_{i=1}^k\mu_i\right).
```

Thus (6) increases every such hinge sum. For `a=log s>0`, integrate
against the positive function `a^2 exp(-a t)`; the identity

```math
a^2\int_{-\infty}^{\infty} e^{-at}(t-\mu)^+dt=e^{-a\mu}
```

gives

```math
\sum_i s^{-\mu_i}\ge\sum_i s^{-\mu_i^{\rm ref}}
\qquad(s>1).                                          \tag{7}
```

At `s=1` both sums are 25. The zero-place partition, physical energies,
type entropy and restoration terms are unchanged, so (7) makes the full
weighted Smith cost no better at **any** scale. This is an exact comparison
to the reference; it does not rely on a numerical pressure minimizer.
The reference's reported ratio near `0.7263479477` remains a numerical
description in its source note, not a newly certified optimum here.

## 5. Replay, scope and next question

The [algebra and orbit implementation](../tools/joint_energy_conjugations.py),
[checker](../tools/check_joint_energy_conjugations.py), and
[evidence](../evidence/joint-energy-conjugations.json) retain the complete
orbit, its 21 profile classes, unit counts, source hashes and exact orders.
The checker reconstructs the binary constraint spaces, common flags and
independent diagonal characters. It regenerates the complete compressed
orbit and verifies every crossing endpoint.

All **1,298** reference/candidate profiles are recomputed with independent
fraction-free polynomial arithmetic at integer weights. There is no
series truncation, floating-point pivot decision, or inference that an
unseen pivot is zero. Full generic rank follows from the strict crossing
endpoint. Every prefix inequality (6) is then checked in integers. The
complete focused replay uses only Python's standard library and takes
several minutes.

```sh
make verify-joint-energy-conjugations
python3 tools/check_joint_energy_conjugations.py --generate
make verify
```

The all-field statements are the preservation theorem, commutant dimensions
and orbit cardinality formula. The profile dominance catalogue covers
`F4`, this exact reference, these natural heights, and conjugations fixing
the complete leading weighted coefficient. New extension-field parameters,
different heights, moving that coefficient while retaining its rank,
nonconstant conjugations, arbitrary coefficient recombinations and other
seeds remain outside it. In particular this is not an actual-rank lower
bound on the coefficient span.

The next useful conjugation must move one of the first two coefficient
matrices, or leave this constant energy-block orbit. The construction
also gives a controlled comparison tool for stronger seeds: one whole
place can be preserved while a small matrix algebra organizes changes
at the other. This work used substantial OpenAI assistant support and
retains the Gaussian reference's source attribution. No new Lean work
or improved factoring exponent is claimed.

## 6. Formal rigidity when the leading rank may move

Let `K` be any field of characteristic two and let `epsilon` be a formal
parameter independent of `t`. Suppose `U(epsilon)` preserves the integer
energy blocks, lies in `GL_25(K[[epsilon]])`, and has constant term `I`.
Require only that the realignment of
`U(epsilon) B_-8 U(epsilon)^-1` has rank one over `K((epsilon))`.
Then there is a factorization

```math
\boxed{U(\epsilon)=
 (D_x(\epsilon)\otimes D_y(\epsilon))V(\epsilon),}        \tag{8}
```

where the `D` matrices are diagonal formal units with constant term `I`,
and `V` has constant term `I` and belongs to the leading centralizer `L`
at every order. These separate physical diagonal changes are constant
in `t` and preserve both weighted Smith profiles. This is a formal local
theorem, not a classification of every constant matrix retaining leading
rank one.

To prove it, let `C=realign(B_-8)`. It has rank one with `C_00=1`.
Write `C_rc=a_r b_c`, with `a_0=b_0=1`. For an infinitesimal energy-block
matrix `X`, put `delta C=realign([X,B_-8])`. The derivatives of the
minors `C_rc C_00-C_r0 C_0c` are

```math
\delta C_{rc}+b_c\delta C_{r0}+a_r\delta C_{0c}
 +a_rb_c\delta C_{00}=0.                               \tag{9}
```

The 69-dimensional energy-block domain has a 16-dimensional kernel in
(9). The subspace

```math
\mathfrak d=\{D_1\otimes I+I\otimes D_2:
 D_1,D_2\text{ diagonal}\}
```

has dimension nine. The leading centralizer has dimension ten, their
intersection has dimension three, and their sum is exactly the kernel:

```math
\ker(9)=\mathfrak d+\mathcal L.                         \tag{10}
```

The [tangent checker](../tools/check_joint_leading_tangent.py) constructs
all 625 equations and verifies this equality by exact binary elimination.
Its [certificate](../evidence/joint-leading-tangent.json) retains the
kernel basis and source hashes. Binary equations and bases extend (10)
to every coefficient field of characteristic two.

For the formal lifting, suppose the current residual is
`I+epsilon^k X+O(epsilon^(k+1))`. Rank one and the unit `(0,0)` pivot
imply (9), so decompose
`X=D_1 tensor I+I tensor D_2+L_k` using (10). Set

```math
D_k=(I+\epsilon^kD_1)\otimes(I+\epsilon^kD_2),\qquad
V_k=I+\epsilon^kL_k.
```

Left multiplication by `D_k^-1` and right multiplication by `V_k^-1`
remove the order-`k` residual. The leading rank-one property survives:
separate physical diagonal conjugation scales realignment rows and
columns, while right multiplication by a centralizer unit does not
change the conjugated leading matrix. The residual is therefore
`I+O(epsilon^(k+1))` with the same property. Induction and formal
convergence give (8); the `k`th correction starts at order `epsilon^k`,
so all coefficients of the infinite products are well-defined.

This is a structural statement, not a new effective factoring family over
an infinite field. Formal centralizer parameters and other constant
components of the rank-one condition are not covered by the F4 profile
catalogue. The precise local conclusion is that no additional formal
directions occur beyond the centralizer and separate physical diagonals.

```sh
make verify-joint-leading-tangent
```
