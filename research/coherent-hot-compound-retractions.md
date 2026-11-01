# Coherent inverse systems with complete occupation-block profiles

## Result and scope

There is an effective class of physical reflection restrictions built from
arbitrary homogeneous spaces closed under bit deletion. Both natural cut
profiles reduce to compressed compound matrices. The two compressions are
different and are computed separately. A Boolean derivative inverse system
gives such a physical space together with a compatible coordinate retraction.
This permits coherent relations at several physical heights.

The resulting 32-letter example strictly improves the natural coordinate
pressure at every scale and has an exact rational certificate
`rho < 1781/2501`. It is also optimal at every scale within the complete
121-map binary occupation-preserving boundary class described below.

The implementation constructs each occupation block directly. It does not
construct a four-leg tensor before compressing it. The exact controls below
reproduce the existing eight-letter coordinate and nonorthogonal profiles,
and two new six- and seven-letter inverse-system examples. This is a
construction and computation result. No factoring exponent above one half
is established.

The energy-boundary specialization is already proved in
[the preceding coherent-boundary note](https://github.com/1337isnot1337/integer-factoring/blob/8d3a2da11eca5e1a7d064b645c698719ac141098/research/coherent-hot-boundary-crossings.md).
That predecessor also has a complete rank-one `F4` search through sixteen
letters. The present contribution is the deletion-invariant space interface,
the inverse-system construction, and direct sparse native block evaluation;
the boundary argument itself is not claimed as new. The literal control on
current main is the [eight-letter retraction](eight-letter-reflection-retraction.md).

## 1. Physical spaces and the two exact cuts

Work over a field of characteristic two. Let the occupation basis on `s`
physical modes be `e_A`, for subsets `A` of `{0,...,s-1}`. Give it height

\[
h(A)=\sum_{a\in A}(2a+1).
\]

For each mode, the deletion operator is

\[
D_a e_A=[a\in A]e_{A\setminus\{a\}}.
\]

Let `P` be an embedding with `q` columns and `L` a retraction with `q`
rows. Assume:

1. `LP=I_q`;
2. each column `P_j` and row `L_j` is supported at one common height `h_j`;
3. the column space of `P` and the row space of `L`, interpreted in the
   occupation basis, are invariant under every `D_a`.

For occupation-block evaluation also require each pair `P_j,L_j` to have
one common occupation number `n_j`. The algebraic reduction before this
last block split does not require that additional grading.

Apply `L,L,P,P` to the four physical roles `x,y,u,v` of a reflection
crossing. The crossing endpoint is preserved because `LP=I`. Choose a
sufficiently large reflection interval for the finite support. In its
deflated variable `t=w^2`, the regular zero cut has one-particle cross block

\[
B(t)=\begin{pmatrix}D(t)&J(t)\\J(t)&D(t)\end{pmatrix},
\quad
D_{ab}(t)=(1+t)t^{2(a\mathbin{|}b)+1},\qquad
J_{ab}(t)=[a=b]t^{2a+1}.
\tag{1}
\]

Write `E_B(R,S)=det B[R,S]` when the two occupation masks have the same
size, and zero otherwise. In the regular Gaussian factorization, the
within-row and within-column factors have entries

\[
U(R,R')=[R'\subseteq R]\mathrm{Pf}A[R\setminus R'],
\quad
V(S',S)=[S'\subseteq S]\mathrm{Pf}F[S\setminus S'].
\tag{2}
\]

These are polynomials in bit deletions with identity diagonal. The deletion
coefficients depend only on the removed bits. Invariance in condition 3
therefore gives exact induced maps

\[
(L\otimes L)U=U_L(L\otimes L),\qquad
V(P\otimes P)=(P\otimes P)V_P.
\tag{3}
\]

Every nonidentity deletion strictly lowers height. The induced factors in
(3) are triangular with identity diagonal, so their determinants are one.
Consequently the compressed zero cut has exactly the Smith profile of

\[
(L\otimes L)E_B(P\otimes P).
\tag{4}
\]

The natural infinity kernel exchanges the entire `y,u` roles. Its regular
compound matrix is instead compressed by

\[
(L\otimes P^{\mathsf T})E_B(L^{\mathsf T}\otimes P).
\tag{5}
\]

The same invariance proof applies to both factors in (5). Its profile need
not equal the profile of (4). All four map spaces retain their own deletion
condition; a retraction identity alone is insufficient.

The ambient scalar monomial shifts at zero and infinity are opposite.
Removing a common scalar content changes those shifts by opposite amounts
again. They cancel in the product of the two counting functions. This
statement is about that product; it does not identify the raw polynomial
degree or scalar content of a compressed tensor.

When every label has occupation `n_j`, (4) and (5) both split by
`n_x+n_y`. If `a_k` counts physical labels of occupation `k`, the block
orders are

\[
N_r=\sum_k a_k a_{r-k},\qquad \sum_rN_r=q^2.
\tag{6}
\]

The implementation orders each column pair oppositely to its row pair.
Since `B(1)` exchanges its two mode sets, every compressed block is the
identity at `t=1`. This proves its full generic rank separately from the
series precision used for Smith elimination.

## 2. Derivative inverse systems give compatible maps

Identify an occupation vector with a multilinear polynomial in variables
`x_0,...,x_(s-1)`, where each variable has height `2a+1`. Let
`f_1,...,f_r` be nonzero forms, each homogeneous in height and ordinary
degree. Different top forms may have different heights and degrees.
Define

\[
\mathcal U=\mathrm{span}\{D_Af_j:A\subseteq\{0,...,s-1\},\ 1\le j\le r\}.
\tag{7}
\]

The deletions commute and square to zero, so `U` is deletion invariant.
Every derivative is homogeneous in both gradings.

Order monomials by their binary occupation masks, with the largest mask
leading. Take a reduced echelon basis of `U`. Let `C` be its leading masks,
and normalize the associated basis vectors `p_c` to have coefficient one
at `c` and zero at every other leading mask.

**Lemma.** `C` is a bit-deletion downset.

If `c` leads a polynomial `p` and `a` belongs to `c`, deleting `a` preserves
the numerical order among all monomials containing `a`: it subtracts the
same `2^a` from each such mask. Terms not containing `a` disappear.
Therefore `D_a p` has leading mask `c` with `a` removed. This derivative
is in `U`, so its leading mask is one of the echelon pivots. This proves
the lemma.

Let the columns of `P` be the `p_c`, and let the rows of `L` be coordinate
evaluation at the masks in `C`. Reduced echelon form gives `LP=I`.
The column space is `U`, while the row space is the coordinate downset
span on `C`; both are deletion invariant. Echelon elimination only combines
vectors sharing a monomial and therefore sharing height and occupation.
Thus each paired label has the two required gradings.

This is an effective construction of every hypothesis in section 1.
It can correlate several energy levels. It is not confined to a single
top boundary of an energy ball.

For example, on four modes the form

\[
f=x_1x_2+x_0x_3
\]

has height eight and degree two. Its inverse system has basis the vacuum,
all four single modes, and `f`, giving a six-letter code. Coordinate
projection onto the leading word of `f` gives its retraction. Adding the
top monomial `x_0x_2` gives a seven-letter example. Both complete physical
cuts are checked against the literal ambient Pfaffian tensor below.

### 2.1 Homogeneous differential equations give another construction

Let `T_j` be finite linear combinations of products of the bit derivatives.
Require each `T_j` to be homogeneous in height and differential degree.
Inside the coordinate space of words with height at most `H`, define

\[
\mathcal U_H=\bigcap_j\ker T_j.
\]

The cutoff space is deletion closed. All bit derivatives commute with
every `T_j`, so their common kernel is also deletion closed. The leading
word argument above again supplies a coordinate retraction. Unlike a
single top form, one differential equation can impose coherent relations
at many successive physical heights while retaining every unconstrained
lower-height state.

The [constructor](../tools/coherent_differential_kernels.py) solves the
finite homogeneous linear systems separately in each height and occupation
cell. It checks the equations on the complete resulting basis, the exact
retraction identity, and both deletion-closure conditions. For example,
`T=D_1 D_2 + D_0 D_3` imposes a height-eight relation on four modes and all
of its higher-height consequences. A separate
[screen](../tools/search_coherent_differential_kernels.py) retains every
completed pair of profiles and its numerical pressure value.

The [kernel checker](../tools/check_coherent_differential_kernels.py)
enumerated 2,355,232 ambient vectors in 272 complete small spaces and
reconstructed every one of their 512,284 kernel vectors. It also compared
both cuts of the ten-letter four-mode kernel with the literal Pfaffian
tensor in both map orientations, checking 20,000 endpoint entries and
20,000 infinity-restoration entries. The
[retained controls](../evidence/coherent-differential-kernel-controls.json)
record the exact counts and source hashes.

## 3. Direct block construction and cost accounting

The [native constructor](../tools/coherent_hot_compounds.cpp) forms each
compressed entry as the complete sum of selected polynomial determinants.
It uses sparse occupation supports for the physical maps. No unproved
low-rank completion is supplied to the evaluator.

The selected determinant degree is bounded by the sum of its row degrees,
and the wrapper also uses the conservative bound `r*max_degree(B)`.
The series window is strictly larger than that complete coefficient bound.
Every block checks its polynomial endpoint and independently checks the
sparse-map endpoint. The Smith engine must then resolve all `N_r` pivots.
An insufficient window causes a retry or an explicit failure, never a
partial profile reported as complete.

The binary implementation supports at most sixteen physical modes, block
order 4,096, and precision 8,192. These are implementation limits. For
fixed finite input data the mathematical construction and complete exact
rank factors remain effective. More generally, enumerating the derivative
spaces, complete matrices, and finite-field eliminations costs `2^O(s)`
operations for an explicitly listed `2^O(s)` collection of top forms.
The polynomial coefficient degrees are `O(s^2)` for the mode interval used
in (1). The bounds charge all matrices and their complete coefficients.

### 3.1 A tighter exact coefficient window

The hot kernel has an additional monomial factorization. For mode indices
`a,b` on either side, write

\[
B_{a,b}(t)=t^{a+b+1}K_{a,b}(t),
\]

where `K_ab=(1+t)t^(a xor b)` on the same side and `K_ab=[a=b]`
on opposite sides. For a selected `k`-minor with row and column heights
`h(R),h(S)`, every determinant term contains the common monomial

\[
t^e,\qquad e=\frac{h(R)+h(S)}2.
\]

Since `a xor b <= a+b`, every determinant term of `K` has degree at most
`e`. Thus the complete original `B`-minor has degree at most
`h(R)+h(S)`. Height homogeneity makes this factor common to every sparse
term of a compressed entry as well.

The [specialized constructor](../tools/coherent_hot_compounds_fast.py)
stores the common monomial as a separate series valuation. Its native
determinant calculation uses the proved degree bound on `K`; it retains
the original retry rule until every generic pivot is found. This improves
the starting precision bound for high-rank height cutoffs. For height-25
physical maps the complete coefficient bound is 100, regardless of the
number of occupied modes. The Smith workspace can still need a larger
window, which is charged and recorded separately.

The [comparison checker](../tools/check_coherent_hot_fast.py) reproduced
985,992 complete polynomial entries and all 4,240 Smith orders in its
four code controls, including both displayed 32-letter maps. These are
exact comparisons with the original generic constructor. The
[retained control](../evidence/coherent-hot-fast-controls.json) also records
both series windows and timings; it does not claim a uniform speed factor.

For profiles `a_i,b_i` in `t`, put

\[
Z_0(u)=\sum_i e^{-ua_i},\qquad Z_\infty(u)=\sum_i e^{-ub_i}.
\]

With one common physical type, a trial `rho` has Gibbs weights proportional
to `e^{-u h_j/rho}`. The complete balanced pressure test is

\[
\log Z_0(u)+\log Z_\infty(u)
 <4\rho\log\sum_j e^{-u h_j/\rho}.
\tag{8}
\]

Both cuts and all physical entropy remain in (8). Numerical root searches
are only discovery tools. A new factoring claim would require a strict
exact certificate below `rho=1/sqrt(2)` and the established effective
crossing transfer. No such certificate is supplied here.

## 4. Independent exact controls

Run

```sh
python3 tools/check_coherent_hot_compounds.py \
  --output evidence/coherent-hot-compound-controls.json
```

The [checker](../tools/check_coherent_hot_compounds.py) records:

| Control | Count |
|---|---:|
| Complete literal Pfaffian cut-profile comparisons | 8 |
| Literal compressed endpoint entries | 11,889 |
| Literal infinity restoration entries | 11,889 |
| Independent complete coefficients by permutation determinants | 280 |
| Independent quotient-ring rank checks | 152 |
| Native block endpoint entries | 7,206 |
| Complete small inverse-system forms | 258 |
| Literal derivative reconstructions | 65,328 |
| Rejected invalid map inputs | 4 |

The small inverse systems include all three nonzero forms in the
four-mode, degree-two, height-eight cell and all 255 forms in the
eight-mode, degree-four, height-32 cell. These are map-construction checks;
they do not claim a pressure search unless the separate complete-profile
artifact records one. The literal six- and seven-letter comparisons
exercise coherent spaces that omit some lower-energy coordinate states.

The [control artifact](../evidence/coherent-hot-compound-controls.json)
retains the maps, both full profiles, and source hashes. A second agent
independently checked the map orientations and inverse-system leading-term
argument, and replayed the six-letter literal cuts in both transposed
orientations. No new Lean work was used. Implementation and analysis used
substantial OpenAI assistant assistance.

## 5. Reproducible larger searches

The [boundary driver](../tools/search_coherent_hot_retractions.py) enumerates
all 120 binary rank-three idempotents on the four two-particle height-16
words, with the four-particle height-16 word retained. It also includes the
other occupation-preserving choice, retaining all four two-particle words
and omitting the four-particle word. Every case has 32 physical letters.
Its separate 96-letter mode makes bounded projector-preserving walks in a
ten-dimensional height-25, three-particle boundary.

The [inverse-system driver](../tools/search_coherent_inverse_systems.py)
enumerates top forms in a specified homogeneous monomial cell. Every saved
completed trial includes both complete profiles. Each driver retains its
in-progress input and distinguishes a time limit from exhaustion of its
declared finite list.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_coherent_hot_retractions.py \
  --family q32 --seconds 1800 \
  --output evidence/coherent-hot-q32-retractions.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_coherent_inverse_systems.py \
  --modes 8 --particles 4 --height 32 --seconds 7200 --max-trials 255 \
  --output evidence/coherent-inverse-8-4-32.json
```

Finite numerical minima are not all-scale optimality results or stronger
factoring bounds. The complete 32-letter catalogue has 121 successful
trials and no unresolved cases. The next section gives an exact result
from its best numerical candidate.

## 6. A 32-letter strict improvement at every scale

Retain all 28 coordinate words with height below 16, and retain the
four-particle word 15. On the remaining two-particle height-16 cell, order
the coordinate words as

\[
(24,36,66,129).
\]

Equivalently, these are the mode pairs `{3,4}`, `{2,5}`, `{1,6}`, and
`{0,7}`. Write their basis vectors as `e_0,e_1,e_2,e_3`. The three new
embedding columns and retraction rows are

\[
P_0=e_0+e_1,\quad P_1=e_0+e_2,\quad P_2=e_0+e_3,
\]

\[
L_0=e_0^{\mathsf T}+e_2^{\mathsf T}+e_3^{\mathsf T},\quad
L_1=e_0^{\mathsf T}+e_1^{\mathsf T}+e_3^{\mathsf T},\quad
L_2=e_3^{\mathsf T}.
\]

They satisfy `LP=I_3` in characteristic two. Every deletion lands in the
retained lower-height space, so section 1 applies. The reference coordinate
map keeps `e_0,e_1,e_2`. Both maps have exactly the same 32 physical heights.

Let `C_0,C_infinity` be the two reference counting polynomials in the
deflated variable `r`, and let `Z,I` be the corresponding coherent counting
polynomials. Recomputing all four complete Smith profiles gives the exact
integer polynomial identity

\[
C_0(r)C_\infty(r)-Z(r)I(r)=r^{32}(1-r)Q(r),
\tag{9}
\]

where `Q` has degree 147, all 148 coefficients are nonnegative,
`Q(0)=1`, and `Q(1)=552960`. The full coefficient list is retained in the
[exact certificate](../evidence/coherent-q32-pressure.json). Thus (9) is
strictly positive for **every** `0<r<1`. Because the physical height
multisets agree, this proves a strict two-cut pressure improvement at
every nonzero finite scale. It is an all-scale comparison of these two
maps, not a global optimum claim over other maps.

There is also an exact rational rate certificate. Set

\[
\rho=\frac{1781}{2501},\qquad
r=\left(\frac{5405}{5406}\right)^{1781},\qquad
A(x)=\sum_{j=1}^{32}x^{h_j}.
\]

The [checker](../tools/check_coherent_q32_pressure.py) uses outward integer
intervals to prove

\[
\bigl(Z(r)I(r)\bigr)^{2501}
<A\!\left(\left(\frac{5405}{5406}\right)^{2501}\right)^{7124}.
\tag{10}
\]

The reference coordinate map fails (10) at this same rate and scale.
Every physical Gibbs weight in this witness is rational. The checker
recomputes 4,096 Smith orders, checks (9) coefficient by coefficient, and
uses 14 independent rational arithmetic controls. Certification uses no
floating arithmetic.

For orientation, the discovery optimization gave the following values:

| 32-letter map | Numerical minimum `rho` | Numerical scale `u=-log(r)` |
|---|---:|---:|
| Reference coordinate map | 0.712119959170746 | 0.329728540764358 |
| Displayed coherent map | 0.712115012500715 | 0.329469478855912 |

The strict certified bound is `rho < 1781/2501 = 0.712115153938...`.
It remains above `1/sqrt(2)`, so the factoring exponent is unchanged.

Reproduce the exact result with

```sh
python3 tools/check_coherent_q32_pressure.py \
  --output evidence/coherent-q32-pressure.json
```

### 6.1 Exact optimality within the 121-map binary boundary class

The displayed map has the smallest product of the two counting polynomials
at every `0<r<1` among all 121 maps in section 5. This includes every
binary rank-three idempotent on the four two-particle boundary states
with word 15 retained, and the remaining occupation-preserving choice
that retains all four two-particle states and omits word 15.

For each candidate, subtract the displayed map's product polynomial from
the candidate's product polynomial. There are two exact zero differences:
the displayed map and its transposed projector, with kernel and covector
exchanged. Every other difference is strictly positive on `0<r<1`.

The [catalogue checker](../tools/check_coherent_q32_catalogue.py) proves
the latter claim by factoring each difference as `r^a (1-r)^b Q(r)`.
In 82 cases `Q` has nonnegative power coefficients. In the other 37 cases
it has an exact expansion

\[
Q(r)=\sum_{i=0}^{d}c_i r^i(1-r)^{d-i},\qquad c_i\ge0,
\]

with at least one positive coefficient. It reconstructs both polynomial
identities exactly. All 121 maps have the same physical height multiset,
so these are comparisons of the complete balanced pressure costs.

The [completed certificate](../evidence/coherent-q32-catalogue-optimum.json)
also recomputes all 242 complete cut profiles, comprising 247,808 Smith
orders, and repeats 61,952 deletion-closure checks. It uses no numerical
pressure fit to decide a comparison. This proves optimality within the
stated binary class for every scale, while leaving extension fields,
other heights, and other physical map classes open. In particular, the
[mixed-occupation extension](coherent-hot-mixed-boundary.md) has 375
additional binary maps outside this 121-map class.

```sh
python3 tools/check_coherent_q32_catalogue.py \
  --output evidence/coherent-q32-catalogue-optimum.json
```
