# Exhaustive auxiliary lattice gauges and finite minimum certificates

## Result and scope

For auxiliary dimension two, rational changes of basis in the
[Smith transfer](auxiliary-gauge-smith-transfer.md) admit a small exact
search space at every fixed depth. Over a finite field with `Q` elements,
the number of distinct gauge classes at positive depth `h` is
`(Q+1)Q^(h-1)`. Over `F_2`, depths at most four, six, and eight contain
46, 190, and 766 classes respectively.

These classes include nonconstant rational changes of auxiliary basis.
They are more general than a constant basis change followed by diagonal
monomial weights. Each class has an explicit polynomial representative,
so its two local spectra can be computed exactly.

There is also a finite certificate that a gauge minimizes one Laurent
sum over **all depths**, for a fixed real `s>=1` and the chosen parameter
field. It suffices to compare the gauge with its `Q+1` neighbors in the
lattice tree. This certificate does not assert that descent terminates,
does not optimize `s`, and does not cover additional gauges introduced
by ramifying the polynomial parameter.

The results improve the finite search interface. They supply no crossing
family with an improved factoring exponent.

## 1. Which changes preserve the local Smith spectra

Let `F` be a finite field, `R=F[[z]]`, and `L=F((z))`. Fix a local
matrix tensor `K(z)` with two auxiliary states. Its entries may be
rational or Laurent series. In one chosen physical cut, both auxiliary
indices lie on the same half. Assume the rectangular flattening has full
generic column rank `n`, as it does in the crossing transfer.

A gauge `G in GL_2(F(z))` acts by

```math
K\longmapsto GKG^{-1}.                              \tag{1}
```

Multiplying `G` by a nonzero scalar does not change (1). Replacing it
by `UG`, with `U in GL_2(R)`, conjugates (1) by `U`. This is an
invertible row operation over `R` on the rectangular flattening: on its
auxiliary matrix entries it is the map `X -> UXU^{-1}`. Therefore every
Smith order is unchanged. The same argument works when the auxiliary
indices are placed on the opposite half and that flattening is
transposed.

Fixed physical monomial weights may be applied before or after this
operation. They commute with auxiliary conjugation. Thus the result also
applies inside the [balanced physical-gauge transfer](balanced-physical-gauge-transfer.md).

## 2. Polynomial representatives for every rational gauge

The rows of `G` generate an `R`-lattice in `L^2`. Scalar multiplication
normalizes it to an integral primitive lattice: it lies in `R^2` and
is not contained in `z R^2`. Row operations over `R` put its basis into
the unique form

```math
H_{a,b,p}(z)=
 \begin{pmatrix}z^a&p(z)\\0&z^b\end{pmatrix},
\quad a,b\ge0,\quad \deg p<b,\quad
\min(a,b,\operatorname{ord}_0p)=0.                  \tag{2}
```

Use `ord(0)=infinity`; when `b=0`, the only allowed `p` is zero.
To obtain (2), normalize the ideal of first coordinates to `z^a R`,
eliminate the other first coordinate, normalize the second diagonal
entry to `z^b`, and reduce the upper-right entry modulo `z^b`.
The last reduction leaves the displayed finite polynomial. Primitivity
is exactly the final minimum condition. The normalized first-coordinate
ideal, the second-coordinate kernel, and the reduced upper-right residue
also prove uniqueness.

All operations can be made over rational functions when `G` is rational:
the eliminated units and their inverses are rational, and reduction
modulo a finite power uses a finite Taylor expansion. Consequently (2)
represents every rational gauge exactly for the purpose of its local
Smith spectra. It is not merely a finite-jet approximation.

The depth is

```math
h=a+b.                                             \tag{3}
```

Indeed the primitive lattice has Smith orders `(0,h)` relative to
`R^2`. Its gauge and inverse have the finite descriptions

```math
H=\begin{pmatrix}z^a&p\\0&z^b\end{pmatrix},\qquad
H^{-1}=z^{-h}\begin{pmatrix}z^b&-p\\0&z^a\end{pmatrix}.
                                                        \tag{4}
```

Thus one can compute the polynomial numerator `H K adj(H)`, then restore
the common shift `-h` after the Smith calculation. At infinity, first
replace `K(z)` by `K(w^{-1})` and use the same construction in `w`.
The polynomial reversal shift must also be restored.

## 3. The exact finite counts

There is one class at depth zero. At positive depth `h`, divide (2)
into three cases:

- `a=0,b=h`: there are `Q^h` polynomials `p`.
- `b=0,a=h`: there is one choice, `p=0`.
- `a,b>0`: `p` must have nonzero constant coefficient. Summing over
  `1<=b<=h-1` gives `sum (Q-1)Q^(b-1)=Q^(h-1)-1` choices.

Their sum is

```math
N_h=(Q+1)Q^{h-1},\qquad
N_{\le H}=1+(Q+1)\frac{Q^H-1}{Q-1}.                \tag{5}
```

For fixed `s`, the two cuts can minimize their Laurent sums independently.
There is no need to enumerate all pairs of gauges: compute the profile
pool for each cut and select its smallest sum. The auxiliary half can
also be chosen independently. The gauges and their rational boundary
costs remain fixed certificate data in the transfer theorem.

## 4. The lattice tree and its neighbors

Vertices are primitive lattices up to scalar homothety. Two vertices
are adjacent if one lattice contains a representative of the other with
quotient of dimension one over `F`. The neighbors of a lattice with
row basis `G` are represented by

```math
\begin{pmatrix}1&c\\0&z\end{pmatrix}G\quad(c\in F),
\qquad
\begin{pmatrix}z&0\\0&1\end{pmatrix}G.              \tag{6}
```

Normalize each product to (2). These are the `Q+1` one-dimensional
quotients of the two-dimensional residue space. Every nonroot vertex
has one neighbor of depth `h-1` and `Q` of depth `h+1`; the root has
`Q+1` neighbors of depth one. Indeed a primitive positive-depth lattice
has one-dimensional image in `R^2/zR^2`. Exactly one of its index-`Q`
sublattices is contained in `zR^2`; division by `z` makes its depth
`h-1`. The others remain primitive and have depth `h+1`. This gives a
connected tree: a cycle would have a maximum-depth vertex with two
neighbors of smaller depth.

Any two lattices have a common diagonal basis by the elementary-divisor
theorem over `R`. The unique path joining them therefore lies in a
family

```math
G_t=\operatorname{diag}(z^t,1)S,\qquad t\in\mathbb Z,       \tag{7}
```

for one fixed invertible rational or Laurent-series matrix `S`, up to
the left units and scalar factors already shown to preserve spectra.
Such a path is the only lattice-tree geometry needed below.

## 5. Convexity of each Laurent sum on a path

Fix a physical cut, an expansion point, its fixed physical weights, and
`s>1`. Let `nu_1(G)<=...<=nu_n(G)` be the local Smith orders after the
auxiliary gauge. Put

```math
\sigma_i(G)=-\nu_i(G),\qquad
F_s(G)=\sum_{i=1}^n s^{\sigma_i(G)}.                 \tag{8}
```

Both the infinity sum `U(s)` and the zero sum `V(s)` in the transfer
are instances of (8).

Along (7), the auxiliary entry with indices `(a,b)` is multiplied by
`z^{t([a=0]-[b=0])}`. Since both auxiliary indices are on one side,
this simply scales rows of the rectangular flattening. Each selected
minor therefore gains **one** monomial factor, whose exponent is affine
in `t`. It follows that the smallest valuation among its nonzero
`k by k` minors is a minimum of affine functions. By the determinantal
description of Smith orders,

```math
\sum_{i=1}^k\sigma_i(G_t)
 =-\min_{\substack{I,J\\|I|=|J|=k}}
       \operatorname{ord}_0\det M_t[I,J]             \tag{9}
```

is convex in `t`, for every `k`. The spectra in (8) are listed in
decreasing order.

These partial-sum inequalities imply convexity of `F_s`. One direct
proof uses, for every real `u`,

```math
\sum_i(\sigma_i-u)_+
 =\max_{0\le k\le n}\left(\sum_{i=1}^k\sigma_i-ku\right).
```

The right side is convex in `t` by (9). With `c=log s>0`, integrate
against the positive weight `c^2 e^{cu}` and use

```math
e^{c x}=c^2\int_{-\infty}^{x}(x-u)e^{cu}\,du.
```

The result is convexity of `sum_i exp(c sigma_i)=F_s` along every
lattice path. When `s=1`, the function is the constant `n`.

## 6. A finite certificate for a global minimum at fixed s

Suppose a vertex `G` satisfies

```math
F_s(G)\le F_s(G')\quad\text{for all }Q+1\text{ neighbors }G'.
                                                        \tag{10}
```

For any other vertex, follow the unique path from `G`. Convexity makes
successive differences nondecreasing. Its first difference is
nonnegative by (10), so its value can never fall below `F_s(G)`.
Thus (10) certifies a global minimum over every rational auxiliary
gauge in the chosen parameter field, at this fixed `s` and physical
weight choice.

For rational `s`, all comparisons in (10) use exact rational arithmetic
on finite integer spectra. A descent routine may inspect the neighbors
and move to a strictly smaller value when one exists. We make no claim
that this routine always terminates or that the infimum is attained.

The certificate has two further limits. It does not by itself minimize
the product over different `s`. Also, substituting `z=w^m` changes the
parameter field used by the search: ramification introduces new vertices
and new branches at fractional distances in the original scale. A
neighbor certificate for `F((z))` does not exclude these additional
gauges over `F((w))`. One can apply the theorem again to each fixed
ramified polynomial tensor.

## 7. A common scalar projection can still make the search ineffective

Suppose a constant linear functional `ell` on auxiliary matrices gives
a scalar polynomial tensor `P=ell(K)` with `P(1)=E_q`. Suppose also that
both selected auxiliary gauges preserve this same projection:

```math
\ell(G_\infty K G_\infty^{-1})=\ell(K),\qquad
\ell(G_0 K G_0^{-1})=\ell(K).                        \tag{11}
```

Then the full-half Smith certificate is at least as large as the
certificate of this scalar pencil with the same physical weights.
Indeed, constant contraction of the auxiliary indices maps the full
rectangular matrix `M` to the scalar flattening `C M`. After clearing a
common Laurent pole, the row module of `CM` is contained in that of
`M`. The finite torsion module `R^n/row(M)` is consequently a quotient
of `R^n/row(CM)`. The maps on every layer `z^(j-1)M/z^j M` are
surjective, so the elementary-divisor lengths of the first module are
componentwise no larger than those of the second. Equivalently,

```math
\nu_i(\text{full})\le\nu_i(\text{scalar})
\quad\text{for every }i.                            \tag{12}
```

The same common pole shift cancels on both sides. Applying (12) at
infinity and zero shows `U_full>=U_scalar` and `V_full>=V_scalar`.
The physical mean corrections are identical, so they do not change the
comparison.

Two useful cases are immediate. A diagonal auxiliary gauge preserves
each fixed diagonal coordinate, including `K_00`. If `K(1)e_0=E_q e_0`,
that coordinate is a scalar crossing pencil. More generally, if
`K(1)=E_q I_D` and the characteristic does not divide `D`, the normalized
trace `ell=Tr/D` is preserved by **every** rational auxiliary gauge and
has endpoint `E_q`.

For external alphabet two, a scalar pencil produces position-invariant
coefficient tensors even after restriction to the middle Hamming layer.
The pinned [all-field binary-word and balanced-layer theorem](https://github.com/1337isnot1337/integer-factoring/blob/bb3bf763ed054f6fdaf882649c037ad810e16b97/research/all-field-binary-word-barrier.md)
then forces `rho>=3/4`. Thus a common projection as in (11) prevents
this certificate from reaching the half-rate threshold. Independent
nonconstant gauges can remove the coordinate-projection restriction;
an auxiliary dimension divisible by the characteristic avoids the
normalized-trace argument. Neither observation proves the existence
of an improving tensor.

## Verification and attribution

Run `python3 tools/check_auxiliary_lattice_gauges.py`. The checker
enumerates the canonical forms, checks their counts and lattice
memberships, verifies all neighbor relations, and compares Smith
profiles before and after nonconstant regular left changes of basis and
scalar homotheties. It checks selected spectra independently through
block-Toeplitz ranks, tests the path-convexity inequalities, and verifies
finite neighbor certificates against a larger exhaustive depth range.
Small additional controls compare the full spectra with their constant
coordinate and normalized-trace contractions.

The proof supplies the general statement; finite examples only check its
transcription. The construction and checker were developed with
substantial OpenAI assistant support. No Lean formalization or stronger
factoring exponent is claimed.
