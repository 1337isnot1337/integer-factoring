# Ordered physical gauges and exact finite rank-sum controls

Matrix-valued physical gauges turn any scalar crossing path into an ordered
auxiliary construction with the exact endpoint preserved. The matrices
can fail to commute, so the resulting word tensor need not have simultaneous
site symmetry. This supplies a uniform construction, not a favorable
uniform rank bound. A complete small binary control rejects a fixed
two-dimensional gauge family under the finite-block rank-sum envelope.
The complete factoring exponent remains `1/2`.

## 1. The endpoint-preserving matrix gauge

Let `T(t;x,y,u,v)` be any scalar polynomial crossing with
`T(1;x,y,u,v)=[x=v][y=u]`. Choose invertible constant `D` by `D` matrices
`F_x,G_y` over its coefficient field, and a boundary `f^T e=1`. Define

```math
\mathcal K(t;x,y,u,v)
 =T(t;x,y,u,v)F_xG_yG_u^{-1}F_v^{-1},                     \tag{1}
```

and retain the physical order in

```math
H_m(t)=f^T\mathcal K(t;x_1,y_1,u_1,v_1)\cdots
                \mathcal K(t;x_m,y_m,u_m,v_m)e.           \tag{2}
```

At the crossing endpoint, `u=y,v=x`, so the four auxiliary factors in
each local matrix cancel to the identity. Otherwise the scalar factor
is zero. Consequently

```math
\mathcal K(1;x,y,u,v)=[x=v][y=u]I_D,\qquad H_m(1)=E_m.   \tag{3}
```

The cancellation happens inside each local factor. It does not require
commuting the factors at different positions or choosing a free boundary
after seeing the physical input. Thus the coefficients of (2) sum to the
full crossing with known scale one, at every length.

For fixed alphabet, auxiliary dimension, field and local degree `L`, the
complete degree is at most `Lm` and there are at most `Lm+1` monomial
pieces. Naively constructing all physical entries, polynomial coefficients,
and exact matrix rank factors takes `2^O(m)` operations over a fixed finite
field. For fixed rational tables, clearing denominators and charging exact
arithmetic still gives `2^O(m)` bit cost; the scalar crossing scale changes
by the known nonzero denominator factor. This fits the existing finite-seed
table-construction interface. Growing auxiliary dimensions require a new
cost calculation and are not covered by that fixed-data statement.

The missing mathematical step is an upper bound on the actual rank sum
of useful coefficient pieces, uniform in length. A finite-field endpoint
also needs its own exact-count or modular factoring transfer; it is not an
automatic stronger finite-alphabet factoring algorithm. The
[auxiliary coefficient-space analysis](auxiliary-coefficient-span.md)
explains that distinction and the finite-block envelope used below.

## 2. A complete binary two-dimensional control

Use the [three-letter Steiner cubic](steiner-cubic-crossings.md) over
`F_2`, auxiliary dimension two, and fix

```math
F_0=I,\quad F_1=\begin{pmatrix}1&1\\0&1\end{pmatrix},
\quad F_2=\begin{pmatrix}1&0\\1&1\end{pmatrix}.           \tag{4}
```

Normalize `G_0=I`, vary `G_1,G_2` over all six members of `GL_2(F_2)`,
and vary all six nonzero binary boundary pairs with `f^T e=1`. This gives
exactly 216 cases. It is complete for that stated normalization and fixed
first table; it is not a classification of all matrix gauges or boundaries
in other dimensions or fields.

For each two-position case, enumerate all 128 binary functionals on the
seven coefficient positions, including dependencies. The optimization
is the total cost

```math
\min\sum_j\sqrt{r(A_j)c(A_j)},\qquad\sum_j A_j=E_2,       \tag{5}
```

over its complete coefficient span. The best retained interval is

```math
\frac{65709677}{1048576}\le\operatorname{opt}
 \le\frac{65709678}{1048576}.                            \tag{6}
```

A positive decomposition attaining `9+24 sqrt(5)` lies in that interval;
it uses coefficient masks 126 and 1, with rank pairs `(36,80)` and `(81,1)`.
The corresponding selected three-position continuation has interval

```math
\frac{558364386}{1048576}\le\operatorname{opt}
 \le\frac{558364387}{1048576},                           \tag{7}
```

and a positive path costing `27+78 sqrt(42)`, about `532.4978`.
The three-position case is selected, not a complete continuation of all
216 cases. Both displayed cases use `G=(I,7,I)` in row-major binary mask
encoding and `e=f=(1,0)^T`.

Every retained lower bound, including all 216 two-position cases, exceeds
`3^(3m/2)`. Hence these fixed-block decompositions fail even the weaker
`rho=3/4` rank-sum envelope, and therefore cannot beat the half-rate target
through that envelope. The much smaller one-position minmax threshold
`16` for the scalar Steiner span is not a rank-sum improvement.

## 3. Integer dual certificates

Let `w` be a coefficient mask and let `A(w)` be the corresponding tensor.
On the XOR group of masks, an edge from `v` to `v+w` represents adding
that tensor. Its weight is `sqrt(r(A(w))c(A(w)))`. Paths to the all-ones
mask are exactly coefficient-span decompositions of the endpoint. Zero
tensors and dependencies are included as zero-weight edges.

With `M=2^20`, integer square roots give rigorous lower weights
`floor(M sqrt(rc))` and upper weights `ceil(M sqrt(rc))`. The evidence
contains an integer potential `phi`, with `phi(0)=0`, such that

```math
\phi(v+w)-\phi(v)\le\lfloor M\sqrt{r(A(w))c(A(w))}\rfloor
\quad\text{for every }v,w.                              \tag{8}
```

Telescoping (8) gives the lower bound in (6) or (7), without trusting the
shortest-path search. The upper path gives its own exact XOR identity and
rounded upper cost. The checker recomputes all rank pairs with the opposite
pivot order and verifies **4,587,520** integer dual inequalities. It replays
57,344 cut-rank entries including the zero tensors. The comparison to
`3^(3m/2)` is made by squaring integer numerators, with no floating-point
decision.

## 4. Replay and next algebra

The [constructor](../tools/ordered_steiner_gauges.py) retains the exact
physical order, validates every polynomial endpoint, and uses the stated
boundary. The [checker](../tools/check_ordered_steiner_gauges.py) and
[evidence](../evidence/ordered-steiner-gauges.json) retain all tables,
boundaries, rank pairs, potentials, positive paths and source hashes.

```sh
make verify-ordered-steiner-gauges
python3 tools/check_ordered_steiner_gauges.py --generate
make verify
```

The finite result leaves other first tables, larger auxiliaries, extension
fields, length-dependent constructions, and cancellations beyond the fixed
block envelope open. A concrete next algebra uses clock and shift matrices
with `YX=zeta XY`. Their ordered word weights have an explicit quadratic
phase and global charge constraints. That representation may support a
uniform rank calculation through finite-field Fourier transforms; no
favorable rank bound for it is asserted here. No new Lean work was used.
