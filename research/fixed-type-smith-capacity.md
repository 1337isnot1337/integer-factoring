# Convex fixed-type Smith pressure and exact flag supports

## Result and purpose

For a fixed polynomial crossing and fixed physical paired types, the
two-place Smith pressure is convex in the positive valuation scale and
all physical gauge energies. Nested nonzero minor flags give explicit
global affine lower bounds. These provide both a numerical cutting-plane
search and an exact rational certificate format for a fixed-type exclusion.

This addresses the gauge optimization inside an actual crossing family.
The physical types remain fixed during each convex problem. Optimizing
the types themselves is a separate problem. No new factoring exponent
is established; the current complete score remains `1/2`.

The crossing transfer and paired-type entropy accounting are the existing
[paired-type theorem](paired-type-crossing-pressure.md). This note gives
the convex representation and global minor-flag supports used to solve
its fixed-type inner problem.

## 1. Weighted determinantal orders

Fix a matrix `M(z)` of generic rank `r` over a field's rational functions,
and fix the local valuation `v`, including the raw valuation at infinity.
Write the parameters as `xi=(tau,a,b)`, where `tau>0` scales valuations
and `a_i,b_j` are row and column energies. Its `k`th determinantal order is

```math
S_k(\xi)=\min_{\substack{|I|=|J|=k\\\det M[I,J]\ne0}}
 \left(\tau v(\det M[I,J])+\sum_{i\in I}a_i+\sum_{j\in J}b_j\right),
\qquad S_0=0.                                               \tag{1}
```

These are finite minima of affine functions. The Smith energies
`kappa_i=S_i-S_(i-1)` are sorted increasingly. For rational parameters,
clearing denominators reduces these facts to ordinary valuation-ring
elimination. The finite-minimum formulas extend them continuously to
real parameters.

Let `F` be any nested row/column flag of nonzero minors through size `r`.
Let `A_k^F(xi)` be the affine cost of its `k`th minor and let
`c_i^F=A_i^F-A_(i-1)^F`. Equation (1) gives `S_k<=A_k^F` for every `k`.
For any decreasing probability vector `p_1>=...>=p_r>=0`, summation by
parts yields

```math
\sum_i p_i\kappa_i
=\sum_{k=1}^r(p_k-p_{k+1})S_k
\le\sum_{k=1}^r(p_k-p_{k+1})A_k^F
=\sum_i p_i c_i^F,\qquad p_{r+1}=0.                       \tag{2}
```

The ordinary Gibbs variational inequality now proves

```math
\boxed{\quad
f(\xi):=\log\sum_i e^{-\kappa_i(\xi)}
\ge H(p)-\sum_i p_i c_i^F(\xi).
\quad}                                                     \tag{3}
```

Every expression on the right is affine in `xi`. At any chosen `xi`,
valuation pivoting supplies an optimal nested flag attaining all `S_k`.
Taking `p_i=e^(-kappa_i)/sum_j e^(-kappa_j)` makes (2)--(3) equalities.
Thus `f` is exactly the supremum of these affine functions, and is
convex. This argument also works for rectangular matrices at fixed
generic rank; the executable exact controls below use square matrices.

For a physical tensor, row and column energies are linear sums of the
individual leg energies. Composition with this linear map preserves
convexity. Exact grade-block decompositions are allowed: merge their
optimal flags in increasing pivot-energy order to obtain a global flag.

## 2. Fixed paired types and the full cost

Fix positive rational types `pi_X,pi_Y`, with entropies `H_X,H_Y` and
`H_X+H_Y>0`. Put `tau=log s>0`, and parameterize each physical energy
directly as `g_(place,leg)(letter)=tau*h_(place,leg)(letter)`.
Use the existing raw infinity convention, so no polynomial degree
is silently removed from the infinity orders.

For a target parameter `rho`, the complete pressure is

```math
\begin{aligned}
\mathcal P_\rho(\xi)
={}& f_\infty(\xi)+f_0(\xi)
 +\langle\pi_X,G_X\rangle+\langle\pi_Y,G_Y\rangle\\
 &-2\rho(H_X+H_Y),                                         \tag{4}
\end{aligned}
```

where `G_X` sums the `x,v` leg energies over both cuts, and `G_Y`
sums the `y,u` energies. Equation (4) is convex for fixed types.
The attained sufficient parameter is

```math
\rho_{\rm attained}=\rho+
              \frac{\mathcal P_\rho}{2(H_X+H_Y)}.            \tag{5}
```

A strict negative value at `rho=1/sqrt(2)`, followed by a rational
replay using the source transfer, would supply a post-half witness.
No such value is recorded here.

Independent paired types matter: the `x,v` distribution may differ
from the `y,u` distribution. A common-type search only covers the
special case `pi_X=pi_Y`.

## 3. Global dual certificates

Choose a minor flag and decreasing probability vector for each cut.
Adding their affine supports (3) to the linear and constant terms of
(4) gives a global affine lower bound on pressure. Convex combinations
of several such supports remain valid.

For a certificate over all physical energies, every energy coefficient
must vanish. Equivalently, the flag probability marginals on each leg
must equal that leg's prescribed type, allowing the stated convex
combination of flags. The remaining scale coefficient must be
nonnegative, since `tau>0`. If the entropy constant is also nonnegative,
the complete fixed-type pressure is nonnegative everywhere.

All these coefficient conditions are rational when the flag probabilities,
types, combination weights, and target `rho` are rational. The entropy
comparison is also exact: its exponential becomes a product of rational
numbers to integer powers after clearing denominators. A floating LP
solution needs this rational reconstruction before it becomes such a
mathematical certificate.

On a compact declared parameter box, affine supports give a lower
epigraph LP. The value at an evaluated point gives a primal upper value.
This supports a finite search with both values and a reported gap.
Floating LP arithmetic is recorded as a numerical estimate.

## 4. An exact sharp control

The [elementary binary pencil](linear-smith-barrier.md) has raw first-cut
infinity orders `(-1,0,0,1)` and second-cut zero orders `(0,1,1,2)`.
For each flag choose

```math
p=(9,3,3,1)/16.
```

The two scale coefficients cancel exactly. Their entropy sum is
`log(256/27)`. Thus (3) gives the global unweighted bound
`B(s)>=256/27`. Equality holds at `s=3`, where these are the actual
Gibbs probabilities. The checker verifies the entropy identity by
rational powers and checks the original crossing endpoint. This
recovers the known sharp bare bound through the flag method; it is
not a new factoring improvement or an all-gauge exclusion.

## 5. Implemented numerical scope

The [searcher](../tools/search_fixed_type_smith_capacity.py) applies (4)
to the [quantum-Pascal seed](quantum-pascal-trace-crossings.md). Its
fixed types are binomial distributions with supplied probabilities.
One energy per leg is fixed to zero; the remaining energies lie in
`[-2,2]`, and the valuation scale lies in `[0.002,1]`.

Each native profile returns its selected rows, columns, and exact
integer prefix valuations. These give the affine pivot increments
and the Gibbs-weighted support gradient. The code checks every
retained support at every newly evaluated point, and checks each new
support at all retained points. Missing promised pivots trigger a
larger series precision, up to the backend's explicit limit; an
unresolved pivot raises an error.

The recorded three-letter example uses `pi_X=pi_Y=Binomial(2,3/4)`.
Its numerical primal and box lower values are both near `rho=0.900932`,
with the saved gap. This is far from the post-half threshold. The
record says exactly which box and fixed types it covers. Other types,
larger energy ranges, and different seeds are separate searches.

## 6. Reproduction and proof status

```sh
make verify-smith-capacity
OPENBLAS_NUM_THREADS=1 python3 tools/search_fixed_type_smith_capacity.py \
  --alphabet 3 --probability-x .75 --probability-y .75 \
  --iterations 500 --output /tmp/fixed-type-capacity.json
```

The default exact checker needs only Python's standard library. It
checks 600 rational log-convexity inequalities, 600 global weighted
flag comparisons, 1,800 prefix-minor inequalities, and 1,800 rational
entropy supports, plus the sharp elementary control. Its inequalities
use fraction-free polynomial minors and rational arithmetic.

The optional numerical search requires NumPy, SciPy, and a C++17 compiler.
Its [record](../evidence/fixed-type-quantum-pascal-capacity.json) is a
bounded numerical experiment. The [exact controls](../evidence/smith-pressure-convexity.json)
support the separate all-size proof in Sections 1--3. Passing either
replay does not provide independent mathematical review.

The variational proof, controls, and implementation used substantial
OpenAI Codex assistance. The existing polynomial crossing, native Smith
backend, paired-type transfer, and elementary binary bound retain
their original attribution. No Lean source is changed.
