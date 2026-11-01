# Iterating exact Smith fibers over F4

## 1. Scope

The complete-entry construction in
[finite-jet Smith fibers](finite-jet-smith-fibers.md) used binary reference
frames. Its output tensors can have coefficients in F4. That restriction
prevented a new output from becoming the reference for the next exact
linear synthesis step.

The constructor
[`f4_two_cut_smith_fibers.py`](../tools/f4_two_cut_smith_fibers.py)
removes this restriction. It constructs regular frames over F4, retains
every polynomial coefficient of every physical entry, and solves the
resulting affine equations over F4. The source can be any explicit strict
crossing with coefficients in F4. The construction supports both row
divisibility constraints and finite Hasse equations.

This is an exact construction method. It does not establish a new
factoring exponent, a better crossing pressure, or completeness of the
chosen regular-frame family.

## 2. Complete polynomial coordinates

Let the four physical alphabets have the same fixed size q. Let T(w) have
degree at most d and the literal crossing endpoint

\[
 T_{xyuv}(1)=\mathbf 1_{x=v}\mathbf 1_{y=u}.
\]

Every degree-at-most-d tensor with this endpoint has a unique expression

\[
 S_e(w)=T_e(w)+(w+1)f_e(w),\qquad \deg f_e<d,
\]

where e denotes one of the q^4 entries. There are q^4 d independent F4
coordinates before the frame equations are imposed. No support,
parity, trace, derivative, or tensor-product constraint is implicit in
these coordinates.

Choose integer physical heights separately for the row-zero and
column-infinity cuts. Add a common scalar power so that each weighted
matrix has polynomial entries. At zero this scalar is

\[
 c_0=-\sum_{\ell=1}^4\min_i h_{\ell,i};
\]

at infinity, after putting t=1/w, it is

\[
 c_\infty=d-\sum_{\ell=1}^4\min_i h_{\ell,i}.
\]

Write a weighted matrix as M_T(t). The regular-frame constructor returns
polynomial matrices U(t),V(t) whose constant matrices are invertible and
such that U M_T V is diagonal through a declared finite precision. The
precision is larger than every row target used below. The complete
independent Smith oracle supplies the source orders used to choose this
precision. See
[`f4_regular_smith_frames.py`](../tools/f4_regular_smith_frames.py).

For a target s_a on row a, impose

\[
 \operatorname{coeff}_{t^j}(U M_S V)_{ab}=0
 \quad\text{for every }b\text{ and }0\le j<s_a.
\]

These are affine equations over F4 in the full coefficient vector of f.
Their coefficients are computed exactly. At zero a coefficient f_{e,k}
contributes a shifted copy of (1+t)U_{a,i}V_{j,b}. At infinity its
contribution is

\[
 t^{c_\infty+h(e)-1-k}(1+t)U_{a,i}V_{j,b}.
\]

The exponent before the frame product is nonnegative because
k<=d-1. Therefore truncating that product at the row target does not
discard a coefficient needed by the equations. This also explains the
one-unit shift in the infinity coefficient extractor.

The two stored bitplanes encode coefficients in F4. They are not two
independent binary systems. All pivot normalization, elimination, affine
right-hand-side updates, and back-substitution use multiplication in F4.

## 3. What the frame constraints prove

The matrices U and V are invertible over the local power-series ring.
Thus they preserve every determinantal ideal of the weighted matrix.
If row a is divisible by t^{s_a}, a minor on row set I is divisible by
t^{sum_{a in I}s_a}. Writing s_(1)<=...<=s_(q^2), every k-th
determinantal order therefore satisfies

\[
 \Gamma_k(M_S)\ge\sum_{a=1}^k s_{(a)}.
\]

Subtracting the common scalar shift from each channel gives the
corresponding bound for the original weighted cut. This is a sufficient
bound, not a claim that the candidate has the same complete Smith
profile as the reference. The checker and search compute the candidate's
complete profile separately.

To reuse a common set of floors, the constructor first solves those
equations and chooses free coordinates in their affine solution space.
It reduces every higher row equation into those same coordinates. A
target increase is then a further exact affine elimination. Sampling
means choosing the remaining free F4 coordinates and back-substituting
through both systems. It does not mean rounding a numerical optimizer.

Finite Hasse equations use the same free coordinates. Since the
variation is (w+1)f, the functional on coefficient f_{e,k} is

\[
 D^{[j]}w^k(a)+D^{[j]}w^{k+1}(a).
\]

The new F4 implementation reduces this functional through its F4
affine base pivots, including the affine right-hand side. Thus it can
combine finite interpolation with the same exact regular-frame bounds.

All steps have finite polynomial cost in the explicit coefficient
window, matrix size, and equation count. Iteration requires no new
mathematical hypothesis: after an output passes the endpoint and degree
checks, its own regular frames can be constructed and used as a new
reference. The resulting finite search still explores only the frames,
heights, and targets actually selected.

## 4. Independent controls and search boundary

The scoped checker is
[`check_f4_two_cut_smith_fibers.py`](../tools/check_f4_two_cut_smith_fibers.py).
It compares every coefficient equation with the established binary
constructor on the original reflection seed, including the common
fiber dimension. It then uses both nonbinary primitive tensors in
[`finite-fiber-basis-input.json`](../evidence/finite-fiber-basis-input.json)
as references. For each one it constructs complete frames, solves at
two target levels, imposes additional finite Hasse conditions, and
replays the full transformed polynomial matrices and complete native
Smith profiles.

```sh
python tools/check_f4_two_cut_smith_fibers.py \
  --output evidence/f4-two-cut-smith-fiber-controls.json
```

The search driver is
[`search_iterated_f4_smith_fibers.py`](../tools/search_iterated_f4_smith_fibers.py).
It records the exact input and source hashes, all retained full
coefficient witnesses, actual trial counts, and its stopping reason.
Its pressure values are numerical. A common scalar factor, when
present, is divided exactly and normalized at the endpoint before a
new pressure is evaluated. The original frame certificate is replayed
on the undivided candidate; the normalized tensor gets a fresh
complete-profile computation. Scalar division is not treated as a
free preservation of the old row targets.

Neither a completed finite search nor failure to find a new pressure
proves an all-family lower bound. Any numerical crossing of the
factoring threshold still needs a separate exact certificate and the
full all-width transfer argument.

## 5. Exact parameter deflation before another iteration

Both primitive reference tensors above have only even polynomial
exponents. They can be written exactly as T(w)=R(w^2), where R has
degree 10 or 11. The helper
[`deflate_crossing_parameter.py`](../tools/deflate_crossing_parameter.py)
checks every coefficient and the literal endpoint before returning R.

This substitution is a computational simplification. It does not by
itself improve the two-place pressure. More generally, if T(w)=R(w^r)
and the physical heights for T are r times those for R, every weighted
minor order at zero or infinity is multiplied by r. Thus every complete
Smith order is multiplied by r. The degree, physical type means, and
inverse local scale change by the same factor. Choosing
s_T=s_R/r makes all terms in the two-place pressure identical. Here s
denotes the logarithmic Chernoff scale, so a channel of order nu
contributes exp(-s nu).

The exact coefficient inputs and this native pressure control are
recorded in
[`deflated-finite-fiber-input.json`](../evidence/deflated-finite-fiber-input.json).
They are reproduced by
[`prepare_deflated_fiber_inputs.py`](../tools/prepare_deflated_fiber_inputs.py):

```sh
python tools/prepare_deflated_fiber_inputs.py \
  evidence/finite-fiber-basis-input.json \
  --output evidence/deflated-finite-fiber-input.json
```

After deflation, an exact frame iteration is free to use every
coefficient in the smaller degree window, including odd powers of the
new parameter. Those changes create new tensors; their costs require
fresh complete profiles. The deflation identity alone supplies no
bound for these new outputs.
