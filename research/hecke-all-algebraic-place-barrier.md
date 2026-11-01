# A determinant budget excludes every two-point Hecke shift

## Result and scope

For each retained [three-row Hecke crossing](three-row-hecke-crossings.md)
over F8, on five or nine letters, no choice of two distinct geometric
points over an algebraic closure improves the **unweighted two-place
supporting-line criterion**, at any positive supporting angle.
The balanced bare-pressure infimum is exactly one.

This includes moving the zero and infinity points by a fractional-linear
change of variable and extending the coefficient field. It does not
exclude physical-height gauges, joint tensor changes, or mixed-place
allocations. In particular, using a whole conjugate orbit simultaneously
in CRT reconstruction is a multiple-place grouping, not the two-geometric-
point problem proved here. No factoring-exponent improvement is claimed.

## 1. Bound unseen roots without factoring a determinant

Let `A(t)` be an `n` by `n` polynomial matrix over a finite field `F_Q`,
with entries of degree at most `D` and nonzero determinant. Homogenize
its determinant to degree `nD`. At a finite point use the usual local
Smith orders. At infinity use the orders of `u^D A(1/u)`.

The sum of all local order sums over the algebraic projective line,
with multiplicity, is exactly `nD`. Suppose all `Q+1` rational points
have been checked, and their order sums total `S`. The unaccounted
multiplicity is `R=nD-S`.

Every nonrational root has at least two distinct Frobenius conjugates,
all with the same multiplicity. Therefore every nonrational point has
local order sum at most `floor(R/2)`. Consequently

```math
U(A)=\max\left(\max_{p\in\mathbb P^1(F_Q)}\sum_j\nu_{p,j},
                    \left\lfloor\frac{R}{2}\right\rfloor\right)
```

bounds the order sum at **every algebraic point**, without a determinant
factorization. This is a bound on the sum of Smith orders at a point,
not on a single rank defect or on the degree of an irreducible factor.

## 2. Apply the bound to both cuts

The primitive Hecke tensors have no common scalar polynomial factor.
Their strict swap endpoint proves generic rank `n=q^2` in both cuts.
A fractional-linear substitution preserves primitiveness of the
homogeneous degree-`D` tensor. Its two local profiles are therefore those
at the selected geometric points, with the same reconstruction degree.

The complete rational-place calculations give:

```text
alphabet q     nD       remaining cut budgets    all-point sum bounds
    5          900          792, 542                  396, 332
    9         6561         6341,4842                 3170,2421
```

Thus `nD-U_0-U_1` is respectively 172 and 970, both positive.
The statement covers either assignment of the two cuts to the two points.

## 3. The entire unweighted supporting-angle family is excluded

Let `mu,nu` be the selected local order lists. For positive `a,b` with
`a^2+b^2=1`, write the uniform-type supporting cost as

```math
\Phi(s)=\frac{Ds+a\log\sum_j e^{-s\mu_j/a}
                    +b\log\sum_j e^{-s\nu_j/b}}{2\log q},
\qquad s\ge0.
```

The derivative numerator is `D-E_s(mu)-E_s(nu)`. Each tilted mean is
nonincreasing in its nonnegative tilt, since its derivative is minus a
variance. Hence

```math
\Phi'(s)\ge\frac{D-(U_0+U_1)/n}{2\log q}>0.
```

At zero, `Phi(0)=a+b>=1`. The supporting cost therefore cannot fall
below one at any positive scale or angle. At the balanced angle, the
corresponding bare `rho` has infimum one, attained as a scale limit.
This argument does not introduce or optimize physical-height parameters.

## 4. Exact replay and independent controls

Run the short five-letter check with:

```sh
python3 tools/check_hecke_all_algebraic_places.py --alphabet 5
```

Use `--alphabet 9` for the larger construction. Each command rebuilds the
literal tensor, checks its hash against the prior confluent-jet replay,
verifies scalar primitiveness, and computes both full profiles at all
eight F8 points and infinity.

A separate bit-plane Gaussian elimination checks every local determinant
order using the kernel dimension of a block-Toeplitz matrix. If `K` exceeds
all returned orders, that kernel dimension is exactly their sum.
This is independent of the scalar-restriction DVR elimination used for
the profiles. The original zero/infinity profiles are also compared with
the prior source-pinned receipts.

The [five-letter receipt](../evidence/hecke-q5-all-algebraic-places.json)
and [nine-letter receipt](../evidence/hecke-q9-all-algebraic-places.json)
retain all profiles, Toeplitz ranks, remaining budgets and source hashes.
[`finite_field_place_budget.py`](../tools/finite_field_place_budget.py)
contains the small degree-budget calculation and the independent rank
routine. The calculation uses exact integers and finite-field arithmetic;
no numerical optimization or finite list of extension fields is substituted
for the all-algebraic-point argument.
