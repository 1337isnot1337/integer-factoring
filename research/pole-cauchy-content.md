# Global factor counts from a Cauchy determinant slope

The [pole trace](interval-pole-trace.md) has a second exact global
interface: two determinant valuations of an implicitly given Cauchy matrix determine
the factor-pair count. Padding handles arbitrary rectangular intervals
without an exponential number of square queries. The padding and known
Vandermonde contents have polynomial binary-input cost. Computing the
two global determinant valuations at the desired cost remains open.
There is no improved factoring bound or new Lean formalization.

## 1. A square matrix with the exact factor incidence as its reduction

Let `A={a_1,...,a_m}`, `B={b_1,...,b_m}` be distinct positive integer
nodes. For odd `N`, choose `h` with `|N-a_i*b_j|<2^h` throughout the
matrix. Put `z_k=N+2^k*w`, where `w^2+w+1=0`, and take `k>h+1`.
Define

```math
(M_k)_{ij}=\frac{2^k w}{z_k-a_i b_j}.                \tag{1}
```

The pole-filter proof shows that every entry is 2-adically integral,
and its reduction modulo two is one exactly when `a_i*b_j=N`.
These ones occupy distinct rows and columns. Therefore the **actual**
rank of `M_k mod 2` over `F_4` is exactly the number `D` of factor pairs.
The matrix also satisfies the exact rank-one displacement relation

```math
\operatorname{diag}(a_i)M_k\operatorname{diag}(b_j)-z_kM_k
=-2^k w\,\mathbf1\mathbf1^T.                        \tag{2}
```

Over `Q_2(w)`, column scaling converts the unscaled kernel to a Cauchy
matrix with nodes `a_i` and `z_k/b_j`. These node families and the
resonant denominators must be retained in any use of a Cauchy algorithm.

## 2. The determinant slope gives the count

Write `V(A)=prod_(i<j)(a_j-a_i)` and similarly for `B`. The exact
Cauchy determinant identity, with its column factors canceled, is

```math
\det\left(\frac1{z-a_i b_j}\right)
=\frac{z^{m(m-1)/2}V(A)V(B)}{\prod_{i,j}(z-a_i b_j)}. \tag{3}
```

One obtains (3) by applying the ordinary Cauchy identity to `a_i,z/b_j`.
The powers of each `b_j` from the column scaling, Vandermonde, and
denominator cancel. This also fixes the positive sign in (3).

Let `E_k=v_2(det M_k)` and
`E_*=sum_(a_i*b_j!=N) v_2(N-a_i*b_j)`. All nonzero differences have
valuation below `h`, while a zero difference becomes `2^k*w` of
valuation `k`. Since `N` is odd, `v_2(z_k)=0`. Thus (3) gives

```math
E_k=k(m-D)+v_2(V(A))+v_2(V(B))-E_*,\qquad
\boxed{D=m-(E_{k+1}-E_k)}.                          \tag{4}
```

All determinants here are nonzero over `Q_2(w)`. Their contents can
nevertheless be enormous. The unknown constant `E_*` cancels from the
slope. Equation (4) is an exact reduction, not a way of evaluating
either content. Substituting a large unknown cross-product valuation
into (4) does not make that evaluation inexpensive.

## 3. Implicit padding preserves arbitrary rectangle queries

First clip both original factor intervals to `[1,N+1)`; this cannot
remove a factor pair. Let their lengths be `r,s` and set `m=max(r,s)`.
Append `m-r` nodes from `[N+1,N+1+m-r)` to the first domain when needed,
and likewise for the second. A new node exceeds `N`, so it participates
in no factor pair with a positive integer. Both domains now have size
`m` and at most two ordered disjoint intervals, and `D` is unchanged.
Zero-size matrices have determinant one and give count zero directly.

The padded endpoints have `O(n)` bits when `N` and the original endpoints
do. One may choose `h,k=O(n)`. Every exact entry of (1) has a binary
rational-coordinate description of polynomial length. The matrix still
has potentially exponentially many entries; it is supplied by its
interval endpoints and kernel, not expanded.

## 4. The known Vandermonde contents really are cheap

Set `H(t)=prod_(j=0)^(t-1) j!` and `S(t)=v_2(H(t))`. For one interval
of length `t`, its Vandermonde equals `H(t)` up to sign. Legendre's
formula gives

```math
S(t)=\sum_{q=2,4,8,\ldots<t}\sum_{j=0}^{t-1}\lfloor j/q\rfloor.
```

If `t=q*u+r`, `0<=r<q`, the inner sum is
`q*u*(u-1)/2+u*r`. There are `O(log(t+1))` terms, with polynomial-size
binary intermediates.

For two ordered disjoint intervals `[L,U)`, `[V,W)`, the cross product
of their differences is

```math
\prod_{a=L}^{U-1}\prod_{b=V}^{W-1}(b-a)
=\frac{H(W-L)H(V-U)}{H(W-U)H(V-L)}.                 \tag{5}
```

This follows by telescoping the factorials for each fixed `a`.
Consequently the cross content is four evaluations of `S`, and the
within-interval contents add two more. All normalization data in (3)
except `E_*` is therefore effective in polynomial bit cost for the
padded domains. The reference computes `S` for length `2^200+12` as
a 399-bit exponent, without constructing its Vandermonde product.

## 5. A normalized trace from two Cauchy powers

There is a direct scalar target that avoids the unknown constant `E_*`.
Use the indeterminate `t` and write `M(t)_(ij)=t*w/(N+t*w-a_i*b_j)`.
Let `M circle M` denote the entrywise square, while the inverse and
matrix multiplication below are ordinary matrix operations. Then

```math
t\,\partial_t M=M-M\mathbin{\circ}M,\qquad
M(2t)=\frac{2M(t)}{1+M(t)}\quad\text{entrywise}.       \tag{6}
```

Jacobi's determinant identity gives
`t*(partial_t det M)/(det M)=m-tr(M^(-1)*(M circle M))`.
Taking the same logarithmic derivative directly in (3) therefore proves
the exact identity

```math
\operatorname{tr}\bigl(M(t)^{-1}(M(t)\mathbin{\circ}M(t))\bigr)
=\sum_{i,j}M(t)_{ij}-\binom m2\frac{t w}{N+t w}.      \tag{7}
```

For `t=2^k`, with `k>=p+h`, the right side is integral and equals `D`
modulo `2^p`. Thus the required normalized scalar trace is integral,
even though the inverse matrix can be nonintegral. Computing a matrix
inverse modulo two would be invalid when `D<m`.

Equations (6)--(7) connect the pole filter directly to the first two
Cauchy Hadamard powers. They require only a short parameter update to
describe the next matrix; evaluating the global scalar remains open.
One can target this trace or the determinant slope directly, without
separately computing the potentially larger background contents.

## 6. Exact block elimination and the cross normalization

There is a closed block rule, but its normalization must be charged.
Partition both ordered node sets into prefixes `A_1,B_1` of the same
size `r` and remaining sets `A_2,B_2`. Set
`C_(a,b)=1/(z-a*b)`. Direct elimination gives

```math
(C_{22}-C_{21}C_{11}^{-1}C_{12})_{a,b}
=\frac{z^r\prod_{u\in A_1}(a-u)\prod_{v\in B_1}(b-v)}
 {\prod_{v\in B_1}(z-av)\prod_{u\in A_1}(z-ub)}
 \frac1{z-ab},\quad a\in A_2,b\in B_2.              \tag{8}
```

For one pivot, subtracting the rank-one update leaves the numerator
`z*(a-u)*(b-v)`. Iterating this identity proves (8); diagonal factors
on eliminated rows and columns cancel from the next Schur update.
For `M=c*C`, its Schur complement is `c` times (8).
Thus the kernel remains Cauchy, with diagonal row and column factors
described by products over the eliminated node intervals.

Let `D_ij` count factor pairs in `A_i cross B_j` and use `z=z_k`.
The determinant valuation of the row and column factors in (8) has
`k`-slope exactly `-D_21-D_12`: the static differences and `z^r` have
constant valuation, while each cross denominator contributes slope one
precisely at a genuine factor pair. Therefore the Schur determinant's
slope is

```math
(m-r)-D_{22}-D_{21}-D_{12}.                         \tag{9}
```

The eliminated block contributes `r-D_11`, recovering (4). These are
global cross-rectangle counts. Plain block elimination has moved their
evaluation into the diagonal normalization factors; dropping those
factors would lose the answer. Formula (8) provides a compact block
description, but no smaller bound on the remaining scalar queries.

The normalization can be nonintegral even when the initial `M` is
integral. Take `N=15`, both node sets `{3,4,5}`, and eliminate `{3,4}`.
The resulting one-by-one Schur complement is

```math
\frac{4z_k^2}{2^k w\,(z_k-25)(z_k-20)^2}.
```

For `k>=3` its valuation is `1-k`, with slope `-1`. The exact checker
obtains valuations `-6,-7` at `k=7,8`. The full determinant has slope
one and the eliminated determinant slope two. This is a control for
the required nonunit normalization, not a lower bound on every possible
shared-work algorithm.

## 7. The exact remaining computational target

Given odd `N`, the two padded interval domains, and `k=O(n)`, compute
the slope `E_(k+1)-E_k`, or the normalized scalar in (7) to the precision
needed to recover `D`. Include normalization, parameter construction,
and bit cost. Computing `E_k,E_(k+1)` separately is sufficient. Each content has `O(n)` output
bits: it is nonnegative and at most `k*m+O(n*m^2)` by (3).
An all-input query algorithm with cost `2^(O(n/(log n)^beta))` for a fixed
`beta>1/2` would transfer through the existing polynomial number of
factor-prefix queries to a stronger score. No such cost bound
is established here.

The [arithmetic Cauchy power formulas](arithmetic-cauchy-powers.md)
provide exact Smith orders in their stated unit-denominator,
arithmetic-node regime. The odd/odd part of (1) has resonant
denominators, and its Cauchy column nodes are reciprocal progressions.
Extending or replacing the normalization in this regime is the concrete
cross-lane question in
[Discussion #53](https://github.com/1337isnot1337/integer-factoring/discussions/53).

## 8. Exact finite controls

A [Discussion #53 contribution](https://github.com/1337isnot1337/integer-factoring/discussions/53#discussioncomment-19324116)
independently supplied the same Cauchy determinant normalization and
reported exact rational checks at sizes one through four. Its section
invariant and counterexample are retained with credit in the pole-trace
note. That check concerns the determinant formula, rather than independent
review of the entire global evaluation target.

Run `make verify-interval-pole-trace`. The
[dense checker](../tools/check_pole_cauchy_slope.py) performs exact
Gaussian elimination in `Q(w)`, compares the whole determinant with
(3), and checks the slope against literal factor pairs after padding.
It also solves for the inverse/Hadamard-square matrix over `Q(w)`,
compares the exact scalar trace in (7), and checks the entrywise flow
in (6). Literal block solves check (8) and the nonintegral example.
Explicit two-pair, square, empty, and padded unequal cases are
included. It independently checks the one- and two-interval Vandermonde
contents by literal differences. The
[receipt](../evidence/pole-cauchy-slope.json) records finite controls;
these support the ordinary proof and do not establish a fast evaluator.
