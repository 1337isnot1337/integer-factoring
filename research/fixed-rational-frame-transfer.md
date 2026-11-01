# Fixed rational frames have polynomial coefficient-transfer cost

## Transfer lemma

Fix a field `F`, a number `d>=1` of variables
`Z=(x,y_1,...,y_(d-1))`, and a polynomial matrix `T(Z)`. Suppose a fixed
rational-function factorization is supplied:

```math
T(Z)=L(Z)\operatorname{diag}(x^{a_1},\ldots,x^{a_q})R(Z),
\qquad a_i\ge0.                                           \tag{1}
```

Every entry of `L,R` is assumed regular at `x=0` over
`F(y_1,...,y_(d-1))`. The frames, dimensions, exponents, and rational
functions are fixed independently of the tensor length `m`. The same
argument permits any fixed finite family of such data.

For every fixed constant `B`, there is a constant `C_B` such that all
integer target exponents with `|beta_j|<=Bm` satisfy

```math
\operatorname{rank}_F[Z^\beta]T(Z)^{\otimes m}
\le C_B(m+1)^d
   \#\left\{(i_1,\ldots,i_m)\in[q]^m:
                         \sum_j a_{i_j}\le\beta_x\right\}. \tag{2}
```

Each separable coefficient factor in this bound has a uniform
polynomial-time entry algorithm in `m`, measured in field operations.
Over a fixed finite field this is a polynomial bit-cost bound. The
matrix itself has exponential dimension; the assertion concerns its
rank certificate and coefficient-entry queries.

Thus finite Laurent frames are sufficient but not necessary for the
coefficient-splitting step in the
[Smith-profile transfer](smith-profile-crossing-transfer.md).
Fixed rational frames also suffice. This lemma supplies no new local
Smith spectrum or factoring exponent by itself. It does not cover
arbitrary formal-series frames.

## 1. One positive Laurent cone for the finite family

Represent each nonzero frame entry as `P/Q`, where `P,Q` are ordinary
polynomials and `Q(0,y)` is nonzero. Such a representation exists after
canceling common powers of `x`, because the entry is `x`-adically regular.

Choose positive integer weights on the `y` variables that assign
different weights to all distinct monomials in the finite list of
`Q(0,y)` parts. For example, use successive powers of an integer larger than
every coordinate exponent in those denominators. Then choose the
positive weight of `x` larger than every `y`-weight occurring in the
denominator list. For every `Q`, its unique minimum-weight monomial
`c Z^v` lies in `Q(0,y)`.

Write

```math
Q=cZ^v(1+H),\qquad
Q^{-1}=c^{-1}Z^{-v}\sum_{j\ge0}(-H)^j.                   \tag{3}
```

Every exponent in `H` has strictly positive weight and a nonnegative
`x` coordinate. Negative coordinates in the other variables are
permitted. Each numerator term contributes a fixed shift `u-v`, whose
`x` coordinate is also nonnegative.

Let `G` be the union of these finitely many denominator increments,
and let `C` be their nonnegative integer span. For fixed finite sets of
shifts `S_L,S_R`, every frame-entry support is contained in
`S_L+C` or `S_R+C`. All nonzero generators have positive weight. Hence
only finitely many generator combinations occur below any weight bound;
the expansions and all their products are well-defined formal Laurent
series. No analytic convergence is used.

## 2. Split only the two complete word factors

If either outside frame is zero, (2) is immediate. Otherwise ignore
zero entries when forming the shift sets. For a channel word
`i=(i_1,...,i_m)`, let `a(i)=sum_j a_(i_j)`.
Its contribution to the tensor power of (1) is

```math
x^{a(i)}\,\ell_i(Z)r_i(Z),
```

where `ell_i` is a column vector and `r_i` is a row vector. Their
supports lie in `mS_L+C` and `mS_R+C`, respectively. Their `x`
exponents are nonnegative. Consequently this word contributes nothing
when `a(i)>beta_x`.

Otherwise put `eta=beta-a(i)e_x`. Its coefficient is

```math
\sum_{\gamma+\delta=\eta}
              [Z^\gamma]\ell_i(Z)\,[Z^\delta]r_i(Z).       \tag{4}
```

Each summand has rank at most one. There are polynomially many possible
`gamma`, even though the frame expansions can be infinite.

To see this, let `b_L,b_R` be the minimum weights of the two finite
shift sets. A contributing split has `w(delta)>=mb_R`, hence

```math
w(\gamma)\le w(\eta)-mb_R
             \le Bm\sum_jw_j-mb_R.
```

Write `gamma` as a sum of `m` left shifts and nonnegative generator
multiples. If `mu` is the least positive generator weight, their total
multiplicity is at most

```math
\frac{m\max\{0,B\sum_jw_j-b_L-b_R\}}\mu.                 \tag{5}
```

If there are no generators, only the `m` shifts remain. All shifts and
generators are fixed, so each coordinate of `gamma` is `O(m)`. There
are at most `C_B(m+1)^d` integer exponent vectors in that box. Summing
the rank-one terms in (4) over eligible words proves (2).

This split is made after forming each entire word factor. Splitting
independently at all `m` positions would give an unnecessary exponential
overhead.

## 3. Uniform coefficient queries

A requested factor entry is a product of `m` functions from the fixed
rational list. Use (3) for each function. The numerator shifts have
total weight bounded below by a fixed constant times `m`. The target
and the coefficient splits in Section 2 have weight `O(m)`. Therefore
only denominator increments of total weight `O(m)` are needed.

Their total multiplicity is `O(m)`, and every coordinate lies in an
`O(m)` box. A dictionary of at most `O(m)^d` coefficients therefore
suffices. Repeated polynomial multiplication, geometric-series
expansion, and truncation by the positive weight compute the requested
coefficient with polynomially many field operations. Fixed leading
coefficients `c` are the only scalars inverted. This construction works
in every characteristic; it does not divide by factorials or by
variable binomial indices.

The same argument permits a fixed-dimensional ordered auxiliary matrix
product on one side of the cut. Multiply its coefficient tables in the
original order. Their supports obey the same bound, and the auxiliary
dimension is fixed. Both auxiliary indices must remain on that side,
as in the existing ordered Smith-profile transfer.

## 4. Application to generic local Smith calculations

For a fixed polynomial matrix, a generic `x`-adic Smith calculation over
`F(y)[[x]]` can be realized by rational frames. One construction takes
ordinary polynomial Smith form over the PID `F(y)[x]`, then writes each
nonzero diagonal polynomial as `x^a u(x)` with `u(0)!=0`. Absorb its
unit factor into an outside frame. The resulting frame entries are
rational functions of all variables and regular in `x`.

This is a finite preprocessing step independent of tensor length. A
supplied rational certificate can also be checked by ordinary rational
identities. The lemma then applies to its local profile. At infinity,
first make the usual reciprocal change of variable and clear the
fixed polynomial degree.

Different cuts may use different positive cones. Their final
coefficient matrices are the same intrinsic coefficients of the
polynomial tensor; their separate rank factorizations need not use the
same intermediate Laurent expansion.

The fixed-family qualification is essential to this proof. It supplies
the uniform denominator supports, generator sizes, positive weights,
and coefficient algorithm. None of these conclusions follows for an
arbitrary formal series or for a new rational family whose description
grows with the tensor length without further bounds. The number `d`
must also be fixed for the overhead to be polynomial.

## 5. Exact controls

The checker uses `q=y-z-x` and the polynomial matrix

```math
T=\begin{pmatrix}q&1\\x&0\end{pmatrix}
 =\begin{pmatrix}1&0\\x/q&1\end{pmatrix}
   \begin{pmatrix}1&0\\0&x\end{pmatrix}
   \begin{pmatrix}q&1\\0&-1/q\end{pmatrix}.                 \tag{6}
```

Its local Smith orders are `(0,1)`. The rational frames have infinite
Laurent expansions. The common weights `(4,1,2)` also handle
`1/(x+yz)`, whose denominator increment is `(1,-1,-1)` and has weight
one. Closed-form coefficients of powers of this latter function give
independent scalar controls.

Over `F_2,F_5,F_7`, the checker reconstructs 150 tensor-power coefficient
matrices from their rank-one coefficient splits and compares 10,500
entries with direct polynomial tensor powers. It includes negative
auxiliary-variable exponents, where the rational frame terms must
cancel to give a zero polynomial coefficient. It also checks 75 scalar
coefficients, the support bounds, the nonnegative `x` increments, and
rejection of a denominator with an uncanceled `x` pole.

Run:

```sh
python3 tools/check_fixed_rational_frame_transfer.py --output /tmp/fixed-rational-frame-transfer.json
cmp /tmp/fixed-rational-frame-transfer.json evidence/fixed-rational-frame-transfer.json
```

The implementation is [rational_frame_cones.py](../tools/rational_frame_cones.py),
with its [checker](../tools/check_fixed_rational_frame_transfer.py) and
[receipt](../evidence/fixed-rational-frame-transfer.json). No Lean
formalization is included.
