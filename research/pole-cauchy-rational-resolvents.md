# Normalized rational contractions of implicit pole Cauchy inverses

## Proved scope

The [polynomial circuit moment theorem](pole-cauchy-circuit-moments.md)
compresses polynomial weights, even when their degrees are exponentially
large. Its final obstruction concerns entrywise polynomial approximations
to the singular factor filter. This note takes a different route: it
computes rational weights directly, with no polynomial interpolation on
the original nodes.

Let `A,B` be equally sized, nonempty sets of distinct positive odd
integers, given as disjoint unions of consecutive odd progressions. Let
`r_A,r_B` be their respective progression counts, `r=r_A+r_B`, and
`m=|A|=|B|`. Put

\[
 z=N+2^k w,\qquad t=2^k w,\qquad w^2+w+1=0,\qquad
 M(a,b)=\frac{t}{z-ab},
 \tag{1}
\]

where `N` is positive and odd, and `k>=1`. Choose positive integer
queries `xi,eta` with `xi` outside `A` and `eta` outside `B`. Define

\[
 \mathcal B(\xi,\eta)
 =\sum_{a\in A,\,b\in B}
       \frac{(M^{-1})_{b,a}}{(\xi-a)(\eta-b)}.
 \tag{2}
\]

If an odd query occurs, require `k` to exceed all valuations of its
integer distances from its domain. The implementation uses the sufficient
condition

\[
 k\ge h:=\max\left(
 \{\operatorname{bitlen}|\xi-a|:a\in A,\ \xi\text{ odd}\}
 \cup
 \{\operatorname{bitlen}|\eta-b|:b\in B,\ \eta\text{ odd}\}
 \cup\{0\}\right).
 \tag{3}
\]

The maximum in (3) needs only progression endpoints. There is no additional
condition on the original differences `N-ab`; the identity and algorithm
also apply when the original rectangle is not separated by `2^k`.

**Theorem.** For every requested absolute precision `p>=1`, (2) can be
computed with error of two-adic valuation at least `p`, in deterministic
bit time polynomial in `r`, the binary input length, `k`, and `p`. The
numerical domain size is never enumerated. The output is a binary scale
`s` and a quadratic residue `U` modulo `2^P`, with

\[
 \mathcal B(\xi,\eta)=2^sU+O(2^p),\qquad
 P=p-s\le p+(r+2)k.
 \tag{4}
\]

The scale is a certified lower bound on valuation, not necessarily its
exact value. Cancellation is allowed. For factoring applications, the
usual choice `k=O(input length)` can satisfy both (3) and separation of
the original rectangle.

This is a complete contraction over both matrix axes, including weights
whose integer denominators are even. A polynomial number of such
contractions and binary-encoded rational coefficients with polynomial
denominator valuations can therefore be combined at polynomial cost.
Constructing weights whose final scalar equals the required
Hadamard-square trace remains open. No factoring exponent improves here.

## 1. An exact four-product identity

Set `x_a=z/a` and

\[
 P_X(Y)=\prod_{a\in A}(Y-x_a),\qquad
 P_B(Y)=\prod_{b\in B}(Y-b).
\]

Let `u=M^{-1}(1/(xi-a))_(a in A)` and
`F(Y)=sum_b u(b)/(Y-b)`. The matrix equation gives

\[
 F(x_a)=\frac{z}{t(\xi x_a-z)}.
\]

The rational function with these values, denominator `P_B`, and
numerator of degree below `m` is

\[
 F(Y)=\frac{z}{t(\xi Y-z)}
 \left(1-\frac{P_X(Y)P_B(z/\xi)}
                  {P_B(Y)P_X(z/\xi)}\right).
 \tag{5}
\]

Indeed, the apparent pole at `Y=z/xi` cancels. After that cancellation,
the numerator has degree at most `m-1`, and substitution at every `x_a`
gives the required value. Uniqueness follows from interpolation at `m`
distinct points. None of the relevant denominators vanishes: `xi` is
outside `A`, `eta` is outside `B`, and `z` is nonrational.

Evaluating (5) at `eta` and cancelling the elementary scale factors in
`P_X` yields

\[
 \boxed{\displaystyle
 \mathcal B(\xi,\eta)
 =\frac{z}{t(\xi\eta-z)}(1-K),\qquad
 K=
 \frac{\prod_{a\in A}(z-a\eta)\prod_{b\in B}(z-\xi b)}
 {z^m\prod_{a\in A}(\xi-a)\prod_{b\in B}(\eta-b)}.
 }
 \tag{6}
\]

Each of the four products is affine along every input progression.
Equation (6) involves no product over the full Cartesian rectangle.

## 2. Why normalization costs only polynomial precision

Let `n_A(c,l)` count nodes of `A` congruent to `c` modulo `2^l`, and
define `n_B` similarly. For odd `eta`, unramifiedness gives

\[
 v_2(z-a\eta)=\min(v_2(N-a\eta),k),
 \qquad
 \sum_{a\in A}v_2(z-a\eta)
 =\sum_{\ell=1}^k n_A(N\eta^{-1},\ell).
 \tag{7}
\]

The convention `v_2(0)=infinity` includes exact zeros of the rational
part. If `eta` is even, all these factors are units, so their valuation
is zero. An analogous statement holds for the other numerator product.

An even `xi` makes every `xi-a` a unit. For odd `xi`, condition (3)
ensures

\[
 \sum_{a\in A}v_2(\xi-a)
 =\sum_{\ell=1}^k n_A(\xi,\ell).
 \tag{8}
\]

Consequently `E=v_2(K)` is exactly the sum, over `1<=l<=k`, of

\[
 \mathbf1_{\eta\text{ odd}}
       [n_A(N\eta^{-1},\ell)-n_B(\eta,\ell)]
 +\mathbf1_{\xi\text{ odd}}
       [n_B(N\xi^{-1},\ell)-n_A(\xi,\ell)].
 \tag{9}
\]

All counts are obtained by binary division and arithmetic progressions.
For one odd progression, counts of two odd residues at a fixed level
differ by at most one. Within one domain, they therefore differ by at
most its number of progressions. Each domain count also differs from
`m/2^(l-1)` by less than its progression count.

If both queries are odd, regroup (9) by domain and use the first bound.
If exactly one query is odd, compare its two domain counts using the
second bound and `|A|=|B|`. If neither is odd, (9) is zero. At level one,
each nonzero comparison is exactly zero. Thus

\[
 \boxed{|E|\le r(k-1).}
 \tag{10}
\]

Large numerator and denominator contents can separately have
exponentially large valuations. Their difference in (10) is small,
and their exact valuations are retained as binary integers. This is
the precision step needed to make the rational identity algorithmic.

Write `K=2^E U` with `U` a unit, and set

\[
 d=v_2(\xi\eta-z)
   =\min(v_2(\xi\eta-N),k),\qquad
 V=2^{-d}(\xi\eta-z).
\]

Here `0<=d<=k`, including `xi*eta=N`, when `d=k`. Then (6) becomes

\[
 \mathcal B=2^s\frac{z}{wV}
 \begin{cases}
 1-2^E U,&E\ge0,\\
 2^{-E}-U,&E<0,
 \end{cases}
 \qquad s=\min(0,E)-k-d.
 \tag{11}
\]

Every division on the right is by a unit. Equations (10) and (11) give
`s>=-(r+2)k`, so computing all units modulo `2^(p-s)` proves (4).
The bracket may vanish to any order; no assumption about cancellation
or about the scalar being nonzero is used.

### The height condition cannot simply be dropped

Let `m=2^s>=2`, take both domains to be `1,3,...,2m-1`, choose both
queries `2m+1`, and set `k=1`. Every numerator factor in (6) has
valuation one. Each integer denominator product has valuation

\[
 v_2(2^m m!)=2m-1.
\]

Thus `E=2-2m<0`, `d=1`, and the exact scalar has valuation `-2m`.
A representation of the form (4) needs exponentially many residue bits
for fixed absolute precision. The scalar identity remains true, but a
polynomial precision promise requires a condition such as (3).

## 3. Affine products without a separation assumption

The earlier [inverse-entry construction](pole-cauchy-conditioning.md)
normalizes affine products when each nonzero rational part has
valuation below the pole displacement. For (6), one must also handle

\[
 \prod_{j=0}^{L-1}(c+dj+2^k w)
 \tag{12}
\]

when some nonzero `c+dj` is divisible by `2^k`.

Put `g=v_2(gcd(c,d))`. If `g>=k`, factor out `2^k` from every term.
The remaining progression has coefficient one on `w`; both its even
and odd rational parts are units. Split an odd progression step into
the two index parity classes, which gives even steps.

If `g<k`, first factor out `2^g`. An even remaining step has odd
constant part and already gives a unit progression. For an odd remaining
step, split its indices into the valuation classes
`v_2(c'+d'j)=v<k-g`, plus the single congruence class
`c'+d'j=0 mod 2^(k-g)`. Each low class is an affine unit progression
after division by its exact power of two. The last class becomes a
unit progression with coefficient one on `w`; split its index parity
if necessary. There are at most `k+2` pieces.

The [modular unit progression algorithm](modular-interval-jet-arithmetic.md)
evaluates every piece in polynomial time in the written endpoints and
the unit precision. Its Newton calculations use precision bounded by a
polynomial, rather than expanding the full product. Integer denominator
products in (6) use the existing nonzero affine-product routine. Powers
such as `z^m` use binary powering.

These constructions, together with (9), form a uniform algorithm for
(2). The proof charges progression decomposition, residue counts,
valuation arithmetic, every unit inversion, and absolute precision.

## 4. Implementation and checks

The implementation is
[pole_cauchy_rational_resolvents.py](../tools/pole_cauchy_rational_resolvents.py).
Run

```sh
python3 tools/check_pole_cauchy_rational_resolvents.py \
  --output evidence/pole-cauchy-rational-resolvents.json
```

The deterministic checker uses exact rational arithmetic in `Q(w)`:

- 320 affine-product comparisons against literal products, including
  458 nonzero rational parts at or above the displacement valuation and
  10 zero rational parts;
- 72 complete inverse contractions against exact matrix inversion,
  with 18 cases in each query parity combination, 14 negative ratio
  valuations, and 10 cases `xi*eta=N`;
- 71 of those 72 matrices deliberately lack separation of the original
  rectangle, which the new theorem does not require;
- a `2^100+3` dimensional control at 16 absolute bits, independently
  reduced to a literal `3` by `3` inverse using complete unit cycles;
- a `2^60+3` dimensional execution with odd rational poles at eight
  absolute bits, explicitly marked as lacking an independent inverse
  reference at that dimension;
- exact missing-height counterexamples at dimensions `2,4,8`, and five
  rejected malformed inputs.

For the independent large control, both queries are even. All factors
are congruent to one modulo two. Such units have order dividing
`2^(P-1)` modulo `2^P`; the odd-node residue cycle also has length
`2^(P-1)`. Since `100>=2P-2`, the leading `2^100` nodes contribute
identity to each normalized product and to `z^m`. The last three nodes
then give the same residue as the small exact inverse. This explanation
is separate from the interval-product implementation.

## 5. Remaining scalar target

The accessible weights now include

\[
 H(a,b)=\sum_{j=1}^R
       \frac{c_j}{(\xi_j-a)(\eta_j-b)}
 \tag{13}
\]

for polynomial `R` and polynomial input and denominator precision.
Each term can be evaluated by (6), and the final values can be summed
with a common certified scale. This works directly for rational
weights; it does not use the polynomial interpolation excluded by the
finite-difference obstruction.

To evaluate the factoring trace, it would suffice to construct (13)
such that

\[
 \sum_{a,b}(M^{-1})_{b,a}H(a,b)
 \equiv \sum_{a,b}(M^{-1})_{b,a}M(a,b)^2
 \pmod {2^p}.
 \tag{14}
\]

Equation (14) asks only for equality of the final scalar. No such
construction is proved here. Requiring equality at every matrix entry
would impose a substantially stronger task and is not needed by this
target. The exact identity (6), its polynomial normalization bound,
and the implementation are the established progress.

## Provenance

This is an AI-assisted ordinary proof and implementation, developed from
the repository's interpolation formula and modular unit-product
arithmetic at base `8c4729454e26b107c7eab0b524c1b55f2cf4a843`, together
with the current circuit-moment checkpoint. It does not add a Lean
formalization or claim independent mathematical review. The factoring
exponent remains `1/2`.
