# One pure correction cannot beat 82/113 at any tensor length

## Theorem and scope

Let `T_5` be the retained five-letter crossing. Let `k>=1`, and give
`T_5` to the `k`th tensor power its additive natural physical heights.
Over any field `K` of characteristic two, consider

\[
 T_5^{\otimes k}+\lambda(z)\,
 a(z)\otimes b(z)\otimes c(z)\otimes d(z),
 \qquad a,b,c,d\in K(z)^{5^k}.
\]

Assume the resulting tensor is a strict polynomial crossing after a
permitted common Laurent monomial translation. There is no degree bound
or entry-floor restriction on the correction. Each of its four vectors
may correlate all `k` copy indices arbitrarily.

The fixed-natural-gauge paired-type rate satisfies

\[
 \boxed{\rho\ge\frac{82}{113}
           =0.7256637168\ldots>\frac1{\sqrt2}.}
\]

This holds for every positive tensor length, including one copy. It
strengthens the earlier [two-copy 5/7 floor](pure-two-copy-flag-floor.md).
It does not cover sums of pure corrections, different physical heights,
or arbitrary changes of the scalar seed. It is an obstruction to this
construction and pressure calculation, not a factoring lower bound.
The complete factoring score remains `alpha=1/2`.

## 1. The reference profiles and the geometric flag

At the two local places the one-copy reference minima are `-8` and
`8`. After subtracting those minima, both complete partitions equal

\[
\begin{aligned}
 P(s)={}&1+s+s^3+s^4+s^5+s^6+s^7+2s^8+2s^9\\
       &+2s^{10}+2s^{11}+3s^{12}+2s^{13}+s^{14}
          +2s^{15}+s^{16}+s^{20}.
\end{aligned}
\]

The physical energy partition is

\[
 Q(s)=1+s^2+s^6+s^8+s^{10}.
\]

Let `U` be the eleven-dimensional constant image flag through relative
order nine, reshaped as five-by-five physical matrices. The two image
flags and two coimage flags agree. The
[geometric rank-one certificate](geometric-low-flag-obstruction.md)
proves that `U` contains no nonzero rank-one matrix over any field of
characteristic two. Its proof consists of 25 ordinary polynomial
identities over `F_2`; it does not infer the assertion over arbitrary
fields from a finite-field enumeration.

The reference profiles and these constant flag identities are preserved
under extension of the binary coefficient field. The displayed nonzero
binary minors remain nonzero in every such extension.

## 2. The obstruction survives every tensor length

Use tensor-product regular frames for the reference. Every relative
one-copy order is nonnegative. Consequently the `k`-copy constant image
flag through relative order nine is contained in

\[
 U^{\otimes k}.
\]

Group all row-copy coordinates together and all column-copy coordinates
together. This larger tensor space still contains no nonzero rank-one
physical matrix. Indeed, suppose it contained `x y^T`. Choose values for
the other `k-1` row coordinates that leave a nonzero slice of `x`, and
independently choose the other column coordinates to leave a nonzero
slice of `y`. Contracting those coordinates leaves a nonzero rank-one
five-by-five matrix in `U`, a contradiction.

The same argument applies to both local image spaces. It permits
arbitrary correlations inside `x` and `y`; no product structure across
the `k` copies is assumed for the correction vectors.

## 3. A partition bound that retains all higher jets

The regular-row argument from the
[pure-correction theorem](pure-two-copy-flag-floor.md) applies at each
tensor length. A pure correction has cut rank at most one. After
saturating and normalizing its image line over the local valuation
ring, its nonzero constant reduction has physical matrix rank one.
It cannot lie in the reference flag through order nine.

A unit coordinate on a higher-order row therefore permits regular row
operations that clear the entire correction from the lower rows.
Those operations leave every original truncated row through precision
ten intact. This uses the whole image vector, including its higher
jets, rather than only the correction's first visible coefficient.

For the relative order lists `alpha` and `beta`, this gives

\[
 \sum_i(x-\beta_i)_+\ge\sum_i(x-\alpha_i)_+
 \qquad(x\le10).
\]

Rank-one valuation interlacing bounds the derivative difference below
by `-1`. Integrating the hinge functions against a positive exponential
therefore gives, for `0<s<1`,

\[
 Z_{\mathrm{new}}(s)\ge P(s)^k-s^{10}.
\]

Arbitrary negative Laurent orders cause no exception: multiply both
matrices by a common local monomial until they are integral, apply the
truncated-rank argument, and translate their order lists back. The
opposite reference minima `-8k` and `8k` cancel in the two-place
partition product. Thus the complete product satisfies

\[
 B_{\mathrm{new}}(s)\ge\bigl(P(s)^k-s^{10}\bigr)^2.
\]

## 4. Reduce every tensor length to one scalar inequality

Set

\[
 \delta=s^{10},\qquad A(s)=P(s)-\delta.
\]

Every coefficient of `A` is nonnegative, and `A(s)>=1` on `[0,1]`.
For every integer `k>=1`,

\[
\begin{aligned}
 P(s)^k-A(s)^k
 &=\delta\sum_{j=0}^{k-1}P(s)^{k-1-j}A(s)^j\\
 &\ge\delta.
\end{aligned}
\]

Hence

\[
 P(s)^k-s^{10}\ge A(s)^k,
 \qquad B_{\mathrm{new}}(s)\ge A(s)^{2k}.
\]

The physical partition on each paired boundary is `Q(s)^k`. At a
proposed floor `rho=m/n`, the paired Gibbs pressure is therefore
nonnegative whenever

\[
 \boxed{A(s)^n\ge Q\bigl(s^{n/(2m)}\bigr)^{2m}
       \quad(0\le s\le1).}                            \tag{1}
\]

The physical Gibbs expression is increasing in the proposed rate:
the derivative of `r log sum(exp(-E/r))` is the corresponding entropy.
Thus nonnegative pressure at `m/n` excludes a smaller optimized rate.
This step charges both complete local profiles and both physical types.

The certificate below proves (1) for `m=82`, `n=113`. The comparison
with the factoring threshold is the exact integer inequality

\[
 2\cdot82^2=13448>12769=113^2.
\]

## 5. Exact certificate for every scale

The certificate uses integer arithmetic throughout. No floating-point
optimizer, logarithm, or approximate sign is used to accept an interval.

### The interval from zero to 1/1024

Put `D=1024` and `u=s^(n/m)`. Since `n>m`, one has `u<=s` on `[0,1]`,
and

\[
 Q\bigl(s^{n/(2m)}\bigr)=1+u+u^3+u^4+u^5\le1+4u.
\]

The checker verifies

\[
 (D+4)^{2m-1}<2D^{2m-1},\qquad
 n^m D^{n-m}>(16m)^m.
\]

For `u<=1/D`, the first inequality bounds the derivative of
`(1+4u)^(2m)` by `16m`. Therefore

\[
 Q\bigl(s^{n/(2m)}\bigr)^{2m}\le1+16m s^{n/m}.
\]

Since `A(s)>=1+s`, Bernoulli's inequality gives `A(s)^n>=1+ns`.
The second integer comparison proves
`ns>=16m s^(n/m)` throughout `0<=s<=1/D`. This proves (1) there,
including equality at zero.

### The interval from 1/1024 to one

The remaining interval is partitioned into 2,807 exact dyadic
intervals, with largest denominator `2^18`. On an interval `[a,b]`,
the record retains a dyadic number `v` with denominator `2^48` and
checks both

\[
 v^{164}\ge b^{113},\qquad A(a)^{113}>Q(v)^{164}.
\]

Both comparisons are between integers after clearing powers of two.
The first proves `v>=b^(113/164)`. Monotonicity of `A` and `Q` then
proves (1) on the whole interval. The checker verifies that the
intervals have no gaps or overlaps and cover exactly `[1/1024,1]`.

The dyadic upper roots are constructed by integer Newton iteration.
Both neighboring integer-power comparisons are checked after that
iteration. A negative control confirms that the same lower partition
does not justify the stronger but false floor `3/4`: an explicit
rational substitution reverses the required inequality.

## 6. Reproduction and evidence

The [checker](../tools/check_pure_all_copy_flag_floor.py) first replays
the geometric polynomial identities and their finite-field controls.
It reconstructs the reference profiles, verifies the small interval,
and constructs and checks every dyadic pressure interval. It also
includes 180 exact rational controls for the tensor-length factorization.

The [complete record](../evidence/pure-all-copy-flag-floor.json) retains
the coefficients, all intervals and upper roots, the geometric
certificate hash, the source hashes, and the negative control.

```sh
python3 tools/check_pure_all_copy_flag_floor.py > /tmp/pure-all-copy-floor.json
cmp /tmp/pure-all-copy-floor.json evidence/pure-all-copy-flag-floor.json
```

The ordinary proof supplies the all-length conclusion; finite controls
alone are not used to infer it. The theorem excludes this one-pure-term
family at fixed natural heights over every characteristic-two field.
The multi-term and independent-height searches remain separate tasks.
No Lean source is added or changed.
