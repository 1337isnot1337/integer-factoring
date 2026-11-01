# Exact separation ranks of finite-precision Cauchy powers

## Result and scope

The [inverse-weight compilers](pole-cauchy-product-power-weights.md)
evaluate many separated contractions without listing either matrix axis.
The factor-count identity still needs a singular joint kernel. This note
determines exactly how many separated terms a direct approximation of that
kernel requires after power-of-two padding.

Every positive entrywise power has the same Smith profile as the original
pole Cauchy matrix, multiplied by the power. Consequently, its minimum
separation rank at any specified absolute precision is an explicit
reciprocal-residue overlap count. At the precision sufficient for the
inverse/Hadamard-square trace, an explicit arithmetic family requires
exponentially many terms. Exponentially many of these terms arise from
nonfactor rows.

This is a bound on representing the **whole joint kernel** as a sum of
separated functions. It does not bound direct scalar aggregation, a
representation with additional operators, or an algorithm that first
changes or prunes its domains. The complete factoring score remains `1/2`.

## 1. The power-of-two domain contract

Let `m=2^s`, let `N` be odd, and let each positive odd domain `A,B` contain
every odd residue modulo `2m` exactly once. The nodes may be arbitrary
integer lifts of those residues. This includes the
[power-of-two interval padding](pole-smith-overlap-spectrum.md#3-power-of-two-padding-gives-direct-diagonal-equivalence).
Work over the unramified quadratic extension of `Q_2`, with integer ring
`O` and `w^2+w+1=0`. Put

\[
 t=2^k w,\qquad z=N+t,\qquad
 M_{a,b}=\frac{t}{z-ab},\qquad k\ge s+1.
 \tag{1}
\]

For each row `a`, let `b(a)` be the unique column satisfying

\[
 ab(a)\equiv N\pmod {2m},
 \qquad
 \sigma_a=k-\min\{k,v_2(N-ab(a))\},
 \tag{2}
\]

where `v_2(0)=+infinity`. The quadratic displacement ensures
`v_2(N-ab+t)=min(k,v_2(N-ab))`, including equality of the two valuations.
In particular, every entry of `M` is integral.

The matched entry in row `a` has valuation `sigma_a<=k-s-1`.
Every other entry has valuation at least `k-s`. After permuting columns
by `b(a)`, the existing diagonal equivalence is

\[
 M=T(I+E),\qquad T_{a,a}=M_{a,b(a)},\qquad E\in2O^{m\times m}.
 \tag{3}
\]

The matching has a short query for one row. Its complete construction or
the histogram of (2) is not assumed to have a short global algorithm.

## 2. Every positive entrywise power has an exact Smith profile

For an integer `j>=1`, write `M^{circ j}` for the entrywise power, rather
than a matrix power. The same column permutation gives

\[
 M^{\circ j}=T^j(I+E_j),\qquad E_j\in2^jO^{m\times m}.
 \tag{4}
\]

Indeed, the diagonal ratios equal one. Every off-diagonal ratio is the
`j`th power of the corresponding ratio in (3), so its valuation is at
least `j`. The matrix `I+E_j` is unimodular. Therefore

\[
 \boxed{\operatorname{SmithOrders}(M^{\circ j})
       =\operatorname{sort}\{j\sigma_a:a\in A\}.}
 \tag{5}
\]

This argument does not use a determinant formula for a squared Cauchy
matrix, a permanent computation, or an assumption about the factor count.
It also proves that every positive entrywise power is nonsingular.

## 3. Exact minimum separation rank at absolute precision

For an integral matrix `X` and `q>=1`, define

\[
 R_q(X)=\min\{\operatorname{rank}_K Y:X-Y\in2^qO^{m\times m}\},
 \tag{6}
\]

where `K` is the quadratic field. Congruence in (6) already forces `Y`
to be integral. This equals the minimum number of arbitrary separated
terms `u_i(a)v_i(b)` in such an approximation. No polynomial-degree or
description restriction on those functions is needed for the lower bound.

If the Smith orders of `X` are `lambda_1<=...<=lambda_m`, then

\[
 R_q(X)=\#\{i:\lambda_i<q\}.
 \tag{7}
\]

For the upper bound, replace the Smith diagonal entries divisible by
`2^q` by zero. For the lower bound, transform a proposed approximation
by the same unimodular matrices. In the principal block with
`lambda_i<q`, divide row `i` by its nonzero Smith entry. Its diagonal
becomes one and its error becomes divisible by two. This block is
invertible, so its order is a lower bound for the rank of the approximation.

Define `H_r=m` for `r<=s+1`, and for `s+1<=r<=k` define

\[
 H_r=\#\{a\in A:ab(a)\equiv N\pmod {2^r}\}.
 \tag{8}
\]

These are exactly the [balanced overlap counts](pole-smith-overlap-spectrum.md).
Equations (5)--(7) give

\[
 \boxed{R_q(M^{\circ j})
   =H_{\,\max\{s+1,\ k-\lfloor(q-1)/j\rfloor\}}.}
 \tag{9}
\]

There is also a direct optimal approximation: retain row `a` precisely
when `j sigma_a<q`, and replace every other row by zero. Every entry
in a removed row is divisible by `2^q`. The retained rows have independent
matched pivots. Constructing this approximation explicitly still requires
finding and emitting all the retained rows.

In particular, if `q>j(k-s-1)`, every row is retained and the minimum
separation rank is the entire dimension `m`.

## 4. The natural trace precision has an exact overlap threshold

Choose `h` with `|N-ab|<2^h` throughout the domains. To recover the
factor count modulo `2^p`, take the standard padding

\[
 k=h+p,\qquad L=k-s-1,\qquad q=p+L,
 \tag{10}
\]

enlarging `h` to at least `s+1` if necessary. The conditioning theorem
gives `2^L M^{-1}` integral. Therefore replacing `M^{circ2}` by any
entrywise approximation modulo `2^q` changes
`Tr(M^{-1} M^{circ2})` by a multiple of `2^p`.
This is a sufficient, uniform error contract for the inverse-weight
compilers. It is not claimed to be necessary for every special scalar
identity that might exploit cancellation.

For the original kernel, the [Cauchy trace identity](pole-cauchy-content.md)
is

\[
 \operatorname{Tr}(M^{-1}M^{\circ2})
 =\sum_{a,b}M_{a,b}-{m\choose2}\frac{t}{z}
 \equiv D\pmod {2^p}.
 \tag{11}
\]

The nonfactor terms vanish at the chosen displacement; the correction
also vanishes. By (9), the minimum separation rank at precision (10) is

\[
 \boxed{R_q(M^{\circ2})=H_{r_*},\qquad
 r_*=\max\left\{s+1,\left\lceil\frac{h+s+2}{2}\right\rceil\right\}.}
 \tag{12}
\]

The requested output precision `p` cancels from this threshold. Raising
`p` while increasing the displacement as in (10) does not lower the
number of separated terms required by this particular route.

### Larger padding reduces the optimal rank to the true factor count

The obstruction at (10) does not apply to every displacement. Keep the
absolute kernel precision `q=p+k-s-1`, but allow a larger `k`. If

\[
 k\ge\max\{h,s+1,2h+p-s-3\},
 \tag{12a}
\]

then

\[
 \boxed{R_{p+k-s-1}(M^{\circ2})=D.}
 \tag{12b}
\]

Indeed, a nonfactor entry has valuation at least `k-h+1`, so its square
is divisible by `2^q` under (12a). Every factor entry equals one.
Thus the entire Hadamard square is congruent to the Boolean factor
incidence matrix at this precision. Its nonzero entries form a matching,
so that incidence matrix has rank `D`, and (9) gives the matching lower
bound. One may enlarge `k` further to `h+p` when also using the direct
factor-count congruence (11).

The larger displacement still has polynomial bit length. What remains
missing is a compressed construction of this optimal approximation:
its support specifies the factor pairs. In particular, the existence
of a rank-`D` approximation does not provide its separated factors or
an evaluator for their contraction. This is a positive representation
target after padding, rather than an exponential rank obstruction.

## 5. An exponential contribution from certified nonfactor rows

Let `s>=8`, `m=2^s`, and take

\[
 A=B=\{1,3,\ldots,2m-1\},\qquad
 N=m^2+(m+1)^2.
 \tag{13}
\]

Here `N` lies strictly inside the product range of the rectangle, and
the maximum absolute difference is `N-1`. Thus its bit length and a
valid difference parameter are both `n=h=2s+2`. Equation (12) becomes

\[
 r_*=\left\lceil\frac{3s+4}{2}\right\rceil,
 \qquad u=\left\lceil\frac{r_*}{2}\right\rceil
          =s+1-\lfloor s/4\rfloor.
 \tag{14}
\]

For every integer `j` with

\[
 |j|\le 2^{s-u}-1,
 \qquad a_j=m+1+2^u j,\quad b_j=m+1-2^u j,
 \tag{15}
\]

both nodes belong to the declared domains and

\[
 N-a_jb_j=m^2+2^{2u}j^2>0,
 \qquad 2^{r_*}\mid N-a_jb_j.
 \tag{16}
\]

The divisibility uses `2u>=r_*` and `2s>=r_*`. These pairs have distinct
rows. Since `r_*>=s+1`, each `b_j` is the unique matched column `b(a_j)`.
Equation (16) proves that none of these rows can contain a true factor
pair anywhere in `B`: any such column would have to be the same matched
column. There are exactly

\[
 2^{s-u+1}-1=2^{\lfloor s/4\rfloor}-1
 \tag{17}
\]

certified nonfactor rows. All true factor rows also contribute to `H_{r_*}`.
Consequently, for **every** requested precision `p>=1` with (10),

\[
 \boxed{R_{p+L}(M^{\circ2})
 \ge D+2^{\lfloor s/4\rfloor}-1
 =D+2^{\lfloor(n-2)/8\rfloor}-1.}
 \tag{18}
\]

This is an exponential lower bound in the binary input length for an
explicit arithmetic family, at the stated sufficient trace precision.
The bound permits arbitrary separated functions, so replacing explicit
polynomials by binary powers or more expressive one-variable circuit
weights cannot remove it.

### A variant with exactly zero factor incidence

If a zero-incidence example is desired, replace (13) by

\[
 N=16m^2+(m+1)^2.
 \tag{19}
\]

Every product in the rectangle is now strictly below `N`, so `D=0`.
Here `n=h=2s+5`, `r_*=ceil((3s+7)/2)`, and
`u=ceil((3s+7)/4)`. The same symmetric pairs give

\[
 H_{r_*}\ge2^{s-u+1}-1\ge2^{s/4-2}-1
 \qquad(s\ge12).
 \tag{20}
\]

The factor count in this second family is already available from the
two product corners. It is included to separate the representation
cost from genuine factor incidence, not as a difficult factoring input.

## 6. Consequences for the remaining extraction engine

Equation (9) says precisely what a separated-kernel compiler must find:
the rows whose canonical modular-inverse matches survive a specified
extra number of low-bit tests. The description of one match is small.
The number of surviving matches can be exponential even when they are
known not to be factors.

This closes a route that first approximates the whole Hadamard square
at (10), emits an explicit sum of separated terms, and invokes the
existing inverse-contraction compiler on each term. It does not exclude
a direct trace identity, cancellation between implicitly described
families of terms, the more strongly padded rank-`D` target (12b), or an
aggregation method that preserves the modular inverse relation without
emitting its rows. Those remain appropriate targets for a faster
extraction engine.

## 7. Exact evidence

The [checker](../tools/check_pole_cauchy_hadamard_rank.py) independently
computes small Smith profiles by exact pivot elimination in `Q(w)`.
It tests positive entrywise powers, the matched-pivot valuations, the
row-pruned optimal approximation, and the overlap threshold (12).
Separate integer checks verify the explicit nonfactor witness families,
including huge encoded domains without enumerating them.

```sh
python3 tools/check_pole_cauchy_hadamard_rank.py --output /tmp/hadamard-rank.json
cmp /tmp/hadamard-rank.json evidence/pole-cauchy-hadamard-rank.json
```

The [receipt](../evidence/pole-cauchy-hadamard-rank.json) records the finite
controls and distinguishes literal counts from theorem-based large-domain
certificates: 72 independent Smith profiles, 278 optimal-approximation
controls, and 952 off-diagonal ratio checks on 18 matrices. Larger
padding adds 18 exact rank-`D` profiles, 297 entrywise projection
checks, and 36 trace checks for both the whole square and its incidence
approximation. At `s=12`,
the literal overlap counts in the two explicit families are respectively
39 and 13, while both factor counts are zero. At `s=200`, the first
family certifies `2^50-1` nonfactor rows in a 402-bit input without
enumerating its `2^200` rows. It supplies no global overlap evaluator
and no improved factoring exponent.
