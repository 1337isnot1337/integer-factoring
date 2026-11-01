# Exact rectangular Bezout denominators from two Mahler variables

## Result

The [balanced denominator theorem](balanced-bezout-denominator.md)
extends to unequal domain sizes. The extension gives every Newton
coefficient valuation individually, not only the largest denominator.
Its proof also gives a short precision-dependent formula for the
normalized unit of each coefficient.

Let `A_i=1+4f(i)`, `0<=i<m`, and `B_j=1+4g(j)`, `0<=j<s`, be positive
integers, where both node maps are dyadic isometries on their index sets:

\[
 v_2(f(i)-f(j))=v_2(i-j),\qquad
 v_2(g(i)-g(j))=v_2(i-j)\quad(i\ne j).
 \tag{1}
\]

For `N=3 mod4`, put

\[
 P(X)=\prod_{i<m}(X-A_i),\qquad
 Q(X)=\prod_{j<s}(N-B_jX).
 \tag{2}
\]

No root of `P` is a root of `Q`, already by reduction modulo four.
For every `0<=k<m`,

\[
 \boxed{v_2\bigl((1/Q)[A_0,\ldots,A_k]\bigr)
 =-s-k+v_2{\binom{s+k-1}{k}}.}
 \tag{3}
\]

Consequently the least exponent `r` for a certificate
`UP+VQ=2^r`, with coefficients in `Z_(2)`, is

\[
 \boxed{r=\max_{0\le k<m}
       \left(s+k-v_2{\binom{s+k-1}{k}}\right).}
 \tag{4}
\]

The ideal is unchanged when the two monic root products are exchanged,
so the answer is symmetric in `m,s`. At `m=s`, (4) agrees with the
earlier closed formula `r=m+v_2((m-1)!)`.
This is a precise denominator and coefficient-construction result for
the separated residue classes in (1). It supplies no singular factor
count evaluator and does not improve the complete factoring score.

## 1. A normalized divided-difference functional

Let `x_0,...,x_a` be elements of `Z_2` with
`v_2(x_i-x_j)=v_2(i-j)`. Define

\[
 C_{a,r}(x)=a!\,[x_0,\ldots,x_a]\binom{X}{r}.
 \tag{5}
\]

Then

\[
 C_{a,r}(x)\in\mathbb Z_2,\qquad
 C_{a,r}(x)=0\ (r<a),\qquad C_{a,a}(x)=1.
 \tag{6}
\]

For integrality, use the Lagrange formula. The denominator of its
`i`th term has valuation

\[
 \sum_{j\ne i}v_2(x_i-x_j)
 =v_2(i!)+v_2((a-i)!)\le v_2(a!).
 \tag{7}
\]

Also `binom(x_i,r)` is integral: this rational polynomial is continuous,
is integer-valued on the dense set of nonnegative integers, and `Z_2`
is closed. Thus multiplication by `a!` makes each Lagrange term integral.
The other assertions in (6) follow from degree and leading coefficient.
The case `a=0` is included.

## 2. The double Mahler expansion of a unit Cauchy kernel

Consider

\[
 K(x,y)=\frac1{1+2(y-x)}\qquad(x,y\in\mathbb Z_2).
 \tag{8}
\]

Its double forward differences at zero are exactly

\[
 A_{r,t}=\Delta_x^r\Delta_y^tK(0,0)
 =\frac{(-1)^t2^{r+t}(r+t)!}
        {\prod_{v=-r}^{t}(1+2v)}.
 \tag{9}
\]

To prove (9), write `K(x,y)=G(y-x)` for `G(z)=1/(1+2z)`.
The shift identity `Delta_x=-E^(-1)Delta_y` reduces the mixed difference
to `(-1)^r Delta^(r+t)G(-r)`. Induction gives

\[
 \Delta^dG(z)=\frac{(-1)^d2^dd!}
                   {\prod_{i=0}^{d}(1+2z+2i)}.
 \tag{10}
\]

Every denominator in (9) is odd, so
`v_2(A_(r,t))=r+t+v_2((r+t)!)`. In particular the expansion

\[
 K(x,y)=\sum_{r,t\ge0}A_{r,t}\binom{x}{r}\binom{y}{t}
 \tag{11}
\]

converges uniformly on `Z_2^2`. On nonnegative integer pairs, (11) is
the finite double binomial-inversion identity. Continuity and density
extend it to the whole domain. This proves (11) without assuming a
separate Mahler expansion theorem.

## 3. Exact parity differences on two isometric prefixes

Take two independent isometric prefixes `x_0,...,x_a` and
`y_0,...,y_b` in `Z_2`. Set

\[
 E_i=2x_i,\qquad O_j=1+2y_j,
 \qquad D=[E_0,\ldots,E_a,O_0,\ldots,O_b]\varepsilon,
 \tag{12}
\]

where `epsilon` is zero on even elements and one on odd elements.
The Lagrange formula, or two successive divided differences of (8), gives

\[
 D=2^{-a-b}[x_0,\ldots,x_a]_x[y_0,\ldots,y_b]_yK(x,y).
 \tag{13}
\]

For example, the first divided difference in (13) is
`2^a/prod_i(1+2y-2x_i)`; the second supplies exactly the Lagrange
denominators between the odd nodes.

Apply the two finite divided-difference functionals to (11). By (5)--(6),

\[
 D=\frac1{2^{a+b}a!b!}
    \sum_{r\ge a,\ t\ge b}A_{r,t}C_{a,r}(x)C_{b,t}(y).
 \tag{14}
\]

The term `(r,t)=(a,b)` equals

\[
 B_{a,b}=\frac{(-1)^b\binom{a+b}{a}}
                    {\prod_{v=-a}^{b}(1+2v)}.
 \tag{15}
\]

For any other term, put `d=r+t-a-b>=1`. Its valuation is at least

\[
 d+v_2((a+b+d)!)-v_2(a!)-v_2(b!)
 \ge d+v_2\binom{a+b}{a}.
 \tag{16}
\]

Thus (15) is the unique term of least valuation, and

\[
 \boxed{v_2(D)=v_2\binom{a+b}{a},\qquad D/B_{a,b}\equiv1\pmod2.}
 \tag{17}
\]

The result holds for every pair of prefixes in (12), including `a=0`
or `b=0`. No ordering relation between the two independent isometries
is required.

## 4. Newton coefficients and the exact certificate exponent

Fix `k`, let `C_j=N/B_j`, and write
`Q=alpha prod_j(X-C_j)` with `alpha=(-1)^s prod_j B_j`, an odd unit.
Under `X=1+2Z`, the row nodes become `2f(i)` and the reciprocal column
nodes become

\[
 1+2y_j,\qquad y_j=(N/B_j-3)/4.
 \tag{18}
\]

The `y_j` are integral and form an isometric prefix, since inversion
and multiplication by odd units preserve the pairwise valuations of
the `B_j`.

Put `P_k=prod_(i<=k)(X-A_i)`. The sum of the residues of
`1/(P_k prod_j(X-C_j))` is zero, as its denominator has degree at least
two. The Lagrange formula therefore gives

\[
 [A_0,\ldots,A_k,C_0,\ldots,C_{s-1}]\varepsilon
 =-\alpha\,(1/Q)[A_0,\ldots,A_k].
 \tag{19}
\]

Here `epsilon` labels the reciprocal nodes by one and the row nodes by
zero. The affine change multiplies a divided difference of order `k+s`
by `2^(-k-s)`. Apply (17) with `a=k,b=s-1`. Since `alpha` is a unit,
this proves (3).

The reduced inverse of `Q` modulo `P` is its interpolation polynomial

\[
 V_0(X)=\sum_{k=0}^{m-1}(1/Q)[A_0,\ldots,A_k]
                  \prod_{i<k}(X-A_i).
 \tag{20}
\]

The monic Newton basis in (20) and the ordinary monomial basis are
related by an integral triangular matrix with integral inverse. Their
minimum coefficient valuations are therefore equal. This gives the
right side of (4) as the exact denominator of `V_0`.

If `2^r V_0` is integral, monic division of `2^r(1-V_0Q)` by `P`
produces an integral `U`. Conversely, any integral Bezout certificate
with right side `2^r`, reduced modulo the monic polynomial `P`, has
second coefficient `2^r V_0`. Thus higher-degree certificates cannot
improve the exponent, proving (4).

### Binary evaluation of the maximum

Write `pop(n)` for the number of ones in the binary expansion of `n`.
Legendre's factorial identity gives

\[
 \kappa_k=v_2\binom{s+k-1}{k}
 =\operatorname{pop}(k)+\operatorname{pop}(s-1)
       -\operatorname{pop}(s+k-1).
 \tag{21}
\]

Let `B=bits(m+s)`. The binomial valuation counts binary carries, so
`0<=kappa_k<=B`. The final index already gives `r>=m+s-1-B`.
Any index more than `B` places before the final index gives a strictly
smaller value than that lower bound. It therefore suffices to check
the last `min(m,B+1)` indices in (4). Every integer is binary encoded;
this computes the exact exponent in polynomial bit time.

Put `delta=m+s-1-r`, so `0<=delta<=B`. After multiplication by `2^r`,
a Newton coefficient `ell` places before the last one has valuation
`ell-delta+kappa_k>=ell-delta`. Modulo `2^p`, only the final `p+delta`
Newton terms can survive.

There is an equally short bound for the entire normalized polynomials.
Use `Y=X-1`. Each root shift `A_i-1` has valuation at least two.
A Newton term `ell` places below top degree contributes to the monomial
coefficient `d` places below top degree with valuation at least
`ell-delta+2(d-ell)>=d-delta`. The same argument for `U` exchanges the
two root sets. Its reciprocal-root shifts `N/B_j-1` have valuation one,
which still gives `ell-delta+(d-ell)=d-delta`. Thus

\[
 \begin{split}
 2^rV_0(1+Y)&\equiv\sum_{0\le d<p+\delta}v_dY^{m-1-d},\\
 2^rU_0(1+Y)&\equiv\sum_{0\le d<p+\delta}u_dY^{s-1-d}
 \pmod {2^p},
 \end{split}
 \tag{21a}
\]

where terms with negative exponents are omitted. Both complete
polynomials therefore have at most `p+delta` possibly nonzero
coefficients, with binary encoded degrees. This bound uses the exact
Newton valuations; it does not rely on a truncated division heuristic.

The [binary scale helper](../tools/rectangular_bezout_scale.py) implements
this maximization and supplies complete coefficient windows. Its
slightly tighter carry cap `bits(m+s-1)` is valid by the same argument.

## 5. A constructive formula with only low surplus indices

The same proof computes normalized units. Use the leading term `B_(a,b)`
from (15). For `i,j>=0`, put

\[
 R_{i,j}=\prod_{\ell=1}^{i}(1-2a-2\ell)
          \prod_{\ell=1}^{j}(1+2b+2\ell).
 \tag{22}
\]

These are odd units. Cancellation of the common odd product in (9)
and (15) yields the finite precision formula

\[
 \boxed{\frac{D}{B_{a,b}}\equiv
  \sum_{i+j<p}\frac{(-1)^j2^{i+j}(a+b+1)_{i+j}}{R_{i,j}}
       C_{a,a+i}(x)C_{b,b+j}(y)\pmod {2^p}.}
 \tag{23}
\]

The notation `(u)_d=u(u+1)...(u+d-1)` denotes a rising product, with
empty product one. There are `p(p+1)/2` terms. Each factorial ratio
in (23) has at most `p-1` factors, even when `a,b` are exponentially
large. Equation (16) justifies discarding every omitted term before
performing any large-index arithmetic.

The remaining normalized functionals also use only low moments. Let
`e_d` and `h_d` denote elementary and complete homogeneous symmetric
polynomials. Then

\[
 C_{a,a+i}(x)=\frac1{(a+1)_i}
   \sum_{\ell=0}^{i}(-1)^{i-\ell}
     e_{i-\ell}(0,1,\ldots,a+i-1)\,
     h_\ell(x_0,\ldots,x_a).
 \tag{24}
\]

Expand the falling-factorial numerator of `binom(X,a+i)` and use
`[x_0,...,x_a]X^(a+ell)=h_ell(x)` to prove (24). This is the same
low-surplus identity used in the balanced construction.
The division in (24) is exact in `Z_2` by (6). Its required extra
precision is only

\[
 v_2((a+1)_i)=O(i\log(a+i+1)),
 \tag{25}
\]

which is polynomial in the surplus and the binary index length.

Thus, given the short moment lists through degree `p-1`, (23)--(25)
compute the normalized difference in polynomial bit time. The unit
part of (15) is computable by normalized factorial products and an
odd arithmetic-progression product, as in the
[balanced modular construction](balanced-bezout-modular.md).
For ordinary row progressions and reciprocal column progressions,
the same short-jet primitives provide the required moments.
For an arbitrary isometry supplied only as an evaluation black box,
this note does not assume that its global moments are inexpensive.

## 6. Larger separation scales

The proof also covers every `h>=1`. Suppose

\[
 A_i=1+2^{h+1}f(i),\qquad B_j=1+2^{h+1}g(j),
 \qquad N\equiv1+2^h\pmod {2^{h+1}},
 \tag{26}
\]

with the isometries in (1). Substitute `X=1+2^h Z` in place of
`X=1+2Z`. The normalized even and odd nodes are unchanged in kind,
and the exact formulas become

\[
 v_2((1/Q)[A_0,\ldots,A_k])
   =-h(s+k)+v_2\binom{s+k-1}{k},
 \qquad
 r=\max_{0\le k<m}\left(h(s+k)-v_2\binom{s+k-1}{k}\right).
 \tag{27}
\]

The binary maximization and the low-surplus unit formula retain
polynomial dependence on the encoded sizes and the requested precision.
The separation hypothesis still excludes actual factor incidence.

## 7. Exact checks and replay

The [checker](../tools/check_rectangular_bezout_denominator.py) tests
the double forward-difference formula by direct alternating sums.
It then compares rational divided differences with the exact valuation
and normalized-unit formulas on independent random binary-tree
isometries, including rational fractional-linear changes of nodes.

The [receipt](../evidence/rectangular-bezout-denominator.json) records:

- 81 direct double-difference identities and 64 finite Mahler reconstructions;
- 128 rectangular parity valuations and 384 normalized-unit truncations;
- 1,536 integral normalized functionals and 1,536 low-moment identities;
- 864 exact inverse Newton coefficients at separation heights one, two, and four;
- 192 independently constructed modular Hensel inverses and rational Bezout certificates;
- 576 shifted-polynomial coefficient tail bounds.

Four large controls evaluate only the binary scale formula, with sizes
through `2^400`. For example, sizes `2^200+13` and `2^137+7` have defect
two, so both normalized polynomials have at most `p+2` possible shifted
coefficients modulo `2^p`. This is a consequence of the proved formula;
the large controls do not construct those polynomials or enumerate nodes.

```sh
python3 tools/check_rectangular_bezout_denominator.py --output /tmp/rectangular-bezout.json
cmp /tmp/rectangular-bezout.json evidence/rectangular-bezout-denominator.json
```

All checks use exact rational or integer modular arithmetic. They are
finite evidence for the ordinary proofs above, not a formalization or
a general singular-GCD implementation.
