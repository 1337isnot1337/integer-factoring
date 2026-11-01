# Primitive Bezout polynomials beyond the first residue split

## Constructive result

There is a polynomial-bit-time modular Bezout constructor for reflected
arithmetic domains with exponentially many matched residue pairs.
Their integer roots and reciprocal roots occupy the same low residue
classes and separate only after `log2(m)` further bits. The construction
computes the actual reciprocal deformation; it does not replace those
roots by a nearby arithmetic progression.

Let `m=2^t>=2`, `h>=2t+1`, `D=2^h`, and let `a,eta` be odd integers
with `a>D(m-1)>0` and `N=a^2+Dm*eta>0`. Define

\[
 A_i=a+Di,\quad B_i=a-Di\quad(0\le i<m),\qquad
 P(X)=\prod_i(X-A_i),\quad Q(X)=\prod_i(N-B_iX).       \tag{1}
\]

The optimal Bezout denominator exponent is exactly

\[
 R=h(2m-1)+2v_2((m-1)!).                            \tag{2}
\]

Let `U_0,V_0` be the unique rational pair of degrees below `m` with
`U_0P+V_0Q=1`. For every numerical precision `p>=1`, one algorithm
computes the complete polynomials `U=2^R U_0`, `V=2^R V_0 modulo 2^p`
as at most `4p` possible highest powers of `X-a` in each polynomial.
Both have degree `m-1` and odd leading coefficient. All omitted
coefficients are zero at the requested modulus. The cost is polynomial
in the binary parameter lengths and `p`; `m,R` are retained in binary.

The [sibling-separated constructor](balanced-bezout-modular.md) uses
parity at the first local residue split. Here every `A_i` and `N/B_i`
agree modulo `2^(h+t)` and differ at the next bit. The deeper matching
uses a separate perturbation expansion with an explicit tail bound.
The old one-bit-per-term parity truncation is not assumed.

This proves a compressed primitive for one structured inverse pairing.
It does not handle a general modular inverse permutation. The complete
factoring score remains `1/2`. Section 6 gives the exact changing
root-incidence rule in the wider reflected family; the primitive
constructor itself assumes odd `eta` and hence coprimality.

Files: [constructor](../tools/matched_progression_bezout.py),
[checker](../tools/check_matched_progression_bezout.py), and
[evidence](../evidence/matched-progression-bezout.json).

## 1. Ordered isometry and the translated reference problem

Normalize by `z=(X-a)/D`, and put

\[
 c=m\eta/a,\qquad\lambda=D/a,\qquad
 g(j)=\frac{j+c}{1-\lambda j}.                    \tag{3}
\]

The first root set becomes `0,...,m-1`; the second becomes
`g(0),...,g(m-1)`. Since `v2(c)=t` and `v2(lambda)=h>t`,

\[
 v_2(g(j)-j)=t,\qquad
 v_2(g(i)-g(j))=v_2(i-j).                        \tag{4}
\]

Indeed `g(j)-j=(c+lambda*j^2)/(1-lambda*j)` and
`g(i)-g(j)=(1+lambda*c)(i-j)/((1-lambda*i)(1-lambda*j))`.
For distinct indices in the range, `v2(i-j)<t`, so every other
cross difference has that smaller valuation. Consequently the order

\[
 0,1,\ldots,m-1,g(0),g(1),\ldots,g(m-1)          \tag{5}
\]

is an ordered dyadic isometry of `0,...,2m-1`. Its second-set
indicator is the `t`th binary digit, because `g(j)=j+m mod 2m`.
All nodes are distinct, establishing coprimality.

Let `d_r(lambda)` be the indicator divided difference on the first
`m+r` nodes, for `1<=r<=m`; its order is `k=m+r-1`.
At `lambda=0`, ordinary finite differences give

\[
 d_r(0)=
 \frac{(-1)^{r-1}\binom{m+r-2}{r-1}}
 {\prod_{j=-m+1}^{r-1}(c+j)}.                    \tag{6}
\]

One proof uses
`Delta^q (z^(overline M))^(-1)
=(-1)^q M^(overline q)/(z^(overline(M+q)))`.
This follows by one subtraction and induction. Apply it to the
Lagrange sum on the `r` translated nodes and divide by `(r-1)!`.

Put `K=2m-1` and `s=K-k=m-r`. Every nonzero integer shift `j`
with `|j|<m` satisfies `v2(c+j)=v2(j)`. Legendre's formula in (6)
therefore gives

\[
 v_2(2^{R-hk}d_r(0))
   =hs+v_2(s!)-v_2(s+1)=:\delta_s.              \tag{7}
\]

This includes `r=1`. For `s=0` the value is zero. For every `s>=0`,
`delta_s>=floor(s/2)`: if `s` is even, `v2(s+1)=0`; if `s=2j-1`,
then `v2(s+1)=1+v2(j)<=j`.

## 2. A charged reciprocal perturbation

Keep `c` fixed and view the Lagrange formula as a formal series in
`lambda`. A difference from a translated second-set node to an integer
node factors as

\[
 g(j)-i=(j+c-i)
 \frac{1+\lambda ij/(j+c-i)}{1-\lambda j}.        \tag{8}
\]

Its constant difference has valuation at most `t`. Every degree-`q`
coefficient in the remaining factor or its inverse therefore has
valuation at least `-tq`. Within the second root set, the difference
formula after (4) has integral unit factors and gives at least as
strong a bound.

The constant barycentric denominator at order `k` has valuation at
most `v2(k!)`, by the isometry (5) at `lambda=0`. It follows that

\[
 v_2([\lambda^q]d_r(\lambda))\ge-v_2(k!)-tq.     \tag{9}
\]

After substituting the actual `lambda` and multiplying by `2^(R-hk)`,
the degree-`q` contribution has valuation at least

\[
 -t+h(K-k)+v_2(K!)-v_2(k!)+(h-t)q.              \tag{10}
\]

Here `2v2((m-1)!)-v2(K!)=-t`. Formula (10) proves convergence and
shows that a truncation after

\[
 Q=\left\lceil\frac{p+t}{h-t}\right\rceil        \tag{11}
\]

terms suffices modulo `2^p`, uniformly over all retained prefixes.
This is the full precision charge for the deformation. As `h>=2t+1`,
one has `Q<=p`.

Every nonconstant perturbation term in (10) has valuation at least
`h-2t>=1`, plus a nonnegative prefix term. Thus (7) and (10) prove
that all `2^(R-hk)d_r(lambda)` are integral and the leading one is
odd. The unique indicator interpolant `F` satisfies `F=U_0P` and
`1-F=V_0Q`. Both input polynomials have unit leading coefficients.
The primitive leading coefficient therefore proves (2), including
optimality among Bezout polynomials of arbitrary degree by reduction
modulo the monic polynomial `P`.

Equations (7) and (10) also show that all Newton terms with
`K-k>=2p` vanish modulo `2^p`. The terminal window survives the
deeper pairing, although its proof differs from the parity case.

## 3. Compute every retained perturbation coefficient

Let `F_m(T)=1/prod_(i=0)^(m-1)(T-i)` and define

\[
 L_b(c)=\frac{\Delta^{r-1}}{(r-1)!}
               [T^bF_m(T)]_{T=c}.              \tag{12}
\]

Write `S(b,v)` for a Stirling number of the second kind. For
`b<m` and `b<r`, the falling-factorial expansion of `T^b`, followed
by the same finite-difference identity as in (6), gives

\[
 L_b(c)=\sum_{v=0}^{b} S(b,v)\,
 \frac{(-1)^{r-1}\binom{m+r-v-2}{r-1}}
      {\prod_{j=-m+1}^{r-v-1}(c+j)}.             \tag{13}
\]

Denote the `v`th rational product ratio in (13) by `R_v(c)`.
Only `b<=2(Q-1)` is needed below. If `m` or `r` is at most this
bound, the retained-prefix condition `m-r<2p` implies `m<4p`.
The implementation then uses literal rational interpolation with a
polynomial-size node list. It rejects a request to use that fallback
on a large, nonterminal prefix.

The exact barycentric formula for the deformed nodes is

\[
 d_r(\lambda)=
 \frac{\prod_{j=0}^{r-1}(1-\lambda j)}{(1+\lambda c)^{r-1}}
 \frac{\Delta^{r-1}}{(r-1)!}
 \left[(1-\lambda j)^{r-2}
       F_m\!\left(\frac{j+c}{1-\lambda j}\right)\right]_{j=0}.
                                                               \tag{14}
\]

This follows directly by multiplying the pair differences in (4).
The prefactor in (14) and its formal `lambda`-inverse both have
integral coefficients and unit constant term. Multiplying by either
therefore preserves the coefficient lower bound (9). The same tail
cutoff (11) is valid when the inner factor is truncated first and
the exact unit prefactor is multiplied afterward, as in the code.
Write `H^[q]` for the divided derivative `H^(q)/q!`.
The coefficient of `lambda^s` in the square-bracket expression is

\[
 \sum_{q=0}^{s}(-1)^{s-q}
     {r-2-q\choose s-q}j^s(j+c)^qF_m^{[q]}(j+c). \tag{15}
\]

The binomial coefficient has its ordinary integer extension to
negative upper arguments. To express the finite difference of (15)
using (13), expand `j^s=((j+c)-c)^s` and use

\[
 T^{\ell+q}F_m^{[q]}(T)
   =\sum_{u=0}^{q}(-1)^{q-u}{\ell+q\choose q-u}
            (T^{\ell+u}F_m(T))^{[u]}.           \tag{16}
\]

Identity (16) follows from the product rule, or induction using
`T H'= (TH)'-H`. Differentiation commutes with the finite difference.
Thus the degree-`s` coefficient of the second factor in (14) is

\[
 \sum_{q=0}^{s}\sum_{\ell=0}^{s}\sum_{u=0}^{q}
 (-1)^{s-u}{r-2-q\choose s-q}{s\choose\ell}
 (-c)^{s-\ell}{\ell+q\choose q-u}
 L_{\ell+u}^{[u]}(c).                            \tag{17}
\]

Equations (13) and (17) reduce the whole computation to short jets
of explicit product ratios. The number of terms depends polynomially
on `Q`, not on `m` or the divided-difference order.

## 4. Binary normalization of the product ratios

The exact valuation of `R_v` is

\[
 E_v=v_2\binom{m+r-v-2}{r-1}
       -t-v_2((m-1)!)-v_2((r-v-1)!).             \tag{18}
\]

Compute the jet of `2^(-E_v)R_v(c+2^t Z)` to degree `Q-1`.
Every denominator factor becomes

\[
 c+j+2^tZ=\frac{m\eta+aj+a2^tZ}{a}.
\]

For `j!=0` in the displayed interval, its constant valuation is
`v2(j)<t`. In the bucket `j=j_0+2^(e+1)b` of valuation `e`, division
by `2^e` leaves an odd constant jet and the even progression step
`2a`. The isolated `j=0` factor leaves the unit jet `eta+aZ`.
There are only `t+1` buckets. The existing
[unit progression algorithm](interval-unit-jets.md) computes their
products at any polynomial precision without visiting their nodes.
The extra odd denominator `a` is handled by modular binary powering.

The binomial numerator in (18) is handled as a ratio of three
normalized factorials. Its enormous integer value is never formed.
Exact factorial valuations use `v2(n!)=n-s_2(n)`; odd factorial
parts use valuation buckets and odd progression products.

If the normalized ratio jet has coefficient `b_(v,u)`, then

\[
 R_v^{[u]}(c)=2^{E_v-tu}b_{v,u}.                 \tag{19}
\]

In (17), after multiplying by `lambda^s` and `2^(R-hk)`, the
corresponding power-of-two exponent is

\[
 e=R-hk+E_v-tu+hs+t(s-\ell).                    \tag{20}
\]

All remaining factors are integer binomial or Stirling coefficients,
an odd power of `eta`, and an odd denominator power of `a`.
Compute the minimum `e` over the finite retained formula and set
`guard=max(0,-min e)`. Work at precision `p+guard`, add every term
after multiplication by `2^guard`, and divide only the complete sum
by that power of two. Its integrality follows from Section 2.
This avoids dividing a nonintegral individual summand in the residue
ring.

The guard has polynomial numerical size. Changing `v` from zero
in (18) alters only `O(v)` nearby factorial factors, each with
`O(t)` bits. Thus `|E_v-E_0|=O(v(t+1))`. Since
`R-hk+E_0=delta_(K-k)>=0`, (20) gives `guard=O(Q(t+1))`.
Every literal shift and allocated coefficient ring is therefore
polynomial in `p` and the encoded input length.

The prefactor in (14) is a unit. Its progression product and its
negative integer power use ordinary modular unit operations at the
same precision. Multiplying it before the final guard division
preserves all requested bits.

## 5. Complete shifted polynomials from the terminal coefficients

The first `m` interpolation values are zero. Therefore every nonzero
Newton term has the factor `P`, and the normalized first Bezout
polynomial has the expansion

\[
 U(X)=\sum_{r=1}^{m}2^{R-h(m+r-1)}d_r(\lambda)
                 \prod_{j<r-1}(X-N/B_j).         \tag{21}
\]

Only `m-r<2p` can survive. Every root shift `N/B_j-a=Dg(j)` is
divisible by `2^h`. A coefficient at reverse distance `j` in one
of the products in (21) has valuation at least `hj`. With
`q=floor((p-1)/h)`, the complete reverse support of `U` is thus
contained in degrees at most `2p-1+q`.

These short prefix products are also effective unit progressions:

\[
 1-(N/B_j-a)T
  =\frac{a-Dj-D(aj+m\eta)T}{a-Dj}.              \tag{22}
\]

Both numerator and denominator products in (22) have odd constant
terms and even progression steps. Removing a final reversed linear
factor uses a unit series inverse. No root list is expanded.

In the coordinate `Y=X-a`, let stars denote reversal at each
polynomial's actual degree. The identity `UP+VQ=2^R` gives

\[
 V^*=-U^*P^*/Q^*\pmod{T^m}.                    \tag{23}
\]

The scalar correction has reverse degree `2m-1`, outside this range.
Every degree-`j` coefficient of `P^*/Q^*` has valuation at least
`hj`: both root-shift products have this bound, and the inverse
series recurrence preserves it. Thus the complete reverse support
of `V` is at most `2p-1+2q`, which gives fewer than `4p` coefficients.
The same bound includes `U`.

All computations in (21)--(23) use `O(p)`-degree jets, polynomial
precision, and binary degrees. This proves the complete modular
constructor claimed at the start.

## 6. A nonconstant root-incidence rule in the wider family

Keep `m=2^t`, `h>=t`, positive reflected nodes, and odd `a`, but now
allow any integer `eta` for which `N>0`. The primitive theorem above
still requires its stronger height bound and odd `eta`; this section
does not assert that extension.

For a possible cross pair,

\[
 N-A_iB_j=D\bigl(m\eta+a(j-i)+Dij\bigr).        \tag{24}
\]

Reduction modulo `m` forces `i=j`, because `D` is divisible by `m`,
`a` is odd, and `|i-j|<m`. Equation (24) then vanishes exactly when

\[
 i^2=-m\eta/D\quad\hbox{is an integer square with }0\le i<m.
                                                               \tag{25}
\]

An integer square-root calculation therefore gives the exact gcd:
it is one if (25) fails, and the single monic factor `X-(a+Di)`
if it holds. This is nonconstant incidence, including `eta=0` and
negative parameters producing a retained root. It uses this family's
identity matching permutation, not a general inverse-graph oracle.

For odd `eta` in the primitive theorem, (25) is impossible because
`h>t`. When `eta` has enough powers of two to allow (25), the
constant pair-depth premise behind the perturbation proof changes.
That is the exact boundary of the present normalized constructor.

## 7. Exact controls and scope

The checker compares 84 normalized Newton coefficients with direct
rational divided-difference tables; 54 comparisons use the actual
perturbation path rather than the polynomial-size fallback. It also
compares all 252 coefficients of 27 complete Bezout pairs with
independent rational Lagrange inversion and exact division, and
checks every polynomial identity. Thirty-seven literal domain checks
verify the changing-incidence rule (25).

The implicit controls use `m=2^32`, `2^64`, and `2^200`, at output
precisions 70, 74, and 12. They compare higher-precision outputs and
the complete modular polynomial identities. The first two controls
have leading coefficients different from the translated-grid value
in (6), verifying that reciprocal deformation is actually retained.
The first has two nonzero shifted powers in each Bezout polynomial.
The third verifies binary degree and scale handling at a larger size.
None of these controls lists a progression node.

```sh
python3 tools/check_matched_progression_bezout.py \
  --output evidence/matched-progression-bezout.json
```

These are algebra and representation controls, not factoring timing
records. The constructor handles a specific deeply matched family;
the general singular factor-count contraction remains unresolved.
The proof and implementation were developed with substantial OpenAI
assistant support, without new Lean formalization or an improved
complete factoring exponent.
