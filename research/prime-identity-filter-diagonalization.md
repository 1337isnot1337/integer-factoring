# Exact finite-difference profiles for prime-state identity filters

## Result and scope

For every prime characteristic, the first identity-block summand of a
truncated symmetric Pascal interval has an exact, fixed-size zero
profile. Its active tensor has two modes, each with \(p\) states.
After a finite-difference transform of the digit weights, every charge
block is diagonalized by explicit polynomial unimodular matrices.

The result applies to arbitrary Laurent digit weights, including
specializations with tied valuations and exact cancellations. It
recovers the original Gaussian spectrum when every digit weight is
one. It supplies an exact profile formula for new searches; it does
not assert an improved factoring exponent.

## 1. The universal active tensor

Work over a field of characteristic \(p\), with valuation normalized
by \(\operatorname{ord}(z)=1\). Rational exponents can be handled in
a fixed Puiseux extension. Equivalently, clear their denominators by
one fixed ramification and then divide all orders by that common
positive factor.

In the algebra where \(Y^p=X^p=U^p=V^p=0\), put

~~~math
E(H)=\sum_{d=0}^{p-1}\frac{H^d}{d!},\qquad
A=z^a,\quad B=z^b,\quad T=AB,\quad c=1-T,
\qquad a+b>0.
\tag{1}
~~~

Every factorial appearing in a denominator is a unit. Consider

~~~math
G=E(AYU)E(BXV),\qquad
Q=(AY+cV)(U-A^{-1}X).
\tag{2}
~~~

For digit weights \(\alpha_0,\ldots,\alpha_{p-1}\), define

~~~math
\mathcal F_\alpha(Q)
=E(-Q)\sum_{d=0}^{p-1}\alpha_d\frac{Q^d}{d!},
\qquad \widetilde G=\mathcal F_\alpha(Q)G.
\tag{3}
~~~

Because \(Q^p=0\), (3) is an ordinary polynomial in \(Q\).
For \(p=3\), it is
\(1+(\alpha_1-1)Q+(2\alpha_2-\alpha_1+2)Q^2\) when
\(\alpha_0=1\).

The [regular concentration proof](qutrit-identity-zero-reduction.md)
uses only signed Pascal identities, elementary matrix shears, and
the truncated symmetric functor. Every step in its Sections 1 and 2
therefore holds in every characteristic, with the same signs. For
an interval of length \(n\), odd ramification \(P\), and aggregate
zero gauge \(J\), it gives (2) with

~~~math
a=P(1-J),\qquad b=P(2n-1+J).
\tag{4}
~~~

The remaining Gaussian spectators have orders

~~~math
P(2i+1-J)\ (1\le i<n),\qquad
P(2i+1+J)\ (0\le i<n-1).
\tag{5}
~~~

Each order \(k\) contributes its \(p\) digit orders
\(0,k,\ldots,(p-1)k\). The finite profile below is consequently
enough to determine the full interval zero profile.

If every \(\alpha_d\) is regular at the crossing endpoint and
equals one there, the multiplier (3) equals one at that endpoint.
This follows from \(E(-Q)E(Q)=1\). Thus it preserves the existing
typed crossing endpoint. Laurent monomials at zero can be cleared
by a common monomial with endpoint value one, as in the other
Laurent crossing constructions. Endpoint preservation and the local
profile formula are separate from any pressure improvement.

## 2. Removing regular internal factors

For the zero cut, rows are indexed by the exponents of \((Y,X)\),
and columns by those of \((U,V)\). Put \(Y'=AY\) and
\(X'=A^{-1}X\). A row monomial of charge

~~~math
q=\deg Y-\deg X=\deg U-\deg V
\tag{6}
~~~

acquires exactly the scalar \(A^q\). We retain this scalar
explicitly; the potentially nonregular substitution is not treated
as a valuation-preserving change of basis. Omit primes temporarily.
The universalized Gaussian is \(E(YU)E(TVX)\), and
\(Q=(Y+cV)(U-X)\).

The Gaussian multiplication identity gives

~~~math
E(YU)E(TVX)E(-Q)=E(YX)E(-cVU)E(VX).
\tag{7}
~~~

One may prove this without dividing by any multiple of \(p\).
Regroup the Gaussian products by the original row variables
\(Y,V\). The addition formula for truncated exponentials applies
to terms sharing one such variable: every omitted term contains
its \(p\)-th power. A rank-one bilinear Gaussian associated with
\(Q\) is exactly \(E(Q)\), so the same identity applies to
\(E(-Q)\).

The first factor on the right of (7) is internal to the row side;
the second is internal to the column side. Their coefficients are
regular, and their constant terms in the nilpotent variables are
one. Multiplication by either factor has a regular inverse. Since
\(Y'X'=YX\), the first factor is also regular in the original
variables. Removing the two factors leaves the exact coefficient
kernel

~~~math
K_\alpha
=\sum_{d=0}^{p-1}\frac{\alpha_d}{d!}
 (Y+cV)^d(U-X)^d E(VX).
\tag{8}
~~~

This removal preserves the charge blocks. It introduces no
unpaid basis change at zero.

## 3. The finite-difference weights

Define

~~~math
\boxed{\quad
\beta_k=\sum_{h=0}^{k}(-1)^{k-h}\binom{k}{h}
 c^h\alpha_{p-1-k+h},\qquad 0\le k<p.
\quad}
\tag{9}
~~~

The inverse transform is

~~~math
\alpha_d=\sum_{k=0}^{p-1-d}(-1)^k
 \binom{p-1-d}{k}c^{p-1-d-k}\beta_k.
\tag{10}
~~~

For example, (10) follows by the ordinary binomial inversion after
reversing the sequence \(\alpha_d\). Both transformations are
triangular with constant unit diagonal and have polynomial entries.

Fix a charge \(-(p-1)\le q\le p-1\). Its row and column index
sets are

~~~math
I_q=\{\max(0,q),\ldots,\min(p-1,p-1+q)\}.
\tag{11}
~~~

Row \(i\) means \((\deg Y,\deg X)=(i,i-q)\); column
\(j\) means \((\deg U,\deg V)=(j,j-q)\). Expanding (8) gives

~~~math
K_q(i,j)=
\sum_{d=\max(i,j)}^{\min(p-1,i+j-q)}
\alpha_d c^{d-i}(-1)^{d-j}
\frac{d!}{i!j!(d-i)!(d-j)!(i+j-q-d)!}.
\tag{12}
~~~

An empty sum is zero. All factorials in denominators in (12) have
arguments between zero and \(p-1\).

## 4. Explicit Pascal factors and their unit determinants

Let

~~~math
J_q=\{\max(0,-q),\ldots,\min(p-1,p-1-q)\}.
\tag{13}
~~~

It has the same cardinality as \(I_q\). Define matrices with row
and column sets \(I_q,J_q\) and \(J_q,I_q\), respectively:

~~~math
L_q(i,k)=
\begin{cases}
\displaystyle
(-1)^{i-q}c^{p-1-i-k}
\frac{(i+k)!}{i!(i-q)!},&i+k<p,\\[4pt]
0,&i+k\ge p,
\end{cases}
\tag{14}
~~~

~~~math
R_q(k,j)=
\begin{cases}
\displaystyle
\frac{(j+k)!}{k!(k+q)!j!(j-q)!},&j+k<p,\\[4pt]
0,&j+k\ge p.
\end{cases}
\tag{15}
~~~

Then the exact polynomial matrix identity is

~~~math
\boxed{\quad
K_q=L_q\operatorname{diag}(\beta_k:k\in J_q)R_q.
\quad}
\tag{16}
~~~

### Proof of the identity

Substitute (10) into (12). If \(i+k\ge p\), the coefficient
of \(\beta_k\) is already zero: its original summation would
require both \(d\ge i\) and \(d\le p-1-k\). Otherwise it
contains the factor \(c^{p-1-i-k}\), with a nonnegative exponent
independent of \(d\). Extract that factor and use

~~~math
\binom{p-1-d}{k}=(-1)^k\binom{d+k}{k}
\quad\text{in characteristic }p.
\tag{17}
~~~

Both sides vanish when \(d+k\ge p\), in the range at issue.
Writing \(N=i-q\), \(t=d-j\), and \(m=i+k\), the remaining
sum reduces to

~~~math
\sum_{t=0}^{N}(-1)^t\binom Nt
 (j-i+t+1)_m
=(-1)^N\frac{m!}{(m-N)!}
 (j-q+1)_{m-N}.
\tag{18}
~~~

Here \((x)_m=x(x+1)\cdots(x+m-1)\). Equation (18) is the
standard finite-difference identity for a rising factorial. If
\(m<N\), the left side is zero. Terms corresponding to an
originally negative factorial denominator vanish because the rising
product contains zero. Terms with \(d\ge p\), added to complete
the finite-difference sum, contain a factor \(p\) in their
numerator, while all remaining factorial denominators are units.

When \(k+q=m-N\ge0\), (18) yields the coefficient

~~~math
(-1)^{i-q}c^{p-1-i-k}
\frac{(i+k)!(j+k)!}
 {k!(k+q)!i!j!(i-q)!(j-q)!},
\tag{19}
~~~

which is exactly (14) times (15). If \(i+k\ge p\) or
\(j+k\ge p\), the coefficient is zero. If \(k\notin J_q\),
either the finite difference vanishes or every row factor vanishes.
These observations prove (16), including all boundary ranges.

### The factors preserve valuations

Order both index sets increasingly. The matrices \(L_q,R_q\)
are triangular after reversing one index order. Their antidiagonals
have \(i+k=p-1\) and \(j+k=p-1\). The power of \(c\) in
(14) is then zero, and every displayed factorial is a unit.
Their determinants are therefore nonzero constants in the prime
field. Thus both matrices are unimodular over \(\mathbb F_p[c]\).

In particular, (16) remains a regular equivalence after substituting
any regular \(c\). It does not require \(c\) to be invertible,
and it makes no genericity assumption about the digit weights.

## 5. The complete active zero spectrum

Put \(b_k=\operatorname{ord}(\beta_k)\), allowing
\(b_k=+\infty\) for an identically zero weight. Restoring the
charge scalar \(A^q\), (16) gives the complete profile

~~~math
\boxed{\quad
\{\,a(r-k)+b_k:0\le r,k<p\,\}.
\quad}
\tag{20}
~~~

Indeed, \(\beta_k\) occurs in exactly the \(p\) charges
\(-k,-k+1,\ldots,p-1-k\). Infinite orders describe zero
directions; the other entries are all finite Smith orders.

For positive \(s\), the active partition function is consequently

~~~math
Z_{0,\rm act}(s)
=\left(\sum_{r=0}^{p-1}s^{-ar}\right)
 \left(\sum_{k=0}^{p-1}s^{ak-b_k}\right).
\tag{21}
~~~

One then multiplies by the explicit spectator factors from (5).
Any input-factorial normalization of the truncated symmetric kernel
changes the cut only by constant diagonal units and does not affect
these orders.

If all \(\alpha_d=1\), equation (9) gives
\(\beta_k=(c-1)^k=(-T)^k\). Formula (20) becomes
\(\{ar+bk:0\le r,k<p\}\), exactly the two unperturbed Gaussian
digit modes. If only the highest digit has a pole of order \(h>0\)
and all other digit weights are regular, every \(\beta_k\) has
order \(-h\), since its coefficient of \(\alpha_{p-1}\) is
\(c^k\), a unit. The active spectrum then becomes
\(\{a(r-k)-h\}\); the high positive Gaussian order is lost.

For monomial or polynomial digit weights, each \(\beta_k\) is
explicit. Equal powers of \(z\) must be combined in the prime
field before taking a valuation. The resulting calculation is exact
at ties and uses no series-precision cutoff. Independent direct
truncated-algebra expansions verify the polynomial identities for
\(p=2,3,5,7,11\); the proof above supplies the uniform result.

The [prime-state profile implementation](../tools/prime_state_filter_profiles.py)
evaluates (9) and (20) by exact exponent aggregation. The independent
[checker](../tools/check_prime_state_filter_profiles.py) expands the
original truncated algebra, verifies the polynomial factors, compares
complete finite Smith profiles, and checks physical crossing endpoints.
Its [recorded evidence](../evidence/prime-state-filter-profiles.json)
also includes complete original cut comparisons, so it checks the
factorial normalization and both valuation places independently of the
closed profile formula.
