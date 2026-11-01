# Implicit inverse contractions with external rational poles

## Constructive result

The complete inverse contraction with weights

\[
 f(a)=\frac1{a-c},\qquad g(b)=\frac1{b-d}
 \tag{1}
\]

can be computed in polynomial bit time on binary-encoded unions of
odd progressions. Here `c,d` are arbitrary integers outside their
respective matrix domains. The denominators may be even on every
odd node. This extends the evaluable weight family beyond the unit
inversions of the [regular circuit compiler](pole-cauchy-analytic-weights.md).

The algorithm uses an exact rational identity and four normalized
affine products. It does not evaluate a formal series outside its
convergence domain. The numerical matrix cardinality is never
enumerated. The required precision is bounded using the
[implicit inverse conditioning certificate](pole-cauchy-domain-precision.md).

This is an algorithm for two separated external poles. The singular
joint kernel depending on `N-ab` is a different input contract; its
factor-count aggregation remains unresolved. There is no new
factoring exponent claim.

## 1. A rational identity for the whole contraction

Let `A,B` be square domains of distinct positive odd nodes, represented
as explicit unions of disjoint step-two progressions. Use

\[
 t=2^kw,\qquad z=N+t,\qquad w^2+w+1=0,
 \qquad M_{a,b}=\frac{t}{z-ab}.
 \tag{2}
\]

The positive odd integer `N` and the displacement `k>=1` satisfy
the existing condition `abs(N-ab)<2^k` for every nonzero difference
on the two matrix domains. Write `m=|A|=|B|`.

For integers `c` outside `A` and `d` outside `B`, put

\[
 Q(c,d)=
 \frac{\displaystyle\prod_{a\in A}(z-ad)
                     \prod_{b\in B}(z-bc)}
      {\displaystyle z^m\prod_{a\in A}(a-c)
                         \prod_{b\in B}(b-d)}.
 \tag{3}
\]

Then the full contraction is

\[
 \boxed{\displaystyle
 F(c,d):=\sum_{a\in A,b\in B}
        \frac{(M^{-1})_{b,a}}{(a-c)(b-d)}
     =\frac{z}{t(cd-z)}\bigl(1-Q(c,d)\bigr).}
 \tag{4}
\]

To prove this identity, use the exact moment-generating identity from
[joint polynomial moments](pole-cauchy-analytic-kernels.md):

\[
 \sum_{a,b}\frac{(M^{-1})_{b,a}}{(1-aX)(1-bY)}
  =\frac zt\,\frac{1-R(zX)T(Y)}{1-zXY},
 \tag{5}
\]

where

\[
 R(U)=\frac{\prod_b(U-b)}{\prod_a(U-z/a)},\qquad
 T(V)=\frac{\prod_a(1-zV/a)}{\prod_b(1-bV)}.
\]

Both sides of (5) are rational functions; their formal expansions
agree at the origin. They therefore agree algebraically. Substitute
`X=1/c,Y=1/d` and simplify to obtain (3)--(4) when `c,d` are nonzero.
Clearing denominators extends the identity to `c=0` or `d=0` as well.
There is no appeal to convergence at these substituted values.

There is also a direct interpolation proof. Put
`P_X(Y)=prod_a(Y-z/a)`, `P_B(Y)=prod_b(Y-b)`. For `c!=0`, define

\[
 H_c(Y)=\frac zt\,
   \frac{1-(P_X(Y)/P_B(Y))R(z/c)}{z-cY}.
\]

The apparent pole at `Y=z/c` is removable. The function is proper
and has only simple poles at the nodes of `B`, so write it as
`sum_b u_b/(Y-b)`. At `Y=z/a`, it equals `a/(t(a-c))`; hence
`Mu` is exactly the row weight in (1). The desired contraction is
`-H_c(d)`. Substituting the products yields (3)--(4). For `c=0`,
use `H_0(Y)=(1-P_X(Y)/P_B(Y))/t` directly.

The factor `cd-z` cannot vanish for integer `c,d`, since `z` has
the nonzero quadratic-basis coefficient `2^k w`. The exclusions
`c not in A` and `d not in B` ensure that all other displayed
denominators are nonzero.

## 2. Exact affine valuations without listing the factors

For an integer affine progression `x_i=u+vi`, the unramified
quadratic extension gives

\[
 v_2(x_i+2^kw)=\min(v_2(x_i),k),
 \tag{6}
\]

including `x_i=0`. At equality of the two valuations, the residue
is `1+w`, which is nonzero in the quadratic residue field. Hence

\[
 \sum_{i=0}^{L-1}v_2(x_i+2^kw)
 =\sum_{\ell=1}^k
       \#\{0\le i<L:u+vi\equiv0\pmod{2^\ell}\}.
 \tag{7}
\]

Each count comes from one linear congruence. Divide by
`gcd(v,2^ell)` when it divides `u`, invert the remaining odd
coefficient, and count the resulting residue progression in
`[0,L)`. A remaining modulus of one is handled directly.
This costs polynomially many bit operations in the
encoded parameters and `k`.

For a nonzero integer product, use the same count through
`bitlength(max(abs(u),abs(u+v(L-1))))`. No integer factor can have
a larger valuation. The external-pole exclusions ensure that these
denominator products contain no zero.

Thus the exact valuation

\[
 q=v_2(Q(c,d))
 \tag{8}
\]

is computed before any product-unit precision is chosen. The
individual product valuations can be proportional to `m`, but
their binary encodings have only polynomial length. Their signed
difference in (8) is exact.

## 3. Normalized affine products at an arbitrary displacement

External poles create a new case: the integer differences `N-ad`
or `N-bc` can have valuation at least `k`, even though `k` separates
the original matrix products. The algorithm handles this case
directly and keeps the original matrix displacement.

First remove the common power of two from the affine coefficients,
clipped at `k`. If the normalized step is even and the normalized
constant is odd, all remaining factors are units with an even
step. The existing [modular interval product](modular-interval-jet-arithmetic.md)
evaluates them with polynomial bit cost.

After removing a common valuation, use the remaining displacement
exponent in place of `k`. For an odd normalized step, partition
the indices into their exact integer valuation classes below it.
The class of valuation `j` is one progression modulo `2^(j+1)`.
After division by `2^j`, its
normalized step is even. Its imaginary coefficient is divisible
by two, and the real coefficient is odd.

The remaining indices have integer valuation at least `k`. After
division by `2^k`, their factors have the form

\[
 a+bj+w.
 \tag{9}
\]

Every one is a unit. If `b` is odd, split `j` by parity; each
subprogression then has an even step. Both possible constant
residues, `w` and `1+w`, are units. This requires only two final
progressions. A constant affine sequence is handled by binary
exponentiation.

The [implementation](../tools/pole_cauchy_rational_poles.py) therefore
evaluates every product in (3) as an exact valuation and a unit
modulo a specified power of two, including external differences
beyond the original displacement and exact integer zeros in the
quadratic numerator factors.

## 4. Conditioning bounds the numerical precision

Let the effective inverse certificate supply `L` with
`2^L M^-1` integral, and define

\[
 u_A=\max_{a\in A}v_2(a-c),\qquad
 u_B=\max_{b\in B}v_2(b-d),\qquad
 s=v_2(cd-z)=\min(v_2(cd-N),k).
 \tag{10}
\]

Use `s=k` when `cd=N`. The two denominator maxima are computed
by binary search on the congruence counts used in (7). They are
bounded by the encoded endpoint and pole lengths.

Every summand of `F(c,d)` has valuation at least `-L-u_A-u_B`.
Summing the whole matrix does not lose 2-adic precision. Since
the prefactor in (4) has valuation `-k-s`,

\[
 v_2(1-Q)\ge k+s-L-u_A-u_B,
 \qquad
 \boxed{q\ge\min(0,k+s-L-u_A-u_B).}
 \tag{11}
\]

This proves that a negative product valuation cannot conceal a
loss proportional to the numerical matrix dimension.

For clarity, the accepted progression representation also bounds
`L` itself uniformly. If the row and column domains contain
`c_A,c_B` explicit pieces, each residue count differs from its
domain's average by less than its piece count. Thus the range
certificate has `epsilon_ell<=c_A+c_B-1`, with `epsilon_1=0`, and

\[
 L\le(k-1)\max(0,2c_A+2c_B-3).
 \tag{12}
\]

For a requested absolute precision `p`, put `T=p+k+s` and
`e=min(0,q)`. The program returns

\[
 F(c,d)=2^{e-k-s}R+O(2^p),
 \tag{13}
\]

with `R` known modulo `2^(T-e)`. To compute it, write

\[
 1-Q=2^e\bigl(2^{-e}-2^{q-e}U_Q\bigr).
 \tag{14}
\]

If `q<T`, the unit `U_Q` is needed only modulo `2^(T-q)`.
Changing its lift at that precision changes the bracket in (14)
by a multiple of `2^(T-e)`. If `q>=T`, the entire `Q` term vanishes
at the requested precision, and all four unit-product evaluations
are skipped. A huge positive valuation is never used as a literal
shift-allocation size.

Equations (11)--(12) bound every required numerical precision by
a polynomial in the explicit input lengths, `k`, and `p`. There
are only polynomially many affine valuation classes, and each
unit progression is evaluated by the existing modular product
algorithm. The factor `z^m` is a unit power computed from the
binary encoding of `m`. This proves the claimed polynomial bit
cost, with no matrix-node enumeration.

## 5. Relation to regular residue selectors

When a pole is odd, its denominator is even at all odd nodes.
The regular circuit compiler rejects that inverse gate. The
algorithm here instead uses the domain-specific fact that the
pole lies outside the chosen nodes and evaluates the exact
rational product identity (4).

Accordingly the [regular selector precision theorem](regular-analytic-selector-precision.md)
continues to hold with its stated unit-inverse gate restriction.
This new algorithm admits a different family of singular rational
weights. It does not evaluate an arbitrary nonregular circuit, or
the two-variable product pole needed for the factor-count trace.

## 6. Exact finite checks and implicit control

The [checker](../tools/check_pole_cauchy_rational_poles.py) compares
normalized affine products with literal exact products, including
integer zeros in quadratic factors and valuations above the
displacement. It verifies the rational identity (4) exactly and
compares complete finite inverse contractions with dense inversion.
Controls include zero poles, negative poles, `cd=N`, progression
unions, and external differences whose valuations exceed `k`.

The large control uses

\[
 N=2^{401}+45,\quad m=2^{200},\quad
 A=\{1,3,\ldots,2m-1\},\quad B=\{3,5,\ldots,2m+1\},
 \qquad c=-1,\quad d=2m+3.
\]

The poles are outside both domains, while both weights have even
denominators everywhere on their domains. Their maximum denominator
valuations are both 201. The computation uses exact valuation
separation and checks a twelve-bit answer against a sixteen-bit
replay. The [evidence](../evidence/pole-cauchy-rational-poles.json)
records the full precision certificate, residues, and source hashes.

```sh
python3 tools/check_pole_cauchy_rational_poles.py \
  --output evidence/pole-cauchy-rational-poles.json
```

Use `--small-only` to run the independent finite controls without
the large implicit example.
