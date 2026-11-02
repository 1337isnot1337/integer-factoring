# Rough prime-power targets obstruct prescribed one-shot smooth descent

## Result and stopped scope

The [near-height enumeration theorem](near-height-smooth-descent.md) does
not acquire a uniform success theorem by increasing its scalar height.
The obstruction occurs at actual prime-power targets on composite inputs
that survive ordinary preprocessing, rather than only at arbitrary classes.

Fix constants `C>0`, `alpha>0`, and for sufficiently large `D` put

```math
F(D)=\sqrt D\exp\left(\frac{C\log D}{(\log\log D)^\alpha}\right).
                                                               \tag{1}
```

Let `Z(D)>=3` be a nondecreasing cutoff satisfying

```math
\log Z(D)\,\log\log Z(D)=o(\log D),                    \tag{2}
```

and let the smoothness cutoff satisfy `3<=Y(D)<=Z(D)`. Then there are
arbitrarily large odd composite, non-perfect-power integers `N` such that:

- `N` has no prime divisor at most `Z(N)` and `N==3 (mod 4)`;
- `3` is an unramified split prime at discriminant `Delta=-N`;
- if `k` is the least integer with `3^k>ceil(F(N))`, then
  `red(g_3^k)` represents **no** `Y(N)`-smooth positive norm below `3^k`,
  even allowing imprimitive representations.

The target orientation is the canonical root one at norm three. Thus the
usual prescribed near-height strict-descent statement is **false**, for
every fixed scalar constant and logarithmic saving exponent in (1).
The complete ellipse remains cheap to enumerate; its points can all lie
on `y=0` and carry the same rough coefficient. No positive lower bound on
primitive smooth candidates follows from its area or the known boundary norm.

Condition (2) includes the quasipolynomial smoothness cutoffs used in the
current density route. It also permits a larger trial-division cutoff,
for example `Z(D)=exp(sqrt(log D))` times a fixed constant. The examples
need not survive arbitrary preprocessing of cost comparable to the full
factoring bound. This is not a factoring lower bound, an obstruction to
whole orbit windows, or an exclusion of adaptive/multiple target choices.
In particular, failing one target need not prevent a useful partial
presentation. The complete proven score is unchanged at **1/2**.

The [construction question and pivot](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359793)
were recorded before the lower-sieve proof. This stops the proposed
prescribed one-target descent mechanism; changing only `C`, `alpha` or
the cutoff inside this scope is not a continuation. Further work needs
collective or adaptive targets, another source of relations, or a different
arithmetic representation.

A subsequent [classical-sieve extension](full-trial-sieve-obstruction.md)
replaces condition (2) by `log Z=o(log D)`, covering every subexponential
trial-division cutoff, including full positive-score trial budgets. The
elementary proof here and its narrower source scope are retained. The
extension still excludes no adaptive/collective targets or arbitrary
preprocessing of the same cost.

## 1. A rigid prime-power target at any near-height ratio

Let `m=3^k` and `1<a<m`. Set

```math
N=4ma-1,\qquad Q=(a,-1,m),\qquad\Delta=-N.
```

The form is primitive and strict reduced. Since `Delta==1 (mod 4m)`,
the oriented norm-`m` form `(m,1,a)` represents the `k`th power of the
unramified norm-three ideal with root one. One reduction swap gives
`Q=red(g_3^k)`. The primitive vector `(0,1)` supplies the known smooth
boundary norm `m`.

Every smaller positive norm has `y=0`. For `|y|=1`,
`a x^2-xy>=0`. For `|y|>=2`, completion of the square gives

```math
Q(x,y)\ge(m-1/(4a))y^2>m.
```

At `y=0` all norms are divisible by `a`, and only `x=+/-1` is primitive.
If `a` has no prime divisor at most `Z>=Y`, no norm below `m` is
`Y`-smooth. Increasing `m/sqrt(N)` can add imprimitive vectors on that
line, but it adds no primitive candidate below the cap.

## 2. An elementary lower-sieve estimate

Here is a sufficient uniform estimate, proved directly to avoid an
unidentified sieve or prime-distribution hypothesis. Fix `m=3^k`, a
cutoff `Z>=3`, a prime `q>Z`, and an arithmetic progression

```math
a=a_0+q^2t,\qquad 4ma_0-1\equiv q\pmod{q^2}.           \tag{3}
```

For each prime `r<=Z`, exclude the residues for which `r` divides
`a(4ma-1)`. There is one forbidden residue at `r=2,3` and two at every
other prime. The residues are distinct and the progression step is a
unit at every such prime. Write these counts as `nu_r`, and put

```math
\mu=\sum_{r<=Z}\nu_r/r,\quad
P=\prod_{r<=Z}(1-\nu_r/r),\quad
U=6\log(1+\log Z).
```

First `mu<=U`. For `s=1+1/log Z`, each `r<=Z` satisfies
`r^(-s)>=e^(-1)/r`. The finite Euler product is at most `zeta(s)`,
and the integral comparison gives `zeta(s)<=1+log Z`. Hence
`sum_(r<=Z)1/r<=e log(1+log Z)`, proving the assertion since `2e<6`.
Also `nu_r/r<=2/3`, and `log(1-u)>=-3u` on that interval, so

```math
P\ge e^{-3\mu}\ge e^{-3U}.                              \tag{4}
```

Choose an even integer `J>=32(U+1)` with `J=O(U+1)`, and let
`L=J-1` be odd. If `e_j` is the degree-`j` elementary symmetric sum
of the numbers `nu_r/r`, Bonferroni gives

```math
P-e_J\le\sum_{j=0}^L(-1)^j e_j\le P,
\quad e_J\le\mu^J/J!\le(e\mu/J)^J\le\tfrac12e^{-3U}.
                                                               \tag{5}
```

If there are fewer than `J` primes, take `e_J=0`; the same statement
holds. The factorial estimate follows from summing `log j` or its integral.
The oversized constant 32 makes the last inequality immediate, uniformly
also for small nonnegative `U`.

On any interval of `X` consecutive integers `t`, simultaneous forbidden
conditions for a squarefree product `d` occupy `nu(d)=prod_(r|d)nu_r`
residues. Their exact count is
`X nu(d)/d+O(nu(d))`, with absolute error at most `nu(d)` by CRT.
The total absolute rounding error through order `L` is at most

```math
\sum_{j=0}^L\binom{\pi(Z)}j2^j\le(2Z)^{L+1}.
```

Applying the odd Bonferroni truncation to the actual excluded events,
(4)--(5) therefore give the explicit lower bound

```math
\#\{t:\gcd(a(4ma-1),\prod_{r<=Z}r)=1\}
\ge\tfrac12 X e^{-3U}-(2Z)^{L+1}.                       \tag{6}
```

The error has logarithm `O(log Z log log Z)`. In particular, if
`X=m^(1-o(1))` and `log Z log log Z=o(log m)`, the right side is
positive for sufficiently large `m`. This is an existence proof, not an
efficient uniform factoring sampler: exhaustive construction in this
long interval is not assigned subexponential cost.

## 3. Aligning the target with the prescribed exponent

For `m=3^k` sufficiently large, let `D_0` be the unique large real solution
of `F(D_0)=m/2`. The logarithmic derivative of (1) tends to `1/2`, so
`F` is eventually increasing and

```math
D_0=m^{2-o(1)},\qquad
F(2D_0)/F(D_0)=\sqrt2+o(1).
```

For every integer `N` in `[D_0,2D_0]`, eventually

```math
m/3<F(N)\le\lceil F(N)\rceil<m.
```

Thus `k` is exactly the prescribed least exponent, with a fixed margin
at both endpoints. Choose `Z_* = ceil(Z(ceil(2D_0)))` and a prime
`q` in `(Z_*,2Z_*]`, using Bertrand's postulate. Since `q!=3`, choose
the residue in (3), which forces `v_q(4ma-1)=1`.

Restrict `a` to the integer interval giving `D_0<=4ma-1<=2D_0`.
It has length `m^(1-o(1))`; after restricting modulo `q^2`, the
corresponding interval of `t` still has length `X=m^(1-o(1))`.
Condition (2) supplies exactly the error condition for (6), and also
`Z_*=m^o(1)`. Hence choose an `a` for which both `a` and `N=4ma-1`
have no prime factor at most `Z_*`. Eventually `1<a<m`.

Now `N>q` and `v_q(N)=1`, so `N` is composite and not a perfect power.
It has no prime factor at most `Z(N)`, and `N==3 (mod 4)`, so the
ordinary factoring discriminant is exactly `-N`. The target in Section 1
has no smaller `Y(N)`-smooth norm, proving the result. The proof does not
need `a` prime, a factorization of `a`, or an assumption about class-group
orders or indices.

## 4. Exact controls and attribution

The [near-height checker](../tools/check_near_height_smooth_descent.py)
and [receipt](../evidence/near-height-smooth-descent.json) retain two
predeclared bounded CRT/sieve constructions, with `m=3^24,3^40` and
cutoffs `Z=1009,10007`. They reject four and 109 CRT candidates,
respectively, before finding the stored examples. Gcd against the full
primorial verifies both roughness conditions. The supplied proper primes
`1013,10009` have valuation one, certifying compositeness and excluding
perfect powers. Both examples survive trial division through
`2^ceil(sqrt(bits(N)))`, as well as through their smoothness cutoffs.
Their complete strict ellipses contain 44 and 181 vectors modulo sign,
but only one primitive vector in each, and no smooth success.

The coordinating worker supplied two additional exact **semiprime,
fundamental-discriminant** controls in
[Discussion 51, comment 19359806](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359806)
and [comment 19359827](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359827):
`7423=13*571`, `m=64`, `a=29`, and `26623=79*337`, `m=512`, `a=13`,
both at cutoff 11. The second has `m>sqrt(N)`. The checker independently
reconstructs the prime-two powers, verifies the supplied prime factors,
and exhausts primitive norms through the cap: only `a` and `m` occur.
The original private receipt hashes and discussion credit are retained
in the combined receipt. These controls are not an infinite semiprime
or fundamental-discriminant theorem; the lower-sieve construction above
only asserts rough composite non-perfect-power inputs.

A separate exact rational panel checks all 256 CRT intersections at
primes through 19, compares with 257 direct residue tests, and verifies
the odd Bonferroni truncations through order seven and their rounding
envelopes. It checks the finite combinatorics, not the asymptotic theorem
by numerics. Run `make verify-near-height-smooth-descent`.

Bonferroni, CRT, the Euler product, Bertrand's postulate and square
completion are classical. Sol worker 1 developed the combined arithmetic
counterfamily and prescribed-target alignment with substantial OpenAI
assistant support. Existing form arithmetic and norm-decoder credit is
preserved, as is the coordinating worker's contributed semiprime panel.
This is an author-checked ordinary proof with finite controls,
without independent mathematical review, new Lean or a better complete
factoring bound.

The exact lower-sieve, CRT and semiprime controls pass in the focused
target at `5fbb3f1`. The earlier
[full default replay](../evidence/near-height-descent-full-verification.json)
at `7863baa` predates those additions; no full replay of the follow-up is
claimed. The [ledger](../docs/REVIEW_LEDGER.md) distinguishes these scopes.
