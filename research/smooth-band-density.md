# Smooth forms from several prime bands

This note strengthens the uniform smooth-form density bound under the
[pinned finite-order Hecke zero-free input](polylog-prime-form-generators.md).
It counts recognizable forms with almost the minimum possible number
of prime factors for their size. The construction gives a sharper
prime cutoff and short-word budget. It does not supply a fast algorithm
to find a smooth norm in a specified form class.

Write `Delta < 0` for a quadratic-order discriminant, `D = abs(Delta)`,
`b = bits(D)`, and `G = C_Delta`. All logarithms in the proof are natural
unless a base is displayed. The notation

$$L_D[1/2,c]=\exp\bigl((c+o(1))\sqrt{\log D\log\log D}\bigr)$$

describes a parameter size, not a proved factoring running time.

The main quantitative consequence is the following. For every
sufficiently large discriminant, the existing biased-prime translate
argument gives a smooth norm in every class with

$$Z=L_D[1/2,\sqrt8],\qquad
 \Omega(n)\le L=(1/\sqrt8+o(1))
       \sqrt{\frac{\log D}{\log\log D}},\qquad
 n\le D^{1+o(1)}.\tag{1}$$

Here all prime factors of `n` are split and coprime to `Delta`.
The same bound on `Omega(n)` applies to its signed prime-form word.
An explicit integer recipe appears in Section 4. The needed norm
selection algorithm remains the missing obligation.

## 1. Uniform split primes in a short multiplicative interval

Let

$$w(t)=\begin{cases}t(1-t)&0\le t\le1,\\0&\text{otherwise},\end{cases}
 \qquad \lambda=16/9.$$

The principal-character estimate in
[the explicit-cutoff note](explicit-prime-form-cutoffs.md), Section 3,
is

$$A_1(x)=x/6+E(x),$$

$$|E(x)|\le100x^{31/32}(b+1)
 +\tfrac12\sqrt x\log x+\tfrac12(b+1)\log x.\tag{2}$$

The sum `A_1` counts both orientations of every split odd prime
coprime to `Delta`, with weight `log(p) w(p/x)`.

Use the nonnegative compactly supported combination

$$\psi(t)=w(t)-\frac{21}{16}w(4t/3)
                  +\frac{27}{64}w(16t/9).\tag{3}$$

Direct polynomial calculation gives

$$\psi(t)=\begin{cases}
0&0\le t\le9/16,\\
\tfrac43t(t-9/16)&9/16\le t\le3/4,\\
t(1-t)&3/4\le t\le1,\\
0&t\ge1.
\end{cases}$$

In particular, `0 <= psi <= 1/4`. Its integral is `259/6144`.
The corresponding combination of (2) therefore gives

$$\sum_{p\ {
m split},\ p\nmid2\Delta}
       2\log p\,\psi(p/x)
 =\frac{259}{6144}x+E_\psi(x),$$

$$|E_\psi(x)|\le300x^{31/32}(b+1)
       +\tfrac32(\sqrt x+b+1)\log x.\tag{4}$$

The sum of the absolute coefficients in (3) is `175/64 < 3`,
which justifies this common error envelope.

Put

$$T=\bigl(10^6(b+1)\bigr)^{32}.$$

For `x = z^32 >= T`, divide the error in (4) by `x` and use
`log z <= z` to get

$$\frac{|E_\psi(x)|}{x}
 \le\frac{300(b+1)}z+\frac{48}{z^{15}}
                  +\frac{48(b+1)}{z^{31}}<\frac1{1000}.$$

Since `259/6144 - 1/1000 > 1/32`, the weighted sum in (4) exceeds
`x/32`. Each prime contributes at most `(log x)/2`, so there are
at least `x/(16 log x)` eligible primes in `(x/lambda,x)`.
We will use the weaker bound

$$\#\{p\in(x/\lambda,x):p\text{ split},\ p\nmid2\Delta\}
       \ge\frac{x}{32\log x}\qquad(x\ge T).\tag{5}$$

There is also a reciprocal-mass consequence. If `B >= T` and
`log(C/B) >= log(lambda)`, tile `(B,C)` by
`J = floor(log(C/B)/log(lambda))` complete intervals of ratio
`lambda`, starting from `C`. In each interval (5) gives reciprocal
prime mass at least `1/(32 log C)`. Since
`J >= log(C/B)/(2 log(lambda))` and `log(lambda) < 1`,

$$\sum_{\substack{B<p<C\\p\ {
m split},\ p\nmid2\Delta}}\frac1p
 \ge\frac{\log(C/B)}{64\log C}.\tag{6}$$

Endpoint primes can be omitted. They have zero weight in (3),
and the complete open intervals used above are disjoint.

## 2. A finite smooth-form density theorem

Suppose `R,Y > 1`, `4R^2 < D`, and put

$$y=\log Y,\qquad u=\frac{\log R}{y}\ge4.$$

Choose `0 < eta <= 1/2`, set `s = ceil((1+eta)u)`, and assume

$$\frac{\eta y}{4s}\ge\log(16/9),\qquad
 Y^{1/4}\ge T,\qquad
 Y^{1/4}\ge\frac{512s^2}{\eta},\qquad
 Y^{5/8}\ge64sy.\tag{7}$$

Let `S` consist of the reduced primitive forms with leading
coefficient `a <= R` such that `a` is squarefree, coprime to
`Delta`, and has exactly `s+1` prime factors, all at most `Y`.
Then

$$|S|\ge
 \frac{R}{32(s+1)s^2\log Y}
       \frac{(\eta/256)^s}{s!}.\tag{8}$$

### Prime bands and their products

Let `ell` index a largest coefficient of

$$P(z)=(1+z+\cdots+z^{s-1})^s.$$

There are at most `s^2` coefficients and `P(1)=s^s`, so

$$[z^{\ell}]P(z)\ge s^{s-2}.\tag{9}$$

Define

$$w_0=\frac{\eta y}{4s},\qquad
 a_0=\frac{\log R-3y/4-w_0\ell}{s}.$$

For each `j = 0,...,s-1`, use eligible primes in the open band

$$\exp(a_0+jw_0)<p<\exp(a_0+(j+1)w_0).\tag{10}$$

Every band is contained in `(Y^(1/4),Y)`. To verify the lower
endpoint, use `0 <= ell <= s(s-1)`, `u >= 4`, and
`s <= 7u/4` to obtain

$$\frac{a_0}{y}
 \ge\frac{u-3/4}{s}-\frac\eta4
 \ge\frac{13}{28}-\frac18>\frac14.$$

For the upper endpoint,

$$\frac{a_0+sw_0}{y}
 \le\frac{u}{s}+\frac\eta4
 \le\frac1{1+\eta}+\frac\eta4<1.$$

By (6), (7), and (10), each band's reciprocal prime mass is
at least

$$\alpha=\frac{w_0}{64y}=\frac\eta{256s}.\tag{11}$$

Consider ordered band-index patterns `(j_1,...,j_s)` whose sum
is `ell`. There are at least `s^(s-2)` patterns by (9).
For each pattern choose distinct primes from the indicated bands.
Deleting at most `s-1` previously chosen primes removes reciprocal
mass at most `s Y^(-1/4) <= alpha/2`, by (7). Consequently
the total reciprocal weight of the ordered prime choices for each
pattern is at least `(alpha/2)^s`.

Each unordered product occurs exactly `s!` times across these
ordered patterns. Thus the resulting squarefree products `m`
satisfy

$$W:=\sum_m\frac1m
 \ge\frac{(\eta/512)^s}{s!s^2}.\tag{12}$$

The prescribed index sum also fixes their size:

$$\frac{R}{Y^{3/4}}<m<\frac{R}{Y^{3/4-\eta/4}}.\tag{13}$$

### One final prime and all ideal orientations

For each `m`, put `t = R/m`. Equation (13) gives

$$Y^{5/8}\le Y^{3/4-\eta/4}<t<Y^{3/4}.$$

By (5), at least `t/(32y)` eligible primes lie in `(t/lambda,t)`.
The last condition in (7) permits deleting all primes dividing `m`
while retaining at least `t/(64y) = R/(64my)` choices.
For every retained prime `p_*`, the product `n = m p_*` is
squarefree, lies in `(R/lambda,R)`, has `s+1` factors, and is
`Y`-smooth and coprime to `Delta`.

A fixed `n` can arise from at most `s+1` choices of the final
prime `p_*`. Therefore (12) produces at least

$$\frac{R}{64(s+1)y}
       \frac{(\eta/512)^s}{s!s^2}$$

different values of `n`.

Every value gives `2^(s+1)` different proper ideal classes by
choosing an orientation at each split prime. More explicitly,
the choices give different residues `B` modulo `2n`, with the
required parity and `B^2 = Delta (mod 4n)`. Center the middle
coefficient `B` modulo `2n`; then `abs(B) <= n` and

$$C=\frac{B^2+D}{4n}>n,$$

because `4n^2 < D`. Thus `(n,B,C)` is reduced. The classes are
distinct: unique reduced representatives identify the orientations,
and different `n` give different leading coefficients. Since all
primes are unramified and odd, the boundary identification
`B = +/- n` cannot occur. Multiplying by `2^(s+1)` proves (8).

The set `S` is recognizable by factoring its leading coefficient
through primes up to `Y`. This recognition statement supplies no
efficient method to locate an element in an arbitrary translate.

## 3. Density with the leading smoothness term retained

The order class-number bound needed here follows from
[Adzaga et al., Lemma 4.3](https://people.math.wisc.edu/~chidambaram3/papers/QCPlusQuotients.pdf#page=9),
which gives `h(Delta) <= sqrt(D)(log(4D)+2)/pi` for quadratic orders,
including nonfundamental discriminants. Since `D < 2^b`, `log 2 < 1`,
`log 4 < 2`, and `pi > 3`, this implies

$$h(\Delta)<\sqrt D\,(b+2).\tag{14}$$

Indeed, the preceding bound is strictly smaller than
`sqrt(D)(b+4)/3`, and `(b+4)/3 < b+2`. This proves the stated estimate for the binary input length.

Take `R` to be a fixed positive fraction of `sqrt(D)`, still
with `4R^2 < D`. Choose

$$\eta\asymp\frac1{\log u},\qquad
 s=\lceil(1+\eta)u\rceil.$$

When (7) holds, (8) and (14), together with `s! <= s^s`, give

$$\log\frac1\delta
 \le u\bigl(\log u+\log\log u+O(1)\bigr)+O(\log b),
 \qquad \delta=|S|/h(\Delta).\tag{15}$$

The `u log u` term has coefficient one. The number of factors
is `s+1 = (1+o(1))u`. Restricting all factors to one short band
would lose a power of the number of bands. Selecting a common
sum of band indices in (9) retains this leading term.

## 4. Parameters computed entirely with integers

The following recipe proves a uniform finite version once its
displayed eligibility checks pass. It does not require the
factorization of `D`, the class number, or the group structure.
Let

$$t=\lceil\log_2 b\rceil,\qquad
 v=\lceil\sqrt{8bt}\rceil,\qquad
 r=\lceil b/2\rceil-3,$$

$$Y=2^v,\qquad R=2^r.$$

Require `r >= 4v`. Let `q` be the smallest integer with
`r <= v 2^q`, and put

$$\eta=1/q,\qquad
 s=\left\lceil\frac{(q+1)r}{qv}\right\rceil.$$

Then `q >= 2`, so the required range for `eta` holds. The
following integer comparisons suffice for (7):

$$v\ge4sq,\qquad
 2^v\ge[10^6(b+1)]^{128},\qquad
 2^v\ge(512s^2q)^4,\qquad
 2^{5v}\ge(64sv)^8.\tag{16}$$

For example, the first comparison gives
`eta log(Y)/(4s) >= log(2) > log(16/9)`. The last comparison
uses `log(Y) <= v`. All four comparisons hold for sufficiently
large `b`. They are deliberately conservative.

Also `4R^2 < D`, and (14) gives

$$\frac R{h(\Delta)}>\frac1{8(b+2)}.$$

Consequently an explicit density lower bound from (8) is

$$\delta\ge\delta_0=\frac1{\mathcal D},\qquad
 \mathcal D=256(b+2)(s+1)s^2v\,(256q)^s s!.\tag{17}$$

Apply the existing integer prime-sample construction with
`epsilon = delta_0/4`. In closed form, define

$$k_\epsilon=\lceil\log_2(4\mathcal D)\rceil,
 \qquad X=[40000(b+k_\epsilon+2)\mathcal D]^{32},
 \qquad Z=\max(X,Y).\tag{18}$$

The sample has nonprincipal character bias below `delta_0/4`.
The positive-correlation argument in the
[generator and sample note](polylog-prime-form-generators.md)
therefore ensures that every class, after multiplication by one
sampled prime class, belongs to `S/S`.

Each element of `S` has a prime-form word of length exactly
`s+1`. Absorbing the sampled prime into the difference of two
such words gives a word of length at most

$$L=2s+3,\qquad K=L+1=2s+4.\tag{19}$$

After canceling opposite orientations and removing integer
content, this word supplies a primitive represented norm

$$n\le XR^2,\qquad \Omega(n)\le L,$$

supported on split primes at most `Z` and coprime to `Delta`.
The direct signed-word conversion and its diagonal-dominance
consequence are proved in
[the short-norm note](smooth-norm-short-words.md) and
[the direct-discriminant note](direct-discriminant-short-norms.md).
For a generator-power target `g_i^K`, finding any such norm
with `Omega(n) <= K-1` gives the required strict relation.
The search may accept any `n <= Z^L`; it need not retain the
smaller witness bound `XR^2`.

### Asymptotic cutoff and multiplicity

Put `ell_D = log D` and
`H_D = sqrt(ell_D log(ell_D))`. Since
`b = ell_D/log(2)+O(1)` and
`t = log(ell_D)/log(2)+O(1)`, the integer choice of `v` gives

$$\log Y=(\sqrt8+o(1))H_D.$$

More generally, if `log Y = (c+o(1))H_D`, then

$$u=(1/(2c)+o(1))\sqrt{\ell_D/\log\ell_D}.$$

Equations (15), (17), and (18) give

$$\log(1/\delta_0)=(1/(4c)+o(1))H_D,
 \qquad \log X=(8/c+o(1))H_D.$$

Balancing `c` and `8/c` gives `c = sqrt(8)`, proving (1).
In particular,

$$L\log Z=(1+o(1))\log D,\qquad
 \log(XR^2)=(1+o(1))\log D.\tag{20}$$

The coefficients describe this proof's parameter sizes. They
are not a global optimum over possible arithmetic algorithms.
Nor does the result make enumerating all `Z`-smooth integers
up to `Z^L` subexponential.

## 5. Reproducible checks and remaining obligation

The [parameter module](../tools/smooth_band_density.py) implements
(16)--(19) without floating-point logarithms. It marks small
parameters as ineligible instead of applying an asymptotic
inequality prematurely. The
[checker](../tools/check_smooth_band_density.py) validates the
kernel identities, the rational shell-error budget, the central
coefficient bound, the distinct-choice weighted counting step,
and the integer parameter calculations. Run:

```sh
python3 tools/check_smooth_band_density.py
```

The [evidence record](../evidence/smooth-band-density.json)
contains compact parameter summaries. The tests cover exact
finite arithmetic; they do not independently establish the
external zero-free theorem.

The resulting task is to select, for each needed form class,
a factored norm with at most `L` prime factors from this cutoff.
Once that integer is supplied,
[the orientation decoder](factored-norm-orientation-search.md)
recovers coordinates in
`2^O(sqrt(log D/log log D))` times a polynomial in the written
bit lengths. This is smaller than the prime cutoff's scale.
Selecting the integer within comparable deterministic time is
still open in this route, so the proved speedrun score does
not change.

This note used substantial OpenAI assistant support. Its external
analytic dependency and source credit are the same as in the
prime-form generator note.
