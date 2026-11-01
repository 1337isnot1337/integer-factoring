# Short smooth norms from a powered prime sample

The prime sample can be powered before it is absorbed into a norm.
This reduces the smoothness cutoff from
`exp(O(sqrt(log D log log D)))` to `exp(O(sqrt(log D)))`, while
keeping the number of prime factors `O(sqrt(log D))` and the norm
at most `D^(1+o(1))`. The sample's product support is used only in
an existence proof, so it is never enumerated by this reduction.

Under the same
[pinned finite-order Hecke zero-free input](polylog-prime-form-generators.md),
write `D = abs(Delta)` and `b = bits(D)`. For every sufficiently
large negative quadratic-order discriminant, every proper ideal
class primitively represents an integer `n` satisfying

$$\gcd(n,\Delta)=1,\qquad
 n\text{ is }Y\text{-smooth},\qquad
 \log_2Y=(1/\sqrt2+o(1))\sqrt b,$$

$$\Omega(n)\le L=(\sqrt2+o(1))\sqrt b,
 \qquad \log_2 n\le b+O(\sqrt b\log b).\tag{1}$$

All prime factors can be taken odd and split. Once a candidate
integer `n` is supplied, factoring it over the bounded base and
recovering its orientations and coordinates takes

$$2^{(1/\sqrt2+o(1))\sqrt b}\tag{2}$$

bit operations. Building the factor base has the same bound.
This is a per-candidate construction and decoding bound, not a
norm-selection or complete factoring bound.

## 1. A small-bias walk costs only its number of steps

Let `mu` be a distribution on prime forms with nonprincipal
character bias at most `epsilon`. Its `t`-fold convolution has
bias at most `epsilon^t`, because each Fourier coefficient is
raised to the `t`-th power. Every point in its support has a
prime-form word of length at most `t`.

Suppose a recognizable smooth set `S` has density at least
`delta_0`. Choose `epsilon^t <= delta_0/4`. The
[positive relation-kernel lemma](positive-relation-kernel.md)
then implies, for every class `g`, that

$$g u=a b^{-1}\quad\text{for some }a,b\in S,
 \quad u\in\mathrm{supp}(\mu^{*t}).\tag{3}$$

If the forms in `S` have squarefree leading coefficients at most
`R`, each with exactly `k` prime factors at most `Y`, and the
prime sample also uses primes at most `Y`, absorbing the walk
into (3) gives a primitive represented norm with

$$\Omega(n)\le2k+t,
 \qquad n\le Y^tR^2,
 \qquad n\le Y^{2k+t}.\tag{4}$$

The primitive ideal is obtained by canceling rational-integer
content, as in the
[single-prime absorption proof](direct-discriminant-short-norms.md).
Cancellation only decreases the norm and multiplicity.

No algorithm is asked to enumerate `supp(mu^{*t})`, whose size
could be as large as the `t`-th power of the factor-base size.
The remaining algorithmic task continues to be a direct norm
search in the prescribed target form. Equation (3) proves that
this search has an answer within the budget (4).

## 2. A simpler squarefree count suffices at this scale

Let

$$r=\lceil b/2\rceil-3,\qquad R=2^r,
 \qquad Y=2^v,
 \qquad k=\lfloor r/v\rfloor\ge1.\tag{5}$$

Then `4R^2 < D` and `Y^k <= R`. Assume

$$Y\ge[10^6(b+1)]^{32},\qquad Y\ge64vk.\tag{6}$$

The split-prime shell bound in
[the multiband note](smooth-band-density.md), Section 1, gives
at least `M >= Y/(64v)` split odd primes below `Y`, coprime
to `Delta`. In particular `M >= k`.

Every product of `k` distinct such primes is at most `Y^k <= R`.
Its `2^k` orientation choices give distinct reduced forms: the
centered middle coefficient has absolute value at most the leading
coefficient, while `4R^2 < D` makes the last coefficient strictly
larger. The unramified coprimality condition excludes the boundary
identification between opposite middle coefficients.

Thus the set `S` of all reduced forms with squarefree leading
coefficient at most `R`, coprime to `Delta`, and exactly `k`
prime factors at most `Y` has

$$|S|\ge2^k\binom Mk\ge\left(\frac{Y}{32vk}\right)^k.\tag{7}$$

The elementary inequality `binom(M,k) >= (M/k)^k` follows by
comparing each ratio `(M-j)/(k-j)` with `M/k`.
The order class-number bound
`h(Delta) < 2^ceil(b/2) (b+2)` now gives

$$\frac{|S|}{h(\Delta)}\ge\delta_0=\frac1{\mathcal D},$$

$$\mathcal D=(b+2)\,2^{\lceil b/2\rceil-kv}(32vk)^k.\tag{8}$$

Only the bit length of this explicit denominator will enter the
walk length. For `v` near `sqrt(b)`,

$$\log_2\mathcal D=O(\sqrt b\log b).\tag{9}$$

This simpler count loses more density than the multiband theorem.
That loss is harmless here: powering the prime sample absorbs it
into only `O(log b)` additional factors.

## 3. An exact integer choice for the prime sample and walk

Take `v = 32w`, with `w` a positive integer. Define

$$a=\left\lceil\log_2\bigl(10000(b+w+2)\bigr)\right\rceil,
 \qquad h=w-a,$$

and require `h >= 1`. Set `epsilon = 2^(-h)` and use the
existing [explicit prime-sample construction](explicit-prime-form-cutoffs.md).
Its cutoff is

$$X=[10000(b+h+2)2^h]^{32}\le2^{32w}=Y.\tag{10}$$

Indeed, `h <= w` and `10000(b+w+2) <= 2^a`.
The dyadic weight construction has all nonprincipal biases
strictly below `epsilon`.

Let

$$d=\lceil\log_2(4\mathcal D)\rceil,
 \qquad t=\lceil d/h\rceil,
 \qquad L=2k+t,
 \qquad K=L+1.\tag{11}$$

Then `epsilon^t <= delta_0/4`. Equations (3)--(4) prove the
existence of a primitive norm in every class with multiplicity
at most `L`, all of whose prime factors lie below the single
cutoff `Y`.

For the relation program, include the prime forms up to `Y` in
the generator list and use target classes `g_i^K`. A successful
norm search in each target gives a word of length at most `L`,
so `K e_i - v_i` is a strictly diagonally dominant relation row.
Full generation follows from the earlier polylogarithmic prime
bound, since (6) places that bound below `Y`. The
[polynomial ambiguous-class extractor](polynomial-ambiguous-extraction.md)
then applies once an annihilator has been constructed from the
relation matrix.

These are conditional reduction steps. They do not estimate the
unknown norm search by the size of its eventual output.

## 4. Balance the root construction and orientation search

Choose the integer `w` as the least positive integer with

$$2048w^2\ge b,\qquad v=32w.\tag{12}$$

Then

$$v=(1/\sqrt2+o(1))\sqrt b,
 \qquad k=(1/\sqrt2+o(1))\sqrt b,$$

$$h=(1/(32\sqrt2)+o(1))\sqrt b,
 \qquad d=O(\sqrt b\log b),
 \qquad t=O(\log b).$$

The inequalities (6) and `h>=1` hold for every sufficiently large
`b`; the implementation checks them exactly. It follows that

$$L=2k+t=(\sqrt2+o(1))\sqrt b,$$

$$Lv=2kv+tv\le2r+tv
       \le b+O(\sqrt b\log b).$$

This proves (1). Since `b log(2) = log D + O(1)`, the last
bound is also `n <= D^(1+o(1))`.

The [root constructor](prime-form-root-construction.md) builds
the factor base in `Y poly(b,log Y)` time. Trial division over
this base factors a supplied candidate `n` within the same bound.
The [orientation decoder](factored-norm-orientation-search.md)
then decides representation and recovers coordinates in
`2^ceil(omega(n)/2) poly(b,log n)` time. Since
`omega(n) <= Omega(n) <= L`, the sum of these two costs proves
(2). In natural logarithms its leading exponent is

$$\sqrt{\frac{\log2}{2}}\,\sqrt{\log D}.
$$

### The cost of processing every relation row

There are at most `O(Y)` oriented generators. Decoding one
already selected norm per generator therefore takes at most

$$2^{v+\lceil L/2\rceil}\mathrm{poly}(b)
 =2^{(\sqrt2+o(1))\sqrt b}.\tag{13}$$

This counts all decoding calls, rather than treating their maximum
as the cost of the batch. The norm-selection calls and determinant
calculation must still be charged separately. Any fixed polynomial
cost in the number of generators remains `exp(O(sqrt(b)))`,
but can change the leading constant in that scale.

### A broader tradeoff

The same construction works for integer widths satisfying
`v/log b -> infinity` and `v=o(b)`, after rounding to a multiple
of 32. In that regime,

$$k\sim\frac{b}{2v},\qquad
 t=O\left(1+\frac{b\log b}{v^2}\right),\qquad
 \frac{t}{k}\longrightarrow0.$$

Thus `L ~ b/v`, and `Lv <= b+o(b)`. The two per-candidate
exponents are `v` for the bounded factor base and `b/(2v)`
for the orientation table. Their balance explains (12).
This optimization concerns the present decoder and setup only.

## 5. Exact parameter replay and the remaining search problem

The [parameter calculator](../tools/powered_prime_norm_parameters.py)
computes (5)--(12) with integer arithmetic and reports every
failed finite eligibility check. It does not construct a walk,
enumerate its support, or select a norm.

Run:

```sh
python3 tools/check_powered_prime_norm_parameters.py
```

The [checker](../tools/check_powered_prime_norm_parameters.py)
checks the squarefree binomial inequality, density denominators,
cutoff domination, exact walk bias, and multiplicity budgets.
It also compares finite convolution powers and bias bounds on
cyclic groups. The
[record](../evidence/powered-prime-norm-parameters.json) includes
large input-bit parameter controls without writing enormous
norms or prime lists.

The new endpoint is a deterministic selector for the target form:
return one integer `n` with at most `L` prime factors from this
base that the orientation test accepts. No such selector with
`exp(O(sqrt(b)))` cost is supplied here. The proved speedrun
score consequently remains unchanged.

This note used substantial OpenAI assistant support. The pinned
analytic source and credit are the same as in the generator note.
