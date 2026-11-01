# A uniform one-half leading coefficient for supplied-norm processing

The [height-only unsigned interface](arbitrary-factored-norm-representation.md)
can approach its limiting processing coefficient with one uniform parameter
choice. The factor base is quasipolynomial in the discriminant bit length.
Under the same pinned analytic input, every required target has a smooth
norm of height `2^((1+o(1))b)`. Processing one supplied unsigned norm for
each target then costs

$$2^{(1/2+o(1))b/\log_2 b}. \tag{1}$$

This is a processing theorem. The deterministic selector is still missing,
so (1) is not a complete factoring bound.

## 1. Integer parameters without writing the density denominator

Let `b=bits(abs(Delta))`, put

$$q=\lceil\log_2 b\rceil,\qquad w=q^2,\qquad v=32q^2,
\qquad Y=2^v,$$

and use the notation

$$r=\lceil b/2\rceil-3,\qquad k=\lfloor r/v\rfloor,$$

$$a=\left\lceil\log_2(10000(b+w+2))\right\rceil,
\qquad h=w-a.$$

Require the same finite conditions as the
[powered-prime theorem](powered-prime-smooth-norms.md):

$$k\ge1,\quad h\ge1,\quad
2^w\ge10^6(b+1),\quad 2^v\ge64vk. \tag{2}$$

All hold for sufficiently large `b`. The last two can be checked by
integer bit lengths; neither `Y` nor a prime list has to be materialized
to validate the parameters.

The density denominator in that theorem is

$$\mathcal D=(b+2)2^{\lceil b/2\rceil-kv}(32vk)^k.$$

For parameter calculation it suffices to use the integer upper bound

$$d_+=2+\lceil\log_2(b+2)\rceil+\lceil b/2\rceil-kv
       +k\lceil\log_2(32vk)\rceil. \tag{3}$$

Every factor in `4 D` is at most its corresponding power of two, so
`4 D<=2^d_+`. Choose

$$t=\lceil d_+/h\rceil,\qquad L=2k+t,
\qquad B=vL,\qquad H=2^B. \tag{4}$$

Then `2^(-ht)<=1/(4 D)`. The existing powered-prime proof therefore
supplies a primitive represented norm in every class, supported on odd
split primes at most `Y`, with norm at most `H`. Increasing the walk
length to this integer upper bound only enlarges the accepted budget.

The prime-form base generates the full class group under the same
analytic input. No unknown prime factors or class orders enter (2)--(4).

## 2. Height and base size

Since `q=log_2 b+O(1)`, direct substitution gives

$$k=\frac{b}{64q^2}+O(1),\qquad
h=q^2-q+O(1).$$

Also `32vk=16b+O(q^2)`, while
`ceil(b/2)-kv=O(q^2)`. Consequently (3) gives

$$d_+=\frac{b}{64q}+O\left(\frac b{q^2}+q^2\right),$$

$$t=\frac{b}{64q^3}\left(1+O(1/q)\right)+O(1).$$

It follows that

$$B=v(2k+t)
 =b+\frac{b}{2q}+O\left(\frac b{q^2}+q^2\right)
 =(1+o(1))b. \tag{5}$$

The number of base primes is at most

$$Y=2^{32q^2}=2^{O((\log b)^2)}. \tag{6}$$

Any fixed polynomial cost in the number of base primes and in `b`
is therefore `2^O((log b)^2)=2^o(b/log b)`. This includes deterministic
base construction, modular roots, trial division of a supplied smooth
integer, and the relation engine's polynomial overhead.

## 3. Unsigned decoding and the complete batch

For each prime-form generator `g_p`, choose the least `K_p` with
`p^K_p>H`. A selector may return any odd `Y`-smooth integer `n<=H`,
coprime to `N`, represented by `g_p^K_p`. It need not return a primitive
representation, coordinates, or a factorization.

The [arbitrary-norm decoder](arbitrary-factored-norm-representation.md)
uses at most

$$\sqrt{(B+1)\tau(n)}\,\mathrm{poly}(b+B)
 =2^{(1/2+o(1))B/\log_2 B}\,\mathrm{poly}(b)$$

time. Equation (5) turns this into (1). Its primitive divisor yields a
logarithm-weighted dominant relation. The number of target rows and every
fixed polynomial relation-processing cost contribute only (6), so the
same leading coefficient holds for the full batch of one supplied norm
per target. The two decoder passes in the current wrapper add costs;
they do not multiply their exponential factors.

If a selector spends `T_sel` time and writes `s` candidates, including
failed attempts, the bound is instead

$$T_{\mathrm{sel}}+
s\,2^{(1/2+o(1))b/\log_2 b}. \tag{7}$$

Thus a selector with total cost `2^O(b/log b)` would still complete a
score-one procedure. Its leading running-time coefficient must include
the candidate count in (7). Equation (1) alone does not bound selection.

A coordinate-supplying selector avoids orientation search. Here its
processing overhead is quasipolynomial, which preserves any fixed
speedrun score. Prime and prime-power preprocessing, the full-generation
requirement, and the scoped meaning of failed partial extraction remain
as in the arithmetic framework.

### The powered prime support also fits this scale

There are at most `2Y` oriented prime forms in one sample. Enumerating
every ordered word of length `t` therefore costs at most

$$(2Y)^t\,\mathrm{poly}(b,t)
 =2^{(1/2+o(1))b/\log_2 b}, \tag{8}$$

because `(v+1)t=b/(2q)+O(b/q^2+q^2)` by the parameters above.
Repeated group elements can be left in the enumeration; no efficient
duplicate removal is needed for this bound.

Thus support enumeration itself is affordable at the score-one scale
with this base. It still does not construct a smooth norm. The positive
relation-kernel theorem guarantees that one enumerated word `u` admits
`g u=a b^(-1)` for two elements of the much larger recognizable smooth
set. Finding those two elements, or an equivalent represented norm, is
still unresolved. Equation (8) must not be substituted for that search.

## 4. Exact checks

[`quasipolynomial_norm_parameters.py`](../tools/quasipolynomial_norm_parameters.py)
computes (2)--(4) without constructing the density denominator, the prime
list, or `H`. It also gives an integer divisor-count upper bound for the
mixed-radix table size. The
[checker](../tools/check_quasipolynomial_norm_parameters.py) compares the
new upper walk bound with independently materialized exact density
denominators on moderate controls, and checks larger symbolic bit-length
parameters with exact integer arithmetic.

The [evidence](../evidence/quasipolynomial-norm-parameters.json) records
finite eligibility and the actual height and table-exponent bounds.
Small controls are not claimed to have reached the asymptotic constants.

This note used substantial OpenAI assistant support. The existence input
is the same pinned all-number-field finite-order Hecke zero-free theorem
as in the powered-prime note; the processing bound is deterministic.
