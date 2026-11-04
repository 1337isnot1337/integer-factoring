# Parity-refined theta programs repair the one-mod-four phase gap

The [partial theta compiler](theta-coordinate-programs.md) proves all-index
odd-coordinate programs at phases 3,7,11,19,23,31 and returns exact
bounded-module exclusions at 5,13,17. In each excluded phase the missing
native coordinate is zero. The obstruction is parity: the unsigned coset
`n = +-a (mod q)` forgets whether the representative is `a` or `q-a`
modulo `2q`. Splitting every coset into its two parity classes modulo `2q`
repairs coordinate zero at all three excluded phases. The refinement is a
different, larger representation, not a submodule of the old one.

This note records the construction, a conservative ordinary modular/Sturm
certificate, the finite results, and the exact remaining limitation. The
degree-four evaluator is still `O(poly(q) N^(3/2))`, so the complete
factoring exponent remains 1/2. The work was developed with local
DeepSeek assistance. The ordinary proof and certificates have not received
independent mathematical review; no new Lean result is claimed.

## 1. Parity variables

Let `q` be an odd prime and `h=(q-1)/2`. For `1<=a<=h` and `s in {0,1}`
define

```math
\theta_{a,0}(t)=\sum_{\substack{n>0\\ n\equiv \pm a\pmod{2q}}}t^{n^2},
\qquad
\theta_{a,1}(t)=\sum_{\substack{n>0\\ n\equiv \pm(q-a)\pmod{2q}}}t^{n^2}.
```

The residue classes are disjoint and their union is the old class
`n = +-a (mod q)`, so

```math
\theta_a=\theta_{a,0}+\theta_{a,1}.
```

The trace series is unchanged:

```math
V_q(t)=\sum_{n>0}(t^{n^2}+t^{qn^2}).
```

The candidate module has `2h=q-1` variables and all monomials of degree at
most four. Its scalar columns are the trace column and
`1+binomial(2h+4,4)` monomial columns. The coefficient algebra is still the
intrinsic Tate algebra `A_q` of degree `h`; no field component is selected.
The old module embeds by `theta_a -> theta_{a,0}+theta_{a,1}`, but the new
module also separates the two parity classes.

## 2. Conservative modular lift

Homogenize each monomial of degree `d` by `H^(4-d)`, where
`H=sum_{n in Z} t^(n^2)`. As in the original construction, every
homogeneous term is a finite sum of rank-four theta series for
`Z^4` with the standard sum-of-four-squares form. The coordinate
congruences are now modulo `2q` instead of `q`.

Replacing `q` by `2q` in the rank-four theta lift gives the conservative
group

```math
\mathcal H''_q=\Gamma_0(64q^2)\cap\Gamma_1(8q).
```

The factors are the same as in the original proof: coordinate modulus `2q`
contributes `Gamma_0(4(2q)^2)=Gamma_0(16q^2)`; odd projection increases the
lower-left divisibility by at most four, giving `Gamma_0(64q^2)`; the
finite theta multiplier is an eighth root of unity, killed by `Gamma_1(8q)`.
This last level is deliberately conservative. A sharper multiplier
calculation may put the forms on `Gamma_0(64q^2) intersect Gamma_1(2q)`,
with bound `16q(q^2-1)`; finite checks also pass at that smaller bound for
`q=5,13,17`, but that sharper level is not used here.

The index is

```math
[\mathrm{SL}_2(\mathbb Z):\mathcal H''_q]
=384q(q^2-1),
```

so the inclusive weight-two Sturm bound is

```math
B''(q)=64q(q^2-1).
```

Checking all odd coefficients through `B''(q)` proves an all-index identity
for the projected forms, under the stated lift. The checker validates
through `B''(q)+128` as a margin. This is an ordinary modular/Sturm
certificate, not a Lean theorem.

## 3. Finite results

The compiler performs exact binary elimination over `F_2` on the parity
module. It fits at `min(B'',max(4096,2*columns))` and validates through
`B''+128`. The original compiler's exclusions are replayed in the same
checker: at `q=5,13,17` its specified module still misses exactly native
coordinate zero.

| q | parity vars | B'' | binary rank | columns | program | max quadratic vars | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 2 | 1,536 | 6 | 16 | 1 | 0 | calibration |
| 5 | 4 | 7,680 | 30 | 71 | 15 | 3 | **repaired** |
| 7 | 6 | 21,504 | 84 | 211 | 51 | 3 | regression |
| 11 | 10 | 84,480 | 318 | 1,002 | 205 | 4 | regression |
| 13 | 12 | 139,776 | 546 | 1,821 | 348 | 4 | **repaired** |
| 17 | 16 | 313,344 | 1,224 | 4,846 | 814 | 4 | **repaired** |
| 19 | 18 | 437,760 | 1,654 | 7,316 | 1,085 | 4 | regression |

Every listed parity program has no missing coordinates and no validation
failure. The three new repairs are `q=5,13,17`. The phases 7,11,19 are
regression controls; the original compiler already succeeded there.
Phases 23 and 31 were not recompiled in the parity module.

The checker also evaluates 63 original-index controls through 511 bits,
using the parity programs and literal odd-divisor planes. The total charged
enumeration is 128,926 steps. It refuses the original partial program at
phase 13, and it records an underdetermined phase-17 fit whose validation
fails, so a training fit alone is not promoted to an all-index identity.

## 4. Cost and scope

A degree-four monomial has at most four positive quadratic variables after
Frobenius substitution. Enumerating all but the last variable and using a
square test for the last gives

```math
O(\mathrm{poly}(q)N^{3/2}\mathrm{poly}(\log N))
```

time and polynomial-in-`q` storage. This is worse than the `O(N)` Tate
composition and remains exponential in the input bit length. The parity
refinement therefore does not improve the complete bound or the score. Its
value is structural: it removes the one-mod-four phase obstruction in the
degree-four theta representation, supplies a uniform candidate module that
succeeds at every phase tested so far, and isolates parity as the missing
coordinate.

The following obligations remain:

- prove uniform spanning of the parity module for every odd prime `q`, or
  find a phase where it fails;
- reduce the maximum number of quadratic variables below four, for example
  by a `q`-nome rescaling or trace elimination, to reach a cheaper query;
- connect the parity variables to the tensor/collective route or to a
  composite-index coefficient engine;
- decide whether the sharper `Gamma_1(2q)` level and bound `B'` can be
  proved, which would reduce certification cost but not query cost.

The first two are the relevant next steps for a stronger complete bound.
More finite fits at larger `q` alone do not resolve them.

## 5. Replay

Run `make verify-theta-parity`. The checker regenerates every program,
validates the conservative Sturm bounds, replays the original exclusions,
checks the original-index evaluator controls, and verifies the deterministic
[compressed program artifact](../evidence/theta-parity-programs.json.gz)
against the [receipt](../evidence/theta-parity-programs.json). The receipt
records ranks, columns, program lengths, validation limits, charged steps,
source hashes, and artifact hashes. Use
`python3 tools/check_theta_parity_program.py --write-programs` only to
regenerate the artifact deliberately after a source change.
