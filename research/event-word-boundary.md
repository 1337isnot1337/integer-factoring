# Row-change events are a one-sided smooth-word boundary

## Exact representation and remaining analytic question

Below the reduced norm budget `R<=sqrt(D)/4`, the truncated row-change
events admit a static arithmetic description. They are exactly the
one-sided boundary of the signed smooth-word body in the target prime
direction. Every word omitting the target prime supplies one event.
Consequently:

- event density on a subgroup containing the target prime is **exactly**
  the density of smooth classes whose leading norm omits that prime;
- each same-row interval is complete, and no integer event row repeats;
- event character sums are explicit annular smooth-ideal sums, with one
  prescribed orientation at the target prime.

These are ordinary consequences of the
[principal-word norm gap](small-norm-event-separation.md), not a uniform
estimate for the nonzero target phases. The covered-subgroup density
theorem remains conditional on its existing principal split-prime input.
Actual target indices, long-range connected correlations, and a stronger
complete factoring bound remain open. The complete score stays **1/2**.
The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19356388)
changes the representation of events, not a cutoff or orbit search.

## 1. A smooth-word body with an injective class map

Use only unramified split primes in the factor base, with a fixed
orientation for each, and let

```math
\phi:\mathbb Z^m\longrightarrow G=\operatorname{Pic}(\mathcal O_\Delta),
\quad \mathcal N(v)=\prod_q q^{|v_q|},\quad
\mathcal A_R=\{v:\mathcal N(v)\le R\}.
```

The injective-ball statement is already recorded in
[orda-aska's graph/lattice Lemma 1](leading-power-window-graphs.md#1-the-graph-and-the-lattice).
The map `phi` is injective on `A_R`. Indeed, two distinct words with the
same class would give a nonzero principal word `v-w` of norm at most
`N(v)N(w)<=R^2<=D/16`, contradicting the lower bound `D/4`.
Each corresponding primitive ideal is already reduced after centering,
because its norm is at most `sqrt(D)/4`. Thus `phi(A_R)` is precisely
the recognizable smooth, discriminant-coprime set with leading norm
at most `R`, not a collection of noncanonical ideal witnesses.

Fix a target prime `p` in this base. For a hit with canonical word `v`,
the preceding exponent is a same-row continuation exactly when its word
is `v-e_p`. By injectivity and the reduced-norm property this happens
if and only if `v-e_p` belongs to `A_R`. Therefore the event set is

```math
\mathcal E_{p,R}=\phi(\mathcal B_{p,R}),\qquad
\mathcal B_{p,R}=\{v\in\mathcal A_R:v-e_p\notin\mathcal A_R\}.
                                                               \tag{1}
```

Equivalently, the event test at a hit is simply

```math
\boxed{v_p\le0\quad\hbox{and}\quad p\mathcal N(v)>R.}     \tag{2}
```

For `v_p>0`, subtracting `e_p` divides the norm by `p`; otherwise it
multiplies the norm by `p`. This proves (2), including the strict boundary
inequality. The observable is not inversion symmetric.

## 2. Fibers, density, and unique integer event rows

For every target-prime-free word `u` with `u_p=0` and `N(u)<=R`, define

```math
a(u)=\max\{a\ge0:p^a\mathcal N(u)\le R\}.
```

Its entire fiber in `A_R` is
`u+z e_p`, `-a(u)<=z<=a(u)`. Equation (1) selects exactly its lower
endpoint `u-a(u)e_p`. These fibers have disjoint class images by
injectivity. Along the target cycle they are complete same-row intervals;
when one crosses phase zero, lift the exponents to integers before
comparing row labels. The periodic event definition retains this drift.

Suppose events at two integer exponents had equal integer rows. Their
words would differ by a positive multiple of `e_p`. They would lie in
the same fiber, whose lower endpoint is unique, a contradiction. This
uniqueness does not imply rational independence of different rows.

Let `U` be any subgroup containing `g_p=phi(e_p)`. A fiber belongs to
`U` precisely when its prime-free word does. Thus the exact event density
on `U` is

```math
\boxed{\frac{|\mathcal E_{p,R}\cap U|}{|U|}
 =\frac{|\{u:u_p=0,\ \mathcal N(u)\le R,\ \phi(u)\in U\}|}{|U|}.}
                                                               \tag{3}
```

The same count holds in every coset of `U`. Also the hit count in that
coset is the sum of `2a(u)+1` over these prime-free words; the previous
polynomial density-loss bound is an immediate consequence. No group
order or discrete logarithm is needed to define the event or prove (3).

For a direct lower bound in the covered subgroup range, repeat PR87's
even-color count after removing the single underlying prime `p` from
the split-prime band. This removes at most `O(log Y)` of principal mass,
whereas total mass is comparable to `Y`; all coloring, normalization,
and one-sided Fourier bounds remain uniform in the choice of `p`.
The constructed words omit `p`. Consequently, if

```math
[G:U]\le D^\kappa,\quad\kappa<1/4,\quad
E=\log(\sqrt D/R)=o(b),\quad b=\operatorname{bits}(D),
```

the norm-budget proof gives **directly for the event density**

```math
\frac{|\mathcal E_{p,R}\cap U|}{|U|}
 \ge\exp\left(-E-O\left(\log Y+\frac b{\log Y}\log b\right)\right).
                                                               \tag{4}
```

At the fixed-alpha cutoff and radius of the separation theorem this is
the same `exp(-O_alpha(b/(log b)^alpha))` scale. The index restriction
remains; (3) does not bound an actual cyclic index.

### A static prime-free flag and bounded phase transport

Let `f_0` mark the smooth, discriminant-coprime reduced classes of norm
at most `R` whose leading coefficient is not divisible by `p`. This is
the image of the prime-free words above. It is inversion symmetric and
recognized from the current reduced form alone. Each such hit has decoded
coordinate `v_p=0`, so its integer relation row has `p` coordinate exactly
the exponent `k`. No two of these hits repeat an integer row.

On the actual target cycle, the fiber proof gives a bijection from a
prime-free hit at phase `j` to its event at `j-a(u)`, with

```math
0\le a(u)\le A_p:=\lfloor\log_pR\rfloor.
```

Injectivity of the pure `p`-word fiber implies `2A_p<h`; otherwise
two words in `[-A_p,A_p]e_p` would have the same class. In every cyclic
window `I` of any integer length `W`, this bounded one-sided transport
therefore gives

```math
\boxed{\left|\sum_{k\in I}e_R(k)-\sum_{k\in I}f_0(k)\right|\le A_p.}
                                                               \tag{4a}
```

Reduce `W` modulo `h`, since full cycles have equal counts by (3).
A point transported into the residual window from outside must start
within `A_p` positions beyond its right end. A point transported out
must start within `A_p` positions of its left end. Each of the two
counts is at most `A_p`, and the net difference is their difference,
proving (4a). Reordering transported points is allowed.

If every length-`W` prime-free window hits, every event window of length
`W+A_p` hits: a prime-free point in its subwindow starting `A_p` steps
after the left end is transported inside it. The reverse direction uses
the inverse rightward transport. In particular, their longest zero gaps
differ by at most `A_p`. Their common window mean is `mu`, and the
triangle inequality in `L^(2q)` gives

```math
\| (\mu-X_{e_R})_+\|_{2q}
 \le\|(\mu-X_{f_0})_+\|_{2q}+A_p.                        \tag{4b}
```

Thus the remaining analysis may use the **static, inversion-symmetric
prime-free flag** directly. The finite-order window criterion applies
to it with the proved covered-subgroup density (4); the raw selector
already accepts its hits and needs no event-transport oracle. If an
event window is desired, the extra `A_p=O(b)` steps do not affect a fixed
saving exponent. Under `4p^T R^2<D`, this prime-free flag is also
`T`-separated: two hits `t<=T` apart would give a nonzero principal word
whose `p` coordinate is `t` and whose norm is at most `p^tR^2`.
This supplies the same exact short-range cumulants without history.
Long-range moment estimates for the prime-free flag remain open.

## 3. Exact annular character sums and their Euler product

For a class character `chi`, put `z=chi(g_p)` and define

```math
F_p(x;\chi)=\sum_{u_p=0,\ \mathcal N(u)\le x}\chi(\phi(u)),
\qquad F_p(x;\chi)=0\quad(x<1).
```

Partition the prime-free words according to `a(u)`. Equation (1) gives
the exact unnormalized character sum

```math
\boxed{\sum_{g\in\mathcal E_{p,R}}\chi(g)
 =\sum_{a\ge0}z^{-a}
       \left(F_p(R/p^a;\chi)-F_p(R/p^{a+1};\chi)\right).} \tag{5}
```

The sum ends at `floor(log_p R)`; real cutoffs mean the corresponding
integer norm inequality. If `chi(g_p)=1`, it telescopes to `F_p(R;chi)`.
This explains (3) also by subgroup character orthogonality. The same
identity holds after restricting the prime-free words and events to
any `p`-stable subgroup or coset, for characters of that subgroup.

There is an explicit multiplicative series behind these sums. For
`Re(s)>0`, finite prime support gives absolute convergence of

```math
\sum_{u_p=0}\frac{\chi(\phi(u))}{\mathcal N(u)^s}
 =\prod_{q\ne p}
 \frac{1-q^{-2s}}
 {(1-\chi(g_q)q^{-s})(1-\chi(g_q)^{-1}q^{-s})}.           \tag{6}
```

Each local factor is just the geometric sum over integer exponents
`e`, weighted by `chi(g_q)^e q^(-|e|s)`. The finite product ranges over
the unramified split factor base, not every prime of the number field.
Allowing only the negative orientation at `p` adds the local factor
`(1-z^(-1)p^(-s))^(-1)`. If `C_p(x;chi)` is the summatory function of
that series, (5) is also exactly
`C_p(R;chi)-C_p(R/p;chi)`. Thus the actual event observable is an
annulus of a multiplicative smooth-ideal count with a specified
target-prime orientation.

Equations (5)--(6) are an analytic interface, not a claimed contour
estimate. A zero-free strip for the full Hecke L-function does not,
without another argument, bound these finite smooth products and their
annular sums uniformly in all the required target phases.

### A contributed finite recurrence and exact Mellin multiplier

The coordination contributor supplied the
[Mellin identity](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19356441)
and its [finite group-algebra proof and control source](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19356480).
Retain their distinction between formal words and actual reduced classes.
Extend `M_R=sum_(v in A_R)[phi(v)]` and the boundary-word measure `B_R`
to every real `R>=1`, counting image collisions with multiplicity. In the
free word algebra, and hence after every group projection,

```math
\boxed{B_R+B_{R/p}=M_R-g_p M_{R/p}.}                     \tag{6a}
```

For one prime-free fiber of half-length `a>=1`, the shifted smaller
interval removes exponents `-a+2,...,a`, leaving `-a,-a+1`, exactly the
two boundary endpoints. For `a=0` only the first endpoint remains. This
proves (6a) before any analytic limiting argument.

Let `F_chi(s)` be the full signed-word Euler product, including `p`,
and write `M_chi(R),B_chi(R)` for character evaluation of these formal
measures. For `Re(s)>0`, absolute convergence permits integration of
each word's step function. Equation (6a) then gives

```math
\int_1^\infty B_\chi(R)R^{-s-1}\,dR
 =\frac{1-zp^{-s}}{1+p^{-s}}\frac{F_\chi(s)}s
 =\frac{1-p^{-s}}{s(1-z^{-1}p^{-s})}
                         F_{\chi,\,\mathrm{without}\ p}(s).
                                                               \tag{6b}
```

These integrals range beyond the injective height. They are transforms
of **formal word measures**, and are not asserted to be transforms of
the actual reduced-class event flag at every height. The event theorem
identifies the measures only inside its stated reduced range.

At real `s=sigma>0`, the raw-to-boundary multiplier tends, as `z` tends
to one, to `(1-p^(-sigma))/(1+p^(-sigma))>0`. Thus this local Euler
factor does not by itself cancel a near-zero target-frequency window
denominator. An asymptotic cluster ratio would need a separate uniform
subgroup saddle-point or regular-variation theorem; none is assumed here.

## 4. The long-range correlation count is now explicit

For `U=<g_p>` of order `h`, and integer offsets `t_1=0,t_2,...,t_r`,
the joint event moment is exactly `1/h` times the number of tuples
`v_1,...,v_r in B_(p,R)` satisfying

```math
\phi(v_1)\in U,\qquad
v_i-v_1-t_i e_p\in\ker\phi\quad(2\le i\le r).           \tag{7}
```

Each allowed start corresponds to exactly one such tuple. Repeated
offsets force the corresponding words to agree; they are not discarded.
The connected correlation is obtained by the usual partition formula.
Unlike the whole smooth-word body, each fiber in (7) has only one
boundary representative, so same-row multiplicities have been removed
exactly, not probabilistically.

For example, if `c_rho` are the normalized Fourier coefficients of the
event flag restricted to `U`, its window count has exact variance

```math
\operatorname{Var}\sum_{j<W}e(x+j)
 =\sum_{\rho\ne1}|c_\rho|^2
                  \left|\sum_{j<W}\rho^j\right|^2.       \tag{8}
```

The annular formula specifies the underlying arithmetic sums. A uniform
phase-weighted energy estimate in (8), and higher connected estimates
through the required orders, remain to be proved. Pairwise control alone
is insufficient for the [finite-order window criterion](connected-event-window-criterion.md).

Enumerating the word body, all its boundary words, or all group characters
is not an affordable algorithm supplied by this representation. Its size
can have logarithm of order `b`; no such enumeration or lattice oracle is
charged as polynomial work. The raw selector's existing cost remains
unchanged. An analytic existence estimate must still yield an affordable
uniform window, or a genuinely efficient witness extractor must be given.

## 5. Exact controls

```sh
make verify-event-word-boundary
```

The [checker](../tools/check_event_word_boundary.py) reuses the fixed radii
of the preceding separation control, without choosing another cutoff.
It checks every saved truncated event against (2), partitions all 235/213
hits into 45/55 disjoint fibers, and verifies unique integer event rows.
The existing windows 3260 and 49712 check (4a) at every phase, without
selecting another window or changing a cutoff.
It checks (5) exactly for character orders 1, 2, 5 and 10 by residue
histograms, using the actual target strides. No floating roots of unity
or new large cycle are used. The
[receipt](../evidence/event-word-boundary.json) pins those event flags.

It also integrates the contributed 600 exact recurrence/projection/subgroup
controls at the original 200 integer heights on 67 words in the fixed
`3,5,7` base. The projected group deliberately has collisions. The
original control source SHA256 and both discussion comment identifiers
are retained in the checker and receipt; this panel is not expanded.

A fixed `Delta=-23`, `p=2`, `R=4` counterexample demonstrates the height
hypothesis: the canonical words are `0,1,-1` on the three-cycle, giving
one event, while an unqualified word-boundary test with that excessive
radius gives none. The word map has already ceased to be injective.

Credit: orda-aska's graph/lattice Lemma 1, PR87's even-product density
and primitive-word injection, PR82's
row-change process, the preceding principal-word separation argument,
and the saved selector/Grok data. The geometric-series and Fourier facts
are classical. The coordinated finite recurrence, Mellin multiplier and
control code are credited to the linked Discussion #51 contribution.
Developed with substantial OpenAI Codex assistance; no
new Lean, independent mathematical review or score improvement claimed.

The [frozen default replay](../evidence/event-word-boundary-full-verification.json)
passes at exact source `ac6cd8c1dea8b6989eeaad14a674fc4f77bc3553` in
1085.08 seconds, with clean tracked state before and after. The
[complete compressed log](../evidence/event-word-boundary-full-verification.stdout.txt.xz)
is archived. The later bounded-transport follow-up `df6c461` has a separate
focused/docs replay; the full suite was not repeated for that follow-up.
All full-suite temporary paths were private to its namespace.

The [separately reported entry-tail sign error](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19357329)
concerns two older moment helpers. This checker imports only their unchanged
cyclic-count routine and uses its own event/fiber identities. The frozen
code replay is not an endorsement of those older tail labels.
