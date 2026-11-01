# Uniform-marginal transport closes both fixed-gauge endpoint escapes

For the full degree-at-most-eleven pencil over every characteristic-two
field, arbitrary physical gauges cannot improve the following two normalized
entry-minimum optima:

```math
\sup_{h_0}\{m_0(h_0)-B_0(h_0)\}=14/5,
\qquad
\sup_{h_\infty}\{m_\infty(h_\infty)-B_\infty(h_\infty)\}=-66/5. \tag{1}
```

Here `m_e` is the first raw local Smith order, equivalently the minimum
weighted entry valuation, and `B_e=(1/5) sum_(leg,letter) h_(e,leg)(letter)`.
The paired physical energies may have arbitrary shapes. These are exact
all-field statements about the complete pencil, not a finite gauge screen.

Combining (1) with determinant identities gives strictly positive first-order
supporting-line margins at both axes, when the gauges and the ratio between
thermal scale and small axis weight are fixed. This strengthens the old
[fixed-height determinant target](two-place-determinant-target.md).
It does **not** exclude all interior angles or height/scale blow-up regimes.
No factoring score above **1/2** is established.

## 1. The fixed pencil and an unaffected coefficient

Use the PR57 pencil, reconstructed and credited in the determinant note:

```math
A_q(t)=A(t)+q(t)(K+t^4L),\qquad \deg q\le11,\quad q(1)=0.
```

The two perturbation bits `K,L` of every tensor entry are retained in the
canonical input. A coefficient of degree `e` in an entry cannot change if
`e` is outside `0..11` when `K=1`, and outside `4..15` when `L=1`.
Choose a nonzero source coefficient outside these ranges. It stays nonzero
over every characteristic-two coefficient field and for every allowed `q`.
Consequently its zero entry order is at most `e`, and its infinity entry
order at most `-e`, regardless of other coefficient cancellations.

## 2. Two small transport certificates

For each place choose a probability distribution on such guaranteed
entry/degree pairs `(x,y,u,v;e)`, with **uniform marginal `1/5` on each
of the four physical legs**. Then, for arbitrary physical heights,

```math
m_0(h_0)\le\mathbb E[e]+B_0(h_0),\qquad
m_\infty(h_\infty)\le-\mathbb E[e]+B_\infty(h_\infty). \tag{2}
```

Indeed the minimum weighted entry is at most the weighted average of
these particular entry orders; the marginal identities evaluate the height
average exactly. No optimizer or fixed field size is involved in (2).

The [witnesses](../evidence/pencil-uniform-transport-witnesses.json) have
eleven and nine positive entries. The zero witness has expected degree
`14/5`; the infinity witness has expected degree `66/5`. The checker
verifies each unaffected coefficient, positivity, total mass one and all
twenty exact marginal equalities. It therefore proves the upper bounds (1).

The inherited natural gauges attain both bounds for **every** member of
the pencil. Their uniform height averages are both `26/5`. Every source
zero entry has order at least 8; every zero perturbation has order at least
8. Every source infinity entry has order at least -8; every infinity
perturbation has order at least -7. The unchanged `(0,0,0,0)` entry is
exactly `t^8`, attaining both minima. Thus
`8-26/5=14/5` and `-8-26/5=-66/5`, proving equality in (1).

## 3. The opposite determinant orders

The existing all-field classification gives opposite-cut polynomial
determinant degree `D>=78`. This statement is independent of physical
heights: they only shift determinant orders by their row/column height sum.

The ordinary row determinant is invariant under the entire pencil. Its
generator satisfies `G^2=0`, `GAG=0` and the exact polynomial identity

```math
GA+AG=t^8(1+t)^4(K+t^4L).
```

Over the rational function field, put `f=q/(t^8(1+t)^4)`. Then

```math
A_q=(I+fG)A(I+fG),\qquad (I+fG)^{-1}=I+fG.
```

Both factors have determinant one, so `det A_q=det A`. The checker freshly
computes the full source determinant by exact fraction-free elimination:

```math
det A=t^{314}(1+t^4+t^6+t^{10}+t^{12}+t^{14}+t^{16}).    \tag{3}
```

Thus its zero order is exactly `V=314` over every characteristic-two field.
The resulting best possible mean gaps over **all** physical heights are

```math
m_0-B_0-D/25\le14/5-78/25=-8/25,
\qquad m_\infty-B_\infty+V/25=-66/5+314/25=-16/25.       \tag{4}
```

The zero-side maximum uses the sharp old bound `D=78`, attained at
`q=t^2+t`. The infinity-side gap is constant across the pencil.

## 4. What the endpoint exclusion means

Keep one tensor and its eight physical gauge lists fixed. Write `nu,mu`
for the raw infinity and zero profiles, and let `Z_inf,Z_0` be their
exponential sums. Let `F_X,F_Y` be the paired physical sums. The critical
supporting-line margin is

```math
M(a,b,t)=a\log Z_\infty(t/a)+b\log Z_0(t/b)
                         -\log F_X(t)-\log F_Y(t),\quad a^2+b^2=1.
```

For fixed `x>0`, take `b->0+`, `a=sqrt(1-b^2)`, `t=b*x`. The uniform
finite-sum expansion and the determinant mean
`mean(nu)=B_inf-D/25` give

```math
\lim M(a,b,bx)/b
 =\log Z_0(x)+x(B_0+D/25)\ge8x/25>0.                  \tag{5}
```

The final bound uses `log Z_0(x)>=-x*m_0` and (4).
At the other axis, `a->0+`, `t=a*x`, the mean
`mean(mu)=B_0+V/25` similarly gives

```math
\lim M(a,b,ax)/a
 =\log Z_\infty(x)+x(B_\infty-V/25)\ge16x/25>0.         \tag{6}
```

For every fixed tensor, gauges and positive `x`, sufficiently close points
along these endpoint scalings therefore fail the strict negative criterion.
The conclusion is not uniform over unbounded changing gauges or changing
scale ratios, and is not an all-angle theorem.

## 5. Replay and relation to the completed numerical search

Run `make verify-pencil-uniform-transport`. The
[checker](../tools/check_pencil_uniform_transport.py) verifies both transport
plans, their all-field coefficient protection, natural-gauge attainment,
nilpotent conjugation identities, the full fresh determinant (3), and the
two rational gaps. It reuses the earlier source-pinned all-field opposite
degree theorem; it does not present that dependency as newly independently
recomputed. The [receipt](../evidence/pencil-uniform-transport-controls.json)
records the actual checks and source hashes.

Separately, all 469 registered free-height asymmetric trials have completed:
82,388 successful complete-profile objective evaluations, zero rejected
queries, and 469 successful optimizer returns in about 189 seconds. The best
numerical supporting value is `1.0000749992800997`, above the threshold one
and at the allowed angle boundary. The two-hour ceiling was not exhausted.
The [complete packet](../evidence/two-place-free-gauge-asymmetric-complete.json.xz)
retains all trial outcomes, full best tensor, parameters, profiles, physical
types, code and source commit/tree. This is finite local optimization;
success flags and lack of a found improvement are not a global exclusion.

The transport theorem supplies the global endpoint conclusion independently
of that finite search. Interior-angle constructions and different pencil
generators remain open. Substantial OpenAI assistant support was used;
no new Lean formalization or improved factoring algorithm is claimed.
