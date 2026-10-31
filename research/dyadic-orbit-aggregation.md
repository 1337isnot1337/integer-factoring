# Exact aggregation of complete dyadic fractional-linear families

## The constructive primitive

The [inverse-section framework](dyadic-inverse-sections.md) reduces the
new arithmetic direction to aggregation of section families. Here is a
family for which this can be done in polynomial bit time, even though
the family and its intervals can be exponentially large.

Let `R=2^r`, `1<=s<=r`, and `K=2^(r-s)`. Define permutations

```math
L_k(z)=\frac{z}{1+2^s k z}\pmod R,\qquad 0\le k<K.
```

Fix odd `alpha,beta`, arbitrary `gamma,delta`, and ordinary half-open
intervals `I,J` within `[0,R)`. There is an exact deterministic
polynomial-bit-time algorithm for

```math
\sum_{k=0}^{K-1}\#\{t\in I:
 \alpha L_k(\beta t+\gamma)+\delta\bmod R\in J\}.     \tag{1}
```

The parameters are represented by their residues modulo `R`. The
algorithm visits `r` valuations and a constant number of Euclidean
floor sums per valuation. It visits no individual `k` or interval
element. The time bound is polynomial in `r`, including coefficient
generation and all integer arithmetic.

This is an exact uniform family average. The inverse sections generated
by factoring have correlated varying affine parameters; (1) alone does
not aggregate them. The proved factoring score remains `1/2`.

## 1. Complete orbits are residue classes

The maps compose as `L_k L_l=L_(k+l)`. For `z` of valuation `v<r`,

```math
L_k(z)-z=\frac{-2^s k z^2}{1+2^s k z}.
```

Its valuation is `s+2v+v_2(k)`, unless that is at least `r`.
Consequently the orbit length is

```math
\ell_v=2^{\max(0,r-s-2v)}.
```

Put `h_v=min(r,s+2v)`. Every orbit value is congruent to `z` modulo
`2^h_v`. This residue class has exactly `ell_v` elements, and the
displayed valuation calculation shows the orbit has that same size.
Therefore the whole class is attained, with each value occurring

```math
\mu_v=K/\ell_v=2^{\min(r-s,2v)}                    \tag{2}
```

times as `k` runs through the group. For `z=0`, all `K` values are zero.

## 2. Affine changes retain a linear congruence

Set `z=beta*t+gamma mod R`. For a fixed valuation `v`, the possible
inputs satisfy

```math
t\equiv a_v:=(2^v-\gamma)\beta^{-1}\pmod{2^{v+1}}.
```

The orbit condition on an output `y in J` is

```math
y\equiv\alpha\beta t+\alpha\gamma+\delta\pmod{2^{h_v}}. \tag{3}
```

Let `D_v=2^(v+1)`, `A=alpha*beta`, `H=alpha*gamma+delta`, and let
`b_v=(A*a_v+H) mod D_v`. Write

```math
t=a_v+D_v u,\qquad y=b_v+D_v w.
```

Conditioning each ordinary interval on these progressions produces
ordinary intervals `U_v,W_v`. Equation (3) becomes

```math
w\equiv Au+\frac{Aa_v+H-b_v}{D_v}
\pmod{2^{h_v-v-1}}.                              \tag{4}
```

Thus the contribution of valuation `v` to (1) is `mu_v` times the
number of interval pairs satisfying the affine congruence (4).
There is also the zero-input contribution: `K` if
`-gamma/beta mod R` belongs to `I` and `delta mod R` belongs to `J`,
and zero otherwise. This proves the algorithm once (4) is counted.

## 3. Counting an affine congruence between two intervals

Let `d` be a power of two and `A` odd. For integer intervals of lengths
`n_U=q_U*d+r_U` and `n_W=q_W*d+r_W`, each residue appears `q_U` or
`q_W` times, with the excess occurrences forming one cyclic interval
of lengths `r_U,r_W` respectively. Hence

```math
\#\{(u,w)\in U\times W:w\equiv Au+c\pmod d\}
=q_Uq_Wd+q_Ur_W+q_Wr_U+E.                        \tag{5}
```

Here `E` counts the affine map on the two excess cyclic intervals.
Each cyclic interval splits into at most two ordinary intervals, so
at most four calls to the affine floor-sum primitive in the companion
note compute `E`. This also covers `d=1`, when both excess intervals
are empty. Equations (2)--(5) prove (1) with `O(r)` floor-sum calls.

All integers used have `O(r)` bits after reduction of the input
parameters. A floor sum takes `O(r)` Euclidean divisions on integers
of `O(r)` bits. School arithmetic therefore gives a polynomial bit
bound and polynomial workspace for the complete family, independently
of `K`, `|I|`, and `|J|` as cardinalities.

## 4. Why a factoring section is not yet this family average

At depth `s`, the section indexed by a low root `a` has the matrix

```math
\begin{pmatrix}-b&(N-ab)/2^s\\2^s&a\end{pmatrix},
\qquad b=Na^{-1}\bmod2^s.
```

After projective normalization it is an affine change of a lower
unipotent map, but the affine scales, offsets, denominator coefficient,
and interval selectors all depend on the same root. The uniform
complete `k` average in (1) does not preserve those correlations.
Replacing the actual family by a complete orbit would change the
answer. A precise next question is whether the root-indexed weights
admit a small exact decomposition into complete cyclic-family sums
with fixed affine parameters and tractable residual selectors.

Restricting the parameter `k` to a short interval is a materially
different problem. Even with a fixed unit input `z`, it restricts the
reciprocal coordinate `1/L_k(z)=1/z+2^s*k`; the interval condition on
the output then returns to an incomplete inverse-graph intersection.
The complete-group identity proves no bound for that restriction.

## 5. Exact evidence and replay

```sh
python3 tools/check_dyadic_orbit_aggregation.py --output /tmp/dyadic-orbits.json
cmp /tmp/dyadic-orbits.json evidence/dyadic-orbit-aggregation.json
```

The [implementation](../tools/dyadic_orbit_aggregation.py) uses the
closed orbit formulas. The [independent checker](../tools/check_dyadic_orbit_aggregation.py)
literally visits every group parameter and input for the selected
small controls. The [receipt](../evidence/dyadic-orbit-aggregation.json)
records 896 affine-residue pair checks, 6,003 literal family/interval
checks (including all interval pairs for three small groups), and six
large permutation controls. The largest complete group has `2^249`
elements; those large checks use the independently known full-interval
permutation identity. They are not a literal large-group enumeration
or a factoring benchmark.
