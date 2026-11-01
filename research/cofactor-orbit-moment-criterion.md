# A lower-tail moment criterion for uniform cofactor windows

## Result and scope

The [composite ideal closure](cofactor-ideal-closure.md) needs one
admissible power in each reached generator's prescribed window. This
note gives an explicit sufficient condition in terms of density and
lower-tail moments on that generator's own cyclic subgroup. It does
not replace subgroup density by density in the whole class group.

The new implication is elementary and unconditional. The density and
moment bounds that would give a better factoring exponent remain
unproved. In particular, no factoring-score improvement is claimed.

The implication uses the fact that sliding-window counts change by at
most one per step. A centered moment of order `2q` gives a period
factor of order `h^(1/(q+1))` in the sufficient window. A direct
single-phase argument would give only `h^(1/q)`.

[check_cofactor_orbit_moments.py](../tools/check_cofactor_orbit_moments.py)
also measures exact densities, maximum gaps, correlations, and moments
on complete small arithmetic orbits. It includes prime generators
and composite residual generators actually reached by the closure.

## 1. The actual admissibility predicate

Let `N>=3` be odd, with the discriminant choice used by the closure,

\[
 \Delta=\begin{cases}-N,&N\equiv3\pmod4,\\-4N,&N\equiv1\pmod4.
 \end{cases}
 \qquad H=\left\lfloor\sqrt{|\Delta|/3}\right\rfloor.
 \tag{1}
\]

Fix the known-prime base at cutoff `Y`, including the ramified-two
generator when appropriate. Let `P` be the product of those primes.
For a reduced form with leading coefficient `A`, define

\[
 S(A)=\gcd\bigl(A,P^{\mathrm{bitlength}(A)}\bigr),
 \qquad R(A)=A/S(A).
 \tag{2}
\]

The modular-power implementation computes the same gcd without
constructing the large ordinary power. Define `f` to be one precisely
when either the leading gcd gives a proper factor of `N`, or
`R(A)<=B`, where `B` is the chosen residual cap. The first event ends
the factoring search. The second gives an admissible partial ideal
identity. Useful target involutions are additional early successes;
they need not be counted by `f`.

This is a function of the reduced class. It contains the identity,
whose leading coefficient is one. It is also invariant under class
inversion, since the reduced inverse has the same leading coefficient.

Fix one allowed generator `g`, and let

\[
 G_g=\langle g\rangle,\qquad h_g=|G_g|,\qquad
 s_g=\sum_{x\in G_g}f(x),\qquad \delta_g=s_g/h_g>0.
 \tag{3}
\]

For a positive integer `W`, put

\[
 A_x(W)=\sum_{j=0}^{W-1}f(xg^j),\qquad
 \mu_g=\delta_g W.
 \tag{4}
\]

Every average below is over `G_g`. The identity
`h_g^-1 sum_x A_x(W)=delta_g W` follows by counting each good class
once at each of the `W` positions. Windows may wrap around the
period, and `W` need not be smaller than `h_g`.

## 2. A pointwise lower deviation forces a moment contribution

The sliding counts satisfy

\[
 A_{xg}(W)-A_x(W)=f(xg^W)-f(x)\in\{-1,0,1\}.
 \tag{5}
\]

Suppose `A_x(W)=a<mu_g`, and write `d=mu_g-a`. For every integer
`j>=0`,

\[
 A_{xg^j}(W)\le a+j.
 \tag{6}
\]

Also `d<h_g`: all `h_g` counts are at most `a+h_g-1`, and their
average is `mu_g`. Thus the phases
`j=0,...,ceil(d)-1` are distinct. For every positive integer `q`,

\[
 \begin{aligned}
 \frac1{h_g}\sum_{y\in G_g}
       (\mu_g-A_y(W))_+^{2q}
 &\ge\frac1{h_g}\sum_{j=0}^{\lceil d\rceil-1}(d-j)^{2q}\\
 &\ge\frac1{h_g}\int_0^d(d-t)^{2q}\,dt\\
 &=\frac{d^{2q+1}}{(2q+1)h_g}.
 \end{aligned}
 \tag{7}
\]

Here `u_+=max(u,0)`. The second inequality is a left-endpoint
Riemann-sum bound for a decreasing nonnegative function. No mixing
or independence assumption enters (7).

In particular, one empty window gives

\[
 \boxed{\displaystyle
 \frac1{h_g}\sum_y(\mu_g-A_y(W))_+^{2q}
 \ge\frac{\mu_g^{2q+1}}{(2q+1)h_g}.}
 \tag{8}
\]

This is stronger than retaining only the empty phase's contribution
`mu_g^(2q)/h_g`.

## 3. Explicit conditional window

Let `C,D,q` be positive integers. Suppose the subgroup satisfies

\[
 \delta_g\ge1/D,
 \qquad
 \frac1{h_g}\sum_x(\mu_g-A_x(W))_+^{2q}
       \le(CqW\delta_g)^q.
 \tag{9}
\]

Only a lower-tail moment is required. The stronger bound with
`abs(mu_g-A_x(W))^(2q)` implies (9), but upper outliers do not need
to be controlled to rule out empty windows.

Combining (8) and (9), an empty window would force

\[
 \mu_g^{q+1}\le(2q+1)(Cq)^q h_g.
 \tag{10}
\]

Given a certified upper bound `Hbar>=h_g`, define the integer

\[
 T=\left\lceil
       \bigl((2q+1)(Cq)^q\overline H\bigr)^{1/(q+1)}
    \right\rceil,
 \qquad
 \boxed{W=D(T+1).}
 \tag{11}
\]

If (9) holds at this particular `W`, then every phase has an
admissible query in its window. Indeed,
`mu_g>=W/D=T+1` contradicts (10).

The algorithm can compute (11) by integer powering and binary search.
It does not need the exact group order, the exact density, the moments,
or an oracle that locates a good phase. It simply examines the `W`
consecutive powers already prescribed by the closure algorithm.
The functions `C,D,q` must be supplied by a uniform theorem; the
finite checker measures them only for its complete recorded orbits.

For comparison, the direct single-phase argument gives
`mu_g<=Cq h_g^(1/q)` for an empty window. Thus
`W=2CqD ceil(Hbar^(1/q))` is another sufficient conditional choice.
Equation (11) uses the additional sliding-window structure.

### An elementary input-only period bound

Every class has a reduced form with `1<=A<=H` and
`-A<B_middle<=A`. There are at most `2A` possible middle coefficients
for a fixed `A`; the last coefficient is determined by the
discriminant. Therefore

\[
 h_g\le|\mathrm{Cl}(\Delta)|\le H(H+1)\le2N.
 \tag{12}
\]

For the final inequality, use `H^2<=4N/3` and `H<=2N/3`, the latter
valid for `N>=3`. Primitivity, parity, and the remaining reduction
conditions can only decrease this count. Consequently
`Hbar=2N` is an effective input-only choice, independent of the
pinned analytic generator theorem.

### Several admissible queries per window

Let `Q>=1`. If `A_x(W)<Q` and `mu_g>=2Q`, then its lower deviation
is at least `mu_g/2`. Equation (7) now yields

\[
 \mu_g^{q+1}
 \le 2^{2q+1}(2q+1)(Cq)^q h_g
 \tag{13}
\]

under (9). Thus the explicit window

\[
 W=D\left(2Q+1+
 \left\lceil\bigl(2^{2q+1}(2q+1)(Cq)^q\overline H\bigr)
                       ^{1/(q+1)}\right\rceil\right)
 \tag{14}
\]

contains at least `Q` admissible powers at every phase, provided
(9) holds at that `W`. This can feed the
[per-generator collision quota](cofactor-collision-quota.md).
The powers are distinct integer exponents; their classes may repeat.
The one-hit closure uses only the simpler bound (11).

## 4. What uniform estimates would imply for the score

Assume the initial prime cutoff generates the intended class group,
under the existing [explicit generator contract](explicit-prime-form-cutoffs.md).
For every initial generator and every allowed residual ideal key
of norm at most `B`, impose (9) with common effective parameters
`C(n),D(n),q(n)` and `Hbar=2N`. It is enough to impose these conditions
on the generators actually reached by the deterministic closure;
the all-allowed-keys version is a stronger, simpler sufficient
hypothesis.

There are at most `K_B=B bitlength(B)` residual keys. The closure
therefore uses at most `(m+K_B)W` norm queries, where `m<=Y+1` is
the initial base size. The integral elimination theorem charges
the bits of every expanded word: definition depth is at most `K_B`,
and selected exponents have `O(log(n+W))` bits.

If `Y` is allowed to grow beyond a polynomial in `n`, charge its
whole cost as well. Prime enumeration, the prime product's bit
length, known-prime trial division, group computations, and relation
extraction all have cost polynomial in `n,Y,K_B,W`. Thus

\[
 \log T_{\mathrm{factor}}
 =O\!\left(
    \log(n+1)+\log(Y+1)+\log(B+1)+\log D+
    \log(Cq+1)+\frac n{q+1}
   \right).
 \tag{15}
\]

The uniform exponent constants in that polynomial affect the
constant in (15), not the power of `log n` in its denominator.
All rejected powers are charged in `W`. No exponentially large
prime product or relation vector is treated as a constant-cost
object.

For any fixed `alpha>0`, suppose

\[
 q=\left\lceil(\log(n+2))^\alpha\right\rceil,
 \qquad C\le n^{O(1)},
 \qquad
 \log Y+\log B+\log D
    =O\!\left(\frac n{(\log n)^\alpha}\right),
 \tag{16}
\]

and suppose (9) holds uniformly at the window (11). Then (15)
would give

\[
 T_{\mathrm{factor}}
 \le\exp\!\left(O\!\left(\frac n{(\log n)^\alpha}\right)\right).
 \tag{17}
\]

For a concrete conditional `alpha=1` instance, let
`L=ceil(log_2(n+2))`, `a=ceil(n/L^2)`, and choose

\[
 q=L,\quad
 Y=\max\{X_{\mathrm{gen}}(n),2^a\},\quad
 B=Y^2,\quad D=Y^4,\quad C=n^c
 \tag{18}
\]

for a fixed effective `c`. Here `X_gen` is the polynomial generator
cutoff already supplied by the repository. The remaining claims
would be the subgroup-specific density bound `delta_g>=Y^-4` and
the lower-tail estimate (9) at (11), uniformly for the specified
generators. Those two claims are not established by this note.
Their formulation in (18) is a sufficient research target, not
evidence of a new score.

## 5. Exact correlations identify part of the missing estimate

For one certified orbit, write `f_i=f(g^i)` and

\[
 c_r=\frac1h\sum_{i=0}^{h-1}f_i f_{i+r},\qquad
 \delta=\frac1h\sum_i f_i.
\]

The complete second centered moment is exactly

\[
 \frac1h\sum_i(A_i(W)-\delta W)^2
 =W\delta(1-\delta)
    +2\sum_{r=1}^{W-1}(W-r)(c_r-\delta^2).
 \tag{19}
\]

This follows by expanding the square and grouping pairs by their
displacement. It does not assert that the covariances are zero.
The checker verifies (19) independently from its direct window
counts at all widths through 32 on every complete orbit.

The higher absolute moment is likewise the sum of the centered
`2q`-point correlations over all `2q` choices of positions in the
window. A suitable bound on that sum would imply the stronger
absolute version of (9). A direct lower-tail estimate may require
less. Neither the pinned prime-class bias theorem nor density in
the full class group supplies these correlations along one fixed
power orbit.

## 6. Finite arithmetic controls

The [evidence](../evidence/cofactor-orbit-moments.json) enumerates ten
complete orbits for three inputs, with no period-cap truncation.
The period is certified by the first return of the exact composed
form to the identity. Every density is computed on that individual
orbit. Inverse symmetry is checked at every phase.

The larger input is `N=1000036000099`, with `Y=19` and `B=10000`.
Its recorded prime and residual generators have the following
exact orbit statistics:

| Generator | Recorded norm | Period | Admissible classes | Maximum inadmissible gap |
|---|---:|---:|---:|---:|
| Prime form | 5 | 198850 | 42053 | 100 |
| Prime form | 7 | 198850 | 42053 | 64 |
| Reached residual ideal | 9169 | 198850 | 42053 | 57 |
| Reached residual ideal | 167 | 198850 | 42053 | 39 |

Equal period and equal density therefore do not determine the
maximum gap. The ordering induced by the particular generator
matters.

There is also a concrete failure of the stronger absolute-moment
guess with `C=1`. At `q=6,W=480`, the prime-five orbit's centered
twelfth moment is approximately `36.0467` times `(q W delta)^6`.
For the prime-seven orbit the ratio is approximately `7.08593`.
The evidence stores the exact integer numerator and denominator,
so these failures do not depend on rounding. Both require at least
`C=2` at that particular window. This is not a disproof of a
uniform polynomial-`C` bound, nor a substitute for checking the
weaker lower-tail condition at (11).

For each orbit the checker also evaluates the lower-tail condition
at (11), for `C=1,2,4,8`, using both its certified actual period and
the bound `2N`. These finite checks use the **measured** denominator
`D=ceil(h/s)`. The fact that a period-bound field is input-only
does not make that measured density bound a uniform input-only
theorem.

The weaker lower-tail hypothesis with `C=1` also fails in a recorded
arithmetic case: for the prime-five generator, `q=6`, `D=5`, and
`Hbar=2N`, formula (11) gives `W=1925`, and the exact lower-tail
moment exceeds the right side of (9). At the window recomputed
with `C=2`, that orbit passes. All ten recorded orbits pass the
lower-tail test with `C=2` at their respective windows (11), for
both period bounds. This is a finite observation only.

The remaining controls check integer ceiling roots, the exact
Lipschitz lower bound on clustered cyclic sets with empty windows,
and a nontrivial abstract orbit whose sufficient window is shorter
than its period. Small experimental prime cutoffs are not asserted
to satisfy the full generating-base contract.

```sh
python3 tools/check_cofactor_orbit_moments.py \
  --output evidence/cofactor-orbit-moments.json
```

