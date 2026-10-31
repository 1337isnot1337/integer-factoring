# Searches with the exact prime-state filter profiles

## Result and scope

The [prime-state diagonalization](prime-identity-filter-diagonalization.md)
turns one polynomial identity direction into an exact finite search
problem. We tested four characteristic-five families: monomial digit
weights, their finite-difference duals, continuous cancellation branches,
and polynomial weights synthesized from degree and order constraints.
The recorded searches found no negative critical pressure and do not
improve the proved factoring score of `1/2`.

The main constructive result is a reproducible search interface with
exact finite-field cancellations and exact polynomial witnesses. A
negative finite search does not prove an obstruction. The degree grids,
charge grids, scale grids, and selected Pareto priorities below are
finite. Correlated filters on several identity directions are a separate
family and are not covered by these experiments.

## 1. Common pressure convention

Write `P` for the odd ramification, `n` for the interval length, `J` for
the aggregate zero-cut charge, and `t=P log(s)`. Let `b_k` be the exact
zero order of the finite difference `beta_k`, and let `i_d` be the
exact infinity order of digit `d` in the selected mode.

The two factors that change from the Gaussian background are

\[
\sum_{k=0}^{p-1}\exp\!\left(-[b_k/P-(1-J)k]t\right),
\qquad
\sum_{d=0}^{p-1}\exp(-i_d t/P).
\]

Their Gaussian counterparts have energies `(2n-1+J)k` and
`(1-2n)d`. Every other factor remains explicit. The pressure code
includes all those spectators and the two independent paired-type
physical partition functions, at `rho=1/sqrt(2)`.

`tools/prime_state_filter_pressure.py` separates ground energies from
bounded thermal terms. This limits cancellation when a common Laurent
monomial clears the weights. Its independent numerical controls compare
`96` complete physical-profile sums and check `96` common scalar shifts
of degree one million. The maximum discrepancies in the saved run are
about `5.4e-15` and `2.2e-11`, respectively. These are numerical controls,
not interval-arithmetic pressure certificates.

## 2. Monomial weights and their duals

For monomial weights, put

\[
\alpha_d(q)=q^{-e_d},\qquad q=-v,\qquad e_0=0.
\]

For each displayed `(n,P)`, the other four exponents range over the
distinct values in

\[
\{-2nP,-nP,-2,-1,0,1,2,nP,2nP\}.
\]

This includes all `52` set partitions describing equal exponents of
the five digit labels. It also deliberately includes differences of
`nP`, where coefficients from different powers of `q^(nP)` can cancel.
The five Laurent finite differences are expanded exactly before their
orders are read. Generic leading-coefficient assumptions are not used.

The monomial runs use `n=1,2,4,8,16,32,64` at `P=1`, and `n=16,32`
at `P=5`. There are `625`, `2401`, or `6561` candidates per run,
depending on how many exponent levels coincide. Each candidate is
evaluated on the recorded finite charge/scale grid; the strongest
sixty candidates receive bounded continuous refinements.

The dual family instead prescribes

\[
\beta_k(q)=(-1)^k q^{f_k},\qquad f_k=nPk+e_k,
\]

with the same deviation levels. The inverse finite-difference transform
produces explicit Laurent-polynomial `alpha_d`. At `q=1`, its parameter
`c=1-q^(nP)` vanishes, so every `alpha_d(1)=1`. Thus these filters have
the required endpoint. The checker recomputes their finite differences
and verifies the prescribed monomials exactly.

Dual runs use `n=8,16,32,64` at `P=1` and `n=32` at `P=5`. They test
`6561` candidates each. Their best recorded filters are the Gaussian
weights. For the longer monomial runs the Gaussian weights are also
best. The approximate optimized Gaussian critical pressures are:

| `n` | Critical pressure |
|---:|---:|
| 8 | 0.304991281508 |
| 16 | 0.191309757166 |
| 32 | 0.113988916176 |
| 64 | 0.066190686428 |

Small intervals also admit charge/scale limits with vanishing physical
type entropy. Very small pressure values in such limits are not a
threshold witness. The recorded refinements have finite bounds, retain
better grid points, and use the complete pressure.

## 3. Exact cancellation branches with continuous exponents

A branch partitions digit labels into groups. Within each group, write

\[
e_d/P=E_g+n m_d.
\]

The first digit in each group has `m_d=0`; every other digit has an
offset in `{-1,0,1}`. With five digits this gives

\[
\sum_{g=1}^{5}S(5,g)3^{5-g}=742
\]

branches, where `S(5,g)` is a Stirling number of the second kind.
For each group, coefficient terms are combined exactly into a Laurent
polynomial `R_(k,g)(u)`, with `u=q^(nP)`. Between distinct groups the
generic order formula is

\[
b_k/P=\min_g\{2n\,\mathrm{ord}_u R_{k,g}-2E_g\}.
\]

At a fixed charge and scale, the pressure is convex in these group
exponents. The implementation introduces the five order variables and
their linear upper constraints, then solves the resulting smooth
convex problem numerically. The group exponents are bounded by
`|E_g|<=4n`. This is a search procedure; numerical optimizer success is
not a certificate of global optimality.

The `n=32` run evaluates all `742` branches at `21` charge/scale pairs.
The strongest fifty branches receive an outer charge/scale refinement.
Their exponents are then rounded at the odd ramification `P=1000001`,
and every finite difference is recomputed in `F5`. This replay combines
any extra cancellations between different groups. Generic branch
limits and these exact integer replays are saved separately.

The Gaussian remains best in both stages. Its refined pressure is
approximately `0.113988916176`. In particular, the finite branch search
does not supply a negative critical pressure.

## 4. Synthesized polynomial weights

The residue-space algorithm in
`tools/polynomial_state_filter_synthesis.py` decides whether there are
polynomials with

\[
\deg\alpha_j\le A_j,\qquad
\alpha_j(1)=1,\qquad
\mathrm{ord}_q\beta_k\ge B_k.
\]

This is exact linear algebra over the prime field. Its dependence on
`q` through `q^(nP)` permits residue compression. Every returned sparse
witness is substituted back into the finite-difference transform and
checked directly. A common scalar clearing permits Laurent filters.

For `p=5,n=32`, the search uses common clearing `L=4nP`, fixes `A_0=L`,
and ranges each other degree bound over `L` plus the nine offsets

\[
\{-nP,-nP/2,-2,-1,0,1,2,nP/2,nP\}.
\]

There are `6561` degree vectors. For each one, four priority orders
greedily raise the five zero-order bounds using exact feasibility
decisions. The returned witness may improve on the requested bounds;
its actual degrees and orders are used in every pressure evaluation.

The completed `P=1` and `P=5` runs each produce `9310` distinct actual
degree/order profiles and construct `26244` witnesses. They make
`979338` and `1265336` feasibility decisions, respectively. The Gaussian
is best on the `48` charge/scale points tested.

This enumerates the degree grid and the four specified priority
choices. It does not enumerate every feasible Pareto profile under
those bounds, and it does not prove Gaussian dominance for polynomial
filters.

## 5. Reproduction and evidence

The aggregate record is
`evidence/prime-state-filter-searches.json`. Exact order formulas and
their independent physical-matrix controls are in
`evidence/prime-state-filter-profiles.json`. Numerical pressure controls
are in `evidence/prime-state-filter-pressure-controls.json`.

Representative portable commands are:

```sh
python3 tools/search_prime_state_filters.py --family monomial \
  --n 32 --ramification 1 --output monomial-replay.json
python3 tools/search_prime_state_filters.py --family dual \
  --n 32 --ramification 5 --output dual-replay.json
python3 tools/search_prime_filter_branches.py \
  --n 32 --output branch-replay.json
python3 tools/search_prime_polynomial_filters.py \
  --n 32 --ramification 5 --output polynomial-replay.json
python3 tools/check_prime_state_filter_pressure.py
```

The search commands require NumPy and SciPy. The exact finite-difference
profile and synthesis routines use only the standard library. A future
negative pressure candidate needs an independent rigorous numerical
certificate before it can support a factoring exponent.
