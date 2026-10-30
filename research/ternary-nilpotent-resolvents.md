# Exact ternary nilpotent resolvents for ordered crossings

## A fixed full crossing endpoint

Work over `F_3`. Let the external alphabet and auxiliary space both have
dimension three. On three tensor factors write

```math
A=\operatorname{SWAP}_{x,a},\qquad
B=\operatorname{SWAP}_{y,a},\qquad
E=\operatorname{SWAP}_{x,y}.
```

The operator `T=A E` is a three-cycle. Therefore `N=I-T` satisfies
`N³=0`. The polynomial

```math
K(z)=(I+zN+z^2N^2)A
     =(1+z+z^2)A+(z^2-z)B+z^2E                         \tag{1}
```

has `K(0)=A` and `K(1)=E tensor I_aux`. Its full operator determinant is
the nonzero constant `det(A)`, equal to `-1` in the displayed dimension.
The determinant is not necessarily one.

The constant term acts only on the first physical factor and the auxiliary
factor, so it has one column channel when both auxiliary indices are put
on that physical half. Ordered products of (1), followed by any normalized
auxiliary boundary, give exact polynomial crossings. The Smith and path
interfaces apply without a cross-half commutation assumption.

The local plain row infinity profile is `(2,2,2,2,2,2,2,2,2)`. The useful
column zero profile is `(0,1,1,1,1,1,1,1,1)`. These profiles give `B=81`
and `rho=1`; they are not competitive. Blocking two copies keeps all row
infinity exponents equal to four and also gives the noncompetitive plain
bound `B=9^4=6561`.

## 1. A finite family with guaranteed nilpotence

Let `P=J_3` be a nilpotent Jordan block. Choose any polynomial

```math
f(X,Y)\in F_3[X,Y]/(X^3,Y^3),\qquad f(0,0)=1.
```

There are exactly `3^8=6561` choices. Replace `A` by

```math
\widetilde A=f(P_x,P_a)A,\qquad T=\widetilde A E.
```

Conjugation by the undeformed three-cycle permutes the three commuting
operators `P_x,P_y,P_a`. Hence `T³` is the product of the three cyclic
copies of `f`. Each such copy has cube one, by the Frobenius identity and
`P_x³=P_y³=P_a³=0`. The copies commute, so `T^9=I` and
`N=I-T` has `N^9=0`.

Thus every choice supplies the explicit polynomial

```math
K(z)=\sum_{k=0}^{8}z^kN^k\widetilde A,\qquad K(1)=E.       \tag{2}
```

Trailing zero coefficients can be removed. An exhaustive local screen
found degrees 2, 4, 5, 8 with multiplicities 9, 18, 702, 5832, respectively.
Every endpoint and nilpotence identity was checked exactly. None improved
the plain Smith-profile bound below `B=81`. This finite screen is not an
exclusion of other gauges, blocked ranks, or other polynomial families.

## 2. Nontriangular rank-one perturbations

The preceding family has a common degree filtration. A different search
uses

```math
\widetilde A=(I+uv^T)A,
\qquad u,v\in F_3^9,\quad v^Tu=0.                         \tag{3}
```

The condition makes the pair-space perturbation invertible. Reshape `u,v`
as three-by-three matrices `U,V` in row-major order. Put

```math
B_0=UV,\qquad C_0=(VU)^T.
```

There is a short exact nilpotence test:

```math
\det(\lambda I_{27}-(I-\widetilde A E))
 =\lambda^{18}
  \det\!\left(\lambda^3I_3+
       (\lambda^2+\lambda+1)B_0+(1-\lambda)C_0\right).     \tag{4}
```

Thus the full operator is nilpotent if and only if the three-by-three
polynomial determinant on the right is `lambda^9`.

To prove (4), let `T_0=A E`, `N_0=I-T_0`, and let `U_inj` insert `u`
in the first and auxiliary positions while retaining the middle position.
Let `V_proj` contract those positions against `v`. The perturbation of
`N_0` has rank at most three. Since `N_0³=0`, the matrix determinant lemma
reduces its characteristic polynomial to a three-by-three determinant.
The three contractions are

```math
V_{\rm proj}T_0U_{\rm inj}=B_0,\qquad
V_{\rm proj}T_0N_0U_{\rm inj}=B_0-C_0,\qquad
V_{\rm proj}T_0N_0^2U_{\rm inj}=B_0+C_0.
```

The last equality uses `v^Tu=0` and `-2=1` in `F_3`. Clearing the
three denominators `lambda^3` gives the power `lambda^18` and (4).
Cayley--Hamilton then bounds the nilpotence index by 27. For every survivor,
we construct `sum_(k=0)^26 z^k N^k tilde(A)` and check both nilpotence and
the complete endpoint directly.

A deterministic sample of 200,000 pairs gave 3,607 survivors, with degrees
5 through 11. Plain local profiles again did not improve `B=81`. Unlike
(1), some profiles contain negative infinity exponents. One fully replayed
example has

```text
U = [1,2,0, 0,2,0, 1,2,1]
V = [1,1,0, 2,2,0, 2,0,0]
```

and the row infinity profile

```math
(11,11,11,2,2,2,2,2,-1).
```

The corresponding useful column zero profile is `(0,1,1,1,1,1,1,1,1)`.
The displayed profile is not competitive; its role is a concrete test of
larger nilpotence and negative infinity orders.

## 3. Reproducible polynomial arithmetic

`ternary_resolvent_probe.cpp` is a standalone C++20 program. It needs no
external arithmetic library. It contains packed exact `F_3` polynomial
arithmetic, truncated local Smith elimination, the three-cycle constructor,
and the rank-one constructor with its small determinant certificate.
Its two fixed probes are the undeformed degree-two model and the degree-11
model above. For each it computes all four extended cuts at zero and at
infinity.

A complete full-rank Smith elimination modulo `z^P` identifies the reported
orders when every pivot has order below `P`. The probe fails if the chosen
precision is exhausted; it does not interpret that event as infinite order.
The fixed probes complete at precision 64.

`check_ternary_resolvent_probe.py` independently rebuilds the operators with
ordinary integer arithmetic modulo three. It checks nilpotence and every
endpoint entry, then forms ordinary block Toeplitz matrices. For each
reported order list `(nu_i)`, it verifies

```math
\operatorname{rank}T_h=\sum_i\max(h-\nu_i,0).
```

It checks precisions 1, 2, 3 and one beyond the largest reported order,
with duplicates removed. The two full probes pass 50 Toeplitz checks.
This verifies both the finite-field engine and the selected profiles
without reusing packed arithmetic or the Smith algorithm.

Replay from a repository root after placing the code under `tools`:

```sh
g++ -O3 -std=c++20 tools/ternary_resolvent_probe.cpp -o /tmp/ternary-resolvent-probe
/tmp/ternary-resolvent-probe evidence/ternary-resolvent-probe.json
python tools/check_ternary_resolvent_probe.py --output evidence/ternary-resolvent-probe-audit.json
```

The endpoint constructions and determinant reduction are ordinary proofs.
The finite searches and profile checks supply no improved factoring exponent.

The two larger screens in Sections 1--2 are reproducible separately. Their
sources share the three `f3_*_core.hpp` headers in `tools`. They check each
constructed model's full endpoint and direct nilpotence. The ordinary dense
Toeplitz checker independently audits the two selected representatives;
it does not independently replay every local profile in the larger screens.

```sh
g++ -O3 -std=c++20 tools/f3_nilpotent_swap_deformations.cpp -o /tmp/f3-deformations
/tmp/f3-deformations evidence/f3-nilpotent-deformations-all.json 0 6561
g++ -O3 -std=c++20 tools/f3_rank_one_resolvent_search.cpp -o /tmp/f3-rank-one
/tmp/f3-rank-one 200000 evidence/f3-rank-one-resolvents-200k.json
```

The endpoint proofs, exact finite screens, and independent representative
checks have the distinct scopes stated above. This work used substantial
OpenAI assistant support.
