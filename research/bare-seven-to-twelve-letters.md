# Box-minor Kronecker budgets close the bare criterion through twelve letters

## Result and scope

For the bare, ungauged two-place Smith criterion of the
[shared-minor](shared-minor-smith-floor.md) and
[leading-rank](leading-rank-entropy-floor.md) notes, every polynomial
crossing on 7 through 12 physical letters has

```math
\rho_{\rm bare}=\frac{\log B(s)}{4\log q}>\frac1{\sqrt2}
```

at every scale s>=1. The statement holds over every field, at every
polynomial degree, and with a fixed auxiliary dimension whose two indices
stay together on one physical side of each cut. With the earlier floor
through six letters, **the bare criterion cannot give a half-rate
improvement on any alphabet with q<=12.** The certified floors are:

| q | exact floor on rho_bare |
| --- | --- |
| 7 | 181639/250000 = 0.726556 |
| 8 | 0.721261 |
| 9 | 0.717261 |
| 10 | 0.713605 |
| 11 | 0.710675 |
| 12 | 0.707627 |

Each floor satisfies 2 rho^2>1 exactly. The new ingredient is a family of
proved necessary conditions coupling the two cuts through rectangular
**box minors** (Theorem 1). They hold for every q. In particular, they
violate the formal seven-letter optimizer of the shared-minor relaxation
directly, at every ramification, without the monomial-determinant hypothesis
of the [inverse-degree theorem](inverse-degree-crossing-compatibility.md).
The proof does not use the unique-leading-channel case split.

The certificates use only the box budgets of Theorem 1. Section 5 also
proves a square-box **adjugate inequality**. Integer profiles that satisfy
every box budget and both adjugate inequalities come within 1.2e-5 of every
floor in the table (their pressures are evaluated in floating point). So the
certificates are essentially sharp for all constraints proved here. At
**q=13 these constraints do not close the criterion**: Section 6 gives an
exact integer profile that satisfies every box budget and every adjugate
inequality and has B(s)<13^(2 sqrt 2). That profile is the precise
thirteen-letter construction target. It is not a tensor, and no realizing
crossing is known.

This is a lower bound on one sufficient certificate. It is not a lower bound
on actual coefficient ranks or on factoring time. Independent physical
gauges and types, multipoint jets, cancellation modules, coefficient
recombination, and alphabets q>=13 remain outside it. The proved factoring
score remains 1/2.

**Relation to the inverse-degree aspect bounds.** The
[aspect-bound note](inverse-degree-aspect-bounds.md) treats q=7,8,9,10. It
takes the scalar crossing obtained by reversing Phi(R), where R(1)=I, R has
a polynomial inverse of degree d>=1, and deg R=M. It shows rho>17/24 under
the bare ungauged pressure when M>=A_q d, or when M>=A_q and d<=2. It also
excludes all 2^31 binary level weightings of the fixed nine-letter flag, with
rho>5/7 or rho>17/24 according to the case. Those paths are scalar polynomial
crossings with endpoint E. Theorem 1 therefore applies to them with no
aspect, inverse-degree or degree hypothesis.

For q=7..10, this note implies every bare-criterion exclusion stated there.
The floors 0.726556, 0.721261, 0.717261 and 0.713605 all exceed
17/24≈0.708333, and the nine-letter floor exceeds 5/7≈0.714286.

That note also proves things this one does not imply:

- the per-crossing inequality log B(s)>=H_q−eta_q k log s, in terms of the
  second infinity order k;
- the bound k<=2d linking that order to the inverse degree;
- the exact first-jet rank data of the nine-letter flag: rank 80, and largest
  direction ranks 22, 22 and 18.

These facts concern the crossings themselves, not only the bare floor. They
remain relevant if the certificate changes, for example to gauged criteria,
to which the box budgets have not been extended.

**Credits.** The work builds directly on
[PR #40](https://github.com/1337isnot1337/integer-factoring/pull/40) at
commit c3e6ed0 (c3e6ed076a6314d67e3692884df2302033774b45). It uses that PR's
[all-degree](all-degree-small-alphabet-smith-floor.md),
[shared-minor](shared-minor-smith-floor.md),
[inverse-degree](inverse-degree-crossing-compatibility.md),
[leading-rank](leading-rank-entropy-floor.md) and
[aspect-bound](inverse-degree-aspect-bounds.md) notes. It also uses the
common-minor idea of
[PR #39's universal binary theorem](https://github.com/1337isnot1337/integer-factoring/blob/0ac85be63d6eff434eddd92a32eb5c7e4ee75cb2/research/universal-binary-smith-floor.md).
The proofs and certificates were developed with substantial Anthropic
Claude assistance. No Lean proof is claimed.

## 1. Setting

Use the normalization of the shared-minor note. K(z) is a polynomial
crossing of actual degree M>0 with K(1)=E=[x=v][y=u], and n=q^2.

- R is the cut (x,y)|(u,v). Its Smith orders at infinity, regularized by
  w^M, are mu_1<=...<=mu_n.
- C is the cut (x,u)|(y,v). Its Smith orders at zero are nu_1<=...<=nu_n.
- After removing a common monomial, which does not change B, both
  minimum orders are zero: mu_1=nu_1=0.

If M=0, both cuts are constant permutation matrices, so B=n^2 and
rho_bare=1.

Exchanging the letters y and u maps E to itself and swaps the two cuts.
So the other assignment of places needs no separate treatment. Set
mu_0=0. The pressure is

```math
B(s)=\Big(\sum_i s^{M-\mu_i}\Big)\Big(\sum_i s^{-\nu_i}\Big).
```

Write H(w)=w^M K(1/w), with flattenings H_R and H_C. The Smith form of
H_R over F[[w]] gives, for every j>=0,

```math
H_R\equiv\sum_{i\le j}w^{\mu_i}a_i(w)b_i(w)^T \pmod{w^{\mu_{j+1}}}. \tag{1}
```

## 2. Kronecker layers on a box

Choose letter sets S_1 and S_2 of sizes alpha<=beta. Let N=alpha beta and
r=alpha^2. The **C-box** is the N-by-N submatrix of C with rows
(x,u) in S_1 x S_2 and columns (y,v) in S_2 x S_1. At z=1 it is the
restriction of E, which is a permutation matrix on these index sets.

**Lemma 1.** Let kappa be the Smith orders at w=0 of the C-box of H. Then
for every j>=0, #{k : kappa_k < mu_{j+1}} <= jr.

*Proof.* In tensor form, (1) says
H(x,y,u,v) ≡ Σ_{i≤j} w^{μ_i} a_i(x,y) b_i(u,v) modulo w^{μ_{j+1}}.
On the box, each term is the Kronecker product of the alpha-by-beta matrix
[a_i(x,y)] and the beta-by-alpha matrix [b_i(u,v)]. Over the principal
ideal domain F[[w]], each factor factors through F[[w]]^alpha. So each
product factors through F[[w]]^r.

Reduce modulo k=mu_{j+1}. The image of the box is contained in a module
generated by jr elements. Over the principal ideal domain F[[w]], a
submodule of an m-generated module is m-generated. On the other hand, the
reduced image is a direct sum of nonzero cyclic modules, one for each
kappa_k<k. Its minimal number of generators is therefore #{kappa_k<k}.
QED.

## 3. The box budgets

**Theorem 1.** For every 1<=alpha<=beta<=q, with N=alpha beta and
r=alpha^2,

```math
\sum_{j\ge0}(\mu_{j+1}-\mu_j)(N-jr)^+ +\sum_{i\le N}\nu_i\le NM,
\tag{C}
```

```math
\sum_{j\ge0}(\nu_{j+1}-\nu_j)(N-jr)^+ +\sum_{i\le N}\mu_i\le NM.
\tag{R}
```

*Proof of (C).* Let Delta(z) be the determinant of the C-box of K. Then
Delta(1)=±1, so Delta is nonzero of degree at most NM. The determinant of
the C-box of H is w^{NM} Delta(1/w), so its order at zero equals
NM−deg Delta. That order is Σκ_k.

Lemma 1 bounds this sum layer by layer. For μ_j<m<=μ_{j+1}, at least
(N−jr)^+ of the κ_k are >=m, so Σκ is at least the first sum in (C). The
box is an N-minor of C, so ord_0 Delta >= Σ_{i≤N} ν_i. Since
ord_0 Delta <= deg Delta, adding the two bounds gives (C).

*Proof of (R).* The same argument applies with places and cuts exchanged.
At zero, C≡Σ_{i≤j} z^{ν_i} c_i d_i^T. The **R-box**, with rows
(x,y) in S_1 x S_2 and columns (u,v) in S_2 x S_1, then has Kronecker
layers. Its determinant has order at infinity at least Σ_{i≤N} μ_i. QED.

Special cases recover the earlier budgets:

- alpha=1 gives r=1, and both inequalities become the shared-minor
  budget Σ_{i≤β}(μ_i+ν_i)<=βM.
- alpha=beta gives only the prefix budgets Σ_{i≤N}ν_i<=NM. At alpha=q
  this is the determinant budget.

The new information comes from 1<alpha<beta. For example, with
alpha=q−1 and beta=q, (C) reads

```math
(q-1)\mu_2+\sum_{i\le q(q-1)}\nu_i\le q(q-1)M .
```

**Fixed auxiliary dimension.** Contract the auxiliary indices with the
boundary vectors that give the crossing at one. The contraction is a
constant linear map on one side of each rectangular cut. It carries the
rank-one layers of (1) to rank-one layers of the scalar contraction. By
Cauchy–Binet, it can only increase determinantal-divisor valuations, so
the bound ord_0 Delta >= Σ_{i≤N}ν_i survives. The contraction has degree at
most M and the same endpoint. Hence (C) and (R) hold for the auxiliary
profiles exactly as stated.

## 4. Entropy duality and exact certificates

Put a_i=mu_i/M, b_i=nu_i/M and t=M log s>=0. Then

```math
\log B=t+\log\sum_i e^{-ta_i}+\log\sum_i e^{-tb_i}.
```

For probability vectors pi and sigma, Gibbs' inequality gives

```math
\log B\ge t\,(1-\pi\cdot a-\sigma\cdot b)+H(\pi)+H(\sigma). \tag{2}
```

Suppose nonnegative multipliers on the budgets prove π·a+σ·b<=h on the
whole relaxed polytope. The certificate is an LP dual in increment
coordinates, and every budget coefficient is nonnegative. The terms with
mu_1=nu_1=0 also give log B>=t. Combining this with (2),

```math
\log B\ge \frac{H(\pi)+H(\sigma)}{\max(1,h)}\qquad\text{for all }t\ge0.
```

**Tightness.** This part is used only for the sharpness remarks, never for
a certified floor. Let P be the relaxed polytope of normalized profiles. It
is compact, since the determinant budgets bound every order. Let
h(π,σ)=max_P(π·a+σ·b). It is convex, as a maximum of linear functions.

- For fixed t, the right side of (2) is linear in (a,b) and concave in
  (π,σ). Gibbs' inequality is an equality at the Gibbs distributions. With
  both sets compact and convex, Sion's theorem gives
  min_P log B=max_{π,σ}[t(1−h(π,σ))+H(π)+H(σ)].
- The bracket is affine in t>=0 and concave in (π,σ), and the simplices are
  compact. Sion's theorem again swaps inf over t with max over (π,σ). Hence
  the relaxation value is max{H(π)+H(σ) : h(π,σ)<=1}.
- If h>1, mix with the point mass δ_1 on the zero-order channel. Put
  π'=π/h+(1−1/h)δ_1, and the same for σ. Since a_1=b_1=0, this gives
  h(π',σ')<=1. By concavity H(π')+H(σ')>=(H(π)+H(σ))/h. So the best
  bound H/max(1,h) equals the relaxation value.

By the symmetry of P under exchanging the profiles, one may take σ=π.

The [record](../evidence/bare-seven-to-twelve-letters.json) stores, for
each q from 7 to 12:

- one rational probability vector k_i/D with D=10^7, used for both
  profiles;
- the positive multipliers.

The [checker](../tools/check_bare_seven_to_twelve_letters.py) recomputes
every tail inequality of the dual in exact rationals and computes h. It
bounds H and log q with rigorous rational atanh series and remainders.
It then proves H>2 sqrt 2 max(1,h) log q by an exact squared comparison
and derives the decimal floors of the table.

Every certificate uses only full-width boxes, beta=q. At q=7 the positive
multipliers are C(alpha,7) for alpha=1..7 and R(alpha,7) for alpha=2..7.
The values h−1 lie between −4e-6 and 2e-6.

Three corrupted certificates are rejected:

- keeping only the budgets that predate this note;
- halving the largest multiplier;
- moving probability mass to the top order.

## 5. What the budgets do to earlier targets

Earlier work identified the formal seven-letter optimizer

```math
M=12,\qquad \mu=\nu=(0,7^{\times6},13^{\times42}).
```

It satisfies the determinant and shared-minor budgets, but violates 22 box
budgets. The worst violation is C(3,7): 308>252. Every ramification is
also excluded, because the budgets are homogeneous.

Without the adjugate inequality below, the strengthened relaxation at q=7
has optimum about 0.7265564 [NUM]. The certificate floor 0.726556 is just
below this value. The stored q=7 witness satisfies every box budget and both
adjugate inequalities, and its floating-point value is 0.7265675. Its
normalized orders form a staircase:

- one order each at 0, 0.329, 0.491 and 0.594;
- three orders at 0.695;
- five blocks of seven between 0.862 and 1.188;
- six orders at 1.275 and one at 1.322.

Two other refinements were examined.

**Adjugate inequality.** This applies to a square box, alpha=beta>=2, with
N=alpha^2. Let P(z) be the C-box of K, with rows (x,u) and columns (y,v)
in S_1 x S_1.

- Every entry of adj P is an (N−1)-minor of C. So adj P is divisible by
  z^S with S=Σ_{i<N}ν_i, and its entries have degree at most (N−1)M.
- At infinity, (1) with j=1 makes the C-box of H congruent to a Kronecker
  product A⊗B modulo w^{μ_2}. Since
  adj(A⊗B)=(det A)^{α−1}(det B)^{α−1} adj A⊗adj B, its adjugate is again
  congruent to a Kronecker product.
- Realign the adjugate to rows (y,x) and columns (v,u). Modulo w^{μ_2} the
  realigned adjugate of the C-box of H then has rank at most one. Its
  determinant therefore vanishes to order at least (N−1)μ_2.
- At z=1 the realigned adjugate of the permutation box is again a
  permutation matrix. So the realigned determinant D(z) of adj P is a
  nonzero polynomial, divisible by z^{NS}.
- The determinant at infinity equals w^{N(N−1)M}D(1/w).

Together these give

```math
(N-1)\mu_2\le N(N-1)M-\deg D\le N\Big((N-1)M-\sum_{i<N}\nu_i\Big).
\tag{3}
```

The R-box mirror realigns to rows (u,x) and columns (v,y), and gives
(N−1)ν_2<=N((N−1)M−Σ_{i<N}μ_i). Both inequalities are hypothesis free. They
hold with fixed auxiliaries by the same contraction, which needs only the
first layer and Cauchy–Binet for (N−1)-minors.

Under the inverse-degree theorem's hypothesis det K_col=c z^{nM}, one has
Σ_{i<n}ν_i=nM−ν_n. Then (3) at alpha=q reduces to μ_2<=n/(n−1)(ν_n−M),
which sharpens that theorem's μ_2<=2(ν_n−M). **The inverse-degree form
itself still needs the hypothesis.** In the literal controls of Section 7,
65 scalar crossings without a monomial column determinant have
μ_2>2(ν_n−M), and all of them satisfy (3).

At the formal seven-letter optimizer, (3) at alpha=7 reads 336<=49 in both
orientations. So (3) alone also excludes it without the monomial
hypothesis.

A floating-point solve [NUM] shows how much (3) and its mirror raise the
relaxation optimum: about 1.1e-5 at q=7, 5.7e-6 at q=8, 3.3e-6 at q=9,
2.0e-6 at q=10, 1.3e-6 at q=11 and 8.5e-7 at q=12. The certificates do not
use them.

**Low/high coefficient split.** This bound is not homogeneous. The low part
has C-rank at most Σ(P−ν_i)^+, and the high part has R-rank at most
Σ(P'−μ_i)^+. These rank bounds grow with the scale, so the split cannot move
the scale-free relaxation.

## 6. Thirteen letters: the remaining target

At q=13 the same relaxation has optimum about 0.705479<1/sqrt(2) [NUM].
The record stores an exact witness: M=2000 with the profile below, used
for both mu and nu.

| order | 0 | 492 | 728 | 867 | 984 | 1039 | 1146 | 1290 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| multiplicity | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 6 |

| order | 1558 | 1768 | 1903 | 1981 | 2087 | 2092 | 2151 | 2227 | 2303 | 2379 | 2454 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| multiplicity | 13 | 13 | 13 | 13 | 13 | 26 | 13 | 13 | 13 | 13 | 13 |

The checker verifies every box budget and both adjugate inequalities (3)
in integers. The smallest adjugate slack is 5731. At s=414548/413547 the
checker computes B exactly and proves B<13^(2 sqrt 2) with rational
logarithm bounds. Hence the thirteen-letter bare criterion is **not** closed
by Theorem 1 together with (3).

Closing q=13 needs a further proved coupling between the cuts. Otherwise,
it needs a crossing whose profiles approach this staircase:

- a single leading channel in each cut;
- second-order channels spaced as above;
- blocks of q equal orders, four of them below M and the upper ones above
  M.

## 7. Literal controls

The checker also tests the lemmas directly on actual crossings over F_p,
with complete local Smith profiles computed by exact elimination over
truncated power series. The samples come from a platform-independent
pseudorandom generator, so the stored counts reproduce exactly. There are
seven families, with 140 accepted crossings in all:

- scalar q=3 crossings over F_5, with a column-rank-one constant and a
  row-rank-one top coefficient;
- scalar crossings with long two-sided rank-one jets (R-cut jets at
  infinity, C-cut jets at zero), for q=3 over F_3 and q=2 over F_2;
- auxiliary dimension two, for q=2 and q=3 over F_3. The auxiliary indices
  sit on either physical side of each cut, chosen independently, and only
  the boundary contraction is required to equal E;
- unimodular products K_C=E_C∏(zI+(1−z)N_k) with nilpotent N_k, for q=3
  over F_3 and q=2 over F_5. Their column determinant is a monomial, so the
  determinant budget is an equality.

For each crossing the checker verifies (C), (R) and both adjugate
inequalities on the full profiles. On a C-box and an R-box for every
(alpha,beta), 1200 boxes in all, it verifies:

- Lemma 1: C-boxes at infinity against mu, and R-boxes at zero against nu;
- determinantal interlacing at the other place;
- the box determinant budget.

On 280 square 4-by-4 boxes it also checks the steps of (3) and its mirror
directly. These are the nonvanishing of the realigned adjugate determinant
at one, the divisibility at zero, and the Kronecker-layer order at infinity.

The controls have power:

- Lemma 1 holds with equality in 400 nontrivial layers, 93 of them with
  1<alpha<beta.
- 24 budgets are equalities, but all of them are the column determinant
  budget C(q,q) of the unimodular families, which is tight by
  construction. They are not evidence of power on the other budgets.
- 41 crossings carry auxiliary indices.
- 65 scalar crossings violate the inverse-degree form without its
  hypothesis.

No budget with alpha<q, and no R budget, became an equality in these
samples; for 1<alpha<beta the smallest normalized slack observed is 7/24.
In an independent referee's mutation tests (outside the checker), these
controls rejected Lemma 1 with `r=alpha` (on C-boxes, and on R-boxes
alone), the box determinant budget tightened by one, interlacing shifted
by one, and both adjugate steps with `N-1` replaced by `N`. All of this
checks consistency of the stated lemmas and is not part of the proof.

## Reproduction

- `make verify-bare-seven-to-twelve-letters` runs the exact checker. It
  needs only the standard library and takes about three seconds.
- `python3 tools/generate_bare_seven_to_twelve_certificates.py OUT.json [q ...]`
  reruns the floating-point search. It needs numpy and scipy, and plays no
  role in the proof. The search depends on the environment: SLSQP and
  HiGHS results vary with the library build and thread count. A rerun
  therefore generally gives different multipliers and witnesses. In one
  independent rerun all regenerated certificates differed from the stored
  ones and still passed the exact checker after the built-in rational
  repair. Nothing guarantees this in every environment, and only the stored
  certificates are claimed. Runs over all q have taken from about one to
  about five minutes. Run single alphabets under a timeout, with
  OMP_NUM_THREADS=1.
