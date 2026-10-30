# A sharper ordinary-proof rate for the existing k6 machine

## Result and status

The same deterministic finite-alphabet machine certified in
[PR #17](https://github.com/1337isnot1337/integer-factoring/pull/17) has the
stronger bound

```math
T(n)=\exp\!\left(O\!\left(\frac{n}{(\log n)^{\beta_{12}}}\right)\right),
\qquad
\beta_{12}=\frac{\log(12/11)}{\log12}=0.0350159540186562\ldots.
```

In particular it qualifies for score `1/29`, and for the exact exponent
`beta_12`, under the site's logarithmic-saving convention. This note is an
ordinary mathematical proof. The existing k6 Lean endpoint remains `1/30`;
no new Lean theorem or independent review is claimed. No program is changed.

The source is the patched reference and machine bridge merged at
`059fdf72a6f86cce0b300e04256a52575a0b6e3e`. The reference before that patch is
Isaac Rudich and Louis-Martin Rousseau's
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`. The argument sharpens the
[baseline audit](baseline-audit.md) for the actual degree-twelve schedule;
it does not replace the original counting, compiler, or simulation proofs.

## 1. Schedule and exact depth

The implemented schedule is

```math
K_0=48,\qquad G_d=5K_d^6,\qquad K_{d+1}=25K_d^{12}.
```

Write `a = log 48 + (log 25)/11`. Direct induction gives

```math
\log K_d=a12^d-\frac{\log25}{11}.
```

For all sufficiently large bit lengths `n`, the implemented depth `ell`
is the greatest integer with `K_ell^40 <= n`. Put `P=K_ell` and `G=5P^6`.
Then

```math
P^{40}\le n<25^{40}P^{480},\qquad
P\le n^{1/40},\qquad
\frac nP<25^{1/12}n^{479/480}.
```

The failed next test also gives `log n < 480 a 12^ell`. For `r=11/12`
and `beta_12=log(1/r)/log12`, this implies

```math
r^\ell=(12^\ell)^{-\beta_{12}}
\le (480a)^{\beta_{12}}(\log n)^{-\beta_{12}}.
```

The finite range below the first positive-depth cutoff is absorbed into a
fixed constant. The already certified machine is correct there as well.

## 2. Summable error retains the limiting contraction

At height `d`, let `lambda` be the incoming label width, and let
`N_(d+1)=G_d^4 lambda` and `N_d=K_d^2 Lnext(K_d,G_d,lambda)`.
The actual definitions retained in `FactoringAux/EngineGroup.lean` are

```math
U_{\max}(g,\lambda)
=\left\lfloor\frac{12g^3\lambda+11g^4\lambda+24g^4}{12}\right\rfloor,
\qquad
L_{\rm next}(k,g,\lambda)
=\left\lfloor\frac{kU_{\max}+g^4\lambda}{k^3}\right\rfloor.
```

Dropping floors and using `G_d >= K_d` yields

```math
N_d\le q_d N_{d+1}+2G_d^4,
\qquad q_d=r+\frac2{K_d}.
```

The patched `Kb48_ge` proves `K_d >= 48*2^d`. Consequently

```math
q_d\le\frac{23}{24},\qquad
\frac{q_d}{r}\le1+\frac1{22\,2^d}.
```

For every finite interval of heights, the sum of these relative errors is
at most `1/11`. Expanding the product and bounding the degree-`j`
elementary symmetric sum by the `j`th power of that total gives

```math
\prod_{i=d}^{j-1}q_i\le\frac{11}{10}r^{j-d}.
```

This argument includes the base `K_0=48`. It does not assume the stronger
constant contraction used only above the bottom in the `1/30` proof.

Every recursive `G_d` is at most the top `G`. Applying the product estimate
to the initial term, and the geometric bound `q_d<=23/24` to the additive
errors, gives the following simultaneous bound, where `N_top` is the top
clique volume:

```math
N_d\le\frac{11}{10}r^{\ell-d}N_{\rm top}+48G^4
\qquad(0\le d\le\ell).
```

## 3. All levels have the same principal cost bound

The patched `EngineCostRate.eLev_le` bounds the exponent of a transition
at height `d` by

```math
\frac{14(N_{d+1}+2G^4)}{K_d}+7+(G^4+9)B,
```

where `B=1+floor(log2(W+mb+2G+64))` and one may use
`W=N_top+2ell G^4`. The modulus bound `mb` is the same one used by the
implemented counting engine. Substitution of the volume bound gives

```math
\frac{14(N_{d+1}+2G^4)}{K_d}
\le\frac7{20}r^\ell N_{\rm top}+15G^4.
```

Here the principal coefficient follows from
`K_d>=48*2^d`, `2r>1`, and
`(14*11/10)/(48r)=7/20`. The additive coefficient follows from
`700/48<15`. The bottom enumeration has exponent at most

```math
\frac{11}{10}r^\ell N_{\rm top}+48G^4+2B_0,
\qquad B_0=1+\lfloor\log_2(mb+48+64)\rfloor.
```

These estimates bound the actual `chainMax`, which is the maximum of all
transition exponents and the bottom exponent. They do not infer a full
runtime bound from the bottom alone.

## 4. The other costs and the existing machine bridge

The unchanged source cost definitions and their patched proofs give the
following contributions to `Zu = zT + zC + zP`:

| Contribution | Bound on its exponent |
| --- | --- |
| Top volume | `N_top <= 720 S^4 + 306 G S^3`, with `S^4=O(n)` |
| Top table and grouping search | `O(n/P + G n^(3/4) + G log n + log n)` |
| Recursive transitions and bottom | `O(r^ell N_top + G^4 log n + log n)` |
| Choice paths, modulus, depth, parameter values | `O(ell G^4 log n + G log n + log n)` |

For the logarithms in this table, `W`, `mb`, `P`, `G`, and the next
stopping-test value are bounded by fixed powers of `n` for sufficiently
large `n`. In particular `mb <= cmbn*S^4`, where `cmbn` is a polynomial
in `ell` and `G`; it is the small modulus bound, not the enormous product
of the moduli. Thus their bit lengths are `O(log n)`.

Using `G=5P^6`, `P<=n^(1/40)`, and `ell<=P`,

```math
G n^{3/4}=O(n^{9/10}),\qquad
\ell G^4\log n=O(n^{5/8}\log n).
```

It follows, with one constant uniform in the input, that

```math
Zu(n)=O\left(nr^\ell+n^{479/480}\right)
=O\left(\frac n{(\log n)^{\beta_{12}}}\right).
```

The existing
[`factoring_Zu_exponential`](../variants/k6/FactoringK6/MachineBridge.lean)
already proves that one finite-alphabet machine computes `Nat.minFac` in
time `C*2^(C*(Zu(n)+1))`. Its schedule hypothesis holds by equality for
the actual k6 program. That theorem includes the source's uniform program,
prefix search, fixed-degree compilation, and finite-alphabet simulation.
Substitution proves the claimed new ordinary-proof bound, including the
same all-input correctness and fixed-machine quantifier order.

## 5. Exact comparisons and reproduction

No floating-point value is used in the proof. The comparisons

```math
12^{28}>11^{29},\qquad 12^{27}<11^{28}
```

give `1/29 < beta_12 < 1/28`, hence a strict improvement over `1/30`.
The script below checks those integers, exercises the floor-sensitive
one-step inequality, and checks finite products against the uniform
`11/10` envelope. Those checks supplement the proof; they are not a
verification of all input lengths or a new Lean certificate.

```sh
python3 tools/check_k6_sharp_rate.py
```

## Source map

- [Reference patch](../variants/k6/reference.patch): actual degree-twelve
  schedule, top construction, and preserved cost definitions.
- [Machine bridge](../variants/k6/FactoringK6/MachineBridge.lean):
  `factoring_Zu_exponential` and its complete machine-cost construction.
- Reference `EngineGroup.lean`: `Umax`, `Lnext`, `Kb`, `Gs`.
- Reference `EngineCostRate.lean`: `Ntop_S`, `eTop_le`, `eLev_le`.
- Reference `EngineCostFinal.lean`: `Wn`, `mbn_le`, `NPof_le` dependencies.
- Reference `UniExpChain.lean`: `Kb48_ge`, `chainMax`, runtime depth.
- Reference `UniExp.lean`: `zT`, `zC`, `zP`, `Zof`, `Zu`.

The original reduction and counting infrastructure are due to Isaac Rudich
and Louis-Martin Rousseau. This note concerns a sharper analysis of their
infrastructure with this repository's already certified k6 schedule.
