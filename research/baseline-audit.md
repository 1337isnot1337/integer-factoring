# Quantitative audit of the deterministic factoring baseline

Date: 2026-10-29 (America/Chicago). Audited source commit:
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).

This note derives rates from the source's explicit inequalities. The additional
derivations in this note have not been checked by Lean. They must not be described
as a new Lean-verified factoring theorem. The baseline's build status also needs
to be recorded separately from this source audit.

## 1. What the source actually states

The advertised theorem is

\[
\exists M\;\forall k\ge1\;\exists C_k\;\forall n,
\qquad T_M(n)\le C_k2^{\lfloor n/k\rfloor}.
\]

Here `n` is the length of the binary encoding of the input integer, the output is
`Nat.minFac`, and the machine has finite stack alphabets. This is a uniform
`2^{o(n)}` bound. It does not specify a quantitative function tending to zero.

The separate `SqrtLogLogConjecture` asks for

\[
T(n)=\exp\!\left(O\!\left(\frac{n}{\sqrt{\log n}}\right)\right).
\]

That proposition is explicitly unproved in the source. It is a stronger target
than either quantitative rate derived below.

Primary source locations:

- [`Challenge.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/Challenge.lean): audited theorem definitions and conjecture.
- [`UniMain.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniMain.lean): final assembly.
- [`UniParams.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniParams.lean): runtime depth.

## 2. Exact size recurrence

Set

\[
K_0=48,\qquad G_d=3K_d^8,\qquad K_{d+1}=G_d^2=9K_d^{16}.
\]

The runtime chooses the greatest depth `ell` satisfying `K_ell^40 <= n`, except
that it returns zero for smaller inputs. This exceptional finite range does not
affect an asymptotic upper bound.

Writing `x_d = log K_d`, `a = log 9 / 15`, and `c = log 48 + a`, the recurrence
solves exactly as

\[
x_d=c16^d-a.
\]

Consequently

\[
\ell(n)=\log_{16}\log n+O(1),\qquad
P:=K_{\ell(n)}\le n^{1/40},\qquad
P\ge 9^{-1/16}n^{1/640}.
\]

The last inequality follows directly from `n < (9 P^16)^40`. It is important to
retain this inequality: `P <= n^(1/40)` alone cannot upper-bound the cost `n/P`.

## 3. A sharper contraction already contained in the definitions

The source's uniform asymptotic proof uses `23/24` as a convenient contraction.
The actual `Umax` and `Lnext` definitions give a better one.

For a transition at height `d`, let `lambda` be the old label width, and put

\[
N_{d+1}=G_d^4\lambda,\qquad
N_d=K_d^2L_{\rm next}(K_d,G_d,\lambda).
\]

By the two floor operations in those definitions,

\[
\begin{aligned}
N_d
&\le U_{\max}(G_d,\lambda)+\frac{G_d^4\lambda}{K_d}\\
&\le\left(\frac{11}{12}+\frac1{G_d}+\frac1{K_d}\right)N_{d+1}+2G_d^4\\
&\le\left(\frac{11}{12}+\frac2{K_d}\right)N_{d+1}+2G_d^4.
\end{aligned}
\]

The last step uses `G_d >= K_d`. Define `r = 11/12` and
`q_d = r + 2/K_d`. The source proves `K_d >= 48 * 2^d`, so

\[
q_d\le r\left(1+\frac1{22\,2^d}\right),\qquad q_d\le\frac{23}{24}.
\]

For any finite interval of heights,

\[
\prod_{i=d}^{j-1}q_i\le\frac{11}{10}r^{j-d}.
\]

To check the constant, the sum of the factors' relative errors is at most
`sum_{i>=0} 1/(22*2^i) = 1/11`. For nonnegative `u_i` with sum `s<1`, expanding
the product and bounding every elementary symmetric polynomial by `s^k` gives
`product(1+u_i) <= 1/(1-s)`. Here this is at most `11/10`.

Let `G = G_ell`, and let `N_top` bound the bits at the initial clique. Unrolling
the recurrence, using the sharper product on the initial term and `q_i <=23/24`
on the additive errors, gives at every height `d`

\[
\boxed{\quad
N_d\le \frac{11}{10}r^{\ell-d}N_{\rm top}+48G^4.
\quad}
\]

This is an explicit inequality with a constant independent of the depth.

Primary source locations:

- [`EngineGroup.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean): `Umax`, `Lnext`, `Kb`, `Gs`.
- [`UniExpChain.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniExpChain.lean): `Kb48_ge`, `chainMax`.
- [`EngineCostRate.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostRate.lean): `N_step`, `eLev_le`, `ePolyL_mono`.

## 4. Controlling every level, not just the bottom

A small bottom brute-force term is insufficient on its own: a higher-level table
could dominate the running time. The source's level-cost lemma charges the
transition at height `d` by

\[
e_d\le\frac{14(N_{d+1}+2G^4)}{K_d}+7+(G^4+9)B,
\]

where

\[
B=1+\log_2(W+m_b+2G+64),\quad
W=N_{\rm top}+2\ell G^4.
\]

Floors can only decrease the displayed upper bound. Substitution of the preceding
inequality, and `K_d >=48*2^d`, gives

\[
\begin{aligned}
\frac{14(N_{d+1}+2G^4)}{K_d}
&\le\frac{77}{5K_d}r^{\ell-d-1}N_{\rm top}+\frac{700G^4}{K_d}\\
&\le\frac7{20}r^\ell N_{\rm top}+15G^4.
\end{aligned}
\]

For the main term, `2r>1` implies `(2r)^(-d) <=1`, and the constant simplifies
exactly: `(77/240)/r =7/20`.

The bottom cost is at most

\[
\frac{11}{10}r^\ell N_{\rm top}+48G^4+2B_0,
\quad B_0=1+\log_2(m_b+48+64).
\]

Thus one simultaneous upper bound on all the level exponents and the bottom is

\[
\boxed{\quad
\frac{11}{10}r^\ell N_{\rm top}+48G^4+(G^4+9)B+2B_0+7.
\quad}
\]

This avoids the original `small_zC` proof's depth cutoff `D=80696*k`, which hides
most of the quantitative information.

## 5. Other exponent terms and the machine cost

The source's top clique bound, with `S=2^(sOf n)` and `S^4=O(n)`, is

\[
N_{\rm top}\le720S^4+306G S^3.
\]

The top table's exponent is

\[
O(n/P+G n^{3/4}+G\log n).
\]

The number of paths contributes exponent `O(ell G^4 log n)`. The modulus bound
is polynomial in `n`, `P`, and `ell`, hence its logarithm is `O(log n)`. The
runtime depth search, the next `K` value, and the level-size bookkeeping also
have logarithms `O(log n)` or smaller than the displayed path bound.

In the unchanged baseline `G=3P^8` and `P<=n^(1/40)`, so all the requirements can
be bounded by `V=2^Z` with

\[
\begin{aligned}
Z&=O\bigl(nr^\ell+n/P+P^8n^{3/4}+\ell P^{32}\log n+\log n\bigr)\\
 &=O\bigl(nr^\ell+n^{639/640}\bigr).
\end{aligned}
\]

The source proves that the uniform counting command costs at most
`c*(V+1)^2408`, with fixed `c` and fixed degree independent of `ell`. The
prefix search adds a polynomial factor in `n`. Compilation then incurs a square,
followed by a cubic simulation bound. These fixed polynomial powers multiply
the exponent by constants; they do not change the exponent of the logarithm.

Therefore the unchanged algorithm, under the source's computational model and
its proved cost bridges, admits the handwritten quantitative bound

\[
\boxed{\quad
T(n)=\exp\!\left(O\!\left(
\frac{n}{(\log n)^{\beta_{16}}}\right)\right),\qquad
\beta_{16}=\frac{\log(12/11)}{\log16}.
\quad}
\]

The deduction uses `r^ell = O((log n)^(-beta_16))`, and
`n^(639/640)=o(n/(log n)^beta_16)`.

This strengthens the advertised *description* of the baseline's running time.
It does not create a different factoring algorithm or establish the square-root
logarithm conjecture.

Primary source locations:

- [`EngineCostFinal.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostFinal.lean): `Ntn`, `Wn`, `mbn_le`, `NPof` bounds.
- [`UniExp.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniExp.lean): top and bookkeeping exponents, `small_zC`.
- [`UniCost.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniCost.lean): `uniC_cost`, degree `2408`.
- [`UniFinal.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniFinal.lean): counting-to-factoring and simulation bridge.

## 6. Propagation for a replacement grouping lemma

Suppose a replacement grouping theorem is established with the same bounds on
group sizes and inter-group label width, but requires only

\[
G(K)=12K^4(1+\lfloor\log_2K\rfloor).
\]

The concrete integer choice avoids any need for exact real logarithms in the
algorithm. It is at least `12K^4 log K` for `K>=2`. Assume the grouping can be
found by the same exhaustive search over `K^G` colorings and all remaining
operations have the same fixed polynomial cost accounting.

The group-count recurrence becomes

\[
K_{d+1}=144K_d^8(1+\lfloor\log_2 K_d\rfloor)^2.
\]

The contraction proof in sections 3 and 4 continues to hold: it uses only
`G_d>=K_d`, `K_d>=48*2^d`, and the same label-width inequality.

### Why the logarithm in the grouping size does not lose a logarithmic factor

Let `x_d=log K_d`. Then

\[
x_{d+1}=8x_d+O(\log x_d+1).
\]

More precisely, `x_{d+1}>=8x_d`. Since `K>=48` implies
`1+floor(log_2 K)<=K` and `12<=K`, we also have `G(K)<=K^6`, so
`x_{d+1}<=12x_d`. The latter crude bound gives `log x_d=O(d+1)`. Substituting
that into the first recurrence and dividing by `8^(d+1)` yields

\[
\frac{x_{d+1}}{8^{d+1}}
\le\frac{x_d}{8^d}+O\!\left(\frac{d+1}{8^{d+1}}\right).
\]

The error series converges. Together with the lower bound this proves
`x_d=Theta(8^d)`, and hence

\[
\ell(n)=\log_8\log n+O(1).
\]

There is no extra `log log n` factor in this depth estimate.

### Cost propagation

Retain the same stopping rule `K_ell^40<=n`. Then `P<=n^(1/40)` and, using
`K_{ell+1}<=P^12`, the failure of the next stopping test implies
`n<P^480`, hence `n/P<n^(479/480)`.

The coarse bound `G<=P^6` suffices for all the lower-order terms:

\[
G n^{3/4}\le n^{9/10},\qquad
\ell G^4\log n=O(n^{5/8}\log n).
\]

The second inequality uses `ell<=P`, which follows already from
`K_ell>=48*2^ell`. The search over `K^G` groupings contributes logarithmic cost
`G log K`, which is smaller still. Thus

\[
Z=O\bigl(nr^\ell+n^{479/480}\bigr),
\]

and the corresponding quantitative factoring bound is

\[
\boxed{\quad
T(n)=\exp\!\left(O\!\left(
\frac{n}{(\log n)^{\beta_8}}\right)\right),\qquad
\beta_8=\frac{\log(12/11)}{\log8}=\frac43\beta_{16}.
\quad}
\]

This section is conditional on the replacement grouping theorem and its
algorithmic use satisfying the stated hypotheses. It is a complete propagation
argument, not a proof of that grouping theorem. Changing `16` to `8` raises the
logarithmic saving exponent by `4/3`, not by `2`.

## 7. Useful exact certificates for later formalization

For the baseline, two convenient natural-number inequalities are:

1. `14 * 24^(d+1) <= K_d * 23^(d+1)`, following from `K_d>=48*2^d`.
   Together with the source's `23^(d+1)N<=24^(d+1)W0`, this gives the simpler
   conservative estimate `14N/K_d<=W0` at every level.
2. At every height `d>=1`, `K_d>=6144`, so
   `11/12+2/K_d<=939/1024`. The exact integer certificate
   `16*939^32 < 1024^32` yields a rational target saving exponent `1/32`
   (the first height contributes only a constant).

The latter target is slightly weaker than `beta_16`, but avoids real exponent
identities in a formal statement: a denominator given by five repeated integer
square roots of `log_2 n` has order `(log n)^(1/32)`. Such a statement still
requires a checked bridge to the original factoring machine; checking the
   numerical certificate alone does not verify a factoring theorem.

## 8. Cost audit for an independent grouping-accuracy parameter

A stronger grouping proposal separates the accuracy parameter `D` from the
number `K` of groups. The proposed replacement lemma has hypotheses

\[
D\ge2,\qquad G\ge12D^2K^2\log K,
\]

and conclusions, for a network with total edge weight `U` and weighted vertex
degree at most `P_local`,

\[
|\chi^{-1}(\gamma)|\le2G/K,\qquad
K^2Y_{\gamma,\gamma'}\le U+GP_{\rm local}/D.
\]

Its proof is separate from this audit. For that lemma, the following is a
concrete deterministic schedule:

\[
\begin{gathered}
K_0=48,\quad D_d=(d+48)^2,\quad b_d=1+\lfloor\log_2K_d\rfloor,\\
G_d=12D_d^2K_d^2b_d,\qquad K_{d+1}=G_d^2.
\end{gathered}
\]

The next label width is the integer

\[
L_{\rm next}=\left\lfloor
\frac{D_d U_{\max}(G_d,\lambda)+G_d^4\lambda}{D_dK_d^2}
\right\rfloor.
\]

It satisfies

\[
N_d\le\left(\frac{11}{12}+\frac1{G_d}+\frac1{D_d}\right)N_{d+1}+2G_d^4.
\]

### Contraction constants remain bounded independently of depth

The schedule has `G_d>=K_d` and `K_d>=48*2^d`. Therefore

\[
\sum_{d\ge0}\left(\frac1{G_d}+\frac1{D_d}\right)
\le\frac1{24}+\frac1{47}.
\]

The second bound follows by comparison with the integral of `x^(-2)` from
`47` to infinity. Dividing by `r=11/12` gives a relative error sum at most
`71/1034 < 1/11`. The same product bound `product q_i <=(11/10)r^length`
therefore applies. Each individual `q_d` is also at most `23/24`, since
`1/G_d+1/D_d<=1/48+1/2304<1/24`.

The all-level bounds from sections 3 and 4 survive with the same convenient
constants.

### Exhaustive search is affordable in the claimed asymptotic model

There are `K^G` colorings of `G` vertices by `K` colors. The enumeration is
`K^G`, not `G^K`.

The elementary bounds

\[
K_d\ge48\,2^d\ge d+48,\quad D_d\le K_d^2,\quad b_d\le K_d,
\]

give

\[
G_d=12D_d^2K_d^2b_d\le K_d^8.
\]

At the largest retained depth, write `P=K_ell` and `G=G_ell`. With the original
stopping rule `P^40<=n`,

\[
\log_2(K_d^{G_d})\le G\log_2P\le n^{1/5}\log_2 n.
\]

Thus the full exhaustive search over groupings is smaller than the claimed
`n/(log n)^beta` exponent. A verification pass over a candidate grouping uses
only polynomially many arithmetic operations in the network size, weight bit
length, `K`, and `G`. Even doing such verification for every coloring and every
visited counting path preserves the bound: the exponential contributions add
in the logarithm of the work.

### Other large values and nested loops

The main possible hidden costs have the following bounds.

| Item | Bound relevant to the logarithm of the work |
| --- | --- |
| Exhaustive grouping enumeration | `G log K <= P^8 log P` |
| Choice count at one level | `G^4 log(lambda+2)` |
| All counting paths | `O(ell G^4 log n)` |
| Top clique internal work | `O(n/P + G n^(3/4) + G log n)` |
| Divisor-count modulus bound | Polynomial in `n`, `P`, and `ell` |
| Next stopping-test value | `K_(ell+1)^40 <= P^640` |
| Label-width arithmetic | Fixed-degree polynomials in `D`, `K`, `G`, and `lambda` |

For the path count, the conservative old bit bound
`W=N_top+2ell G^4` remains valid because `q_d<=23/24`. The multiplier exponent
is at most `2ell G^4(W+1)`, as in the source. Therefore the modulus bound remains
polynomial, even though the multiplier itself need not be represented as an
ordinary integer: the algorithm evaluates it modulo the small moduli.

The runtime can compute `b_d` by repeated division by two in a bounded loop, and
compute `D_d` by additions and multiplication. A slower implementation with at
most `K_d` halving iterations is already polynomial in `P`. It therefore does
not need a unit-cost logarithm instruction. All these operations can be written
with a fixed number of nested bounded loops and a fixed number of registers.
The program must keep the uniform source's runtime loop over levels; generating
a syntactically deeper program separately for every `ell` would not justify a
fixed machine bound.

The exact source proof with polynomial degree `2408` does not apply verbatim
after these changes. For example, the label-width denominator changes from
`K^3` to `DK^2`, and the numerator acquires a factor `D`; the corresponding
value-bound lemmas need revision. These are fixed polynomial changes because
`D<=K^2`. They do not force a polynomial degree depending on the input length or
on the number of levels. The source's verified compiler and simulation apply to
a fixed modified source program once its cost/correctness lemmas are supplied.

### Depth and resulting conditional rate

For `x_d=log K_d`,

\[
x_{d+1}=4x_d+8\log(d+48)+2\log b_d+\log144.
\]

The lower bound is `x_(d+1)>=4x_d`. The coarse estimate `G_d<=K_d^8` gives
`x_(d+1)<=16x_d`, so `log x_d=O(d+1)`. Dividing the exact recurrence by
`4^(d+1)` then shows that its accumulated error is bounded by a convergent
series `sum O(d+1)/4^(d+1)`. Hence

\[
x_d=\Theta(4^d),\qquad \ell(n)=\log_4\log n+O(1).
\]

As `n<K_(ell+1)^40<=P^640`, we have `n/P<=n^(639/640)`. Also
`G n^(3/4)<=n^(19/20)` and
`ell G^4 log n<=n^(33/40)log n`. Every remaining contribution is therefore
smaller than the main contraction term, giving

\[
Z=O\bigl(n(11/12)^\ell+n^{639/640}\bigr).
\]

Subject to the stated replacement grouping theorem and its insertion in the
uniform counting program, the rate is

\[
\boxed{\quad
T(n)=\exp\!\left(O\!\left(
\frac{n}{(\log n)^{\beta_4}}\right)\right),\qquad
\beta_4=\frac{\log(12/11)}{\log4}=2\beta_{16}
\approx0.06276544104192941.
\quad}
\]

This is an asymptotic cost audit. The grouping theorem and the changed
mathematical algorithm require their own proofs, and a formal proof of the
changed machine requires the additional cost and correctness bridges described
above.

## 9. A limitation of the arbitrary-network grouping interface

The exponent `2` in the new grouping size cannot be reduced below `2` by a
theorem with the same guarantees for every weighted network.

Fix `K>2`, `D>=2`, and a cycle on `g>=3` vertices, with every edge of weight
one. This network has `P=2` and ordered total `U=2g`. If

\[
2g(1+1/D)<K^2,
\]

the pair-width guarantee would require every integer `Y_ab` to be less than
one, hence zero. Every cycle edge would then have both endpoints in the same
color class. Connectivity forces every vertex into one class, contradicting
`K s_a<=2g` because `K>2`.

Therefore any universal sufficient node threshold for this interface must be
at least of order `K^2`, even for simple unweighted graphs. The explicit
Chernoff construction is within the additional accuracy and logarithmic
factors of this necessary growth exponent.

With the retained relation `K_(d+1)=G_d^2`, the limiting growth base is
therefore at least `4` for a proof that treats all networks through this same
grouping lemma. This is a limitation of that interface, not a lower bound for
factoring. A smaller growth base could still be possible by using additional
properties of the particular networks produced by the factoring construction.
