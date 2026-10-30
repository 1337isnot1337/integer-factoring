# A quantitative contract for further factoring improvements

## Result and scope

This note proves an ordinary mathematical transfer theorem. Any effective
crossing family with fixed rational parameter `rho` in the Appendix A
interface gives factoring exponent `alpha = log2(1/rho)`.

The theorem closes the algorithm and bit-cost connection once that crossing
family exists. It does not construct a family with `rho < 1/sqrt(2)`.
The known half-pieces family gives exponent `1/2`.
The exponent `1/2` is also certified by the [complete Lean endpoint](../formal-appendix/HalfTheorem.lean). The general-parameter transfer in this note remains an ordinary proof.

This work follows the ordinary proof in
[PR #19 at c19044b](https://github.com/1337isnot1337/integer-factoring/tree/c19044bd399e439a86b13d403c4b1c4d28855f01).
That proof supplies the written multiplication-network constructor.
This note generalizes its final composition to rational crossing parameters.
No new Lean proof or compilation accompanies this result.

The sources are pinned:

- Claude, *A 2^o(n) Time Algorithm for 3-SAT*, Appendices A and B,
  [16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
- Rudich and Rousseau, *Subexponential Deterministic Factoring*,
  [c91cc3f](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
- The written constructor and codec in
  [half-top-constructor.md](https://github.com/1337isnot1337/integer-factoring/blob/c19044bd399e439a86b13d403c4b1c4d28855f01/research/half-top-constructor.md).

## 1. The precise hypothesis

Fix integers `a,b` with `b/2 <= a < b`, and set `rho=a/b`.
Assume an effective family satisfies `AppRunningTime.Exist rho c` for
one fixed natural number `c`.
The definition appears in the pinned
[AppRunningTimeStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean).

This hypothesis includes more than matrix rank estimates. For every label
width, it requires the crossing identity, two explicit factorizations,
bounded coefficients, and uniform programs with the stated costs.
After the source's permitted trimming, each factor width is at most
`2*lambda`. Their sum satisfies

```math
c^{\rm row}_i+c^{\rm col}_i
\le 4\rho\lambda+c\log_2(2\lambda).
\tag{1}
```

The number of pieces is polynomial in `lambda`, and the family has a
nonzero common integer scale. The source interface also charges the
programs that enumerate and evaluate the pieces.
An isolated finite tensor decomposition does not establish this hypothesis.

**Transfer theorem.** Under this hypothesis, one deterministic finite-alphabet
machine computes `Nat.minFac` on all natural-number inputs, with

```math
T(n)\le C\exp\!\left(\frac{Cn}{(\log n)^\alpha}\right),
\qquad \alpha=\log_2(b/a),\qquad n\ge2.
\tag{2}
```

Here `n` is binary input length. The constant `C` depends on the fixed
family and its programs. It does not depend on the input or its length.

## 2. Fix one evaluator and one schedule

The assembled
[thmRunningTime_b](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/ThmRunningTime.lean)
provides a constant `C0` and one polynomial-time machine fragment `run`.
Choose a fixed integer `x0 >= 18` with `2^x0 >= c+5`.
Use the source schedule

```math
K_0=2^{x_0},\quad g_d=2^{3d+x_0}K_d,\quad K_{d+1}=g_d^2.
\tag{3}
```

For positive widths `lambda,B`, the fragment evaluates a written
`K_ell`-vertex clique when its input pad has at least `2^need` units.
The sufficient exponent satisfies

```math
\operatorname{need}\le e_\ell K_\ell^2\lambda
 +C_0K_\ell^2\log_2(K_\ell(\lambda+B)),
\qquad e_\ell\le0.85\rho^\ell.
\tag{4}
```

All tables occur explicitly in the input. The fragment takes a fixed
polynomial number of steps in the complete encoded input length.
The source theorem pays for every recursive level and signed intermediate
integer. Its polynomial degree is independent of `ell`.

Choose the largest `ell` with `P^40 <= n`, where `P=K_ell`.
For `n<K_0^40`, use trial division. This finite initial range contributes
only to the fixed constant in (2).
The source schedule gives

```math
P^{40}\le n<P^{400},\qquad
2^\ell\le P,\qquad
\log_2 n<80(3x_0+6)2^\ell.
\tag{5}
```

Successive integer powers and comparisons compute this level.
The next candidate has `O(log n)` bits, by `K_(ell+1)<=P^10`.

## 3. Construct a written divisor-prefix instance

The pinned top-constructor proof gives a deterministic procedure for each
prefix query. It writes a clique `I` with exact integer value

```math
\operatorname{zval}(I)
=\#\{d\mid N:\lfloor d/2^k\rfloor=x\},
\tag{6}
```

and positive `lambda,B` such that

```math
P^2\lambda\le4717n,\qquad PB\le6769n.
\tag{7}
```

Its bit-operation cost is `2^O(n/P)`. This bound includes the exhaustive
grouping search, gate-state enumeration, complete tables, and binary entries.
It also includes orientation and padding rules that preserve (6).
It does not assume free access to a table oracle.

The constructor depends on `P`, not on `rho`. Thus the same construction
serves every fixed family in Section 1.

## 4. An explicit pad for rational rho

The public existential theorem does not state that `need` is computable.
Use the following computable integer instead:

```math
q=\left\lfloor\frac{a^\ell P^2\lambda}{b^\ell}\right\rfloor+1
 +C_0P^2\operatorname{size}(P(\lambda+B)).
\tag{8}
```

Here `size` is binary bit length. Equation (4), `0.85<1`, and
`log2(t)<=size(t)` imply `need<=q`.
The floor plus one also covers an integral quotient safely.
The algorithm writes exactly `2^q` pad units.
It never evaluates a real logarithm or calls the existential function `need`.

The integers `a^ell` and `b^ell` have `O(ell)` bits because `a,b` are fixed.
They are multiplied by integers of polynomial bit length in `n`.
Their computation therefore takes polynomial bit time.

To bound the pad, note that `rho>=1/2` and (5) imply
`rho^ell>=1/P`. Equations (5)--(7) give

```math
\begin{aligned}
P(\lambda+B)&\le9128n<2^{14}n,\\
\operatorname{size}(P(\lambda+B))
 &\le\operatorname{size}(n)+14\le414P,\\
P^2\operatorname{size}(P(\lambda+B))
 &\le414P^3\le414n/P\le414n\rho^\ell.
\end{aligned}
\tag{9}
```

The middle inequality uses `n<P^400`, `size(P)<=P`, and `P>=2`.
The last line uses `P^4<=n`, which follows from the stronger cutoff.
Also `1<=n*rho^ell`. Consequently

```math
q\le(4718+414C_0)n\rho^\ell.
\tag{10}
```

## 5. Complete input, output, and machine cost

The four-symbol source encoding has total input length at most

```math
14+4(2^q+\ell+\lambda+P)
 +(3B+10)P\,2^{(P-1)\lambda}.
\tag{11}
```

The input tables contain signed integers, even though this constructor
produces nonnegative entries. Equation (7) bounds the table exponent
by `4717n/P`. Equations (9)--(10) therefore make (11)
`2^O(n*rho^ell)`.

Writing the pad, encoding the tables, running the fixed fragment, and
reading its output take `2^O(n*rho^ell)` bit operations.
The constructor fits the same bound because `1/P<=rho^ell`.
A fixed polynomial simulation on ordinary finite-alphabet tapes changes
only the constant in this exponent.

Use `n` prefix queries to recover the smallest divisor greater than one.
Subtract the divisor-one contribution before choosing each next bit.
The invariant keeps the least prefix that extends to a nontrivial divisor.
The resulting divisor is prime, since a composite result has a smaller
nontrivial divisor. Handle `N=0,1` separately as `2,1`.

The polynomial factor for prefix search is absorbed by the same bound.
Indeed, `n*rho^ell>=n/P>=n^(39/40)`.
Finally, (5) and `rho=2^(-alpha)` yield

```math
\rho^\ell=(2^\ell)^{-\alpha}
\le\frac{[80(3x_0+6)]^\alpha}{(\log_2n)^\alpha}.
\tag{12}
```

This proves (2) for one fixed machine and all inputs.

## 6. Exact targets for a new crossing construction

The following rows are conditional targets. No new family in these rows
is asserted by this note.

| Required rho | Resulting alpha | Meaning |
| --- | --- | --- |
| `7/10` | `log2(10/7) = 0.514573...` | A small strict advance beyond one half |
| `2/3` | `log2(3/2) = 0.584963...` | A larger crossing improvement |
| `3/5` | `log2(5/3) = 0.736966...` | A substantially smaller rank-width sum |
| `1/2` | `1` | The lower endpoint of this source interface |

For a rational `rho=a/b`, the strict half-exponent test is the exact
integer inequality `2*a*a < b*b`.
For example, `2*7^2=98<100=10^2` certifies the first target.

An advance requires all of Section 1, including uniform coefficient and
decoder costs. A smaller numerical rank at finitely many widths is not
enough. Likewise, changing the partition count without preserving the
operator identities does not establish a new family.

This transfer theorem and the subsequent frontier notes used substantial
OpenAI assistant support. Their evidence is ordinary proof and exact
computation, with the source dependencies identified above.
