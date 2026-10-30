# What must change to pass the exponent one half

## Scope

This note gives three restrictions on direct modifications of the Appendix B
method. It proves no lower bound on factoring or on an evaluator's runtime.
The restrictions concern its stated coefficient recurrence, its grading,
and its explicit operator identities.

The evaluator source is `anthropics/formal-math` at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The [factoring composition](half-rate-factoring.md) uses its exponent
`1/2`, which is also [Lean-certified](../formal-appendix/HalfTheorem.lean).
This note does not change that bound.

The arguments below are ordinary proofs. The accompanying Python program
checks finite examples of the operator obstruction with exact integers.
No new Lean formalization accompanies this work.

## 1. The coefficient recurrence has no further logarithmic saving

The source defines

```math
K_0=2^{x_0},\quad
g_d=2^{3d+x_0}K_d,\quad K_{d+1}=g_d^2,\quad
e_0=\tfrac12,\quad
e_{d+1}=(\rho+3\zeta_d)e_d+3/K_d,
\qquad \zeta_d=2^{-d-5}.
```

Its proved estimates include both sides:

```math
\tfrac12\rho^d\le e_d\le0.85\rho^d,
\qquad
\log_2K_d=(3x_0+6)2^d-6d-2x_0-6.
```

See [LevelsStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/LevelsStatement.lean)
and the assembled arithmetic in
[Num/PartA.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Num/PartA.lean).
For `rho=1/sqrt(2)`, these estimates give

```math
e_d=\Theta\!\left((\log K_d)^{-1/2}\right).
```

Changing the constants in the upper estimate cannot produce a stronger
power of the logarithm for this sequence.

There is also a schedule-independent restriction within the same square-grid
layout. Suppose a schedule has

```math
K_0\ge2,\qquad g_d\ge K_d,\qquad K_{d+1}=g_d^2,
\qquad e_0\ge c_0>0,\qquad e_{d+1}\ge\rho e_d.
```

For any fixed `a>0`, select a level with `K_ell<=n^a`. Induction gives

```math
K_\ell\ge K_0^{2^\ell}\ge2^{2^\ell},
\qquad 2^\ell\le a\log_2n.
```

Writing `beta=log2(1/rho)>0`, we obtain

```math
e_\ell\ge c_0\rho^\ell
=c_0(2^\ell)^{-\beta}
\ge c_0(a\log_2n)^{-\beta}.
```

Thus this coefficient cannot be `O((log n)^(-beta'))` for `beta'>beta`.
At `rho=1/sqrt(2)`, this excludes every `beta'>1/2` in this template.
It still applies when the initial group count depends on `n` but remains
at least two, provided the positive lower bound for `e_0` remains uniform.

This is a restriction on that upper-bound template. An upper bound that
contains `e_ell` does not imply that a machine requires that much time.
A different analysis can use structure that this scalar recurrence omits.
New pieces, different intermediate spaces, or a different network layout
can also change its assumptions.

## 2. Nonlinear mode weights do not preserve the operator grading

One possible route to a better coefficient is slower growth of the spaces
through which the two factorizations pass. The source uses partition spaces:
`dim V_E=p(E)` and `log p(E)=Theta(sqrt(E))`. Its paper observes that growth
like `E^gamma` would give `rho=2^(gamma-1)` in the same dimension calculation.

That observation needs compatible operators and identities. Replacing
ordinary partitions by partitions into squares does not supply them.

Here is an exact restriction. Give each variable `t_k` a positive weight
`w(k)`, with weights additive on monomials. Keep the source's operators
`a[r]`. Suppose that, for each positive `r`, every operator `a[r]` has one
weight shift `s(r)` independent of `a`. This is the homogeneity property
used to bound the dimensions of the intermediate spaces.

Apply two operators to the vacuum `1`. The definitions give

```math
t_1[r]1=t_r,
\qquad (t_1^r)[r]1=t_1^r.
```

For the second equality, only creation operators can survive on the vacuum.
Their positive mode numbers sum to `r`, and there are `r` of them.
Each mode number therefore equals one. Its coefficient is one.

The two outputs must have the same weight shift. Hence

```math
w(r)=s(r)=r\,w(1)\qquad(r\ge1).
```

All compatible positive gradings are therefore scalar multiples of the
original grading. Nonlinear weights such as `w(k)=k^2` violate this
property. For example, at `r=2` the two output weights are four and two.

The source definitions and its homogeneity statement `O1.maps` are in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The proof files
[Ring/Op.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/Op.lean)
and [Ring/OpFacts.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/OpFacts.lean)
assemble the operator properties.

This argument concerns the existing operators on the full polynomial ring.
It does not exclude a different operator family or a new invariant subspace.
Either route needs a replacement homogeneity and crossing proof.

## 3. A finite mode cutoff breaks the required identities

Another possible route retains only finitely many variables. Define truncated
operators by setting both creation and differentiation in omitted modes to
zero. This changes the identities before any runtime estimate enters.

Use the source's notation

```math
\nabla_{-k}f=t_kf,\qquad
\nabla_kf=k\,\partial f/\partial t_k\quad(k\ge1),\qquad \nabla_0=0.
```

Let `Q_s(a,b)` be its signed sum, and let `Q'_s(a,b)` reverse the two
operators in every summand. The source requires `Q_s=Q'_s`.

**Proposition.** Fix any integer `h>=1`. Keep only modes `1,...,h`.
For `a=t_1`, `b=t_h`, and `s=h`, the truncated operators satisfy

```math
Q_h(a,b)-Q'_h(a,b)
=(h+1)\binom{2h}{h-1}\,\mathrm{Id}.
```

These inputs meet the original lemma's conditions: `a` has weight at most
`h`, and `b` has weight `h`. The displayed coefficient is positive.

**Proof.** For a linear field,

```math
t_j[r]=\binom{r-1}{j-1}\nabla_{-r}.
```

At `s=h`, the `i`th summand uses modes `i` and `-i`. The factor from
`a` is `binom(h+1,i)`, which vanishes for `i>h+1`. Before truncation,
the order-swap difference is the scalar

```math
\sum_{i=0}^{h+1}(-1)^i\binom{h+1}{i}
\binom{-i-1}{h-1}(-i)=0.
```

The summand after the binomial sign is a polynomial in `i` of degree `h`.
Its finite difference of order `h+1` vanishes. This also follows directly
from the commutator `nabla_-i nabla_i-nabla_i nabla_-i=-i Id`.

The cutoff deletes only the last nonzero summand. Its value is

```math
(-1)^{h+1}\binom{-h-2}{h-1}(-h-1)
=-(h+1)\binom{2h}{h-1}.
```

The retained sum is its negative, as asserted. Adding mode `h+1` repairs
this particular identity. It does not establish all other source identities.

The failure also reaches the actual crossing construction. Take `h=1`
and its unique partition `(1)`. Then `a=b=t_1`, and the tested state is
`t_1`. For `s=-2,-1,0,1,2`, the truncated `Q_s` coefficients are

| Mode cutoff | Five coefficients | Their sum |
| --- | --- | --- |
| Only mode 1 | `2, -4, 3, -2, 1` | `0` |
| All modes that occur | `1, -2, 3, -2, 1` | `1` |

The required crossing coefficient is `h!=1`. The finite cutoff gives zero.
This is a negative control for the algebraic family, not a purported
positive-bit-length factoring instance.

## 4. Reproduce the finite checks

Run:

```sh
python3 tools/check_boson_truncation.py
```

The checker uses sparse integer polynomials and literal multiplication and
differentiation. It does not evaluate the closed-form defect to construct
the operator output.

It verifies 96 cases of the proposition and 96 cases with the omitted
mode restored. These cover `h=1,...,32` and three input polynomials.
It also checks the source's actual families for one-part partitions at
`h=1,...,12`. For each, the full operator sum gives `h!`; the cutoff at
`h` gives a different answer. Those additional cases are finite evidence.
The all-`h` proposition above has a separate proof.

The recorded output is
[boson-truncation.json](../evidence/boson-truncation.json).
These checks do not prove an asymptotic factoring bound. Their purpose is
to detect invalid operator substitutions before those substitutions become
claimed rate improvements.

## Consequence for further research

The useful next target is a replacement algebraic construction or a smaller
factorization of the existing pieces. A slower dimension growth alone is
insufficient. The construction must also preserve the crossing sum, both
factorizations, a polynomial number of pieces, and affordable exact factors.
This note rules out two direct substitutions while leaving those targets open.
