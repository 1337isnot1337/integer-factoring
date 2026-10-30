# Modular crossing pieces with an invertible scale

## Result and scope

The Appendix B operator construction has a useful version over every prime
field, including characteristics smaller than the partition weight. Its
crossing multiplier can be made exactly one. A stable ideal removes the modes
whose indices are divisible by the characteristic and bounds every remaining
multiplicity. The resulting factors therefore use fewer intermediate channels.

This is a positive construction, with an ordinary proof and exact finite
checks. In characteristic two, Section 7 proves the uniform width coefficient
`rho=1/sqrt(2)`, so the construction retains the full half-rate while using
only one fixed characteristic and an invertible scale. It does **not**
establish a coefficient below `1/sqrt(2)` or a factoring exponent above
`1/2`. Finite modular rank drops do not establish an asymptotic rank bound.
Section 6 gives the ring-arithmetic transfer, including its dependence on
the coefficient algebra's representation length.

The polynomial ring, normally ordered operators, dual polynomials, and
crossing identity are from Claude (Anthropic), *A 2^o(n) Time Algorithm for
3-SAT*, at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The source definitions are in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
No new Lean formalization is claimed.

## 1. Remove the factorial without dividing by zero

For a partition `mu`, put

```math
z_\mu=\prod_{k\ge1} k^{m_k(\mu)}m_k(\mu)!.
```

The source's integer dual polynomial is

```math
a^\mu=\prod_{n\in\mu}(-1)^{n-1}
                 \left(\sum_{j=1}^n\binom njt_j\right).
```

Its stated duality is
`a^mu[-h] ((h!/z_nu)t^nu)=h! [mu=nu]`. Over the rationals, divide
by the nonzero integer `h!/z_nu`. The resulting identity is

```math
a^\mu[-h]t^\nu=z_\nu[\mu=\nu],\qquad |\mu|=|\nu|=h.
\tag{1}
```

Every coefficient on both sides of (1) is an integer. Consequently (1)
can be reduced modulo any prime `p`, even when `p<=h`. This derivation
does not divide by `h!` inside the prime field.

Call a partition *p-regular here* when no part is divisible by `p` and
every part has multiplicity less than `p`. This local definition includes
both conditions; it is stronger than some other uses of that term.
Exactly these partitions have `p` not dividing `z_mu`.
For a regular label set, define

```math
b_p^\nu=z_\nu^{-1}t^\nu\quad\text{over }\mathbb F_p.
```

Equation (1) becomes

```math
a^\mu[-h]b_p^\nu=[\mu=\nu].
\tag{2}
```

The source's operator identities have integer coefficients on each fixed
finite input. Reduce them modulo `p` and use (2). For any set of regular
partitions of weight `h`, the original formula

```math
A_s(x,y,u,v)
=[t^{\mu(u)}]Q_s(a^{\mu(x)},b_p^{\mu(v)})t^{\mu(y)}
```

therefore satisfies

```math
\sum_{s=-2h}^{2h}A_s(x,y,u,v)=[x=v][y=u].
\tag{3}
```

The scale in (3) is one in every characteristic.

## 2. A stable quotient of the state space

Work over `F_p` and set

```math
J_p=(t_{kp}:k\ge1)+(t_k^p:k\ge1),\qquad
\overline V_p=\mathbb F_p[t_1,t_2,\ldots]/J_p.
\tag{4}
```

The regular monomials form a basis of this quotient. It is an infinite
graded algebra; every weight space is finite dimensional. It is not a
finite cutoff in the mode index.

Every source generator `nabla_r` preserves `J_p`. Multiplication preserves
any ideal. For positive modes, `nabla_k=k partial/partial t_k` kills
`t_{jp}` when `k=jp`, since its coefficient is zero modulo `p`; it also
kills `t_j^p` when differentiating that factor. The derivation rule proves
stability on every multiple of each generator of the ideal.

Every normally ordered source operator consequently descends to the
quotient. An output coordinate indexed by a regular partition vanishes on
`J_p`. Hence projecting to the quotient **between** the two operators
does preserve each matrix coefficient in (3). This is precisely the
invariance property that failed for the characteristic-zero finite cutoff.

Let `d_p(E)` count regular partitions of `E`, and set it to zero for
negative `E`. The two source factorizations now have channel counts

```math
R^{\rm row}_{p,h,s}=\sum_{i=0}^{2h}d_p(s+i),\qquad
R^{\rm col}_{p,h,s}=\sum_{i=0}^{2h}d_p(2h-s-i).
\tag{5}
```

They are exact factorizations over `F_p`, and thus over every commutative
`F_p`-algebra. In particular, zero divisors in such an algebra cause no
problem: the only inverses used to define the factors are the nonzero
prime-field constants `z_nu`.

For the selected binary domains both channel counts are positive. If
the lower end of an interval in (5) is nonpositive, that interval contains
weight zero. Otherwise it contains weight `2h`. That weight has a regular
partition: use the single part `2h` when `p` does not divide it, and
`(2h-1,1)` otherwise. The exceptional pair `p=2,h=1` cannot encode even
one bit, so never arises. Thus the integer logarithms used below are defined.

## 3. Reducing the operator labels is valid as well

An optimized implementation also reduces `a^mu` modulo `J_p` before
evaluating its operators. This requires a second check: invariance of
the state ideal alone would not justify changing the operator label.

For an operator factor `t_{kp}`, every creation coefficient is
`binom(l-1,kp-1)`. This is zero modulo `p` unless `p` divides `l`, by
the least-digit case of the binomial congruence. In the latter case the
created variable `t_l` is zero in the quotient. An annihilation
coefficient is

```math
(-1)^{kp-1}\binom{l+kp-1}{kp-1}l.
```

It vanishes when `p` divides `l` because of the last factor, and otherwise
its binomial coefficient vanishes by the same least-digit calculation.
Thus any normally ordered monomial containing an operator factor
`t_{kp}` acts as zero on the quotient.

If an operator label contains `p` identical factors `t_n`, group their
mode assignments into orbits under a cyclic permutation of those `p`
positions. Every nonconstant orbit has size `p`; its normally ordered
terms have equal coefficients and cancel. Constant assignments have
either `p` identical creators or `p` identical annihilators. The former
multiply by `t_l^p`, which is zero in the quotient. The latter contain
the `p`th derivative, which is zero in characteristic `p`. The other
operator factors do not affect this argument because normal ordering
places all derivatives before all multiplications.

By linearity, any operator label in `J_p` acts as zero on the quotient.
Both the state and label reductions used by the checker are therefore
valid.

## 4. Uniform finite construction and its costs

The exact generating series for the channel dimensions is

```math
\sum_{E\ge0}d_p(E)q^E
=\prod_{\substack{k\ge1\\p\nmid k}}
       (1+q^k+\cdots+q^{(p-1)k})
=\frac{P(q)P(q^{p^2})}{P(q^p)^2},
\qquad P(q)=\prod_{k\ge1}(1-q^k)^{-1}.
\tag{6}
```

No partition asymptotic is needed to obtain enough labels at polynomial
weight. For a requested binary width `lambda>=1`, take
`L=4lambda+4` distinct positive parts not divisible by `p`.
Their `j`th smallest member is at most `2j`, since `p>=2`.
The `2^L` subsets are regular partitions, all of weights at most
`S=L(L+1)`. Therefore one weight `h<=S` has

```math
d_p(h)\ge\frac{2^L}{S+1}\ge2^\lambda.
\tag{7}
```

The last inequality holds at `lambda=1` and persists by induction:
the numerator relative to `2^lambda` multiplies by eight each step,
whereas `(4lambda+4)(4lambda+5)+1` grows by a factor below four.
The maximizing weight cannot be zero because `d_p(0)=1<2^lambda`.

Compute the counts through `S` by integer dynamic programming and choose
the least maximizing weight. Retain the first `2^lambda` regular
partitions in a fixed lexicographic order. This specifies a uniform
construction with `h=O(lambda^2)` and `4h+1` pieces.

The source's ordinary partition upper bound gives
`sum_(E<=4h) p(E)=2^O(lambda)`. Thus the state and operator enumeration,
normal-order coefficient recursion, and both factors can be computed
in `2^O(lambda)` field operations, with fixed constants. For a fixed
prime `p`, the field coefficients have constant bit length. The
integer count tables and all indices have polynomial bit length.
The same algorithms work over `F_p` before any larger coefficient
algebra is introduced.

This proves effectiveness, polynomially many pieces, and the exact
channel formulas (5). It does **not** prove that this particular
weight-selection rule gives the sharp half-rate coefficient, or a
coefficient below it. Formula (6) and smaller finite ranks are not a
substitute for the missing uniform width inequality.

## 5. What the finite rank checks establish

The checker constructs the pieces over the prime field from the source
normal-order recurrence, with the quotient just proved. It verifies the
duality and the complete summed crossing, and independently forms the
two operator orders for its recorded cases.

Some characteristic-two examples are:

| Weight `h` | Regular labels | Largest row-rank times column-rank |
| ---: | ---: | ---: |
| 8 | 2 | 4 |
| 12 | 3 | 36 |
| 16 | 5 | 210 |

At weight 16 the central row and column ranks are 15 and 14, compared
with the full matrix dimension 25. These are ranks over `F_2`; a modular
rank drop is not a rational rank drop. The small two-label example does
not supply a large-label family with coefficient one half.

The characteristic-two checker also asks whether the crossing lies in
the span of the lower-rank pieces together with their exchanges
`y<->u`. The exchange fixes the crossing and swaps the two ranks.
For the three displayed cases, the least allowed rank-product cutoff
that spans the crossing remains respectively 4, 36, and 210. Thus this
particular finite recombination screen does not eliminate the costly
pieces. It is not a search over arbitrary tensors or decompositions.

A larger symmetry screen does give a finite improvement at weight 12.
Allow every permutation of the four roles that preserves the two crossing
pairs `{x,v}` and `{y,u}` as an unordered pair of sets. There are eight
such permutations. They preserve the crossing and either retain or exchange
the two flattening ranks.

On the three labels `(7,5)`, `(9,3)`, and `(11,1)`, a sum of 17 permuted
pieces is exactly the crossing over `F_2`, and each has rank product at
most 27. This improves the original maximum 36 on this finite domain.
The checker contains the explicit list of piece indices and permutations,
recomputes all their entries and ranks, and verifies the 81 crossing
entries. Omitting one specified term makes the identity fail.
The evidence records the full certificate. There is no extrapolation from
these three labels to a growing family, and the integer ceiling costs of
a complete recursive evaluator are not inferred from the rank product.

### Why naive type aggregation does not amplify this certificate past half

Take the tensor power of this 17-term identity on the full `3^t`-label
domain. Group the ordered choices by their 17 multiplicities. There are
`binom(t+16,16)` types, which is polynomial in `t`, so this grouping
does repair the exponential number of separate words.

However, the certificate contains
`T(x,y,u,v)=A_(-16)(x,u,y,v)` with row rank nine and column rank three.
The type that chooses `T` in every position contains exactly one word,
namely `T` to the `t`th tensor power. Matrix rank is multiplicative
under tensor products over `F_2`, so this type has ranks `9^t` and
`3^t`. Its rank product is exactly `27^t`. There is no other summand
inside that type to cancel in characteristic two.

Measured against the full domain's label length `log2(3^t)`, the
leading width coefficient of this type family is therefore at least

```math
\frac{\log_2(27^t)}{4\log_2(3^t)}=\frac34>
\frac1{\sqrt2}.
```

Polynomial multiplicities, ceilings, and an `O(log t)` width allowance
cannot change this leading obstruction. Thus merely collecting tensor
words into types cannot turn this finite certificate into an improvement
over the current half-rate family. This statement concerns the full
ternary-word domain and these type pieces. Arbitrary restrictions to
binary-label subsets, or new linear recombinations between types, may
change ranks and require a separate argument; they are not excluded here.

Run:

```sh
python3 tools/check_modular_pieces.py
python3 tools/check_modular_pieces.py --full
```

The larger command records weight 16. The evidence file lists every
case, the rank profile, check counts, and source hashes. Negative
controls retain the old factorial normalization when it is zero and
try to use the quotient over the rationals.

## 6. Conditional transfer to ring arithmetic

Here is the interface needed from a future improvement. Fix a prime
`p`, constants `1/2<=rho<1` and `c`, and a uniform modular family for
every `lambda>=1`. Assume:

1. Its sum is `d_lambda` times the crossing, with
   `d_lambda` a nonzero element of `F_p`.
2. It has at most `(2lambda)^c` pieces, two explicitly computable
   factorizations, and widths whose sum is at most
   `4rho lambda+c log2(2lambda)`.
3. The widths, piece count, multiplier, and every factor or piece entry
   are uniformly computable in `2^(c lambda)` prime-field operations
   and control bit operations, with polynomially bounded output index
   lengths. The count and widths have their stated polynomial size bounds.

These are the algebraic and effective parts of Appendix A's Definition 7.
They do not require integer coefficient-height bounds because coefficient
arithmetic is counted separately. Widths above `2lambda` can be replaced
by the trivial full-input factorization exactly as in that definition.
Unlike the stronger particular Lean program interface, the paper's
Definition 7(c) permits the same exponential computation bound for the
count and widths. This ordinary transfer uses that paper interface.

**Ring-arithmetic transfer.** For any commutative `F_p`-algebra `R`,
there is one deterministic algorithm evaluating an `R`-valued
`(K_l,lambda)` clique using at most

```math
2^{e_l K_l^2\lambda}
\bigl(K_l(\lambda+1)\bigr)^{C_0K_l^2}
\tag{8}
```

ring operations, input-entry queries, and elementary control operations.
Here the schedule and recurrence are those of Appendix A:

```math
K_0=2^{x_0},\quad g_l=2^{3l+x_0}K_l,\quad K_{l+1}=g_l^2,
\quad e_0=1/2,\quad
e_{l+1}=(\rho+3\,2^{-l-5})e_l+3/K_l,
```

where `x_0>=18` and `2^x_0>=c+5`. The constant `C_0` depends only
on the family constants and the fixed prime, not on `R` or its
representation length.

**Proof.** Use the constructions in Appendix A, Lemma 8 and Theorem A(b),
of the pinned [paper](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf).
The grid identity, grouping, and local table definitions use finite sums
and products, so their identities hold over `R`. The grid's final division
is multiplication by `d_lambda^(-number_of_rectangles)`. This is an
invertible prime-field scalar. No inverse of a general ring element is
needed.

All grouping tests, orientation choices, and loop limits depend on widths
and piece indices. They do not inspect the coefficient values. Recount
the proof with every ring addition or multiplication costing one and
with `Xi=g^2(lambda+1)` in place of `g^2(lambda+B)`.
The source's integer height estimates are then unnecessary.

At level zero, enumerate the `2^(lambda K_0(K_0-1)/2)` assignments and
sum their products. At a recursive step put `N=g^4 lambda`.
There are at most `Xi^(c g^4/4)` piece patterns. Each pattern gives
both a row and a column network. Choose the smaller total and apply
the source's balanced grouping. The resulting width `L` satisfies

```math
K^2L\le(\rho+3\zeta)N+(c/4)g^4\log_2(2\lambda)+K^2,
\qquad \zeta=2^{-l-5}.
\tag{9}
```

Equation (9) is precisely its displayed equation (42). A grouped entry
is a sum of at most `2^(2N/K)` products of polynomially many entries.
The source factor computations and control arithmetic cost
`2^(c lambda) Xi^O(1)` per such term. Also `L<=N+1`, so
`K(L+1)<=Xi^(c+5)` after increasing harmless constants.

Insert the induction hypothesis into the source's equations (43)--(44).
The grouping search, local enumeration, and factor computations add at
most `3N/K` to the exponential part. All remaining powers of `Xi`
are absorbed by `Xi^(C_0 g^4)` with the same constant-choice argument
as the source, since `(c+5)K^2<=g^4/2`. This is (8) at level `l+1`.
The term `(2|I|+2)2^(c lambda)` in the source's equation (44) also pays
for computing the piece count, multiplier, and the two widths of every
piece. Allowing those computations to take `2^(c lambda)` therefore
does not change the induction. All constructions are uniform.

For `R=F_p[Z]/(Z^q-1)`, store an element as `q` coefficients. Addition
costs `O(q)` and schoolbook cyclic multiplication costs `O(q^2)` for
fixed `p`. The complete cost is therefore (8) times a fixed polynomial
in `q` and the control bit lengths. A naive deterministic tape
simulation can add a fixed polynomial in this entire operation bound.
It changes the leading absolute constant and keeps the dependence
on `q` polynomial with a fixed degree. It does not introduce the
factor `q^(C_0K_l^2)` that would result from putting `B=q` into the
integer coefficient-height bound.

This is a new ordinary ring-arithmetic argument, not an invocation of
the existing integer-valued Lean theorem on ring inputs. To use it for
a stronger factoring rate, a modular family must still satisfy its
uniform width hypothesis with `rho<1/sqrt(2)`. Sections 1--5 do not
discharge that obligation.

## 7. The full half-rate coefficient over characteristic two

For `p=2`, the regular partitions are exactly the partitions into distinct
odd parts. Let `d(E)=d_2(E)`. The classical asymptotic is

```math
d(E)\sim
\frac{\exp(\pi\sqrt{E/6})}{2\,24^{1/4}E^{3/4}}
\qquad(E\longrightarrow\infty).
\tag{10}
```

It is stated explicitly for distinct odd parts in Remark 3.3 of
Walter Bridges, Johann Franke, and Taylor Garnowski,
[*Asymptotics for the twisted eta-product and applications to sign changes
in partitions*](https://doi.org/10.1007/s40687-022-00355-x), 2022.
An independent derivation for the equinumerous self-conjugate partitions
appears in Section 5.1 of Mark Wildon,
[*Counting partitions on the abacus*](https://arxiv.org/html/math/0609175v1),
2006. The equality of the two counts also follows from the diagonal-hook
bijection: a self-conjugate Young diagram has distinct odd diagonal-hook
lengths, and those hooks reconstruct the diagram.

Replace the generic weight-selection rule in Section 4 by

```math
h(\lambda)=\min\{h\ge1:d(h)\ge2^\lambda\}.
\tag{11}
```

Existence and the polynomial bound follow already from (7):
`h(lambda)<=(4lambda+4)(4lambda+5)`. The count table computes (11)
in polynomial time. No monotonicity assertion about `d(h)` is needed.

Put `a=pi/(sqrt(6) ln2)`. Formula (10) gives, as `h` tends to infinity,

```math
\log_2d(h)=a\sqrt h-\tfrac34\log_2h+O(1).
\tag{12}
```

The least weight in (11) tends to infinity with `lambda`. For large
`lambda`, apply (12) at `h-1` and use `d(h-1)<2^lambda`:

```math
a\sqrt h\le\lambda+O(\log h)
=\lambda+O(\log\lambda).
\tag{13}
```

The change from `sqrt(h-1)` to `sqrt h` is bounded by one. This proves
(13) even though the sequence of small counts is not monotone.

Formula (10), after increasing one constant to cover finitely many small
weights, also gives the global upper bound

```math
d(E)\le C\,2^{a\sqrt E}\qquad(E\ge0).
```

Use it on every term in (5). Uniformly for `-2h<=s<=2h`, the two
integer bit widths satisfy

```math
\begin{aligned}
\lceil\log_2 R^{\rm row}_{2,h,s}\rceil
 &\le a\sqrt{2h+s}+O(\log(h+1)),\\
\lceil\log_2 R^{\rm col}_{2,h,s}\rceil
 &\le a\sqrt{2h-s}+O(\log(h+1)).
\end{aligned}
\tag{14}
```

Each sum has `2h+1` terms, and all its nonzero energies are at most
the energy displayed on the right. The bound thus includes the endpoints
where one of those energies is zero. Concavity of the square root and
(13) now give

```math
c^{\rm row}_s+c^{\rm col}_s
\le2a\sqrt{2h}+O(\log(h+1))
\le2\sqrt2\,\lambda+O(\log(2\lambda)).
\tag{15}
```

Choose one fixed integer `c` large enough to dominate the constant in
(15), the effective construction costs from Section 4, and the finitely
many remaining small widths. Since `log2(2lambda)>=1`, these finite
exceptions are absorbed into `c log2(2lambda)`. All three hypotheses
of Section 6 then hold with `rho=1/sqrt(2)`, common scale one, and
the algorithmically specified selection (11).

This gives a uniform half-rate clique evaluator over every commutative
`F_2`-algebra, with the algebra dimension charged by a fixed polynomial
as proved in Section 6. It is an alternative algebraic foundation for
the current exponent. Obtaining a strictly larger factoring exponent
still requires a smaller uniform width coefficient or another evaluator.

## 8. Composition with the weighted divisor transfer

Sections 6--7 and the [weighted constructor](weighted-prefix-constructor.md)
give an alternative complete ordinary proof of the existing factoring rate.
It uses characteristic two throughout the coefficient algebra.

Choose the same largest level with `P=K_l` and `P^40<=n` as the
[half-rate proof](half-rate-factoring.md). Its established schedule bounds give

```math
n<P^{400},\qquad e_l=O((\log n)^{-1/2}),\qquad P^2\lambda\le4717n.
```

Consequently the logarithm of (8) is

```math
O\!\left(\frac{n}{\sqrt{\log n}}+P^2\log n\right)
=O\!\left(\frac{n}{\sqrt{\log n}}\right).
```

Constructing and reading the weighted tables costs
`2^{O(n/P)}poly(n,q)`. Even a linear scan for each table query only
multiplies the evaluator's operation count by that bound. Ring operations
cost `poly(q)`, and a fixed polynomial tape simulation changes only
absolute constants in the resulting exponent. Since
`n/P<n^{399/400}` and `P^2<=n^{1/20}`, all these costs together give

```math
2^{O(n/\sqrt{\log n})}\operatorname{poly}(n,q)
```

for each weighted query. Apply the
[fixed-characteristic transfer](fixed-characteristic-transfer.md).
Its extra `2^{O(n/log n)}` factor is absorbed, so the final machine
again computes the smallest prime factor within
`2^{O(n/sqrt(log n))}` steps on every input, with the fixed finite fallback.

More generally, any effective family satisfying Section 6 with
`1/2<=rho<1` gives exponent `log_2(1/rho)` through the same argument.
This exponent is at most one, so the sparse-polynomial overhead fits.
This closes the transfer from the proposed modular pieces to factoring.
A family with `rho<1/sqrt(2)` is still the missing improvement.
