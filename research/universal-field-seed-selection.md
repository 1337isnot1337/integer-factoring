# One factoring machine can select its finite field with its seed

## Result and scope

The finite field in the [gap-free seed theorem](gap-free-seed-score-transfer.md)
need not be fixed across input lengths. One uniform factoring machine can
search a stream of explicit seeds over different finite fields, select one
field at the start of an input, and use that field throughout its recursion
and recovery queries. This attains the infimum of the certificate parameter
over the entire stream.

Changing fields at different levels of a single tensor contraction would
not preserve its value. The construction here makes one selection for the
whole input and proves its cost directly; it does not apply a fixed-field
theorem outside its hypotheses.

For a valid seed of binary width `a`, with `M` nonzero pieces and row/column
ranks `R_i,C_i`, retain the exact rational quality

```math
Q(S)=\frac{\max_i(\lceil\log_2R_i\rceil+
                         \lceil\log_2C_i\rceil)
                 +2\lceil\log_2M\rceil}{4a}.
\tag{1}
```

Let a fixed computable stream print explicit finite-field seed records.
Allow invalid records, which are rejected by exact validation. Include the
one-piece crossing over `F_2` as a fallback, and put

```math
\rho_*=\inf_S Q(S).
\tag{2}
```

If `rho_*<1`, one deterministic finite-alphabet machine computes the
smallest prime factor of every binary integer and has score at least
`log2(1/rho_*)`. There is no seed-gap, generation-speed, or field-growth
hypothesis. Every fixed record eventually becomes eligible; that statement
is sufficient for the limiting score.

If the stream is closed under finite type amplification, the infimum in
(2) equals

```math
\inf_S\frac{\log_2\bigl(\sum_i\sqrt{R_iC_i}\bigr)}{2a}.
\tag{3}
```

In particular, a universal enumeration of all explicit crossing
decompositions over all finite fields gives one machine attaining the best
finite-certificate score available to this rank-sum method. The
[finite-field specialization lemma](finite-field-rank-specialization.md)
shows that this infimum is unchanged if one allows **all coefficient
fields**, including characteristic zero and fields without a computable
presentation. Every finite rank profile has an exact finite-field
counterpart. The infimum can also be attained by a search restricted to
prime fields after finite type amplification and linear projection.

These statements give no new numerical bound for (2) or (3). The
unconditional factoring exponent proved in the repository remains `1/2`.

The ingredients are the [vanishing-slack recursion](vanishing-slack-score-transfer.md),
the [finite type construction](gap-free-seed-score-transfer.md), and the
same pinned weighted factoring and sparse-nonzeroness reductions cited
there. This note supplies the uniform selection and cost analysis when
the coefficient field changes with input length.

For the existential endpoint alone, [universal verified factoring
search](universal-verified-factoring-search.md) is a shorter alternative:
every fixed finite seed gives a program, and polynomially verified
universal simulation preserves its score. The explicit selection below
also describes the field, local decomposition, and recursive resource
bounds of a machine within the tensor method.

## 1. Finite fields and exact records

A field record consists of a prime `p`, a positive integer `e`, and a monic
irreducible polynomial `f` of degree `e` over `F_p`. It represents

```math
\mathbb F=\mathbb F_p[Z]/(f),\qquad q=p^e.
```

Elements are coefficient lists of length `e`, with coefficients in
`0,...,p-1`. A seed record contains its width and every entry of every
piece tensor. Validation checks the field, dimensions, entries, the exact
crossing sum, and the two ranks by field elimination. Zero pieces are
discarded before calculating (1).

All validation can be elementary. Test primality by trial division and
irreducibility by testing all monic possible factors of degrees at most
`e/2`. Find field inverses by enumerating nonzero field elements. These
procedures are polynomial in `q` and the explicit table size; no
unproved number-theoretic or rank oracle is used. Record headers are
checked against the actual available record length before allocating a
claimed tensor.

The exact crossing identity is scale one over the selected field. It
therefore extends to every commutative algebra over that same field,
including the cyclic polynomial rings used for divisor recovery.

## 2. Select once from the input bit length

For input bit length `n>=2`, let `L(n)=floor(log2 n)+1` be its integer bit
length, and define

```math
B(n)=\max\{1,\lfloor\log_2 L(n)\rfloor\}.
\tag{4}
```

Run the seed enumerator for `n` steps on one fixed universal machine.
Among the completely printed valid records, accept only those satisfying

```math
a\le B(n),\qquad M\le B(n),\qquad
\lceil\log_2 q\rceil\le B(n).
\tag{5}
```

Apply cheap size and length filters **before** primality, irreducibility,
inversion, tensor allocation, or rank validation. First require `p>=2` and
`e>=1`, and reject `p>2^B` or `e>B`. Compute `p^e` with arithmetic capped
at `2^B`, rejecting it if the cap is exceeded. Compare the claimed width
with `B`, and every remaining dimension and coefficient-list length with
the actual available record length, before allocating the corresponding
objects. This prevents even an invalid short record from requesting an
exponential field search or an allocation based on a huge header. Once
these preliminary checks pass, validate the field and entries, discard
zero pieces, and enforce `M<=B` on the remaining pieces.

Choose minimum rational quality, with a fixed tie rule, also considering
the fallback. Denote the chosen seed, field, and quality by `S_n,F_n,Q_n`.
Keep this choice for every prefix query and every recursive level on this
input. Complete factorization can either pad all cofactor queries to the
original bit length, retaining its choice and schedule, or make a fresh
selection and schedule for each smaller residual integer. The polynomial
number of calls does not change the score. Keeping the original seed but
silently replacing the original length in the cutoff estimates would not
be the same algorithm.

The cap satisfies `B(n)=O(log log(n+2))` and `q<=2^B<=2L(n)`.
After discarding zero pieces, every accepted seed table has size at most

```math
M\,2^{4a}\lceil\log_2q\rceil
\le B(n)^2 2^{4B(n)}=(\log(n+2))^{O(1)}.
\tag{6}
```

The raw bounded output has at most `O(n)` symbols, including any zero
pieces or invalid records. There are at most `O(n)` completed records. The
elementary validation from Section 1, including both matrix ranks, is
therefore polynomial in `n`. The same applies to bounded simulation,
parsing, and exact rational comparisons. This is an ordinary algorithm,
not an unbounded preliminary search.

Every fixed valid seed is eventually completely printed and satisfies
(5), since its output time, width, piece count, and field size are finite
constants. The eligible set grows with `n`. Consequently

```math
\tfrac12\le\rho_*\le Q_n\le1,
\qquad Q_n\longrightarrow\rho_*.
\tag{7}
```

The lower bound is the rank-subadditivity argument in the gap-free note.
No convergence modulus or approximation to `rho_*` is an instruction.

## 3. Use the selected seed at all expanded widths

Use the same fixed schedule and original-input volume cutoff as the
vanishing-slack proof. Write `P=K_ell`, where `P^40<=n` is maximal.
Every expanded call has

```math
\lambda\ge n^{9/20},\qquad V=K_j^2\lambda>\sqrt n.
\tag{8}
```

At that call, amplify `S_n` by type over
`k=floor(lambda/a)` blocks and use the ordinary crossing on the fewer
than `a` remaining bits. The finite type identity gives scale one over
the same `F_n` at every level. Exact factorization and padding give

```math
r_i,c_i\le2\lambda,
\qquad r_i+c_i\le4Q_n\lambda+4B(n)+2,
\tag{9}
```

and piece entropy

```math
f_n(\lambda)\le B(n)\log_2(\lambda+1).
\tag{10}
```

If a finite small-input call has `lambda<a`, its whole crossing is the
remainder; the construction is still exact. Those finite exceptions do
not require a different machine.

Coefficient dynamic programming, complete matrix construction, and exact
rank factorization have `2^O(lambda)` field operations, as in the gap-free
proof. The polynomial factors in `n` and in `log q` are charged separately.
Their constants are uniform: all seeds use the same table, polynomial,
and elimination routines, rather than seed-specific programs. At expanded
widths, the type-state exponent `B(n)log(lambda+1)` is `o(lambda)`.

## 4. Uniform cost of the input-dependent recursion

For the asymptotic analysis, put

```math
\beta_n=Q_n+
       \frac{B(n)+1}{n^{9/20}}+n^{-19/40}.
\tag{11}
```

The source grouping estimate with (8)--(9) yields

```math
V_{j-1}\le(\beta_n+3\zeta_{j-1})V_j.
\tag{12}
```

Here `beta_n -> rho_*<1`. Eventually it is bounded by one fixed constant
strictly below one, while remaining at least `1/2`. The same finite product
estimate as in the vanishing-slack proof gives

```math
V_j\le D\beta_n^{\ell-j}V_\ell,
\qquad \sum_jV_j=O(n),\qquad \lambda_j=O(n),
\tag{13}
```

with constants independent of the selected field and seed. The last
estimate follows from the first two; no estimate of that constant is used
to run the algorithm.

Equations (8), (10), and (13) bound the total piece-assignment exponent by

```math
\sum_j\frac{V_j}{4\lambda_j}f_n(\lambda_j)
\le O\!\left(n^{11/20}B(n)\log(n+2)\right).
\tag{14}
```

As before, this is a bound along a worst recursion path obtained from
each node's sum of child costs; it bounds the whole tree. Repeated parent
queries and local factors contribute `O(n beta_n^ell)` to the exponent.
The top constructor contributes `O(n/P)`, and the cutoff leaves and
structural searches contribute `O(sqrt n)` and smaller terms. All uniform
seed processing and arithmetic-polynomial overheads fit a polynomial
factor in `n` or the already charged exponential factors.

## 5. Coefficient rings and final score

Trial division through `2^ceil(sqrt n)` either returns a factor or proves
the same `2^O(sqrt n)` divisor-support cap as in the vanishing-slack
transfer. Use its sparse cyclic tests over

```math
F_n[X]/(X^l-1).
```

The number of tests and their largest ring dimension are `2^O(sqrt n)`.
A coefficient of `F_n` has at most `B(n)` bits; its arithmetic and the
represented-ring arithmetic therefore remain within `2^O(sqrt n)` per
recursive arithmetic operation, up to fixed polynomial factors. The
isolation argument uses only that distinct divisor monomials have
coefficient one. It is valid over every selected field. Any auxiliary
prime required to differ from the characteristic can simply omit that
one prime, without changing the bound.

Thus the complete bit complexity satisfies

```math
\log_2T(n)\le
O\!\left(n\beta_n^\ell+
 n^{11/20}B(n)\log(n+2)+\sqrt n+n/P\right).
\tag{15}
```

Because `2^ell=Theta(log n)`, `beta_n -> rho_*`, and
`P>n^(1/400)`, equation (15) gives

```math
T(n)\le2^{\,n/(\log n)^{\log_2(1/\rho_*)-o(1)}}.
\tag{16}
```

Polynomial input scanning and finitely many initial inputs are absorbed
in the same score interpretation. The program is independent of the
eventual score tolerance.

Finally, close the seed stream under finite type amplification within
each seed's own field. Section 6 of the gap-free note proves equality
of the rational-quality and rank-sum infima. Its argument is fieldwise
and therefore also proves (3) for the union of these streams. This
supplies one uniform machine for the best enumerable finite-field
certificate score, without treating a change of coefficient field
inside one contraction as a valid algebraic operation.

## Evidence and remaining requirement

The proof reuses the exact finite type identities checked by
`tools/check_gap_free_seed_transfer.py`. The accompanying
`tools/check_seed_field_selection.py` checks small field presentations,
delayed bounded selection, rank invariance under a field extension,
and the integer resource caps. These are finite validation of the
selection mechanism; they do not exhibit a better rank-sum seed or
execute the asymptotic factoring machine. No Lean changes are involved.
