# A uniform score transfer from arbitrarily sparse finite seeds

## Result and scope

The score-preserving vanishing-slack theorem removes the gap condition
from sparse-width interpolation. The good seed widths may be arbitrarily
far apart, and their generator may take arbitrarily long between outputs.
A bounded search followed by type amplification still gives one uniform
machine attaining their limiting score.

This is a refinement of two existing ordinary results:

- [Finite-seed amplification at 1476525](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/finite-seed-amplification.md)
  already turns one sufficiently strong finite rank-sum certificate into
  a strict factoring improvement.
- [Sparse-width interpolation at 05d9c1b](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/sparse-width-interpolation.md)
  gives quantitative constant-factor bounds under polynomial gaps.

The new conclusion here is the uniform endpoint score without a gap or
generation-rate hypothesis. It does not supply a better seed or prove
a factoring exponent above `1/2`. A finite list of unsuccessful or
near-threshold seeds does not prove the required limiting hypothesis.

## 1. A rational quality attached to a finite certificate

Fix a finite prime field `F_p`. A seed of binary width `a>=1` is a
finite list of explicit tensors on `d=2^a` labels such that

```math
E_d(x,y,u,v)=[x=v][y=u]=\sum_{i=1}^{M}A_i(x,y,u,v).
```

Discard zero tensors. Let `R_i,C_i` be their exact positive ranks in
the two cuts `(x,y)|(u,v)` and `(x,u)|(y,v)`. Define

```math
\begin{aligned}
w&=\max_i\bigl(\lceil\log_2R_i\rceil+
                    \lceil\log_2C_i\rceil\bigr),\\
c&=\lceil\log_2M\rceil,\qquad
Q=\frac{w+2c}{4a}.
\end{aligned}
\tag{1}
```

All quantities in (1) are computed exactly with field elimination,
integer bit lengths, and rational comparison. `Q` is a sufficient
quality bound; the rank-sum quantity in the cited finite-seed theorem
can be better for an individual certificate. Section 6 recovers the
infimum of the rank-sum objective while retaining this rational selector.

Let a fixed program enumerate complete finite seed records. It may
print invalid records too: the algorithm verifies every candidate
before using it. Define

```math
\rho_* = \inf\bigl(\{1\}\cup
       \{Q(S):S\text{ is a valid seed printed by the program}\}\bigr).
\tag{2}
```

**Gap-free transfer theorem.** If `rho_*<1`, one effective all-width
scale-one family over `F_p` has separate widths at most `2lambda`,
total widths

```math
4\rho_*\lambda+o(\lambda),
```

piece entropy `O(sqrt(lambda) log(lambda+1))`, and uniform
`2^O(lambda)` construction and factor-query cost. Consequently
[vanishing-slack-score-transfer.md](vanishing-slack-score-transfer.md)
gives one deterministic factoring machine with score at least

```math
\log_2(1/\rho_*).
\tag{3}
```

The infimum in (2) is an analytical parameter. The machine never
computes it or receives a convergence modulus.

The quality of any valid seed is at least `1/2`. Indeed rank
subadditivity in the row cut gives some `R_i>=d^2/M`, so
`w>=2a-log2 M`. Since `2c>=2log2 M`, equation (1) implies
`Q>=1/2`. Thus (3) lies in the existing crossing interface's
permitted range.

## 2. A bounded, effective seed selector

On input label length `lambda`, run the fixed enumeration program for
exactly `lambda` steps on a fixed universal machine. Its output length
is at most a constant times `lambda`; changing this fixed constant
has no effect on the argument. Consider only complete records printed
within this budget. A record is an explicit tensor list, not an oracle
or a compressed assertion of rank.

For each complete record:

1. Check its field, dimensions, entries, and crossing sum.
2. Discard zero pieces and compute both exact matrix ranks.
3. Require `a^2<=lambda` and `M^2<=lambda`.
4. Compute its rational quality (1).

Include the one-piece width-one seed `A_1=E_2`, of quality one,
as a permanent fallback. Select a seed with minimum quality, breaking
ties by any fixed rule. Write its parameters as `a_lambda,M_lambda`
and `Q_lambda`. In particular,

```math
a_\lambda,M_\lambda\le\sqrt\lambda,
\qquad \rho_*\le Q_\lambda\le1.
\tag{4}
```

The fallback satisfies the size filters even at `lambda=1`.
Verification is polynomial in the printed record length: the tensors
are explicit, and field rank is computed by ordinary exact elimination.
A header claiming more entries than the complete printed record
contains is rejected without allocating the claimed object. Thus
bounded simulation, parsing, validation, and selection all have
polynomial cost in `lambda`.

Every particular valid seed is eventually included. Its finite output
time, record length, `a^2`, and `M^2` are constants for that seed.
For every sufficiently large `lambda`, the simulation has printed it
completely and both filters accept it. Therefore, for every `eta>0`,
choose one printed seed with `Q<rho_*+eta`; eventually

```math
\rho_*\le Q_\lambda<\rho_*+\eta.
```

This proves `Q_lambda -> rho_*`. The argument does not require a
bound on how large that eventual threshold is. Invalid records cannot
affect correctness or lower the chosen quality because they are
independently rejected.

## 3. Amplify the selected seed by piece type

Set

```math
k=\lfloor\lambda/a\rfloor,\qquad b=\lambda-ka,
\qquad0\le b<a,
```

where `a=a_lambda`, `M=M_lambda`. Split each external label into
`k` consecutive `a`-bit blocks and one `b`-bit remainder.

For each vector of nonnegative integers
`nu=(nu_1,...,nu_M)` with `sum nu_i=k`, define

```math
B_\nu=
\sum_{\substack{j_1,\ldots,j_k\in\{1,\ldots,M\}\\
                 |\{t:j_t=i\}|=\nu_i\ \forall i}}
A_{j_1}\otimes\cdots\otimes A_{j_k}.
\tag{5}
```

Take the new piece to be `B_nu tensor E_(2^b)`. Every piece-index
word appears in exactly one type. Finite distributivity gives an exact
scale-one crossing sum on all `lambda` bits. No averaging or division
by a factorial occurs, so the identity is valid in every characteristic.

There are at most

```math
(k+1)^M
\tag{6}
```

types. By (4), their entropy is at most
`sqrt(lambda) log2(lambda+1)`. This is smaller than
`lambda/(log lambda)^delta` for every fixed `delta>0`, eventually.

Let `D_nu` be the multinomial coefficient counting the words in (5).
Tensor-product rank and rank subadditivity give

```math
\begin{aligned}
\operatorname{rank}_{row} B_\nu
 &\le D_\nu\prod_iR_i^{\nu_i},\\
\operatorname{rank}_{col} B_\nu
 &\le D_\nu\prod_iC_i^{\nu_i}.
\end{aligned}
```

As `D_nu<=M^k`, the sum of the two rounded rank widths, including
the remainder crossing, is at most

```math
\begin{aligned}
2\log_2D_\nu+
 \sum_i\nu_i(\log_2R_i+\log_2C_i)+4b+2
&\le k(w+2c)+4b+2\\
&=4Q_\lambda\lambda+4(1-Q_\lambda)b+2\\
&\le4Q_\lambda\lambda+4\sqrt\lambda+2.
\end{aligned}
\tag{7}
```

The additive two pays for rounding the two final ranks. Zero new
pieces are omitted. Exact Gaussian elimination gives factorizations
through their actual ranks; pad to the next power of two. Each width
is at most `2lambda`, since its flattening has `2^(2lambda)` rows
and columns.

Combining (7) with `Q_lambda -> rho_*` proves the required
`4rho_* lambda+o(lambda)` total-width bound.

## 4. Complete construction and query cost

An entry of (5) is the coefficient of `X^nu` in

```math
\prod_{t=1}^{k}
 \left(\sum_{i=1}^{M}
 A_i(x_t,y_t,u_t,v_t)X_i\right).
\tag{8}
```

A coefficient dynamic program has at most `(k+1)^M` states.
By (4),

```math
\log_2 (k+1)^M\le\sqrt\lambda\log_2(\lambda+1)=O(\lambda).
```

Its transition count adds only a polynomial factor in `lambda`.
The selected seed tables were explicitly read and verified, so a seed
entry is obtained by finite table access with its storage cost charged.

There are at most `2^O(lambda)` new pieces and `2^(4lambda)`
entries per piece. Building every entry using (8), flattening in both
cuts, computing exact ranks, and constructing exact rank
factorizations all take `2^O(lambda)` bit operations. For instance,
ordinary elimination on order `2^(2lambda)` uses at most
`2^(6lambda)` field operations; multiplying by the piece count and
the coefficient-construction cost still gives a fixed exponential
constant. The field is fixed.

Piece indices are their type vectors, using
`O(M log(lambda+1))=O(lambda)` bits. The complete index set can
be enumerated within the same bound. A factor-entry query may
regenerate the entire construction, so it too has `2^O(lambda)`
time. No local-space bound is asserted. This verifies the effective
interface required by the score transfer, including all construction
costs rather than only abstract ranks.

## 5. Corollary for an unbounded sequence of good widths

Suppose a computable enumerator prints valid seeds at an unbounded
sequence of widths `a_j`, with

```math
w_j\le4\rho a_j+o(a_j),
\qquad \log_2M_j=o(a_j),
\qquad 1/2\le\rho<1.
\tag{9}
```

Then the qualities in (1) satisfy `limsup Q_j<=rho`, so
`rho_*<=rho`. The theorem yields score at least
`log2(1/rho)`. No bound on `a_(j+1)/a_j`, on
`log a_(j+1)/log a_j`, or on seed-generation time is needed.
Polynomially many pieces are a special case of (9).

This corollary concerns effective finite certificates. An uncomputable
choice of good tensors, a conjectural rank profile, or finite evidence
at a bounded list of widths does not meet the hypothesis. Conversely,
even a single seed of quality strictly below `1/sqrt(2)` already
gives a strict improvement by the earlier finite-seed theorem; the
gap-free endpoint result does not remove that underlying mathematical
search problem.

The constant-factor runtime `exp(O(n/(log n)^alpha))` is not
asserted at the infimum when its convergence is arbitrarily slow.
The conclusion is exactly the challenge's score interpretation,
`2^(n/(log n)^(alpha-o(1)))`, for one machine.

## 6. The rank-sum infimum, without a real-number oracle

For a valid width-`a` seed `S`, put

```math
\mathcal S(S)=\sum_{i=1}^{M}\sqrt{R_iC_i},\qquad
\rho(S)=\frac{\log_2\mathcal S(S)}{2a}.
\tag{10}
```

Starting from any computable stream of explicit finite seeds, form a
new stream which also prints every finite type amplification of every
printed valid seed. Dovetail all construction jobs, including jobs
started from earlier amplified seeds. This produces an effective
amplification-closed stream: every particular finite construction is
eventually finished and printed. No bound on that completion time is
needed. A table and its exact ranks can always be constructed by finite
enumeration and elimination; only the selector's bounded prefix is
examined at any target width.

For a fixed seed and a tensor power `k`, let `S^(k)` be its type seed,
with zero type tensors removed. A type `nu` has multiplicities summing
to `k`, and its tensor is the sum of the corresponding ordered words.
With `B_nu=k!/prod_i nu_i!`, rank subadditivity and tensor-product rank
give

```math
\begin{aligned}
\sqrt{R_\nu C_\nu}
&\le B_\nu\prod_i(\sqrt{R_iC_i})^{\nu_i}
 \le \mathcal S(S)^k,\\
\lceil\log_2R_\nu\rceil+\lceil\log_2C_\nu\rceil
&\le2k\log_2\mathcal S(S)+2,\\
M_k&\le(k+1)^M.
\end{aligned}
\tag{11}
```

The second inequality in the first line follows because the preceding
quantity is one nonnegative term in the multinomial expansion of
`mathcal S(S)^k`. Therefore

```math
Q(S^{(k)})\le\rho(S)
  +\frac{2+2\lceil M\log_2(k+1)\rceil}{4ak}
  =\rho(S)+o(1).
\tag{12}
```

For every finite seed, the reverse comparison at the same width is

```math
\mathcal S(S)\le M\,2^{w/2}\le2^{c+w/2},
\qquad \rho(S)\le Q(S).
\tag{13}
```

Equations (12)--(13) show that on any amplification-closed stream,

```math
\inf_S Q(S)=\inf_S\rho(S).
\tag{14}
```

Include the trivial seed in both infima. The gap-free theorem thus
attains the score associated with the exact rank-sum infimum of this
stream. For an arbitrary initial stream, closure may introduce seeds
with a better rank-sum value, so the guaranteed score is **at least**
the supremum of its initial finite-seed scores. Equality with the
initial infimum is not required or asserted.

One can in particular enumerate every explicit crossing decomposition
over a fixed finite field, in increasing finite encoding length, and
verify each record. The collection is already closed under type
amplification. There is consequently one uniform machine attaining
the infimum of the rank-sum parameter over all finite certificates
over that field, in the challenge's limiting-score sense. This is an
existence and uniformity result: the machine need not know the
infimum, and its seed search may be impractically slow. It supplies
no new numerical upper bound for that infimum.

Equation (12), together with `Q>=1/2` for every amplified seed, also
proves `rho(S)>=1/2` for every finite seed. Thus this extension stays
within the same permitted range as the rational theorem.

## Verification status

```sh
python3 tools/check_gap_free_seed_transfer.py
python3 tools/check_gap_free_seed_transfer.py --output evidence/gap-free-seed-transfer.json
```

The checker compares coefficient dynamic programming with independent
word enumeration in five amplification cases over `F_2` and `F_3`.
It checks 393,216 coefficient entries, computes both exact ranks for
21 nonzero type tensors, and verifies 65,536 integer entropy envelopes.
Its rank checks include type multiplicities, discarded zero types,
and a width-two seed with a one-bit remainder.

The construction is an ordinary proof. It reuses the finite type-sum
identity and exact rank-factor construction from the pinned finite-seed
note, but makes the chosen seed depend on the target width under an
explicit bounded search. The cost argument above pays for that growing
seed, growing type count, and regeneration of all factors. No new Lean
proof or numerical factoring benchmark is claimed.
