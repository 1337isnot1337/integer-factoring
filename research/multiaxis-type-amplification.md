# Type amplification removes the piece-entropy cap from the limiting score

## Result and scope

The [fixed-axis cube transfer](multiaxis-cube-transfer.md) has a stronger
limiting-score consequence than its direct accounting of piece assignments.
Fix an integer `d>=2`, put `m=2^(d-1)`, and fix a finite coefficient field
of characteristic two. Suppose a uniform local family has, at every
positive binary label length `t`:

- An exact decomposition of the antipodal equality tensor on `2^d` legs
  into `M_t` common pieces, with exact displayed factorizations at all
  `d` coordinate cuts.
- Separate channel widths `0<=r_i^a<=m t` and a total-width bound

  \[
  w_t:=\max_i\sum_{a=1}^d r_i^a\le dm\rho t+o(t),
  \qquad 1/d\le\rho<1.
  \tag{1}
  \]

- Subexponential piece count, in the precise sense

  \[
  \log_2 M_t=o(t).
  \tag{2}
  \]

- Fixed uniform algorithms for complete enumeration, metadata, and
  scalar factor entries in `2^O(t)` bit time. The field, algorithms,
  dimension, and constants are independent of the factoring input.

Then **one** deterministic finite-alphabet machine factors every input.
Writing `alpha=log_d(1/rho)`, for every fixed `epsilon>0` its time and
working space on all sufficiently large `n`-bit inputs satisfy

\[
T(n)\le2^{\,n/(\log_2 n)^{\alpha-\epsilon}},
\qquad
S(n)\le\frac{n}{(\log_2 n)^{\alpha-\epsilon}}.
\tag{3}
\]

Input is read-only and output is one-way write-only. The conclusion
includes the least prime factor, complete prime factorization, sorted
divisors, and exact prefix or interval counts with read-only endpoints.
There is no assumption on an effective convergence rate in (1)--(2).
No new numerical factoring exponent is obtained: no improving common
local family is supplied here.

For the relaxed hypothesis in Section 5 of the fixed-axis transfer,
direct recursion proves the stronger constant-form time estimate with
score `min(alpha,gamma)`. Amplifying small blocks instead attains score
`alpha`, including when `gamma<alpha`. It need not give the constant-form
bound `exp(O(n/(log n)^alpha))`; the distinction is the `epsilon` in
the leaderboard's definition. In particular, the failure of bare
`o(t)` piece entropy to guarantee a positive saving in the **direct
assignment accounting** is not an obstruction to this amplification.

The method is the fixed-axis extension of the repository's
[two-cut finite-seed endpoint](simultaneous-seed-endpoint.md) and
[type-grouped time/space transfer](simultaneous-time-space-transfer.md).
It retains every block permutation and charges it in every cut.

## 1. The quality of one actual finite seed

For a seed at length `t`, define

\[
c_t=\lceil\log_2 M_t\rceil,
\qquad
Q_t=\frac{w_t+d c_t}{dm t}.
\tag{4}
\]

All quantities in (4) are integers or rational numbers that the
enumeration computes. No real logarithm is required: `c_t` is the
smallest integer with `M_t<=2^c_t`. The endpoint is nonzero, so `M_t>=1`.
Zero pieces may be retained.

At any coordinate cut the endpoint has rank `2^(m t)`. Rank
subadditivity gives a piece of width at least `m t-log_2 M_t` at
that cut. Therefore

\[
Q_t\ge\frac{m t-\log_2 M_t+d c_t}{dm t}\ge\frac1d.
\tag{5}
\]

This argument uses the actual endpoint and its displayed factorizations.
An abstract collection of favorable widths without those identities is
not a seed. The one-piece endpoint at `t=1` is a permanent fallback
with `Q=1`.

Equations (1)--(2) imply `limsup Q_t<=rho`. A finite seed with quality
`Q<1` already gives a conditional exponent `log_d(1/Q)` after the
amplification below. Slowly increasing the admitted seed size attains
the best limiting quality with one machine.

## 2. A bounded selector and materialized small tables

A seed has `M_t<=2^O(t)` pieces, since its complete enumeration already
has that time bound. Each side of a cut has `2^(m t)` labels, and each
factor has at most `2^(m t)` channels. Thus both complete factor tables
at every cut contain `2^O_d(t)` entries per piece. Complete-piece
entries can be computed from any one displayed factorization.

There is a fixed integer `C_1>=1` such that generating all these tables,
checking their cut identities and common endpoint, and reading the
metadata takes at most `2^(C_1 t)` bit time and working bits. The field
is fixed, so its element size and arithmetic cost are constants.
Increasing `C_1` covers the finite encodings and small sizes. This is a
proved bound for the supplied fixed algorithms, not an attempt to infer
their running time from source code.

Set `B=64 C_1` and

\[
t_{\max}(n)=\max\left(1,
 \left\lfloor\frac{\lfloor\log_2(n+2)\rfloor}{B}\right\rfloor
 \right).
\tag{6}
\]

Generate and validate the seeds at `1<=t<=t_max(n)` sequentially.
Retain a minimum-quality seed among them and the fallback, with a
fixed tie rule. Only the current and best seed need be retained.
The geometric sum of the generation bounds is

\[
O(n^{1/64}+1)
\tag{7}
\]

in both time and space, up to family-dependent constants. Bit lengths
and integer division implement (6). There is no input-dependent advice.

Write `t_n,M_n,Q_n` for the selected seed's length, piece count, and
quality. In particular,

\[
\frac1d\le Q_n\le1,\qquad
M_n\le O(n^{1/64}),\qquad t_n=O(\log n).
\tag{8}
\]

If `Q_*` is the infimum of (4) over all seeds and the fallback, the
bounded selector gives `Q_n -> Q_*`: every fixed seed eventually
becomes available, and all selections have quality at least the
infimum. Here `Q_*<=rho<1`. The algorithm never computes this infimum
or a convergence modulus. The argument actually proves (3) with the
possibly larger exponent `log_d(1/Q_*)`.

The same selected seed and coefficient field are used at every level
and for every recovery query on the input.

## 3. Common pieces after type grouping

Fix the selected seed for one input. At a recursive binary label width
`lambda`, write

\[
\lambda=b t_n+u,\qquad0\le u<t_n.
\]

For each vector `nu=(nu_1,...,nu_(M_n))` of nonnegative integers with
sum `b`, define one tensor

\[
B_\nu=
 \left(\sum_{\mathrm{type}(v)=\nu}
 T_{v_1}\otimes\cdots\otimes T_{v_b}\right)\otimes E_u.
\tag{9}
\]

The zero-length endpoint `E_0` is the scalar one. The same tensor in
(9) is used at every axis cut. Summing over all types recovers
`E_lambda` exactly, because it is simply the full tensor product of
the seed crossing sums and the tail endpoint. No factorial is divided
in the coefficient field.

Let `D_nu=b!/prod_i nu_i!`. Across axis `a`, a channel consists of a
word of type `nu`, one original channel at each block, and one full
tail channel. Its cardinality before binary padding is

\[
R_\nu^a=D_\nu\,2^{\sum_i\nu_i r_i^a+m u}.
\tag{10}
\]

Products of the corresponding original left and right factors give
an exact factorization of (9). The tail factorization is the endpoint's
permutation matrix, with `m u` channel bits. Since `D_nu<=M_n^b`,
padding each whole channel set once gives

\[
\begin{aligned}
\sum_{a=1}^d\lceil\log_2R_\nu^a\rceil
 &\le b(w_{t_n}+d c_{t_n})+dm u+d\\
 &\le dm Q_n\lambda+dm t_n+d.
\end{aligned}
\tag{11}
\]

If one padded cut width exceeds `m lambda`, replace that cut by the
identity matrix times the complete piece matrix. This preserves the
tensor and only decreases its width. The complete-entry algorithm
below makes this cap effective; it is not an unspecified rank oracle.

The number of types is at most `(b+1)^(M_n)`, so a type descriptor and
a piece ordinal both require at most

\[
O(M_n\log(\lambda+2))=O(n^{1/64}\log(n+2))
\tag{12}
\]

bits at visited widths `lambda=O(n)`. Thus the original piece entropy
has moved into the small difference between `Q_n` and its limit.
The permutations have not been dropped from either (10) or (11).

### Uniform queries and their space

The fixed-axis recursion selects `P<=n^delta`, where
`delta=1/(100d)`, and stops before processing volume at most `sqrt(n)`.
Every active crossing query therefore has

\[
\lambda>n^\eta,\qquad\eta=1/2-2\delta.
\tag{13}
\]

The retained seed tables in (7) and the descriptors in (12) fit in
`O(lambda)` working bits. In particular, (12) is a storage bound on
the actual count vector, not merely on the number of pieces.

For a factor query, unrank the type word using its first-letter
counts `D_nu nu_i/b`. The remaining mixed-radix channel code decodes
the original channels. Their product cardinality is independent of
the order of the type word, so this numbering is bijective. Multiply
the requested seed factors sequentially. For a complete-piece query,
enumerate the words in `[M_n]^b`, keep the required type, and sum
their products of complete seed entries. This also supplies the
capped factor entries.

The word counter has at most `b c_(t_n)=O(lambda)` bits. Indeed,
`Q_n<=1` in (4) gives `c_(t_n)<=m t_n`. Channel integers likewise
have `O_d(lambda)` bits. Do not retain `b!`, whose length could be
larger. Compute `D_nu` as a product of binomial coefficients, using
the smaller lower argument and the exact recurrence
`z_(j+1)=z_j(h-j)/(j+1)`. Those partial binomial coefficients are
increasing and at most their final values. Every intermediate product
is at most `D_nu<=M_n^b<=2^(m lambda)`; a numerator multiplication
adds only `O(log(lambda+2))` bits.

Only a constant number of field accumulators and one seed lookup are
active at a time. Factor, complete-entry, metadata, and streamed type
enumeration routines therefore use `O_d(lambda)` auxiliary bits and
`2^O_d(lambda)` bit time, with constants uniform in `n` and the selected
seed. These bounds are asserted at the visited widths (13), rather
than at unrelated tiny `lambda` with a large seed already selected.

## 4. The variable quality has a uniform volume envelope

Use exactly the integer schedule, equitable grouping, and geometric
cube partition of the fixed-axis transfer:

\[
K_0=2^{x_0},\quad g_j=2^{a(j+1)+x_0}K_j,\quad
K_{j+1}=g_j^d,\quad \zeta_j=2^{-j-8},\quad 2^a>4d.
\tag{14}
\]

The fixed `x_0` covers the uniform local algorithms above. Choose the
largest level `ell` with `K_ell^(100d)<=n`. The source proof gives
`d^ell=Theta_d(log n)` and
`P=K_ell>n^(delta/(3d))`. The finitely many smaller lengths use a
fixed terminating fallback.

Substitute (11) for the width estimate in that transfer. Its complete
grouping bound, including rounding, gives

\[
V_j\le\left(Q_n+O_d(\zeta_j+1/g_j)
                 +O_d((t_n+1)n^{-\eta})\right)V_{j+1}
\tag{15}
\]

at every active step. The separate caps first give the crude bound
`V_j=O(n)`, exactly as in the fixed-axis proof, so the descriptor bounds
have no circular dependence on the refined estimate.

Since `Q_n>=1/d`, dividing the errors in (15) by `Q_n` is harmless.
The two level-dependent error series converge, and
`(t_n+1)n^(-eta) ell -> 0`. Also `Q_n` is eventually bounded above
by some fixed number below one. Consequently, uniformly along all
visited paths,

\[
V_j\le D Q_n^{\ell-j}n,\qquad\sum_j V_j=O(n)
\tag{16}
\]

for one constant `D`. Neither `Q_n` nor any piece choice changes the
machine's schedule or algorithms.

At a call of order at most `P`, the complete cube assignment has at
most `P^2` piece ordinals. Equations (8) and (12) bound all such
assignments along a path by

\[
O(n^{2\delta+1/64}\log n\log\log n)
\tag{17}
\]

bits and logarithmic enumeration time. This is smaller than `n/log n`.
Grouping maps, width lists, and their exhaustive searches are smaller
as well. The new piece count is therefore fully charged, even though
it need not be polynomial in `lambda` with an input-independent degree.

## 5. Packed workspace, actual queries, and factor output

The fixed-axis proof gives one table-query composition with logarithmic
time cost `O_d(V/g_j+V/K_j^2+lambda)`, apart from logarithmic terms.
An implementation retains the current plane-label assignment, the
internal grouped channels, and only the arguments for the current
old-table request. For a conservative workspace bound also charge one
complete child boundary of length `O_d(V/K_j)`. A local query therefore
uses at most

\[
O_d\!\left(V/g_j+V/K_j+\lambda\right)
\tag{18}
\]

working bits, plus the piece and grouping metadata already charged.
Factor products and their scalar queries are sequential.

There is at most one suspended evaluation frame per level. A table
query then moves toward older parent tables and never calls the
evaluator recursively. It adds at most one table-query frame at each
level. Thus (18) is summed once per level up to a fixed factor, not
once for every old vertex simultaneously. Using (16) and `Q_n>=1/d`,

\[
\begin{aligned}
\sum_j O_d\!\left(V_{j+1}/g_j+V_{j+1}/K_j+
                         V_{j+1}/g_j^{2d}\right)
 &\le O_d(n Q_n^\ell)
 \sum_{j\ge0}d^{j+1}
   (g_j^{-1}+K_j^{-1}+g_j^{-2d})\\
 &=O_d(n Q_n^\ell).
\end{aligned}
\tag{19}
\]

The series converges because the schedule grows doubly exponentially.
At level zero the terminal edge assignment has `O(n Q_n^ell)` bits;
an early terminal call has `O(sqrt(n))` bits.

Use the [streamed top constructor](streamed-top-space.md), rather
than materializing its tables. Its scratch and virtual addresses use
`O(n/P+n^(3/4)+log n)` bits. All top-input bits are regenerated from
the original read-only input. The exact costs of that regeneration
are included in the existing constructor's polynomial query time.

For factor recovery use the paired short descriptors from
[sublinear-space-half-rate.md](sublinear-space-half-rate.md). Its
uniform divisor bound gives descriptor length `O(n/log n)`, and its
scalar bit queries use a fixed commutative extension ring of the
fixed characteristic-two field `F`. One explicit choice is
`F[X,epsilon]/(X^2+X+1,epsilon^2)`, a four-dimensional `F`-algebra;
the quadratic need not be irreducible over `F`. Its image of `X`
satisfies `X^3=1` and `1+X+X^2=0`, which are the exact character
identities used in the recovery. The supplied crossing and cut
factorization identities remain exact in this fixed ring. Its size
does not grow with `n`; the growing ternary hash field is only a
control representation, with `O(n/log n)` stored bits.

Every accepted descriptor is checked by streamed multiplication of
its two virtual factors. Comparisons, duplicate removal, sorted
enumeration, exact counting, and prime-power extraction retain only
a constant number of descriptors. They add `O(n/log n)` workspace
and `exp(O(n/log n))` time, as proved by the cited recovery theorem.
They do not store the returned factor, an ordinary binary prefix of
length `n`, or a mutable cofactor. All output bits are regenerated.

Since `Q_n>=1/d` and `d^ell=Theta(log n)`,
`n Q_n^ell=Omega(n/log n)`. Equations (7), (17)--(19), the top-query
cost, and recovery therefore give the uniform bound

\[
S(n)=O(n Q_n^\ell).
\tag{20}
\]

All loops are finite, the recursion decreases its integer level, and
successor scans strictly increase a represented divisor. After all
oracles are inlined, the only read-only inputs have `O(n)` cells.
A halting fixed-alphabet machine using `S` work cells has at most
`poly(n,S) c^S` internal configurations. Output contents and its
one-way head position cannot affect a transition. Repeating an
internal configuration would prevent termination. Consequently

\[
T(n)\le\mathrm{poly}(n,S)c^S=\exp(O(n Q_n^\ell)).
\tag{21}
\]

The explicit operation accounting gives the same bound: (17) pays
for all piece assignments, the fixed-axis grouping cost is smaller,
its table-query series is dominated by (19), and recovery adds only
`exp(O(n/log n))`. The configuration argument does not treat an
oracle as free or use an exponentially long read-only table tape.

Finally `Q_n -> Q_*<=rho`, and `d^ell=Theta(log n)`. For every fixed
`epsilon>0`, eventually
`log_d(1/Q_n)>=alpha-epsilon/2`. Fixed constants are absorbed by the
remaining positive power of `log n`. Equations (20)--(21) prove (3)
for the same single machine and every epsilon.

## 6. A sharper test for one finite seed

The rational quality (4) is convenient for bounded selection but can
overestimate an individual seed's useful rate. The exact multi-cut
rank-sum quantity is

\[
\mathcal C(S)=\sum_i\left(\prod_{a=1}^d R_i^a\right)^{1/d},
\qquad
\rho(S)=\frac{\log_2\mathcal C(S)}{m t},
\tag{22}
\]

where `R_i^a` are the exact positive ranks of the same nonzero pieces.
Zero tensors are first discarded. This is the fixed-axis counterpart
of the two-cut rank-sum criterion in
[finite-seed amplification](finite-seed-amplification.md) and Section 6
of [the gap-free endpoint](gap-free-seed-score-transfer.md).

One actual finite seed with `rho(S)<1` gives the constant-form bounds

\[
T(n)=\exp\!\left(O\!\left(
 \frac n{(\log n)^{\log_d(1/\rho(S))}}\right)\right),
\qquad
S(n)=O\!\left(
 \frac n{(\log n)^{\log_d(1/\rho(S))}}\right).
\tag{23}
\]

Here the seed is fixed, so constants may depend on its finite tables.
No real-number computation is required by the algorithm.

To prove this, use the exact-rank channels in (10), with cardinalities
`R_i^a` rather than padding each small block separately. For every type,
the multinomial theorem over the nonnegative reals gives

\[
D_\nu\prod_i
 \left(\prod_aR_i^a\right)^{\nu_i/d}
 \le\mathcal C(S)^b.
\tag{24}
\]

The left side is one term of that expansion. The total rounded width
of the `d` type-word channels is consequently at most

\[
d b\log_2\mathcal C(S)+dm u+d
\le dm\rho(S)\lambda+O_S(1).
\tag{25}
\]

The separate caps still apply. With the seed fixed, there are
polynomially many types, their channel numbers are exact mixed-radix
integers of `O_S(lambda)` bits, and the same streamed algorithms have
`2^O_S(lambda)` time and `O_S(lambda)` space. The fixed-axis time proof
and the packed accounting in Section 5 now apply with a constant
contraction `rho(S)` and an `O_S(1)` width error, proving (23).

The contraction is in the required range. For any chosen axis `a`,
rank subadditivity gives `sum_i R_i^a>=2^(m t)`, and positivity of all
other cut ranks gives

\[
\mathcal C(S)^d\ge\sum_i\prod_aR_i^a
 \ge\sum_iR_i^a\ge2^{m t}.
\]

Thus `rho(S)>=1/d`. At the same seed,
`mathcal C(S)<=M_t 2^(w_t/d)`, so `rho(S)<=Q_t` when the old widths
are those of its displayed rank factorizations.

For an exact finite certificate, it is enough to give rational upper
bounds `z_i` with `z_i^d>=prod_a R_i^a` and a positive rational
`u/v` satisfying

\[
\left(\sum_i z_i\right)^v<2^{u m t},
\qquad d u^2<v^2.
\tag{26}
\]

After clearing denominators these are integer comparisons. They prove
`rho(S)<u/v<1/sqrt(d)`, hence an actual factoring exponent above one
half. The certificate must still contain the tensors and all their
exact factorizations. Equation (26) is not a way to certify abstract
rank profiles without a realizing seed.

## 7. Exact controls and remaining construction problem

Run:

```sh
python3 tools/check_multiaxis_type_amplification.py \
  --output evidence/multiaxis-type-amplification.json
```

The checker starts from actual binary tensors and independently
constructs rank factorizations of every required cut. It compares
literal type sums with explicit type-word channels, checks endpoint
sums, includes a nonempty tail, exercises the identity-factor cap,
and checks exact multinomial and selector inequalities. None of the
test seeds is claimed to give an improving quality. Finite controls
verify the construction; the numbered arguments supply the uniform
asymptotic conclusion.

The missing research result remains an actual family whose limiting
quality is below `1/sqrt(d)` for some fixed `d`, or an actual finite
seed below that threshold. Favorable Smith profiles without common
pieces and full factorizations do not meet this requirement.

This ordinary proof adapts the cited repository type and endpoint
arguments to the new multi-axis geometry. It used substantial OpenAI
assistant support. No Lean formalization is asserted.
