# A uniform transfer from any fixed number of cube axes

The [three-axis transfer](cube-crossing-transfer.md) extends to every
fixed integer dimension `d >= 2`. This supplies the algorithmic link that
the [multi-axis Gaussian profile budgets](multiaxis-gaussian-profile-floor.md)
deliberately left open. The dimension, coefficient field, and local
algorithms below are fixed independently of the factoring input.

**Conditional transfer.** Suppose the antipodal equality tensor on the
`2^d` vertices of a cube has a polynomial-piece family whose same pieces
have exact factorizations at all `d` axis cuts. Suppose each separate cut
has width at most `2^(d-1) lambda`, the sum of the `d` widths is at most
`d 2^(d-1) rho lambda + O(log(lambda+2))`, and all piece and factor entries
can be computed in `2^O(lambda)` bit operations. Then deterministic
factoring has time

\[
 \exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right),
 \qquad \alpha=\log_d(1/\rho).
\]

The relevant range is `1/d <= rho < 1`. Improving the proved exponent
`1/2` requires `rho < 1/sqrt(d)`. In particular, a realizing nine-axis
family with `rho < 1/3` would improve the score. The existing favorable
nine-axis **abstract profiles do not supply this family**. This note
constructs the transfer, not an improved local family or factoring bound.

## 1. Exact local hypothesis

Fix a finite field of characteristic two, and let `Q={0,1}^lambda`.
For `v in {0,1}^d`, put `bar v=v xor (1,...,1)`. The endpoint is

\[
 E_\lambda((x_v)_v)
   =\prod_{\{v,\bar v\}}[x_v=x_{\bar v}].                 \tag{1}
\]

There are `2^(d-1)` factors. Assume one fixed constant `C` supplies, for
each positive integer `lambda`, an exact decomposition
`E_lambda=sum_i T_i` with at most `(lambda+1)^C` pieces. Each *same* piece
has a factorization across every coordinate cut, with integer widths
`r_i^a` satisfying

\[
 0\le r_i^a\le 2^{d-1}\lambda,
 \qquad
 \sum_{a=1}^d r_i^a
   \le d2^{d-1}\rho\lambda+C\log_2(\lambda+2).            \tag{2}
\]

A width `r` denotes a channel indexed by `r` bits; unused states may be
zero-padded. Uniform total algorithms enumerate the piece metadata and
evaluate every factor entry in `2^(C lambda)` bit operations. Increasing
`C` covers the fixed number of cuts and legs.

The separate cap in (2) follows from an exact rank factorization when
one is available: each side of an axis cut has `2^(d-1)` physical
letters. It is needed in the geometric cost proof below. Rank
subadditivity at any one cut forces `rho >= 1/d` for a family with
polynomially many pieces. This is an interface restriction, not a
factoring lower bound.

## 2. Partition all generic edges into antipodal cubes

Place a clique of order `P=g^d` on the vertices `[g]^d`. An edge is
generic when its endpoints differ in every coordinate. Its two values
in each coordinate determine one coordinate cube. Each cube has exactly
`2^(d-1)` antipodal edges, and the cubes partition the generic edges.
Consequently

\[
 B=\binom g2^d,
 \qquad E_{\rm generic}=2^{d-1}B,
 \qquad E_{\rm other}=\binom{g^d}{2}-2^{d-1}B.            \tag{3}
\]

Replace each cube's antipodal wires by (1), then enumerate one piece
index for every cube. For one complete assignment, sum its channel
widths along each axis and select an axis with smallest sum. This choice
may depend on the assignment and does not alter its value.

Merge the `g^(d-1)` old vertices in each selected coordinate plane into
one vertex. A cube factor is local to one of its two planes and its
channel joins those planes. Every other old edge is either internal to
a plane or still joins two planes. The result is an ordinary weighted
clique on `g` vertices, not a hypergraph contraction oracle.

For each pair of planes, the number of meeting cubes and the number of
other edges between them are

\[
 \binom g2^{d-1},
 \qquad g^{2d-2}-[g(g-1)]^{d-1}.                         \tag{4}
\]

Give every original edge `lambda` bits and put `V=g^(2d) lambda`.
The separate cap in (2) gives the exact maximum child edge width

\[
 w_{uv}\le 2^{d-1}\lambda\binom g2^{d-1}
 +\bigl(g^{2d-2}-[g(g-1)]^{d-1}\bigr)\lambda
 =\frac V{g^2}.                                        \tag{5}
\]

Let `U=2 sum_(u<v) w_uv` be total oriented width. The smallest of the
`d` cube-width sums is at most their average. Hence

\[
\begin{aligned}
 U&\le 2^d\rho\lambda B+\frac{2C}{d}B\log_2(\lambda+2)\\
   &\quad+g(g-1)\bigl(g^{2d-2}-[g(g-1)]^{d-1}\bigr)\lambda,\\
 \frac UV
 &\le\left(1-\frac1g\right)
 -(1-\rho)\left(1-\frac1g\right)^d
 +\frac{C}{d2^{d-1}}\frac{\log_2(\lambda+2)}\lambda,\\
 &\le\rho+\frac{d-1}{g}
 +\frac{C}{d2^{d-1}}\frac{\log_2(\lambda+2)}\lambda.
                                                               \tag{6}
\end{aligned}
\]

The last line uses `(1-u)^d >= 1-du` and `0 <= rho <= 1`.
Independently, (5) implies `U <= (1-1/g)V`. No asymptotic expansion in
the fixed dimension is used.

## 3. Regrouping and actual table queries

The equitable grouping lemma from the three-axis proof depends only on
the weighted `g`-vertex graph, its total oriented width, and (5). If `K`
divides `g` and

\[
 \frac gK\ge16\zeta^{-2}(2\ln K+1),                     \tag{7}
\]

it gives a partition into `K` equal groups for which every grouped edge
has width at most `((g/(g-1))U+zeta V)/K^2`. Exhaustively scanning
balanced maps finds one in `2^O(g log K)` time, times a polynomial in
the explicit width data. The search uses integer arithmetic. Padding
every grouped edge to its common integer upper bound `L` gives

\[
 K^2L\le\left(\rho+\zeta+O_d(1/g)
 +O_C\!\left(\frac{\log(\lambda+2)}\lambda\right)\right)V+K^2.
                                                               \tag{8}
\]

One plane contains `g^(d-1)` old vertices, each with at most `g^d-1`
incident labels. Summing all labels for one plane-table query therefore
costs

\[
 2^{O_d(V/g+C\lambda+\log(g+\lambda))}                   \tag{9}
\]

old entry queries and bit operations. This includes internal edges and
every cube factor. For fixed `d` their number is polynomial in `g`.

A grouped vertex has `m=g/K` planes. Enumerate the channels between its
planes, whose total width is at most

\[
 \binom m2 V/g^2\le V/(2K^2).                            \tag{10}
\]

For each assignment, query the plane factors separately and multiply.
One entry-oracle composition thus has multiplicative cost

\[
 2^{O_{d,C}(V/g+V/K^2+\lambda+\log(g+\lambda))}.          \tag{11}
\]

This is a concrete enumeration algorithm with charged old table
queries. Generating the piece list, and computing the width metadata
for a complete assignment, takes `2^O(lambda)` time times a polynomial
in the number of cubes. It is included in the `lambda` and logarithmic
terms of (11). No inexpensive plane-table or grouped-table oracle is
assumed.

## 4. A uniform recursion with growth base d

Choose an integer `a` with `2^a > 4d`, then a sufficiently large fixed
integer `x0 >= a`. Use

\[
 K_0=2^{x_0},\qquad g_j=2^{a(j+1)+x_0}K_j,
 \qquad K_{j+1}=g_j^d,\qquad\zeta_j=2^{-j-8}.             \tag{12}
\]

All sizes are powers of two and `K_j` divides `g_j`. Writing
`L_j=log_2 K_j` gives

\[
 L_j=d^j\left[x_0+
 d\sum_{i=0}^{j-1}\frac{a(i+1)+x_0}{d^{i+1}}\right]
 =\Theta_d(d^j).                                        \tag{13}
\]

The convergent positive series proves both bounds. Moreover
`a(j+1)+x0 <= 2L_j`, hence `K_(j+1) <= K_j^(3d)`. The left side of
(7) grows proportionally to `2^(aj)`; its right side is
`O_d((4d)^j)`. Enlarging the one fixed `x0` makes (7) true at every
level. These comparisons prove the schedule and are not real-number
instructions in the machine.

Put `delta=1/(100d)`. For an `n`-bit input choose the largest `ell` with
`P=K_ell <= n^delta`, using the equivalent integer test
`K_ell^(100d) <= n`. The finitely many smaller inputs use trial division.
The growth bound implies

\[
 n^{\delta/(3d)}<P\le n^\delta,
 \quad\ell=O_d(\log\log n),
 \quad\rho^\ell=\Theta_{d,\rho}((\log n)^{-\log_d(1/\rho)}).
                                                               \tag{14}
\]

The existing explicit or polynomial-query top constructor supplies a
`P`-vertex factoring network of volume `O(n)`. It applies at these
smaller values of `P`; its costs remain at most `2^O(n/P)`. Use the
fixed-characteristic weighted witness recovery as in the three-axis
proof.

Stop a branch at volume at most `ceil(sqrt(n))`, or at level zero, and
contract by direct edge-label enumeration. Every other call has order
at most `P`, so

\[
 \lambda\ge V/P^2>n^{1/2-2\delta}.                      \tag{15}
\]

There is also a crude volume bound before any logarithmic-error
estimate is used. The separate cap gives `U <= (1-1/g)V`; the grouping
bound therefore gives `Vchild <= (1+zeta_j)V+K_j^2`. At active calls,
`K_j^2/V <= n^(2 delta-1/2)`. The sum of `zeta_j` is finite and there
are `O(log log n)` levels, so multiplying these crude inequalities
first gives `V_j=O(n)` on every visited path. Consequently
`lambda=O(n)` and `log(lambda+2)=O(log(n+2))` without using the refined
contraction estimate circularly.

Using (8), (12), and (15), along any active path

\[
 V_j\le(\rho+O_{d,C}(\zeta_j+1/g_j)+\varepsilon_n)V_{j+1},
 \quad\varepsilon_n=O_{d,C}(n^{-1/2+2\delta}\log n).      \tag{16}
\]

The rounding cost `K_j^2/V_(j+1)` is included. Both error series in `j`
are summable and `epsilon_n log log n` tends to zero. Multiplying the
bounds gives `V_j=O_(d,C,rho)(rho^(ell-j) n)` until stopping.

The complete remaining costs are:

| Operation | Logarithm of cost |
| --- | --- |
| Explicit top input | `O(n/P)` |
| Piece assignments at one call | `O_C(P^2 log(n+2))` |
| Exhaustive equitable grouping | `O_d(P^(1/d) log(P+2))` |
| One entry-oracle composition | The exponent in (11) |
| Terminal contraction | `O(V_0)` or `O(sqrt(n))` |

The sum of piece and grouping exponents on a path is
`O_(d,C)(n^(2 delta) log n log log n)`, smaller than
`n/(log n)^alpha` for every fixed `alpha` here. In (11) use
`lambda=V/g_j^(2d)` and the convergent series

\[
 \sum_{j\ge0}\rho^{-j}
 (g_j^{-1}+K_j^{-2}+g_j^{-2d})<\infty.                    \tag{17}
\]

All entry-query exponents on a path sum to `O(n rho^ell)`. Branching
multiplies work, so its logarithm adds the charged enumeration terms.
The top-input and terminal costs are smaller. Recursion metadata have
polynomial bit length in the input and displayed table descriptions.
This proves a uniform clique evaluator in `exp(O(n rho^ell))` bit time.

Existing weighted witness recovery adds `2^O(sqrt(n)) poly(n)`
coefficient and factor-recovery cost. Identities over the fixed field
remain valid over its explicit commutative coefficient algebras; their
bit arithmetic is included in this overhead. Equation (14) proves the
asserted factoring bound. One machine contains the fixed local and
constructor algorithms; the input only selects an integer level.

## 5. More than polynomially many pieces

The [two-axis relaxed transfer](subexponential-piece-transfer.md)
already permits a larger piece set and a larger width slack. The same
relaxation holds at every fixed `d`. This section extends that existing
principle; it does not introduce the larger-piece idea itself.

Keep the exact endpoint, all separate width caps, and the uniform local
algorithms. For fixed `gamma>0`, replace the piece and total-width
hypotheses by

\[
\begin{aligned}
 \log_2\max(2,|I_\lambda|)
   &\le\frac{C\lambda}{(\log_2(2\lambda))^\gamma},\\
 \sum_{a=1}^d r_i^a
   &\le d2^{d-1}\rho\lambda+
      \frac{C\lambda}{\log_2\log_2(\lambda+4)}.
\end{aligned}                                                   \tag{18}
\]

Indices use at most `C lambda` bits; complete enumeration of the piece
set and all requested local computations are included in the
`2^(C lambda)` bound. Then the factoring time is

\[
 \exp\!\left(O\!\left(
   \frac n{(\log n)^\alpha}+\frac n{(\log n)^\gamma}
 \right)\right),\qquad \alpha=\log_d(1/\rho),             \tag{19}
\]

and its score is `min(alpha,gamma)`. There is no additional
`log log n` loss when `gamma=alpha`.

Indeed, (15) makes the new per-step relative width error
`O(1/log log n)`. The number of levels is `O(log log n)`, so the
multiplicative error product in (16) is still bounded by a constant.
In particular, on every visited path

\[
 V_j\le D_0\rho^{\ell-j}n,
 \qquad\sum_j V_j\le\frac{D_0}{1-\rho}n                 \tag{20}
\]

for one fixed `D0`. Early levels may have one-step bounds greater than
one; the bounded product, not pointwise contraction, proves (20).

Since `B <= g^(2d)/2^d=V/(2^d lambda)`, the logarithm of the number of
piece assignments at a call is at most

\[
 B\log_2|I_\lambda|
   \le\frac{CV}{2^d(\log_2(2\lambda))^\gamma}
   =O_{d,C,\gamma}(V/(\log n)^\gamma).                  \tag{21}
\]

Sum (21) using (20), rather than multiplying a top-volume bound by
the depth. This gives `O(n/(log n)^gamma)`. Complete local generation
is still covered by the `C lambda` term of (11). The query series,
top construction, termination, and factor recovery are unchanged.
This proves (19).

For example, a hypothetical nine-axis family with `rho=8/25`, piece
entropy `O(lambda/(log lambda)^(3/5))`, and the width slack in (18)
would give score `log_9(25/8)>1/2`. The inequality is exact because
`25/8>3`. Such a family is not constructed here. It illustrates a
strictly broader sufficient regime than polynomially many pieces at
the higher-axis target.

Every fixed bound `log |I_lambda|=O(lambda^beta)`, `beta<1`, fits (18)
for every fixed `gamma`. The same is true of an additive width slack
`O(lambda^beta)`. Thus stretched-exponential piece counts preserve the
entire `rho`-derived rate in every fixed dimension. A bare assertion
`log |I_lambda|=o(lambda)` does not by itself bound the cost of literal
piece-assignment enumeration by a fixed logarithmic saving; the
`lambda/log log lambda` entropy example in the two-axis note applies
to that direct accounting here as well. This is not a restriction on
type regrouping or other algorithms that avoid the literal expansion.

## 6. Why many axes alone give no improvement

Expanding (1) into supported basis tensors makes every cut width zero,
but has `2^(2^(d-1) lambda)` pieces. Its first assignment exponent is

\[
 2^{d-1}\lambda B=\frac12(1-1/g)^dV.                     \tag{22}
\]

This is linear in parent volume for fixed `d` and large `g`. The
polynomial-piece hypothesis fails. Likewise, abstract Smith orders do
not supply common coefficient pieces, cut factorizations, or entry
algorithms.

The new result is a complete transfer at arbitrary fixed dimension.
The proved factoring score remains `1/2`. The next constructive target
is an actual family satisfying (1)--(2) below `1/sqrt(d)`.

## 7. Exact finite controls

```sh
python3 tools/check_multiaxis_cube_transfer.py \
  --output evidence/multiaxis-cube-transfer.json
```

The checker enumerates literal clique edges, partitions generic edges
by coordinate cubes, verifies all plane-pair counts, and checks the
width and volume formulas on varied abstract channel assignments. It
also verifies integer schedule identities and the exponential-piece
negative control. These controls concern geometry and cost arithmetic;
no unconstructed local tensor is treated as a certificate.

This ordinary proof extends the repository's three-axis construction.
It was developed with substantial OpenAI assistant support and has not
been Lean-formalized or independently reviewed.
