# An explicit polynomial divisor oracle with computed benefit certificates

## Result and executable status

For `n>=1`, let

\[
D(n)=\max_{1\le m<2^n}\tau(m),\qquad
M(n)=\min\{m<2^n:\tau(m)=D(n)\},\qquad
q(n)=\lceil\log_2D(n)\rceil.
\]

**Theorem.** The explicit adaptive algorithm below computes `D(n)`,
`M(n)`, and `q(n)` for every input. Its bit runtime is `n^O(1)`.
A version using uniform arithmetic circuits has `O(log(n+2))` working
bits and streamed output.

The algorithm contains **no unspecified numerical Nicolas constants**.
Its stopping test is an exact, computed certificate. The asymptotic
polynomial bound uses two established facts about highly composite
numbers: Nicolas's small-benefit theorem and the uniform gap estimate
proved on the way to his Theorem 5. Those facts are needed for the
runtime analysis, while the algorithm's correctness and termination use
elementary inequalities and its numerical certificate.

The practical implementation
[`adaptive_divisor_preprocessing.py`](../tools/adaptive_divisor_preprocessing.py)
uses ordinary integers, cached primes, and a recursive enumeration.
It satisfies the polynomial-time theorem with polynomial workspace.
The logarithmic-workspace version uses a compact enumeration counter,
recomputed entries, and the uniform iterated-multiplication theorem;
those resource bounds are not attributed to the Python recursion.

This strengthens the existence result in
[polynomial-divisor-preprocessing.md](polynomial-divisor-preprocessing.md)
by eliminating uninstantiated cutoff constants from the algorithm.
The earlier exact
[next-prime exchange oracle](exchange-divisor-preprocessing.md)
remains an independent comparator.

The factoring exponent remains `1/2`. These are auxiliary computations
from an input-length bound, not a factorization of a supplied integer.
No novelty claim is made for classical benefit methods or highly
composite-number algorithms.

## 1. A separable global upper bound

Fix a positive rational `epsilon=a/Q`. For each prime `p`, choose the
largest nonnegative integer `e_p` maximizing

\[
\log(e+1)-\epsilon e\log p.
\]

The marginal gains are strictly decreasing in `e`. Thus the chosen
exponent is the largest `e>=1` satisfying

\[
\left(\frac{e+1}{e}\right)^Q\ge p^a,
\tag{1}
\]

or zero when the first increment is unprofitable. Only primes at most
`2^(1/epsilon)` occur. Define the superior integer and its objective by

\[
N_\epsilon=\prod_p p^{e_p},\qquad
G_\epsilon=\log\tau(N_\epsilon)-\epsilon\log N_\epsilon.
\]

Separability over primes proves, for every positive integer `m`,

\[
\operatorname{ben}_\epsilon(m)
 :=G_\epsilon-\log\tau(m)+\epsilon\log m\ge0.
\tag{2}
\]

This is the classical benefit function. Ties in (1) are included, so
`N_epsilon` is the largest maximizing integer. It is highly composite:
if `m<N_epsilon`, (2) gives
`tau(m)<=tau(N_epsilon)(m/N_epsilon)^epsilon<tau(N_epsilon)`.

Let `X=2^n`. If `N_epsilon<X` and a known feasible integer has divisor
count `D0`, then every integer `m<X` with `tau(m)>=D0` satisfies

\[
\operatorname{ben}_\epsilon(m)
 \le B_\epsilon(D_0):=
 \epsilon\log X-\log D_0+
 \log\tau(N_\epsilon)-\epsilon\log N_\epsilon.
\tag{3}
\]

The right side is computable. It supplies the final stopping certificate.

## 2. A rational superior parameter close to the input bound

For `n>=2`, put

\[
B=4n^2,\qquad L=\lceil\log_2 B\rceil,
\]

and let `Q` be the least power of two at least
`max((n+2)^2,L(L+1))`. Therefore `Q=Theta(n^2)`.
Choose the lower numerator

\[
a_0=\left\lceil\frac Q{L+1}\right\rceil.
\]

Because `Q>=L(L+1)`,
`1/(L+1)<=a_0/Q<=1/L`. The elementary prime bound
`p_j<=4j^2` shows that `N_(a0/Q)` contains the first `n` primes;
hence it is at least `2^n`. At the upper numerator `a=Q`, the largest
superior maximizer is `N_1=2<2^n`.

The product `N_(a/Q)` is nonincreasing in `a`. Binary search on the
integer numerators finds adjacent values with

\[
\boxed{\quad N_\epsilon<X\le N_{\epsilon-1/Q}.\quad}
\tag{4}
\]

The program can stop a product computation as soon as its partial product
reaches `X`. This gives the correct comparison without constructing an
unneeded oversized superior integer.

All parameters stay in a polynomial range. Since
`epsilon>=1/(L+1)`, occurring primes are at most
`2^(L+1)<=16n^2`, and their exponents are less than `2(L+1)`.
Equation (1) is an exact integer-power comparison with polynomially
many output bits. The code normally resolves it by certified logarithmic
intervals and uses the integer-power comparison only if the intervals
overlap. Neither floating point nor a conjectural separation of
transcendental numbers decides a comparison.

## 3. A benefit radius determines explicit height windows

A least maximizer uses consecutive primes with nonincreasing exponents.
Write its conjugate heights as

\[
h_j=\#\{p:v_p(m)\ge j\}.
\]

For a prime `p` and positive level `j`, define the marginal cost

\[
w_j(p)=\epsilon\log p-\log(1+1/j).
\tag{5}
\]

If `v_p(m)>=j` and `w_j(p)>0`, the benefit contains at least the cost
`w_j(p)`: the maximizing exponent omits this unfavorable increment.
If `v_p(m)<j` and `w_j(p)<0`, the benefit contains at least
`-w_j(p)`: the maximizing exponent includes this favorable increment.
Every other prime contributes a nonnegative amount to (2). Thus

\[
\begin{array}{lll}
w_j(p)<-b&\Longrightarrow&v_p(m)\ge j,\\
w_j(p)> b&\Longrightarrow&v_p(m)<j
\end{array}
\qquad\text{whenever }\operatorname{ben}_\epsilon(m)\le b.
\tag{6}
\]

This implication holds for every integer; it does not require a
squarefree exchange or an exponent-asymptotic estimate.

The preceding exchange lemma implies that every highly composite
`m<X` has prime support below `B` and at most

\[
E=\min(n-1,2L-2)
\tag{7}
\]

nonzero height levels. For each `1<=j<=E`, count the primes at most `B`
which (6) forces to be present, and the primes which it still allows.
Use outward-rounded logarithmic intervals: a prime is mandatory only
when its entire interval is below `-b`, and it is excluded only when
its entire interval is above `b`. This gives integer windows

\[
\ell_j\le h_j\le u_j,\qquad0\le\ell_j,u_j\le n-1.
\tag{8}
\]

Since (5) is increasing in `p`, a certified exclusion also excludes
all subsequent primes. The program may stop that prime scan early.
All undecided comparisons are included in the window. Equality at a
threshold never causes a candidate to be omitted.

Enumerate every nonincreasing height profile in (8), rejecting its
reconstructed product when it reaches `X`. A feasible product also
satisfies `sum h_j<=n-1`. The implementation also uses the independently proved
next-prime exchange caps. If `h_1=k` and `q=p_(k+1)`, then

\[
h_j\le\#\{i\le k:p_i^{\lceil j/2\rceil}<q\}.
\tag{9}
\]

These caps hold for all strict divisor records, including the least
maximizer, so they are safe. They are not needed for the asymptotic
candidate-count proof below.

## 4. The adaptive algorithm and its complete correctness proof

Start with the feasible superior integer from (4), and retain it as the
current best. Let `b0` be the reciprocal of the least power of two
strictly greater than `n+1`, so `b0=Theta(1/n)`.

At stage `s=0,1,...`, use the radius `b=2^s b0`.

1. Compute the certified height windows (8).
2. Enumerate the profiles, reject infeasible products, and retain the
   greatest divisor count with the least-integer tie rule.
3. Compute a certified upper bound for `B_epsilon(D0)` in (3).
   Stop if that upper bound is at most `b`; otherwise double `b`.

**Correctness.** Suppose the algorithm stops but its best answer is not
the least maximizer `A`. Then `tau(A)>=D0`, so (3) and the stopping
certificate give `ben_epsilon(A)<=b`. Equations (6)–(9) place `A` among
that stage's enumerated feasible profiles. The objective and tie rule
would have retained `A`, a contradiction.

**Termination without Nicolas's theorem.** At every stage,
`D0>=tau(N_epsilon)`. Hence

\[
B_\epsilon(D_0)
 \le\epsilon\log(X/N_\epsilon)<n\log2<n.
\tag{10}
\]

The fixed rounding error in the upper bound is less than one. Therefore
some stage with `b>=n+1` must terminate. There are `O(log n)` stages.
No unproved estimate about record distribution is required for either
correctness or termination.

The case `n=1` is the direct answer `M=D=1`, `q=0`.

## 5. Why the explicit search is polynomial

The following two established results are used only in this section.
Both are proved in Jean-Louis Nicolas, *Répartition des nombres hautement
composés de Ramanujan*, Canadian Journal of Mathematics **23** (1971),
116–130,
[DOI](https://doi.org/10.4153/CJM-1971-012-6),
[primary PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/4E3B60B2BC31651C55D01FFC9ECB4497/S0008414X00053037a.pdf/repartition_des_nombres_hautement_composes_de_ramanujan.pdf).

- **Small benefit.** If `A` is highly composite and `N_epsilon*` is its
  preceding superior record, then
  `ben_epsilon*(A)<=C p^(-gamma)` for fixed positive constants, where
  `p` is the actual largest prime of `A`. This is Theorem 1, together
  with `x=2^(1/epsilon*)~p` from Proposition 4 and its corollary.
- **Uniform next-record gap.** For fixed `c,C'>0`, every sufficiently
  large highly composite `A` has a later integer `A'` with
  `tau(A')>tau(A)` and
  `A'/A<=1+C'/(log A)^c`. This is the construction in the proof of
  Theorem 5 on printed pages 128–129. Merely knowing the theorem's
  lower bound on the *number* of records would not imply this uniform
  gap; the proof's per-record construction is what is used here.

Let `A=M(n)`. A superior parameter `epsilon*` for its preceding superior
record can be chosen in the mesh interval
`[epsilon-1/Q,epsilon]` from (4). Indeed, the nonincreasing superior
maximizers at the two mesh endpoints lie on opposite sides of `A`.
If `A` is itself superior, choose a parameter at which it maximizes;
its benefit is then zero. Simultaneous prime-exponent ties cause no
problem: all tied maximizers have the same objective `G`.

Put `delta=epsilon-epsilon*`, so `0<=delta<=1/Q`. Since
`G_epsilon*` dominates the objective of `N_epsilon` at `epsilon*`,

\[
\begin{aligned}
\operatorname{ben}_\epsilon(A)
&\le\operatorname{ben}_{\epsilon_*}(A)
  +\delta\log(A/N_\epsilon)\\
&\le Cp^{-\gamma}+\frac{n\log2}{Q}
 =O(p^{-\gamma})+O(1/n).
\end{aligned}
\tag{11}
\]

The next highly composite number is at least `X`; otherwise `A` would
not be the last strict record below `X`. The uniform gap therefore gives

\[
\log(X/A)=O((\log A)^{-c}).
\tag{12}
\]

Only elementary comparisons are needed to express the errors in terms
of `p`. If `k=pi(p)`, then `p<=4k^2` and `A>=2^k`, so
`log A>=c0 sqrt(p)`. Also `p<=4n^2`, so `1/n<=2/sqrt(p)`.
Combining (3), (11), and (12) yields, for some fixed `eta>0`,

\[
B_\epsilon(\tau(A))
=\operatorname{ben}_\epsilon(A)+\epsilon\log(X/A)
=O(p^{-\eta}).
\tag{13}
\]

We may decrease `eta` below `1/4`.

Once `b` exceeds (13) by the small certified rounding margin, the
stage includes `A` and terminates. Since the radius doubles, the final
radius is at most a constant times
`max(b0,B_epsilon(tau(A)),rounding_error)`. Thus

\[
b_{\rm final}=O(p^{-\eta}).
\tag{14}
\]

This statement is a runtime proof, not a bound the executable needs to
know.

To count the windows, let `x=2^(1/epsilon)`. Since
`epsilon>=epsilon*`, `x<=2^(1/epsilon*)=O(p)`. The real center at level
`j` is

\[
\theta_j=x^{\log(1+1/j)/\log2}.
\]

The permitted primes lie in a multiplicative interval with endpoints
`theta_j exp(+-b/epsilon)`, enlarged by the outward arithmetic error.
Choose a fixed `K` with
`log(1+1/(K+1))/log2<eta/2`. For all `j>K`,
`theta_j=O(p^(eta/2))`.

We also have `1/epsilon=O(log n)=O(log p)`. The second equality follows
without a new analytic prime estimate: (12) gives `log A=Theta(n)`;
canonical factorization and the exchange bound give
`log A=O(p(log p)^2)`. Hence `log n=O(log p)`.
Consequently every high-level prime window has length

\[
O\!\left(p^{\eta/2}p^{-\eta}\log p\right)=o(1).
\tag{15}
\]

The arithmetic enlargement is smaller still, by the precision choice in
Section 6. Therefore there are at most a fixed number of possible heights
at each level `j>K`. The first `K` heights have at most `n` possibilities
each, and there are only `E=O(log n)` levels in total. The final
independent window box has size

\[
(n+1)^K C_0^E=n^{O(1)}.
\tag{16}
\]

Earlier stages have narrower windows. Their number is `O(log n)`, and
all arithmetic per profile is polynomial. The Python prefix recursion
visits at most `O(E)` times the full box size; pruning only reduces work.
This proves the polynomial bit-runtime claim for the executable.

The estimates apply for all sufficiently large `p`. Bounded `p` gives
only a bounded set of strict divisor records, by the exchange exponent
cap. Since `D(n)>=n`, such cases occur at only finitely many input lengths.
Their finite running times are absorbed in the fixed polynomial bound.
The algorithm handles them using its computed certificate, with no
knowledge of the exceptional threshold.

## 6. Certified arithmetic and the logarithmic-space version

The code uses precision

\[
F=32\bigl(\operatorname{bitlen}(n+2)+10\bigr).
\]

The integer-only logarithm subroutine is based on the positive atanh
series proved in
[the threshold-center note](polynomial-divisor-preprocessing.md).
For an integer `v`, write `v=2^s r`, `1<=r<2`. Retain only the leading
`F+3` bits of `v`, giving a lower approximation to `r` with error at
most `2^(-F-2)`. With every fixed-point operation rounded down, the
result `ell(v)` satisfies

\[
0\le\log v-\ell(v)\le
(s+1)(8F+5)2^{-F}.
\tag{17}
\]

For `log((j+1)/j)`, use the reduced series directly; its error is at most
`(8F+5)2^(-F)`. Formula (17) remains valid when `v` is a long virtual
product because only its bit length and leading `O(log n)` bits are read.

Every local window comparison and the final upper bound use the signs
of explicit integer interval endpoints. For example, in (3), positive
logarithms are replaced by their upper endpoints and negative logarithms
by their lower endpoints. The resulting upper bound is therefore valid
without a numerical tolerance assumption. Its absolute excess is
`n^(-Omega(1))`, much smaller than `b0`. Local marginal intervals have
smaller errors, sufficient for (15). An interval comparison that remains
undecided is kept in the candidate cover; it is never rounded into an
exclusion.

For (1), undecided intervals are resolved by comparing
`(e+1)^Q` and `e^Q p^a` exactly. All these integers have polynomially
many bits. The special exact tie at `epsilon=1,p=2,e=1` is treated in
this way by the implementation.

A logarithmic-workspace version replaces the Python arrays and recursive
prefix products by the following operations.

- Enumerate the Cartesian window box with a packed mixed-radix counter.
  Equation (16) gives `O(log n)` bits for the entire counter; rounding
  each coordinate width adds only `O(E)=O(log n)` bits.
- Recompute window endpoints, primes, and an exponent from that counter
  when needed. Retain the best descriptor together with its stage index,
  so a descriptor from an earlier, narrower box remains reconstructible.
- Compare products by their regenerated bits, rather than retaining the
  products. Use the same method for exact power comparisons, bit lengths,
  and the leading-bit logarithm routine.

The last step uses Hesse, Allender, and Barrington's uniform iterated
multiplication theorem, JCSS **65** (2002), 695–716,
[primary paper](https://people.cs.rutgers.edu/~allender/papers/division.pdf),
Corollary 6.3. Polynomial-sized constant-depth threshold circuits can be
evaluated with recomputation in polynomial time and logarithmic workspace.
Their virtual factor inputs are generated in the same resources. The
counter, a retained descriptor, the stage and mesh indices, and the
fixed-depth circuit stack together use `O(log n)` bits. Prime scanning
and trial division also fit that bound.

After certification, regenerate the winning product and divisor count
and stream their bits. Computing `q(n)` adds only a bit-length and
power-of-two test.

## 7. Exact verification and finite performance

Run the standard audit with

```sh
python3 tools/check_adaptive_divisor_preprocessing.py
```

Add `--large` to include the recorded 512-, 1024-, and 2048-bit cases.
The audit uses several independent controls:

- An ordinary divisor sieve covers every integer below `2^16`.
- Unrestricted exponent enumeration and the prior exchange oracle check
  overlapping ranges of input lengths.
- Exact rational logarithm intervals independently verify the superior
  exponent inequalities, the adjacent mesh bracket, and the terminal
  benefit certificate.
- A second enumeration proceeds in prime-exponent order, reconstructing
  the terminal height windows with rational intervals. The main oracle
  instead proceeds in height order. Certified partial-benefit lower
  bounds exclude only branches that cannot improve the terminal answer.
  Both enumerations must return the same least witness and divisor count.
- A finite grid compares the certified logarithmic increment test against
  literal integer-power comparisons, including an exact tie.

The full recorded audit checks 41 input lengths, 1,402 superior marginal
inequalities, and 686 literal power comparisons. It independently replays
the terminal benefit search through the 2,048-bit bound.

The JSON evidence records exact counts and outputs. Timings are exploratory
and do not certify the asymptotic theorem. Some finite sizes have more
visited profiles than the exchange-only oracle; no uniform finite-size
speedup is claimed.

At 128 bits the adaptive search visits 4,466 feasible profiles across its
stages, compared with 14,956 in the exchange search, and returns the same
maximum divisor count `318,504,960`. At 256 bits the counts are 3,321
versus 137,450, with exact maximum `121,762,322,841,600`. The distinction
between feasible profiles, recursive nodes, and the unpruned Cartesian
box is preserved in the evidence.
