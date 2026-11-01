# Complete singular pole traces for one-variable circuits

## Constructive theorem

Let `q(X)` be a one-variable integer arithmetic circuit. Its formal degree
can be exponential in its written size. More generally, allow inverse
gates whose inputs are odd at both integer parities. These describe
globally regular rational functions on the dyadic integers.

Let `I=[L,U)` be a binary-described integer interval, possibly with
negative endpoints. Given an integer target `N`, a scale `K >= 0` and
precision `p >= 1`, the complete normalized trace

\[
 \mathcal T_{q,I}(N,K)
 =2^K\sum_{i\in I}\frac1{N+2^K w-q(i)},
 \qquad w^2+w+1=0,
\tag{1}
\]

can be computed modulo `2^p` in deterministic polynomial bit time in the
written circuit, the parameter and endpoint lengths, `K`, and `p`.
Neither the interval nor a polynomial of the original formal degree is
expanded.

The algorithm also returns the **exact number** of interval indices
satisfying

\[
 q(i)\equiv N\pmod{2^K},
\]

represented by polynomially many disjoint dyadic residue balls with
exact interval counts and an explicit representative in each nonempty
ball.

This completes the singular normalization for one-variable circuit
images. It is stronger than the unit-center trace in
[circuit-box-products.md](circuit-box-products.md). It does not evaluate
the two-variable rectangle with `q(a,b)=ab`, so it does not improve the
complete factoring exponent.

## Existing ingredients and the new composition

The repository already has a degree-budgeted valuation tree for
explicitly written one-variable denominators in
[pole-cauchy-rational-weights.md](pole-cauchy-rational-weights.md).
General univariate root counting modulo prime powers is also established
in the literature; see Dwivedi, Mittal and Saxena,
[CCC 2019](https://doi.org/10.4230/LIPIcs.CCC.2019.15).

The contribution here is a different input contract and a complete
normalized trace. Finite-precision function compression first replaces
an unrestricted-degree circuit by a polynomial of degree `O(K+p)`.
The valuation tree then has polynomial size even though the original
degree and interval size can be enormous. Finally, all leaf resolvents
have globally odd norms and can be summed by the
[regular-circuit average engine](circuit-box-averages.md).

No unknown denominator valuation is assumed. The quadratic displacement
in (1) gives the required valuation bound automatically.

## 1. Pointwise precision needed by a normalized pole

Work in the unramified quadratic extension with residue field `F_4`.
For every dyadic integer `u`,

\[
 v_2(2^K w-u)=\min(K,v_2(u)).
\tag{2}
\]

If the two valuations agree, the remaining denominator has residue
`w-u/2^K`, which is nonzero in `F_4`. Thus the normalized summand

\[
 F_K(u)=\frac{2^K}{2^K w-u}
\]

is always integral. Moreover, if `u=v mod 2^(K+p)`, then

\[
 F_K(u)-F_K(v)
 =\frac{2^K(u-v)}{(2^K w-u)(2^K w-v)}
 \equiv0\pmod{2^p}.
\tag{3}
\]

The denominator valuation is at most `2K`, while the numerator valuation
is at least `2K+p`. Summing the errors over an arbitrarily long finite
interval does not lose dyadic precision.

The `K+p` precision is sharp for a uniform pointwise replacement. The
values zero and `2^(K+p-1)` agree one bit earlier, but their normalized
poles differ by a quantity of exact valuation `p-1`.

## 2. Remove the formal degree before partitioning

Set `P=K+p`, and choose the least `D` with

\[
 v_2(D!)\ge P.
\]

Then `D <= 2P`, and the falling factorial

\[
 (X)_D=X(X-1)\cdots(X-D+1)
\]

vanishes at every integer modulo `2^P`. Evaluate the written circuit for
`q(X)-N` in

\[
 (\mathbb Z/2^P\mathbb Z)[X]/((X)_D).
\]

This produces an explicitly written polynomial `f`, of degree below
`D`, such that

\[
 f(i)\equiv q(i)-N\pmod{2^P}
 \quad\text{for every integer }i.
\tag{4}
\]

All coefficients have at most `P` bits. The reduction uses polynomial
work in the written circuit size and `P`, irrespective of the circuit's
formal degree.

For a globally odd inverse input `u`, Euler powering
`u^(2^(P-1)-1)` gives the correct inverse value modulo `2^P` at every
integer. Binary powering uses `O(P)` quotient multiplications. This is
a pointwise replacement, which is exactly what (3) requires. It need
not be an inverse in the coefficient quotient itself.

After this step, the rest of the computation concerns the exact
polynomial represented by these coefficient residues. Equation (3)
has already justified that replacement for the full sum. In particular,
the later use of normalized averages does not divide a pointwise error
by the interval cardinality.

## 3. A finite-precision valuation tree

A state consists of a residue ball `i=a+2^t y` and a polynomial identity
at the available precision,

\[
 f(a+2^t y)\equiv2^s g(y)\pmod{2^P}.
\tag{5}
\]

The coefficients of `g` are represented modulo `2^(P-s)`. Remove their
common power of two. If the resulting content makes `s >= K`, stop at a
high-valuation leaf, written with scale exactly `2^K`:

\[
 f(a+2^t y)\equiv2^K u(y)\pmod{2^P}.
\tag{6}
\]

Otherwise `g` is primitive modulo two. Inspect its two possible input
parities.

- If `g(epsilon)` is odd, the class `y=epsilon+2z` is a unit leaf of
  exact valuation `s`.
- If `g(epsilon)` is even, substitute `y=epsilon+2z`, remove the new
  coefficient content, and continue.

If both parity values are odd, the current state itself is a globally
unit leaf. The average engine accepts this condition directly, so the
extra analytic refinement required by the older Taylor-series
implementation is unnecessary here.

The interval intersection is computed exactly. For a residue ball
`i=a mod 2^t`, its local interval is

\[
 \left[
  \left\lceil\frac{L-a}{2^t}\right\rceil,
  \left\lceil\frac{U-a}{2^t}\right\rceil
 \right).
\tag{7}
\]

Empty intersections are discarded. There is no search through the
represented integers.

### Degree budget

Let `epsilon` have multiplicity `m` as a root of `g mod 2`. In

\[
 g(\epsilon+Y)=\sum_j c_jY^j,
\]

the coefficient `c_m` is odd. The coefficient of degree `m` in
`g(epsilon+2Y)` has valuation exactly `m`. Therefore the removed content
is at most `m`. After that removal, every coefficient of degree larger
than `m` is still even. The child's degree modulo two is at most `m`.

The multiplicities of the roots zero and one sum to at most the parent's
degree modulo two. Thus active degree does not increase. A terminal
constant state can retain a unit of its parent's budget without
producing any children.

Each continuing edge increases the removed valuation by at least one.
Its depth is at most `K`. If the compressed polynomial has degree `d`,
there are at most

\[
 \max(1,d)(K+1)
\]

processed primitive states and only polynomially many unit leaves.
The high-valuation leaves number at most `max(1,d)`. They represent
exactly the roots modulo `2^K` and give the promised compact root count.

Every affine substitution is reduced modulo the currently available
coefficient precision. All state polynomial degrees stay at most `d`,
and every coefficient stays within `P` bits. This avoids expansion of
large exact shifted coefficients.

## 4. The complete resolvent on each leaf

On a unit leaf of valuation `v<K`, write `f=2^v u` and set
`c=2^(K-v)`. On a high leaf, write `f=2^K u` and set `c=1`. In both cases,
the required summand becomes

\[
 \frac{c}{cw-u}
 =-\frac{c(u+c)}{u^2+cu+c^2}
  -\frac{c^2}{u^2+cu+c^2}\,w.
\tag{8}
\]

The scalar denominator is globally odd on the local integer variable.
For a unit leaf, `u` is odd and `c` is even. For a high leaf, `c=1` and
`u^2+u+1` is odd for either input parity.

Thus both coordinates in (8) are globally regular rational circuits.
Their written sizes are polynomial in the compressed degree and
precision. Their complete interval sums can be evaluated by normalized
averages. For a local interval of length `m`,

\[
 \sum h(i)=\frac{m}{2}\Lambda_m(h),\qquad
 \Lambda_m(h)=\frac2m\sum h(i).
\]

Request the integral normalized average modulo `2^(p+1)`, multiply by
the known `m`, and divide the resulting even integer by two. This gives
the ordinary sum modulo `2^p` even when `m` has enormous dyadic valuation.

If `K-v >= p`, the entire unit-leaf contribution vanishes modulo `2^p`
and is skipped. Otherwise (8) includes every term on that leaf, not a
selection of low-order jets at the original singular center.

Combining all leaf sums proves (1). The algorithm uses only unit
inverses after the valuation partition. Its scalar output has two
`p`-bit coordinates in the basis `1,w`.

## 5. Exact zero counts under a height bound

For an integer-valued input circuit, suppose an explicit bound guarantees

\[
 |q(i)-N|<2^B\quad(i\in I).
\]

Taking `K >= B` makes the returned congruence-root count the exact number
of zeros `q(i)=N`. Every nonempty root ball provides such an integer
index directly. Taking `K=B+p` also gives

\[
 \mathcal T_{q,I}(N,B+p)
 \equiv \frac{\#\{i\in I:q(i)=N\}}w\pmod{2^p}.
\]

For regular rational input, an analogous statement needs a bound on the
nonzero values' dyadic valuations, or on their integer numerators. A real
absolute-value bound alone does not control dyadic valuation when odd
denominators can be large.

As a restricted factoring application, let a one-variable circuit curve
give bounded candidates `a(i),b(i)` and `q(i)=a(i)b(i)`. The oracle finds
all exact hits on that curve under the stated height bound. For example,
`q(i)=(i+2)(2i+1)` on `[0,21)` hits 77 only at `i=5`, yielding 7 and 11.
This does not provide a compact curve covering every possible factor
pair.

## 6. Why the two-variable factor box remains open

The degree budget is specifically univariate. A root in one variable
has a multiplicity that bounds its descendants. A positive-dimensional
root set in two variables has no corresponding conserved budget.

For example, `XY-N` with odd `N` has `2^(K-1)` roots modulo `2^K` among
odd residue pairs. After the first parity restriction, its smooth
congruence curve has two continuing parity branches at each level.
A residue-ball partition can therefore have exponentially many
branches although the polynomial has degree two.

The theorem does not convert the independent rectangle coordinates into
one ordinary integer-circuit variable. Binary packing would require
digit extraction or integer quotients, which are outside the supplied
regular-circuit contract. Nor does an inner one-variable trace
automatically have a small admissible circuit in the outer variable.

Consequently, the complete joint factor count for `N-ab`, the unknown
singular product normalization on arbitrary rectangles, and a better
factoring exponent remain open.

## 7. Implementation and evidence

Files:

- [Compiler](../tools/circuit_box_pole_trace.py).
- [Independent checker](../tools/check_circuit_box_pole_trace.py).
- [Deterministic receipt](../evidence/circuit-box-pole-trace.json).

Replay:

```sh
python3 tools/check_circuit_box_pole_trace.py \
  --output evidence/circuit-box-pole-trace.json
```

The receipt records 1,200 direct circuit-compression comparisons,
400 literal valuation-partition identities, 360 complete trace
comparisons against exact rational quadratic arithmetic, and 360 exact
congruence-root counts. Inputs include signed intervals, empty intervals,
constant circuits, repeated roots and globally regular inverse gates.
Nine controls distinguish `K+p` pointwise precision from one fewer bit.

Three large controls use `K=80`, `p=12`, and compression precision 92.
The universal function quotient has degree 96.

| Input and interval | Compressed degree | Processed states | Valuation leaves | Trace in basis `1,w` |
| --- | ---: | ---: | ---: | --- |
| `q(i)=i`, `I=[0,2^160+13)`, target 0 | 1 | 81 | 81 | `(4095,4095)` |
| `q(i)=i^(2^500)`, `I=[7,7+3*2^300+17)`, target 0 | 95 | 8 | 5 | `(4088,4088)` |
| `q(i)=1/(1+2i^(2^500))`, same interval, target 1 | 95 | 8 | 5 | `(4088,4088)` |

The linear control has exactly `2^80+1` congruence roots. Its complete
period blocks have multiplicity divisible by `2^12`, so only the final
13 indices affect the trace; this gives an independent closed-form
reference without expanding a period.

In the two large-degree controls, repeated squaring makes the polynomial
zero on even inputs and one on odd inputs modulo `2^92`. Their singular
trace therefore counts the even indices. The rational control has the
same high-valuation class after subtracting target one. These references
use only parity counts, independently of the compressed polynomial or
valuation tree.

The finite checks supplement the proof. They do not extend the theorem
to two-variable factor rectangles or change the complete factoring
score.
