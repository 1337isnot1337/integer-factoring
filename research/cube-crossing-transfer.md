# A three-axis cube interface for deterministic factoring

## Result and scope

An actual three-dimensional grid reduction can replace the hypothetical
orientation count in [the earlier directions note](post-half-directions.md).
Generic edges of a clique on `[g]^3` split into disjoint groups of four
antipodal edges of a coordinate cube. Factoring one eight-leg tensor across
one of three axis cuts then produces a clique on the `g` parallel planes.

The construction below includes equitable grouping, table-query costs,
piece enumeration, and a uniform recursion. Its conclusion is conditional:

> A uniform eight-leg crossing family with polynomially many pieces,
> separate widths at most `4 lambda`, and total three-cut width at most
> `12 rho lambda + O(log lambda)` gives deterministic factoring time
> `exp(O(n/(log n)^alpha))`, where `alpha = log_3(1/rho)`.

The desired local target is therefore `rho < 1/sqrt(3)`. For example,
`rho=1/2` would give `alpha=log_3(2)>1/2`. **No such eight-leg family is
constructed here.** This is a new sufficient interface, not a factoring
improvement or a claim that three directions alone beat the half frontier.

The [checker](../tools/check_cube_crossing_transfer.py) verifies exact
edge coverage, all three axis widths, the equitable-partition mean and
swap estimates on finite examples, and the recursion arithmetic. It also
retains a deliberately unsuitable pure-tensor decomposition whose piece
enumeration is exponential in `lambda`.

## 1. The required local family

Let `Q={0,1}^lambda`. Index eight physical legs by
`v in {0,1}^3`, and write `bar v=v xor (1,1,1)`. The endpoint is

\[
E_\lambda((x_v)_v)
 =\prod_{\{v,\bar v\}}[x_v=x_{\bar v}],
\tag{1}
\]

where the product has four factors. For each axis `a=1,2,3`, its two
halves are the four legs with `v_a=0` and the four with `v_a=1`.

Fix a finite coefficient field of characteristic two. Suppose fixed
constants `C` and `1/3 <= rho < 1` give, for every positive `lambda`:

1. An index set `I_lambda` of size at most `(lambda+1)^C` and an exact
   tensor identity `E_lambda=sum_i T_i`.
2. Three exact factorizations of each *same* tensor `T_i`, one across
   each axis cut, with channel widths `r_i^a` satisfying

   \[
   0\le r_i^a\le4\lambda,
   \qquad
   r_i^1+r_i^2+r_i^3
     \le12\rho\lambda+C\log_2(\lambda+2).
   \tag{2}
   \]

3. Uniform total algorithms that enumerate the indices and widths and
   evaluate every factor entry in `2^(C lambda)` bit operations.
   The coefficient field and the algorithms are fixed independently of
   the factoring input.

A width `r` means an auxiliary index in `{0,1}^r`. Padding unused states
with zeros is allowed. Factoring different tensors at the different
cuts does not satisfy this interface. Nor does an unevaluated existence
statement for the pieces.

The lower limit `rho>=1/3` is necessary for this particular interface.
At any axis, (1) has matrix rank `2^(4 lambda)`. Rank subadditivity and
the polynomial piece count force some piece to have width at least
`4 lambda-O(log lambda)` at that axis. Equation (2) then gives
`rho>=1/3` in the limit. This is only an interface restriction.

## 2. Exact decomposition of the clique edges

Put `P=g^3` and identify the vertices with triples in `[g]^3`. Call an
edge generic when its endpoints differ in every coordinate. Every
generic edge determines exactly one coordinate cube: for each
coordinate, retain its two endpoint values. The eight cube vertices
have four antipodal pairs, and all four pairs are generic edges.

Consequently the generic edges are partitioned by

\[
B=\binom g2^3
\tag{3}
\]

cubes, with four edges per cube. The remaining number of edges is

\[
\binom{g^3}{2}-4B=\frac32g^4(g-1).
\tag{4}
\]

Replace the four equality wires of each cube by (1). Enumerate one
piece index for every cube. For a fixed assignment, the three available
factorizations give three different exact contractions of that same
summand. Choose the axis with the smallest sum of channel widths.
The selected axis may depend on the complete piece assignment.

Merge the `g^2` old vertices in each parallel plane into one vertex.
Each cube factor is now local to one of two planes, joined by its
auxiliary channel. An original nongeneric edge either stays within one
plane or remains a wire between two planes. Thus the result is an
ordinary weighted clique on `g` vertices, not a hypergraph oracle.

### Exact width and volume bounds

Let `V=g^6 lambda` be the parent's information-volume convention.
For a fixed pair of planes, exactly `binom(g,2)^2` cubes meet both.
The number of original nongeneric edges between those two planes is

\[
g^4-[g(g-1)]^2=2g^3-g^2.
\tag{5}
\]

The separate width cap in (2) gives the particularly useful exact bound

\[
w_{ab}\le
4\lambda\binom g2^2+(2g^3-g^2)\lambda
=g^4\lambda=\frac V{g^2}.
\tag{6}
\]

Write `U=2 sum_(a<b) w_ab` for the resulting oriented volume. Taking
the smallest of the three total cube widths and using (2) gives

\[
\begin{aligned}
U
&\le8\rho\lambda B
 +\frac{2C}{3}B\log_2(\lambda+2)
 +g(g-1)(2g^3-g^2)\lambda,\\
\frac UV
&\le\rho(1-1/g)^3
 +\frac2g-\frac3{g^2}+\frac1{g^3}
 +\frac C{12}(1-1/g)^3\frac{\log_2(\lambda+2)}\lambda,\\
&\le\rho+\frac2g+
 \frac C{12}\frac{\log_2(\lambda+2)}\lambda.
\tag{7}
\end{aligned}
\]

Independently, (6) implies `U <= (1-1/g)V`. All these statements hold
for each piece assignment. No averaging of output values is used.

## 3. Equitable grouping with a deterministic search

The following elementary lemma suffices to regroup the `g` planes into
`K` vertices while keeping the recursion's growth base equal to three.

**Lemma.** Suppose `K` divides `g`, `K>=2`, and a nonnegative weighted
graph on `g` vertices has total oriented weight `U` and maximum edge
weight at most `V/g^2`. If

\[
\frac gK\ge16\zeta^{-2}(2\ln K+1),
\tag{8}
\]

then there is a partition into `K` groups of exactly `g/K` vertices
such that every weight between two distinct groups is at most

\[
\frac1{K^2}\left(\frac g{g-1}U+\zeta V\right).
\tag{9}
\]

An exhaustive scan of balanced maps finds such a partition in
`2^O(g log K)` bit operations, times a polynomial in the graph encoding.

**Proof.** Set `m=g/K` and choose a uniformly random equitable labeled
partition. For fixed distinct groups `A,B`, every original edge meets
these two groups with probability `2m^2/(g(g-1))`. Hence their cut
weight `W(A,B)` has mean `Ug/((g-1)K^2)`.

Expose the first `2m` positions of a random permutation, using the first
`m` for `A` and the next `m` for `B`. Changing one exposed choice can
be coupled by swapping its two possible values in the remaining
permutation. The two final cuts differ by at most one exchanged vertex
in each group. Since all edge weights are nonnegative, their weights
differ by at most `2m` times the maximum edge weight. The Doob martingale
therefore has at most `2m` increments of absolute value at most
`2mV/g^2`.

The elementary bounded-differences estimate gives

\[
\Pr\!\left(W(A,B)-\mathbb EW(A,B)>
       \frac{\zeta V}{K^2}\right)
\le\exp\!\left(-\frac{\zeta^2g}{16K}\right).
\tag{10}
\]

For completeness, this estimate follows by applying
`E exp(t Delta) <= exp(t^2 c^2/2)` to each centered increment with
`|Delta|<=c`, multiplying the conditional estimates, and minimizing
the resulting exponential Markov bound. No independence of increments
is required. Under (8), a union bound over fewer than `K^2/2` cuts has
probability below one. This proves existence. Enumerating at most
`K^g` maps and checking sizes and cut weights gives the deterministic
algorithm. \(\square\)

Pad all grouped edges to the integer ceiling `L` of (9). Equations
(7)--(9) imply

\[
K^2L\le
\left(\rho+\zeta+O(1/g)
 +O\!\left(\frac{\log(\lambda+2)}\lambda\right)\right)V
 +K^2.
\tag{11}
\]

All tests in the exhaustive grouping algorithm can use integer or
rational arithmetic. The logarithms in (8) are only used to prove a
sufficient fixed schedule; they need not be evaluated to choose a
partition.

## 4. Cost of a grouped table query

The grouping argument must not hide contraction of an arbitrary
`g^2`-vertex table. Each new plane table is evaluated by direct summation
over its old incident labels. It contains `g^2` old vertices, each with
at most `g^3-1` labels of `lambda` bits. Thus one plane-table query uses

\[
2^{O(V/g+C\lambda)}
\tag{12}
\]

operations and calls to the old entry oracle, with polynomial factors
in `g` and `lambda`. This deliberately includes both internal labels
and the old labels entering the cube factors.

For one grouped table, first enumerate the channels between its
`m=g/K` constituent planes. By (6), their total width is at most

\[
\binom m2\frac V{g^2}\le\frac V{2K^2}.
\tag{13}
\]

For each assignment, the plane queries are independent and their
values are multiplied. They are evaluated separately. It would be an
unnecessary larger cost to sum jointly over every old label of every
plane. The grouped entry oracle therefore has multiplicative cost

\[
2^{O(V/g+V/K^2+C\lambda+\log(g+\lambda))}
\tag{14}
\]

relative to one call of the old oracle. Zero padding is tested directly.
Thus a query implementation exists at every stage; it is not supplied
as advice or treated as free.

## 5. Uniform recursion and the factoring bound

Fix a sufficiently large integer `x0` and use

\[
K_0=2^{x_0},\qquad
g_j=2^{6j+x_0}K_j,\qquad
K_{j+1}=g_j^3,\qquad
\zeta_j=2^{-j-8}.
\tag{15}
\]

All these numbers are powers of two. In particular `K_j` divides
`g_j`. Their exact logarithms are

\[
\log_2K_j=
\frac{(5x_0+9)3^j-18j-3x_0-9}{2}.
\tag{16}
\]

For `x0>=2`, `K_(j+1)<=K_j^6`. Enlarging the fixed `x0` makes (8)
hold for every `j`: its left side grows as `64^j`, whereas the
right side, apart from a fixed factor, grows at most as `12^j`.

For an `n`-bit factoring input choose the largest level `ell` with
`P=K_ell<=n^(1/100)`. Finitely many smaller inputs use trial division.
The [explicit top constructor](half-top-constructor.md), or its
[polynomial-query version](polynomial-top-queries.md), supplies a
weighted `P`-vertex clique of volume `O(n)`. The same geometry is
available over the fixed characteristic-two coefficient algebras in
[the fixed-characteristic transfer](fixed-characteristic-transfer.md).
The original prefix or weighted witness recovery is unchanged.

Since the next level is at most the sixth power, the selected top
level satisfies

\[
n^{1/600}<P\le n^{1/100},\qquad
\ell=O(\log\log n),\qquad
\rho^\ell=\Theta((\log n)^{-\log_3(1/\rho)}).
\tag{17}
\]

At a call with volume at most `ceil(sqrt(n))`, evaluate by exhaustive
edge-label summation. Also use exhaustive summation at level zero.
Every other call has

\[
\lambda>\frac{\sqrt n}{P^2}\ge n^{12/25}.
\tag{18}
\]

This cutoff is tested by integer arithmetic against the original input
length. It controls every branch, including unexpectedly small child
volumes.

Equation (11), its rounding term, and (18) give, on every active path,

\[
V_j\le(\rho+O(\zeta_j)+\varepsilon_n)V_{j+1},
\quad
\varepsilon_n=O(n^{-12/25}\log n),
\tag{19}
\]

after increasing constants to cover the term `K_j^2/V_(j+1)`.
The sum of the `zeta_j` is finite, and
`epsilon_n log log n` tends to zero. Hence
`V_j=O(rho^(ell-j)n)` until the cutoff is reached.

The complete computation has the following costs.

| Operation | Charged exponent |
| --- | --- |
| Top construction or printed input | `O(n/P)`, or polynomial query setup |
| Piece assignments at one active step | `B log_2 |I_lambda| = O(V log(lambda+2)/lambda)` |
| Deterministic equitable grouping | `O(g log K)` |
| One level of entry-oracle composition | `O(V/g+V/K^2+C lambda+log(g+lambda))` |
| Exhaustive terminal contraction | `O(V_0)` or `O(sqrt(n))` |

Piece-assignment exponents, summed over a path, are
`O(n^(13/25) log(n) log log(n))`, by (18). Grouping-search exponents
are smaller polynomial powers of `n`. Both are negligible compared
with `n/(log n)^alpha` for every fixed `alpha` here. Computing all
piece metadata in `2^(C lambda)` time is also covered by (14).

For the entry-oracle costs, the relevant series is

\[
\sum_{j\ge0}\rho^{-j}
 \left(K_j^{-2}+g_j^{-1}+g_j^{-6}\right)<\infty.
\tag{20}
\]

The last term covers `lambda=V/g_j^6`. Equations (19)--(20) bound
the sum of these exponents along any active path by `O(n rho^ell)`.
Branch enumeration multiplies work, so its logarithm adds to these
path costs. The terminal volume has the same bound, unless the
`sqrt(n)` cutoff is smaller. This gives a uniform weighted-clique
evaluator in `exp(O(n rho^ell))` bit operations.

Finally apply the existing deterministic weighted recovery from
[bit-tagged witness recovery](bit-tagged-witness-recovery.md), or the
fixed-characteristic transfer above. Its additional
`2^O(sqrt(n)) poly(n)` coefficient and recovery cost is smaller than
the claimed exponent. The field operations and these coefficient
algebras have explicit bit implementations. No random bits used in
the existence proof of (9) enter the machine; it uses the exhaustive
grouping search. Together with (17), this proves the conditional
factoring bound stated at the beginning.

## 6. What remains to be constructed

The missing ingredient is precisely the local family (1)--(2), including
its effective factor queries and polynomial piece count. Three
orientations alone do not supply it. A hypothetical partition-style
profile giving `rho=1/sqrt(3)` still yields score exactly one half.

There is a useful negative control. Expanding (1) into its `2^(4 lambda)`
individual supported basis tensors makes every cut rank one. All three
widths are then zero. But enumerating those pieces over the cubes costs

\[
4\lambda B=\frac12(1-1/g)^3V
\tag{21}
\]

bits in the exponent at the first step. The polynomial piece hypothesis
fails, so this decomposition gives no subexponential factoring bound.

The exact geometric interface leaves room for genuinely eight-leg
families. A tensor product of known four-leg crossings is not asserted
to satisfy the required new total-width threshold.

## Reproduction

```sh
python3 tools/check_cube_crossing_transfer.py \
  --output evidence/cube-crossing-transfer.json
```

The finite [evidence](../evidence/cube-crossing-transfer.json) distinguishes
the geometry and cost checks from the unconstructed local hypothesis.
