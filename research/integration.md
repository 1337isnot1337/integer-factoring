# Integration of sharper grouping and crossing compression

## Status

This document specifies a uniform counting algorithm and gives its cost
accounting at the level of explicit table operations. It combines the new
grouping lemma with the rank-56 crossing factorization. It is a candidate
improvement to deterministic factoring, pending implementation of the
generalized evaluator and a complete machine-level proof.

The finite algebraic certificates, grouping theorem, and recurrence analysis
are separate results. Their validity does not mean that the final factoring
machine has already been constructed or that its Lean theorem compiles.

The source is pinned to
`IsaacRudich/SubexponentialDeterministicFactoring` commit
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`. File and theorem names below refer
to that revision.

## 1. Exact integer parameters

For levels numbered from the bottom, compute

\[
\begin{aligned}
K_0&=48,\\
D_d&=(d+48)^2,\\
G_d&=12D_d^2K_d^2(\lfloor\log_2K_d\rfloor+1),\\
K_{d+1}&=G_d^2.
\end{aligned}
\tag{1}
\]

For an `n`-bit input, choose the largest `ell` with `K_ell^40 <= n`, and set
`P=K_ell`. Inputs below the fixed threshold `48^40` can use trial division;
this finite initial range is absorbed into the asymptotic constant. No
practical speed claim is attached to this fallback or to the new algorithm.

The inequality

\[
7^{26}<2^{73}
\]

provides a rational version of the new contraction constant:

\[
\rho_0=\frac{281}{312},
\qquad
\alpha_0=\frac{\log(312/281)}{\log4}
\approx0.07548794932066111.
\tag{2}
\]

Every runtime parameter is computed with integer arithmetic. The logarithms
in (2) only describe the asymptotic bound.

## 2. New local factorization

Let `R,C` be the source's eight-square row and column tensors, and let `E`
be the crossing equality tensor. The source proves `R+C=2E`.

Put

\[
S=8E-R=3R+4C.
\]

In the row flattening, `R` has an explicit integer factorization with 8
components and `S` has one with 56 components. The latter decomposes into
eight XOR classes, with seven components per class. No division is needed.

For `t=ceil(lambda/3)`, zero-pad the strings to `t` three-bit blocks. If
`A_i` is the source's coefficient of `z^i` in the product of `(R+zC)`, use
the scaled piece

\[
B_i=4^t A_i
=[z^i]\prod_{s<t}((4-3z)R_s+zS_s).
\tag{3}
\]

Then

\[
\sum_{i=0}^t B_i=8^t E^{\otimes t}.
\tag{4}
\]

The row factorization uses at most `16^t 7^i` components. An explicit index
is a subset `Q` of the block positions, one of eight components at every
position, and `i` base-7 digits. For `j=|Q|<=i`, use `j` of those digits
to choose the rank-56 components, require the remaining digits to be zero,
and multiply one side by

\[
\binom{t-j}{i-j}4^{t-i}(-3)^{i-j}.
\]

Indices with `j>i` contribute zero. This is an effective factorization:
the index can be decoded and its value computed in polynomial time in
`t` and the arithmetic bit lengths. Binomial coefficients can be computed
by Pascal recursion modulo the working modulus, so this step works for
composite moduli too. Column symmetry gives `16^t 7^(t-i)` components.

Use the integer widths

\[
\begin{aligned}
c_r(\lambda,i)&=\min\left(2\lambda,
  4t+\left\lceil\frac{73i}{26}\right\rceil\right),\\
c_c(\lambda,i)&=\min\left(2\lambda,
  4t+\left\lceil\frac{73(t-i)}{26}\right\rceil\right).
\end{aligned}
\tag{5}
\]

Unused binary labels contribute zero. The two independent ceiling errors
and the final padded block give

\[
c_r+c_c\le 4\rho_0\lambda+10.
\tag{6}
\]

The source's `Level.lean` counting of straight edges and rectangle pairs
therefore gives, for the smaller of the row and column networks,

\[
U\le (\rho_0+1/G)G^4\lambda+3G^4,
\qquad P_{\rm node}\le G^3\lambda.
\tag{7}
\]

Here `U` is the ordered edge total. The `3G^4` term safely covers the
rounding error from (6); retaining the source's `2G^4` unchanged would be
incorrect.

## 3. One recursive step

A level takes a written-out clique on `K_(d+1)=G_d^2` vertices, each pair
having a `lambda`-bit label.

1. Arrange its vertices on the `G_d` by `G_d` grid.
2. Enumerate all assignments of piece indices `i in {0,...,t}` to grid
   rectangles. The source's redundant encoding by all `G_d^4` ordered
   quadruples may be retained.
3. For one assignment, compute the row and column network edge widths
   from (5), and choose the network with smaller ordered total.
4. Enumerate the `K_d^G_d` groupings until finding one with at most
   `2G_d/K_d` vertices per group and
   `D_d K_d^2 Y_ab <= D_d U + G_d P_node` for every group pair.
   The new grouping theorem guarantees success.
5. Write out the coarsened clique on `K_d` vertices, padded to the common
   width

\[
\lambda_d=
\left\lfloor
\frac{D_d U_d^{\max}+G_d^4\lambda}{D_dK_d^2}
\right\rfloor,
\quad
U_d^{\max}=
\left\lfloor
  \frac{281G_d^4\lambda}{312}+G_d^3\lambda+3G_d^4
\right\rfloor.
\tag{8}
\]

The width in (8) depends on `d` and `lambda`, not on the piece assignment.
That is essential: all branches have the same later widths and the same
positive multiplier.

With `N_(d+1)=G_d^4 lambda` and `N_d=K_d^2 lambda_d`, (8) gives

\[
N_d\le(\rho_0+1/G_d+1/D_d)N_{d+1}+3G_d^4.
\tag{9}
\]

Every arithmetic step in (8) is an integer addition, multiplication, or
division with remainder. The rounding directions in (8) and (9) agree.

## 4. Complete algorithm skeleton

The following pseudocode describes the missing generalized evaluator.
`top_clique` is the source's multiplication-network construction, retaining
its original top grouping `g_top=3P^8`. `coarsen_and_tabulate` is the direct
sum over local state assignments from the source's network semantics.

```text
factor(N):
    if N == 0: return 2
    if N == 1: return 1
    n = bit_length(N)
    if n < 48^40: return trial_division_min_factor(N)
    compute K_d, D_d, G_d until K_(ell+1)^40 > n
    prefix = 0
    for k from n-1 down to 0:
        if positive_prefix_count(N, k, 2*prefix, schedule):
            prefix = 2*prefix
        else:
            prefix = 2*prefix + 1
    return prefix

positive_prefix_count(N, k, x, schedule):
    derive all padded widths lambda_d from the initial top width
    M = n + 4*sum(G_d^4 * ceil(lambda_(d+1)/3), d=0,...,ell-1) + 4
    for m from 2 through 2*M:
        T = top_clique(N, k, x, P, g_top=3*P^8, modulus=m)
        value = evaluate_chain(T, ell, modulus=m)
        Hm = product_of_level_multipliers(modulus=m)
        if value != Hm * indicator(1 // 2^k == x) mod m:
            return true
    return false

evaluate_chain(T, d, modulus):
    if d == 0:
        return brute_force_clique_value(T, modulus)
    G = G_(d-1)
    t = ceil(T.width / 3)
    total = 0
    for choice in {0,...,t}^(G^4):
        create row/column networks using the scaled pieces B_i
        W = smaller_total_network(row, column)
        chi = lexicographically_first_valid_grouping(W, K_(d-1), D_(d-1))
        Tnext = coarsen_and_tabulate(W, chi, fixed_width=lambda_(d-1))
        total += evaluate_chain(Tnext, d-1, modulus)
    return total mod modulus
```

This recursive presentation is convenient for the proof. It can be made a
single fixed program by enumerating the mixed-radix choice paths and
iterating down one path at a time, as `UniChain.lean` does. The number of
levels is runtime data, not the number of separately compiled subprograms.

For the redundant `G^4` choice encoding, one step multiplies the clique
value by

\[
(t+1)^{G^4-r_G}8^{t r_G},
\qquad r_G=\binom G2^2.
\tag{10}
\]

The product `H` of these positive integer factors is independent of the
choice path. It must be computed modulo each `m`. One must not divide by
`H` modulo `m`, because `H` need not be invertible. The contribution of
the divisor `1` is subtracted after multiplying by the same `H`.

The modular zero test reuses the logic of `Moduli.lean`: the least common
multiple of `1,...,2M` is large enough to exclude a nonzero integer in the
proved range. For the explicit `M` in the pseudocode,
`log_2 H <= 4 sum_d G_d^4 ceil(lambda_(d+1)/3)` follows from
`log_2(t+1)<=t` and (10). Hence the nonnegative integer
`V=H*cntPre(N,k,x)` is less than `2^(M-4)`, and
`(2M+1)V<4^M` holds for `M>=4`, as required by the source theorem.
This computes `M` from the known width schedule; it requires neither the
unknown divisor count nor construction of the enormous integer `H`.
Since the moduli
extend only to `2M`, this adds a polynomial factor, not `2^M` iterations.

## 5. Cost accounting for the abstract evaluator

The following estimates include the grouping search and choice enumeration.
Their omission would invalidate the proposed bound.

### Parameter sizes and the top construction

For (1), `log K_d=Theta(4^d)`, so

\[
\ell=\frac{\log\log n}{\log4}+O(1).
\]

Also, `K_(d+1)>=K_d^4` and `K_(d+1)<=K_d^16` give

\[
n^{1/640}<P\le n^{1/40}
\]

on the nontrivial input range. The original top construction has volume
`N_ell=O(n+P^8 n^(3/4))=O(n)`. Its leading exponent is bounded by

\[
O(n/P+P^8n^{3/4}+\operatorname{polylog}(n)P^c),
\]

with the fixed polynomial terms accounted for by the source's explicit
top-evaluator bounds. The already proved conservative estimates bound
these extra terms by `O(P^38)`. Thus all top terms are at most
`O(n^(639/640)+n^(19/20))`, apart from fixed constants.

### Chain volumes and table creation

The sums of `1/G_d` and `1/D_d` converge. Moreover, for recursive levels,

\[
G_d^4=K_{d+1}^2\le P^2.
\]

Unrolling (9) therefore gives, uniformly over `0<=d<=ell`,

\[
N_d\le C\rho_0^{\ell-d}N_\ell+C\ell P^2
\tag{11}
\]

for a fixed constant `C`. Since `K_d>=48*2^d` and `2rho_0>1`,

\[
\max_{d<\ell}\frac{N_{d+1}+N_d}{K_d}
=O(\rho_0^\ell n+\ell P^2).
\tag{12}
\]

A group contains at most `2G/K` row or column nodes. Each node has a state
of `G^3 lambda` bits. Enumerating the local states for one group costs at
most `2^(2N/K)`. Writing its entire new table adds at most `2^(N'/K)`
entries.

There is another exponential loop in the source that must be counted:
`CoarsenEx.W` sums over the factor channels of every edge internal to a
group. A row/column pair has at most `G^2 lambda=N/G^2` channel bits, so
an internal pair can require `2^(N/G^2)` operations. There are polynomially
many such pairs, and their sums are computed separately and multiplied;
one does not need to enumerate their joint Cartesian product. Since
`G>=K`, this extra exponent is at most `N/K`.

For one fixed channel index, the new factor decoder requires polynomial
work in the string lengths: scan the subset and its mixed-radix digits,
compute the binomial coefficient by a polynomial dynamic program, and
multiply the block factors modulo `m`. In the trivial-width fallback,
evaluate `4^t A_i` by polynomial coefficient dynamic programming. Calling
the finite checker's `expanded_coefficient` here would instead enumerate
`2^t` subsets; that routine is a tiny-instance verifier, not the proposed
polynomial decoder. Even such a loop has exponent `O(lambda)<=O(N/K)`,
but a formal implementation must charge it if retained.

Thus a conservative tabulation bound includes the exponent
`(2N+N')/K + N/G^2 + O(lambda)` and has cost

\[
\operatorname{poly}(n)2^{O((N+N')/K)}.
\tag{13}
\]

This bound also allows a polynomial simulation of table lookups by
sequential memory: changing a fixed polynomial overhead multiplies the
exponent by a fixed constant and does not change `alpha_0`.

### Grouping search

At each recursive level, the deterministic grouping search costs

\[
K_d^{G_d}\operatorname{poly}(n).
\]

Its binary logarithm is at most

\[
G_d\log_2K_d\le P^{1/2}\log_2P.
\tag{14}
\]

The graph weights are explicit integers with polynomial bit length. The
search is finite and uniform, despite being impractical.

### Number of choice paths

At most `(t_d+1)^(G_d^4)` choices occur at a level. Since every width is
`O(n)`,

\[
\log_2(\#\text{paths})
\le\sum_{d<\ell}G_d^4\log_2(t_d+1)
=O(\ell P^2\log n).
\tag{15}
\]

The recursion can therefore be bounded by the number of paths times
the maximum cost of one path. In particular, it does not require a
factor of `rho_0^(-ell)` beyond the volume bound in (12).

### Moduli, values, and output search

Equation (10) gives `log_2 H=O(sum_d G_d^4(t_d+1))=poly(n)`.
The number and bit lengths of the moduli are polynomial. Modular arithmetic
prevents the enormous exact count from appearing as a unit-cost machine
integer. Prefix recovery needs at most `n` such positivity tests. The
cases `N=0,1` are treated separately, matching `Nat.minFac`.

Combining (11)--(15), the abstract evaluator's cost is bounded by

\[
\operatorname{poly}(n)\,
2^{O(\rho_0^\ell n+\ell P^2\log n+n^{639/640}+n^{19/20})}
=2^{O(n/(\log n)^{\alpha_0})}.
\tag{16}
\]

The sharper constant `rho=(8+log_2 7)/12` and exponent
`0.07551046810477115` also admit an integer-only width schedule. Define
`b(j) = bit_length(7^j - 1)`, so `b(j)=ceil(j log_2 7)` for every
nonnegative integer `j`, including `b(0)=0`. Replace (5) by

\[
c_r=\min(2\lambda,4t+b(i)),\qquad
c_c=\min(2\lambda,4t+b(t-i)).
\]

The common, piece-independent integer upper bound in (8) can then be
replaced by

\[
 U^{\max}_{\rm exact}=G^3\lambda+
 \left\lfloor\frac{G^4(8t+b(t)+1)}4\right\rfloor.
\]

Indeed, `b(i)+b(t-i)<=b(t)+1`, including the endpoint cases. The grid
edge count bounds the smaller ordered total by
`G^3 lambda + G^4 (c_r+c_c)/4`; taking its floor is valid because that
total is an integer. Further, `b(t)<=t log_2 7+1` and
`t<=lambda/3+2/3`, so

\[
 U^{\max}_{\rm exact}
 \le (\rho+1/G)G^4\lambda+3G^4.
\]

The same width formula as (8), now using `U_exact^max`, and the same
cost argument yield the sharper *candidate* exponent under the stated
integration obligations. Computing `7^t-1` and its bit length uses
integers of `O(t)` bits and polynomial time in `t`; the real logarithm
appears only in the analysis. The reference implementation and finite
boundary checks are `tools/exact_widths.py` and `tests/test_exact_widths.py`.
They do not implement the counting engine or prove its machine-level cost.

## 6. What is still missing

The mathematics above fixes the algorithms behind the asymptotic
interfaces. A final certified factoring claim still requires:

- an implementation of the rank-56 decoder in the source expression
  language, with its polynomial evaluation cost proved;
- a generalized `nextT` using (8), the new grouping theorem, and the new
  multiplier (10), together with its value identity;
- a uniform expression program for the schedule (1) and common widths;
- revised cost and value bounds for that actual program, followed by the
  compiler and finite-alphabet Turing-machine simulation;
- a final theorem with the intended input/output encoding and a checked
  axiom list.

The upstream compiler theorems are generic infrastructure, but their
successful application to the old program does not certify a changed
program. Its numerical cost exponents, such as the old fixed polynomial
degree, cannot simply be copied into a new formal claim.

This candidate is weaker than the benchmark's unproved square-root-log
conjecture: both candidate exponents are below `1/2`. It is a quantitative
strengthening of the published `2^{o(n)}` type of bound if the stated
integration obligations are completed.

## Sources and local evidence

- `research/recurrence.md`: parameter growth, contraction, and all-level costs.
- `tools/grouping.py` and `tests/test_grouping.py`: exact finite grouping checks.
- `tools/check_crossing.py`: finite crossing and complement certificates.
- [Upstream grid/value identities](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Level.lean).
- [Upstream piece construction](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Pieces.lean).
- [Upstream cost requirements](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostReq.lean).
- [Upstream modular zero test](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Moduli.lean).
- [Upstream recovery of the smallest prime factor](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Search.lean).
