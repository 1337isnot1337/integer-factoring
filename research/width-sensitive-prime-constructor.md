# A width-sensitive evaluator and a faster prime-constructor schedule

## Result and scope

For any fixed rational `0<epsilon<=1`, the prime/CVP construction can be
scheduled, for all sufficiently large `n`, with

\[
 P=\Theta_\epsilon\!\left(\sqrt{n/\log n}\right),\qquad
 \lambda=\Theta_\epsilon(\log n),\qquad
 P^2\lambda<(12+\epsilon)n,\qquad B_{\rm tab}=1.           \tag{1}
\]

All explicitly written tables, including exact entry computation, cost

\[
 \exp\!\left(O_\epsilon(\sqrt{n\log n})\right).           \tag{2}
\]

The prime search is deterministic polynomial time with an effective
exponent depending on the fixed `epsilon`. It is adaptive and does not need the numerical value
of the analytic prime-distribution constant.

To compose this larger initial clique with the recursive evaluator, we
prove the strengthened circuit bound

\[
 2^{Ck^2\lambda/(\log(k+2))^\alpha}
                  (\lambda+2)^{D_0k^2}.                  \tag{3}
\]

It replaces the factor `(k(lambda+1)+2)^(D_0 k^2)` in the
[asymmetric transfer](asymmetric-modular-transfer.md). The hypotheses on
the crossing pieces and the admissible exponent `alpha` are unchanged.
This is a quantitative construction and overhead improvement. It does
not supply a new crossing exponent or a complete score above `1/2`.

## 1. Remove the unnecessary vertex-count factor from the induction

Use the same pieces, orientation choice, grouping, and geometric margin
as Sections 1--4 of the asymmetric transfer. Write

\[
 F(k,\lambda)=\frac{k^2\lambda}{(\log(k+2))^\alpha},\qquad
 L(\lambda)=\log_2(\lambda+2).
\]

For the selected child `(K,mu)`, that proof supplies fixed `eta<1` and
effective constants such that

\[
 F(K,\mu)\le\eta F(k,\lambda)+C_1k^2 L(\lambda),\qquad
 \mu\le C_2k^2\lambda.                                    \tag{4}
\]

Both possible child counts grow as fixed positive powers of `k` divided
by a logarithm. Thus there is a fixed `delta>0` such that, after a finite
threshold,

\[
 k^\delta/\log k\ll K\ll k^{1-\delta}.                   \tag{5}
\]

All constants are independent of `lambda>=1`. The piece-assignment count
has logarithm at most `C_3 k^2 L(lambda)`. The grouping search, construction
of indices, and replacement of one child query together have logarithm

\[
 o(F(k,\lambda))+O(\log(k(\lambda+1)+2)),                  \tag{6}
\]

uniformly for `lambda>=1`; these are exactly the bounds established in
the existing proof. In particular, no term of order `k^2 log k` occurs
outside the child's induction factor.

The logarithm of the proposed child polynomial factor is bounded by

\[
 D_0K^2L(\mu)
 \le D_0K^2L(\lambda)+C_4D_0K^2\log(k+2).                 \tag{7}
\]

The first term is at most `D_0 k^2 L(lambda)/4` for large `k`. By (5),

\[
 \frac{K^2\log(k+2)}{F(k,\lambda)}
 \ll \frac{K^2(\log(k+2))^{\alpha+1}}{k^2}
 \longrightarrow0                                        \tag{8}
\]

uniformly for `lambda>=1`. Hence the second term belongs in the geometric
exponential margin, rather than in the parent polynomial factor.

Here is a noncircular choice of constants. First choose a fixed `d_0`
large enough to absorb the width error `C C_1 k^2L(lambda)` from (4)
and the piece count, and set `D_0=d_0 C`. Take `C>=1`. Now choose the
finite induction threshold so that the first term of (7) uses at most
one quarter of the parent polynomial exponent, and so that (6) and the
second term of (7) use at most one half of the margin
`C(1-eta)F(k,lambda)`. This threshold depends on `d_0=D_0/C`, not on
the eventual value of `C`. The width and piece terms use another fixed
fraction of the parent polynomial exponent by the choice of `d_0`.
Finally choose `C` large enough for direct enumeration at every smaller
vertex count. There are finitely many such counts, and their costs are
bounded uniformly in `lambda` by `2^O(lambda)` times a polynomial in
`lambda`. Increasing `C` while retaining `D_0=d_0C` preserves every
previous inequality. This closes the induction and proves (3).

The argument changes no query semantics or algebraic operation. The
uniform circuit remains valid over every commutative algebra admitted
by the original transfer, including algebras with zero divisors.

## 2. An effective dyadic prime for a fixed large power of the modulus

There are effective absolute constants `C_*` and `q_*` such that, for
every `q>=q_*`, every reduced residue class `a mod q`, and every
`x>=q^(C_*)`, the interval `(x/2,x]` contains a prime congruent to `a`.

This is a consequence of Thorner--Zaman, *Refinements to the prime number
theorem for arithmetic progressions*, Theorem 1.1 and equations
(1.2), (1.7), with the effective constants explicitly stated there:
[primary paper, arXiv:2108.10878v2](https://arxiv.org/pdf/2108.10878).
We record the substitution because the fixed-power range here is
stronger than the squared-logarithm sufficient condition used in the
[fully explicit BMOR companion](power-two-progression-dyadic-prime.md).

Set their interval length to `h=x/2` and fix their error parameter to
`1/74`. Even in the exceptional-zero case, their effective bounds give
`Lambda>=c q^(-1/2)` for their density parameter, an effective positive constant, and sufficiently
large `x`. Therefore

\[
 \frac{\Lambda h}{\varphi(q)}\gg xq^{-3/2}.               \tag{9}
\]

For a sufficiently large fixed `C_*`, this exceeds
`x^(32/37+1/74)` whenever `x>=q^(C_*)` and `q` is sufficiently large.
The relative error in their theorem is at most an effective constant
times

\[
 \exp\!\left(-c\frac{\log x}
 {\log q+O(1)+(\log x)^{2/5}(\log\log x)^{1/5}}\right).
                                                               \tag{10}
\]

Uniformly for `x>=q^(C_*)`, the denominator is at most
`2 log x/C_*` once `q` exceeds another effective threshold. Choose
`C_*` large enough that this relative error is below `1/2`. The main
term is positive by (9), proving the dyadic-prime assertion. No
ineffective Siegel bound or unproved zero-free hypothesis is used.

## 3. Adaptive exact parameters

Fix `0<epsilon<=1`, let `K` be the least power of two at least
`128/epsilon`, and set `delta=epsilon/48`. For
`u=ceil(log_2(n+2))`, let `A_0` be the largest power of two satisfying

\[
 A_0^2\le\frac{\delta n}{u}.                              \tag{11}
\]

At stage `j`, use `A=A_0/2^j`, as long as `A>=2K`, and put

\[
 \begin{split}
 P&=A/K,&T&=A^2,&D&=A^2/2,&t&=\log_2D,\\
 w&=2\left\lceil\frac{n+1}{2(T-A-1)}\right\rceil,
 &b&=2^w,&s&=2w+t,&X&=2^s.\end{split}                     \tag{12}
\]

Search every integer `Q` in `(X,2X)` congruent to `1 mod T`, certifying
primality by trial division. If none is prime, halve `A` and repeat.
If a prime is found, use it in the
[prime/CVP constructor](prime-transform-boundary-constructor.md), with
residue width `f=s+1`, carry width `c=w+t`, and exactly `K` Fourier
indices in each contiguous group.

If the loop reaches `A<2K` without finding a prime, execute the previous
constructor or direct divisor enumeration. This is a computable branch
of the algorithm and requires no knowledge of an analytic cutoff. The
argument below proves that it occurs on only finitely many input lengths
for each fixed `epsilon`.

For all sufficiently large `n`, this loop succeeds after a bounded
number of stages. Indeed, `A_0` is within a factor two of the square
root in (11). At stage `j`,

\[
 \frac{\log_2(2X)}{\log_2T}
 \ge\frac{2w}{\log_2T}
 \ge\frac{2\cdot4^j}{\delta}.                            \tag{13}
\]

Choose a fixed `j_*` so that the final quantity is at least `C_*`.
As `n` grows, all stages through `j_*` have `T>=q_*` and `A>=2K`.
Section 2 then supplies a prime in the required interval by that stage.
The endpoint `2X` is a power of two and is not prime. The loop does not
need `C_*`, `q_*`, or `j_*` as input: it only performs finite searches
and checks whether one succeeded.

It follows that

\[
 A=\Theta_\epsilon(\sqrt{n/\log n}),\qquad
 w=\Theta_\epsilon(\log n),\qquad
 s=\Theta_\epsilon(\log n).                              \tag{14}
\]

Each complete prime search costs `2^O(s)`, hence polynomial time in `n`.
Finding a quadratic nonresidue and a primitive `T`th root uses the exact
procedure in the prime/CVP note and also costs `2^O(s)`. There is no
factorization of `Q-1`. Finitely many small inputs may use the existing
constructor or direct divisor enumeration; this does not affect (1)--(3).

## 4. Volume and all written tables

All combinatorial and arithmetic hypotheses of the prime/CVP constructor
hold: `A>=2P`, the common blocks have length `K<=A/2`, there is no cyclic
convolution aliasing, `Q>Db^2`, the exterior observation lattice has
minimum at least `sqrt(Q/2)>bA/2`, and every compatible bounded input is
the unique closest lattice point to its box center. Each table entry is
Boolean before weights; exact CVP computes it in `2^O(A) poly(n)` time.
Invalid residue codes, padding, digit restrictions, and prefixes are
postchecked exactly as in that construction.

The exact channel width and volume are

\[
 \begin{split}
 \lambda&=6K^2f+2Kf+4Kc,\\
 \frac{P^2\lambda}{n}
 &=\frac{A^2}{n}\left[
     (12+8/K)w+(6+6/K)t+6+2/K\right].                     \tag{15}
 \end{split}
\]

For every attempted stage with `A>=2K` and `n>=2^40/epsilon`, the
following rational estimates hold:

\[
 \frac{A^2w}{n}<1+\frac{\epsilon}{100},\qquad
 \frac{A^2t}{n}\le\frac{\epsilon}{48},\qquad
 \frac{A^2}{n}\le\frac{\epsilon}{48u},\quad u\ge41.       \tag{16}
\]

For the first, use
`A^2/(A^2-A-1)<=1+2/A`, `A>=256/epsilon`, and the rounding error below
two in `w`; the excess is at most
`epsilon/128 + epsilon/984 + 2^(-40)epsilon(1+epsilon/128)`,
which is less than `epsilon/100`. The other two follow from (11) and
`t<log_2 n<=u`.

Since `K>=128/epsilon`, substitution in (15) gives

\[
 \frac{P^2\lambda}{n}
 <12+\epsilon\left(
    \frac{12}{100}+\frac1{16}+\frac1{1600}
    +\frac{6+3/64}{48}+\frac{6+1/64}{48\cdot41}\right)
 <12+\epsilon.                                           \tag{17}
\]

There are `P 2^((P-1)lambda)` explicitly written table entries. Charging
every CVP call, all arithmetic, the polynomial number of operand-length
and prefix instances, and the prime search gives

\[
 \operatorname{poly}(n)\,2^{(P-1)\lambda+O(A)}
      +\operatorname{poly}(n)
   =\exp(O_\epsilon(\sqrt{n\log n})).                     \tag{18}
\]

This is a full written-table bound, not merely a query bound. If the
tables are queried through an oracle instead, each original query costs
`2^O(A) poly(n)` and that cost must still be included.

## 5. Composition and the scale attained

In (3), the overhead at the parameters (1) is

\[
 D_0P^2\log_2(\lambda+2)
   =O_\epsilon\!\left(\frac{n\log\log n}{\log n}\right)
   =o\!\left(\frac{n}{(\log n)^\alpha}\right)
 \quad(\alpha<1).                                        \tag{19}
\]

The construction and query-replacement costs in (18) are also lower
order. Since

\[
 \log P=\tfrac12\log n-\tfrac12\log\log n+O_\epsilon(1), \tag{20}
\]

the evaluator's leading exponent is at most
`((12+epsilon) C 2^alpha+o(1)) n/(log n)^alpha`, before the separate
coefficient-field and sparse-recovery costs already charged by the
transfer. This remains conditional on the same crossing certificate;
it does not alter `alpha`.

The scale in (1) is maximal up to a constant for this particular equal
block, explicit-residue constructor with `P^2lambda=O(n)`. Indeed,
`Q=1 mod A^2` gives `f>=log_2(A^2)`, so

\[
 P^2\lambda\ge6A^2\log_2(A^2).                           \tag{21}
\]

Thus `A=O(sqrt(n/log n))` and `P<=A` has the same upper bound. This is
a constraint of the stated message representation; it is not a lower
bound on arbitrary factoring networks or algorithms.

## 6. Exact parameter replay

The companion checker
[`check_width_sensitive_prime_constructor.py`](../tools/check_width_sensitive_prime_constructor.py)
tests the exact ceilings, volume identity, adaptive halving ratio, and
rational inequalities on finite boundary grids. It does not numerically
certify the unspecified analytic constants in Section 2, whose proof is
the effective theorem substitution above. The exact command and evidence
counts are recorded in the [evidence](../evidence/width-sensitive-prime-constructor.json).
Run:

```sh
python tools/check_width_sensitive_prime_constructor.py --output evidence/width-sensitive-prime-constructor.json
```

The replay passed 10,752 exact stage checks on 1,698 input boundaries,
six complete small prime/root constructions, and one negative control
where a dyadic interval below the hypotheses contains no required prime.
The uniform margin in (17) is exactly `245731/787200<1`. The largest
observed `(P^2lambda/n-12)/epsilon` was `20005/109568`; this observation
is separate from the proved uniform bound.
