# Exact divisor-count preprocessing from next-prime exchanges

## Result and evidence

The later [polynomial preprocessing theorem](polynomial-divisor-preprocessing.md)
strengthens the uniform existence bound to polynomial time and logarithmic
workspace using Nicolas's classical exponent estimates. Its fixed global
constants are not numerically instantiated. The concrete exchange oracle and
its separate resource bounds remain the subject of this note.

Let

\[
D(n)=\max_{1\le N<2^n}\tau(N),\qquad
q(n)=\lceil\log_2D(n)\rceil,
\]

and let `M(n)` be the least integer attaining `D(n)`. One uniform
deterministic algorithm computes these three quantities in

\[
\boxed{n^{O(\log\log(n+2))}
       =2^{O(\log(n+2)\log\log(n+3))}}
\tag{1}
\]

bit time, using

\[
\boxed{O(\log(n+2)\log\log(n+3))}
\tag{2}
\]

working bits with read-only input and streamed output. The small-space
version uses the established uniform iterated-integer-multiplication
theorem of Hesse, Allender, and Barrington. An elementary implementation
using only the previously checked carry recurrence instead has time
`2^{O(log^2(n+2))}` and workspace `O(log^2(n+2))`.

The practical implementation below uses ordinary integers. It meets the
time estimate (1) with polynomial workspace. Its separate `--streamed`
mode implements the elementary, weaker resource bound. The classical
uniform arithmetic circuits are a proved dependency of (1)--(2), not
an unimplemented primitive silently attributed to either Python mode.

This improves the auxiliary
[exact divisor-count preprocessing](https://github.com/1337isnot1337/integer-factoring/blob/dd28d375ea65863356c67b051f5d4f217562475d/research/exact-divisor-count-parameter.md).
The [primorial-weighted count](primorial-weighted-divisor-preprocessing.md)
shows that there are `2^{Theta(sqrt(n/log n))}` canonical exponent
patterns. The algorithm here avoids most of them by a necessary
condition on a *least maximizer*.

The factoring score remains `1/2`. Computing a maximum divisor count
as a function of the bit bound does not factor a supplied integer.
This note makes no claim that its elementary exchange inequalities
are new in the literature on highly composite numbers.

## 1. A next-prime exchange bounds every exponent

For `n>=2`, canonicalize the least maximizer as

\[
M=\prod_{i=1}^k p_i^{e_i}<2^n,
\qquad e_1\ge\cdots\ge e_k\ge1,
\tag{3}
\]

where `p_i` are the consecutive primes. Replacing prime bases by
smaller unused ones or moving a larger exponent to a smaller prime
can only decrease the integer while preserving its divisor count.
The least maximizer therefore has exactly this canonical form.
The case `n=1` has `D=M=1` and is handled separately.

Let `q=p_(k+1)` be the next unused prime. For every `i`, put

\[
r_i=\left\lceil\frac{e_i}{2}\right\rceil.
\]

**Exchange lemma.**

\[
\boxed{p_i^{r_i}<q.}
\tag{4}
\]

If this failed, replace `r_i` copies of `p_i` by one copy of `q`:

\[
M'=M\frac{q}{p_i^{r_i}}.
\]

The two numbers `q` and `p_i^(r_i)` cannot be equal because `q`
is a distinct prime. Hence the alleged failure gives `M'<M`.
The exact divisor-count ratio is

\[
\frac{\tau(M')}{\tau(M)}
 =\frac{2(e_i-r_i+1)}{e_i+1}\ge1,
\tag{5}
\]

since `2r_i<=e_i+1`. If the ratio is greater than one, maximality
fails. If it is one, the least-witness tie rule fails. Both are
contradictions. This proves (4). Zeroing the exponent when `e_i=1`
causes no problem in the divisor formula, although (4) is already
automatic in that case.

One consequence is the simpler largest-exponent bound

\[
e_1\le2\lceil\log_2q\rceil-2.
\tag{6}
\]

For an odd prime `q`, (4) at `p_1=2` implies
`ceil(e_1/2)<log_2 q`, giving (6). Equivalently, replacing
`2^h` by `q`, where `h=ceil(log_2q)`, contradicts least
maximality whenever `e_1>=2h-1`.

The predecessor's elementary prime estimate is

\[
p_j\le4j^2.
\tag{7}
\]

For reference, the central binomial coefficient gives
`2^m <= binom(2m,m) <= (2m)^(pi(2m))`. With `m=2j^2` and
`log_2 j<=j-1`, this gives `pi(4j^2)>=j`, proving (7).
Also `k<=n-1`, since `2^k<=M<2^n`. Thus

\[
q\le4n^2=:B.
\tag{8}
\]

In particular, only `O(log n)` different exponent levels occur.
Using that observation alone already gives a
`2^{O(log^2 n)}` candidate search. The full inequalities (4)
give the sharper bound below.

## 2. Column heights have a harmonic bit budget

Transpose the exponent diagram by defining

\[
h_a=\#\{i:e_i\ge a\},\qquad a\ge1.
\]

Then

\[
h_1\ge h_2\ge\cdots\ge0,
\quad e_i=\#\{a:h_a\ge i\},
\quad \sum_a h_a=\sum_i e_i\le n-1.
\tag{9}
\]

For any nonzero `h_a`, apply (4) to `i=h_a`. As `e_i>=a`,

\[
(h_a+1)^{\lceil a/2\rceil}
\le p_{h_a}^{\lceil a/2\rceil}
\le p_{h_a}^{\lceil e_i/2\rceil}
<q\le B.
\tag{10}
\]

Put `L=ceil(log_2 B)` and `E=min(n-1,2L-2)`. Equation (6)
shows that `h_a=0` for `a>E`. For `1<=a<=E`, set

\[
J_a=\min\left\{n-1,
       \left\lfloor(B-1)^{1/\lceil a/2\rceil}\right\rfloor-1
       \right\}.
\tag{11}
\]

The strict inequality in (10) is why the integer root uses `B-1`.
Every least maximizer lies in the finite box

\[
0\le h_a\le J_a\qquad(1\le a\le E).
\tag{12}
\]

Its size is at most

\[
\begin{split}
\prod_{a=1}^E(J_a+1)
&\le\prod_{a=1}^E B^{1/\lceil a/2\rceil},\\
\log_2\prod_{a=1}^E(J_a+1)
&\le2L\sum_{r=1}^{\lceil E/2\rceil}\frac1r
=O(\log n\log\log n).
\end{split}
\tag{13}
\]

The harmonic estimate follows from
`sum_(r=1)^R 1/r <= 1+log R`. Rounding each coordinate width
up adds at most `E=O(log n)` bits. Thus a complete height vector
can be stored in `O(log n loglog n)` bits, using
`ceil(log_2(J_a+1))` bits at coordinate `a`.

The box includes nonmonotone profiles, excessive total exponent
budgets, and products at least `2^n`. All are explicitly rejected.
The count (13) already pays for these rejected cases.

## 3. The uniform search

Compute `L,E` and each required `J_a` using integer bit lengths
and integer roots. A root is found by binary search, with repeated
multiplication stopped once it reaches `B`. All intermediate
integers are `O(log(n+2))` bits. No real logarithm or real root
is an algorithmic operation.

Enumerate the box (12) using a mixed-radix counter in the packed
coordinate widths above. A coordinate bound or its offset can be
recomputed when needed. There is no need to retain a recursion
stack with an `O(log n)`-bit remaining budget at every level.
For each vector:

1. Check the monotonicity and the sum bound in (9).
2. Regenerate the canonical integer `M(h)` and reject it if
   `M(h)>=2^n`.
3. Compare its exact divisor count with the retained best profile.
   If larger, retain it. If equal, retain the smaller integer.

The empty profile gives the permanent fallback `M=D=1`.
The exchange lemma ensures that at least one least maximizer is
among these candidates; canonicalization ensures it is the
global least witness. Hence the search computes the claimed
quantities exactly, with no heuristic pruning.

The practical program uses the stronger, input-dependent caps

\[
h_a\le\#\{i\le k:p_i^{\lceil a/2\rceil}<p_{k+1}\},
\qquad k=h_1,
\tag{14}
\]

and visits only nonincreasing profiles with affordable exact
prefix products. Every retained candidate still lies in (12),
and every least maximizer satisfies (14). These improvements
reduce finite work but are not needed for (13).

## 4. Exact products and output within the stated workspace

From a stored height vector, regenerate any exponent by
`e_i = #{a:h_a>=i}`. The candidate integer is a product of at
most `n-1` prime factors with repetition, each at most `4n^2`.
Its total bit length is `O(n log(n+2))`, even before feasibility
is known. The divisor count is the product of the values `e_i+1`
and has at most `n` bits, since `e_i+1<=2^(e_i)` and (9).

Any requested prime `p_i` can be regenerated by trial scanning
integers through `4i^2` and trial divisors through `2i`. This
takes polynomial time in `n` and `O(log(n+2))` working bits.
Therefore any requested bit of the list of small factors is
computable from the packed descriptor in the same resources.

The only additional classical result used for the strongest
simultaneous bounds is uniform iterated integer multiplication:
Hesse, Allender, and Barrington, *Uniform Constant-Depth Threshold
Circuits for Division and Iterated Multiplication* (2002),
[Corollary 6.3](https://people.cs.rutgers.edu/~allender/papers/division.pdf),
[DOI 10.1016/S0022-0000(02)00025-9](https://doi.org/10.1016/S0022-0000(02)00025-9).
It gives polynomial-size, constant-depth majority circuits with
DLOGTIME-uniform wiring for each output bit of an integer product.

Evaluate such a circuit recursively, counting the true inputs of
each majority gate. Gate indices and counts use logarithmic
space, and circuit depth is fixed. Repeated evaluation therefore
takes polynomial time and logarithmic extra space. Here its
virtual input has polynomial size in `n`; each input-bit query
is answered by the factor-regeneration routine above. This
composition takes polynomial time in `n` and only `O(log n)`
extra space. The complete factor list is never written.

Check feasibility and compare two integers by regenerating their
bits from the common bound downwards. Retain the best *height
descriptor*, not its `n`-bit integer or divisor count. After the
search, find the highest nonzero product bit and stream the
output. To obtain `q(n)`, also check whether any lower bit of
`D(n)` is one, distinguishing exact powers of two correctly.

Thus each candidate costs polynomial time and
`O(log n loglog n)` total workspace, including the two packed
descriptors. Combining this with (13) proves (1)--(2). Output
length is polynomial in `n` and is charged to time.

### The entirely elementary alternative

The balanced carry recurrence in the
[primorial-weighted note](primorial-weighted-divisor-preprocessing.md)
computes a requested product bit in `2^{O(log^2 n)}` time and
`O(log^2 n)` extra space. Substituting it in the same candidate
search gives those bounds for the entire computation. This is
the method implemented by `--streamed`; it does not claim the
faster arithmetic consequence of the cited theorem.

## 5. Reproduction and scope

Run:

```sh
python3 tools/check_exchange_divisor_preprocessing.py \
  --output evidence/exchange-divisor-preprocessing.json
python3 tools/exchange_divisor_preprocessing.py 16 32 64 128
python3 tools/exchange_divisor_preprocessing.py --streamed 1 2 3 4 5 6
```

The checker compares the answer with a divisor sieve and an
independent unrestricted exponent search. It checks the exchange
inequality on every record divisor maximum in the sieve, verifies
literal replacement witnesses when a candidate violates it,
checks exact height caps and descriptor reconstruction, and
compares the elementary streamed search with the ordinary
implementation. The 128-bit maximum is an optional larger check.

These finite checks validate the implemented pruning, algebra,
and arithmetic interfaces. The asymptotic search bound is proved
in (13), and the sharper simultaneous time/space theorem uses the
stated classical arithmetic result. No Lean claim or factoring
exponent improvement follows from these checks.

This note and its implementation used OpenAI assistant support.
