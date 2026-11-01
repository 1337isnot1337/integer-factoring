# Exact contraction of a Frobenius rational star

## Status and scope

This note gives a deterministic exact algorithm for one coupled finite-field
sum. A product of rational Frobenius binomials can be summed over all
`2^w` values of its shared variable using a number of field operations
polynomial in the number of factors and in `w`. Repeated poles and zero
numerator coefficients are allowed. The algorithm does not enumerate the
field being summed over.

This is a positive local contraction result. Its output need not retain
the same factor form in a second variable. Section 6 states that remaining
closure issue. Neither this note nor the finite replay changes the complete
factoring exponent.

## 1. The rational-star theorem

Let `b=2^w`, let `E=F_(2^d)` contain `F_b`, and suppose an explicit binary
basis and multiplication algorithm for `E` are given. Thus `w` divides `d`.
Consider

\[
 R(U)=c\,
 \frac{\prod_{i=1}^{r_+}(a_i+b_iU^{2^{e_i}})}
      {\prod_{j=1}^{r_-}(c_j+d_jU^{2^{f_j}})},
 \qquad 0\le e_i,f_j<w.                                      \tag{1}
\]

Here every `d_j` is nonzero, and every root of a denominator factor lies
outside `F_b`. Assume the formal degree inequality

\[
 \sum_i2^{e_i}\le\sum_j2^{f_j}.                              \tag{2}
\]

The numerator slopes `b_i`, numerator constants `a_i`, and scalar `c` may
vanish. Put `r=r_++r_-+1`. There is a deterministic algorithm computing

\[
 \sum_{u\in F_b}R(u)
 \quad\hbox{and}\quad
 \sum_{u\in F_b^\times}R(u)                                 \tag{3}
\]

in `poly(r,w,d)` bit operations when field operations have polynomial bit
cost. One explicit implementation uses at most

\[
 O\!\left(r^2(r+\log(r+1))^2
                 (w+\log(r+1))\right)                       \tag{4}
\]

field additions and multiplications, apart from a polynomial number of
inversions and Frobenius operations. The latter also have polynomial bit
cost: inversion is exponentiation to `2^d-2`, and any Frobenius power needs
at most `d` squarings after its exponent is reduced modulo `d`.

The explicit field is an input to this theorem. A surrounding construction
must separately charge field construction; no field of exponential degree
or unit-cost formal root is implicit here.

If `r_-=0`, condition (2) makes the rational function constant, so its
full-field sum is zero and its nonzero-field sum is that constant. We now
assume at least one denominator factor.

## 2. Pole traces reduce to short coefficient queries

Each denominator factor has the form

\[
 c_j+d_jU^{2^{f_j}}
   =d_j(U+\alpha_j)^{2^{f_j}},\qquad
 \alpha_j=(c_j/d_j)^{2^{-f_j}},                              \tag{5}
\]

where the inverse Frobenius power is taken in `E`. Group equal roots.
For a distinct root `alpha`, let

\[
 m_\alpha=\sum_{j:\alpha_j=\alpha}2^{f_j},\qquad
 H_\alpha(Z)=Z^{m_\alpha}R(\alpha+Z),\qquad
 \Delta_\alpha=\alpha^b+\alpha\ne0.                        \tag{6}
\]

The function `H_alpha` is regular at zero. It is again a scalar times a
product of Frobenius binomials divided by such binomials with nonzero
constant terms. Denominator factors at `alpha` cancel exactly against the
displayed power of `Z`. Cancellations in the numerator simply produce
initial zero coefficients of `H_alpha`; they cause no difficulty.

For every `alpha` outside `F_b`, logarithmic differentiation of
`prod_(u in F_b)(Y+u)=Y^b+Y` gives the formal-series identity

\[
 \sum_{u\in F_b}\frac1{\alpha+u+Z}
       =\frac1{\Delta_\alpha+Z+Z^b}.                       \tag{7}
\]

Consequently the contribution of the principal part at `alpha` to (3) is

\[
 [Z^{m_\alpha-1}]
      \frac{H_\alpha(Z)}{\Delta_\alpha+Z+Z^b}.             \tag{8}
\]

To see the coefficient indexing explicitly, a term
`c_t/(U+alpha)^t` contributes `c_t Z^(m_alpha-t)` to `H_alpha` below
degree `m_alpha`. Its contribution to (8) is `c_t` times the coefficient
of `Z^(t-1)` in (7), exactly its sum over `u`.

By (2), the polynomial part of the partial-fraction expansion has degree
at most zero. Its full-field sum vanishes because `b=0` in characteristic
two. Adding (8) over the distinct roots therefore gives the full sum.
Adding `R(0)` then gives the sum over `F_b^times`.

Since `m_alpha <= r_- b/2`, only `O(r)` terms of

\[
 \frac1{\Delta+Z+Z^b}
   =\sum_{v\ge0}\frac{Z^{bv}}{(\Delta+Z)^{v+1}}           \tag{9}
\]

can affect (8). For each relevant `v`, write `v+1` in binary. Then

\[
 (\Delta+Z)^{-(v+1)}
   =\prod_{j:\operatorname{bit}_j(v+1)=1}
                 (\Delta^{2^j}+Z^{2^j})^{-1}.              \tag{10}
\]

Thus (8) reduces to `O(r)` coefficient queries, each involving `O(r+log r)`
Frobenius binomials or inverse binomials. The queried degree is at most
`rb`, so its binary length is `O(w+log r)`. Section 3 computes each query
without expanding a polynomial of that degree.

## 3. Binary carries compute the coefficient

After dividing denominator factors by their nonzero constant terms, the
required form is

\[
 [Z^K]\prod_i(A_i+B_iZ^{2^{e_i}})
                \prod_j(1+C_jZ^{2^{f_j}})^{-1}.            \tag{11}
\]

Exponents here may exceed `w`; only their maximum and the binary length
of `K` matter. Constants `A_i=0` are permitted. In particular, we never
divide by a numerator constant.

In characteristic two, as formal power series,

\[
 (1+CZ^{2^f})^{-1}
     =\prod_{t\ge0}(1+C^{2^t}Z^{2^{f+t}}).                \tag{12}
\]

At binary position `ell`, introduce one binary choice for each inverse
factor with `f_j<=ell`, with weights `1` and `C_j^(2^(ell-f_j))`. Also
introduce one choice for each numerator factor with `e_i=ell`, with
weights `A_i` and `B_i`. The polynomial in a counting variable `Y`

\[
 E_\ell(Y)=
 \prod_{j:f_j\le\ell}(1+C_j^{2^{\ell-f_j}}Y)
 \prod_{i:e_i=\ell}(A_i+B_iY)                              \tag{13}
\]

records the total weight of choosing exactly `h` ones at that position
as its coefficient of `Y^h`.

A dynamic program stores the incoming integer carry `q`. For each
coefficient indexed by `h`, retain the transition exactly when

\[
 q+h\equiv \operatorname{bit}_\ell(K)\pmod2,
 \qquad q'=(q+h-\operatorname{bit}_\ell(K))/2.              \tag{14}
\]

Multiply its weight by `[Y^h]E_ell`. Start with carry zero and weight one;
after the last bit of `K`, accept only carry zero. Numerator factors with
exponent above that last bit contribute only their constants. Inverse
choices above it are forced to zero because all degrees are nonnegative.
Use one bit even when `K=0`.

If there are `s` factors in (11), every carry lies in `[0,s]`. Constructing
(13) and updating all states costs `O(s^2)` field operations per bit.
The coefficient therefore costs `O(s^2(log(K+1)+1))` operations. In (8)
there are at most `r_-` distinct poles and `O(r)` values of `v`, giving
(4). This proves the theorem.

For pairwise distinct denominator roots, every `m_alpha=2^f<b`, so only
`v=0` occurs in (9). Moreover all branches grouping equal roots disappear.
On that promised domain the formula is a polynomial-size rational
arithmetic circuit in the factor coefficients, with inversions only of
quantities proved nonzero. Over a fixed explicit `E`, powering to
`2^d-2` replaces each inversion by a polynomial-size division-free circuit.

## 4. A radix box on complementary multiplicative groups

The theorem is relevant to the geometric input-box factors that occur
after additive Fourier expansion. Let `u in F_b^times` and let `v` satisfy
`v^(b+1)=1` in `F_(b^2)`. Put

\[
 K(u,v)=\sum_{d=0}^{b-1}(uv)^d.                             \tag{15}
\]

Since `u^b=u` and `v^b=v^(-1)`, away from `uv=1` one has

\[
 K(u,v)=\frac{u+v}{v(1+uv)}.                               \tag{16}
\]

At `uv=1`, which within these two groups means `u=v=1`, the geometric sum
is `b=0`; the rational ratio in (16) is undefined. It must be treated as
that separate value. A nonzero external tag can remove such poles when
its field of definition is suitably disjoint from `F_(b^2)`; checking
that promise is part of the input to Theorem 1.

The exact two-factor contractions, including the exceptional point, are

\[
 \begin{split}
 \sum_{u\in F_b^\times}K(u,v)&=1+v^{-2},\\
 \sum_{v^{b+1}=1}K(u,v)&=1,\\
 \sum_{u\in F_b^\times}K(u,v_1)K(u,v_2)
  &=1+(v_1v_2)^{-2}
       +(v_1v_2)^{-1}(1+\delta_{v_1,v_2}),\\
 \sum_{v^{b+1}=1}K(u_1,v)K(u_2,v)
  &=1+u_1u_2+\delta_{u_1,u_2}u_1^2.
 \end{split}                                               \tag{17}
\]

For the first and third identities, multiplicative character
orthogonality keeps the digit exponents whose sum is a multiple of
`b-1`. The possible totals are `0,b-1,2b-2`. The middle geometric sum is
zero when `v_1=v_2`, and equals `(v_1v_2)^(-1)` otherwise. For the second
and fourth identities the possible totals divisible by `b+1` are `0`
and `b+1`; the latter gives `u_1u_2` for distinct `u_1,u_2` and vanishes
when they coincide. These arguments include `b=2`, where the final two
nonconstant terms in the fourth identity cancel.

More generally, for a nonzero tag `gamma` and `0<=e<w`,

\[
 \sum_{d=0}^{b-1}(\gamma u^{2^e}v)^d
 =\frac{1+\gamma^b u^{2^e}v^{-1}}
        {1+\gamma u^{2^e}v},                               \tag{18}
\]

whenever the denominator is nonzero. Products of these factors have
equal formal numerator and denominator degrees and fit (1), if every
pole lies outside `F_b`. The theorem therefore contracts an entire star
of such tagged input boxes jointly; it is stronger than evaluating a
single box at one fixed frequency.

## 5. Field and representation costs

The binary-carry algorithm uses an explicit field basis and exact field
equality. It needs neither factorization of an integer group order nor
a discrete logarithm. Its only root extraction is inverse Frobenius,
which is repeated squaring in that basis.

If an application constructs `E` by exhaustive search for an irreducible
binary polynomial of degree `d`, that preparation has a separately
charged upper bound `2^O(d)`. It must not be replaced by polynomial work
in `d` without an appropriate construction. Once `E` is available, the
theorem itself has polynomial bit complexity. The sum is represented as
one field element or by the circuit described after (14), not as a
table containing its values for all remaining variables.

## 6. What is still needed for a whole contraction algorithm

The output includes inverses of pole differences and of
`alpha^b+alpha`. Suppose two poles depend on a second frequency variable
`V` as `a V^(2^e)` and `c V^(2^f)`. Their difference is

\[
 aV^{2^e}+cV^{2^f}.                                       \tag{19}
\]

For `e!=f`, factoring out a monomial leaves an exponent
`2^e-2^f`, generally not a power of two. Thus the output does not
automatically have the form (1) as a function of `V`. A large number of
such denominators can also become coupled through the coefficients of
the carry dynamic program.

An iteration must prove closure under these new denominators, or provide
another charged algorithm for their joint sums. Evaluating the circuit
at each remaining assignment would restore the original large table
cost. The local contraction theorem supplies no permission to omit that
cost, and does not establish a factoring exponent above `1/2`.

## 7. Exact replay

The companion checker is
[`check_frobenius_rational_trace.py`](../tools/check_frobenius_rational_trace.py).
It compares the carry coefficient algorithm with direct truncated power
series, and the rational-star sum with full finite-field enumeration. It
also checks (16)--(18), repeated poles with multiplicity at least `b`,
zero numerator constants, cancellation at a pole, proper rational
functions, the nonzero-field correction, and rejection of a pole inside
the summation field. Run:

```sh
python tools/check_frobenius_rational_trace.py --output evidence/frobenius-rational-trace.json
```

The committed [evidence](../evidence/frobenius-rational-trace.json) records
1,920 coefficient comparisons, 400 complete rational traces checked
against 3,000 individual field evaluations, 188 poles with multiplicity
at least `b`, and 194 nonzero higher terms from (9). It also records
404 and 284 checks of the two Gram identities, 332 Cauchy-box identities,
137 tagged-box identities, 58 pole-cancellation controls, and 30 rejected
denominator poles inside the summation field. All comparisons passed.

The mathematical proof above is self-contained. The finite replay checks
the formulas and implementations; it does not replace the proof for
arbitrary dimensions.
