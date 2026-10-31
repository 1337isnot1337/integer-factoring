# A witness for the tapered Hermite precision scale

## Result and boundary

The polynomial-separation evaluator in
[OpenAI's October 25, 2026 paper](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`, recovers a polynomial's
value at zero using a complete two-adic residue grid and tapered Hermite
data. Its added precision has order `MD/u`, where `M` is the wire scale,
`u` is the retained number of valuation levels, and `D` depends on the
chosen center.

That scale is necessary for the polynomial class assumed in its
Hermite lemma. The witness below has value one at zero, yet every
coefficient around every retained center is divisible by
`2^Omega(MD/u)`. Therefore replacing the interpolation routine while
using only those lower-precision residue data cannot recover the value
at zero, even if the full polynomial at every center is supplied.

This result does **not** prove that a contraction polynomial actually
realizes the witness. It does not exclude an algorithm that exploits
more structure of those polynomials, changes the deformation or grid,
uses a different set of data, or computes through another reduction.
It gives no lower bound for factoring. The factoring score is unchanged.

## 1. The admissible polynomial

Let

\[
L=2^a,\quad S=\{-L,\ldots,L-1\},\quad
\alpha=2^{a+1},\quad 1\le u\le a,\quad M\ge u.
\]

Write

\[
h=\lfloor M/u\rfloor,\qquad
m_i=h\max(u-i,0).
\]

Define a rational polynomial

\[
\boxed{F(T)=\prod_{q\in S\setminus\{0\}}
           \left(1-\frac{T}{q}\right)^{m_{v_2(q)}}.}
\tag{1}
\]

It has `F(0)=1`. The grid contains `L/2^i` nonzero elements of
valuation `i`, so

\[
\begin{split}
\deg F
 &=Lh\sum_{i=0}^{u-1}\frac{u-i}{2^i}\\
 &=Lh(2u-2+2^{1-u})\\
 &\le L(M/u)(2u-1)\\
 &\le M(2L-1).
\end{split}
\tag{2}
\]

The last step uses `L>=u`. This is the degree bound required in
the source Hermite lemma. Formula (1) is a single explicit polynomial;
its existence does not rely on a dimension count or an asymptotic limit.

## 2. Exact local coefficient content

For a nonzero rational polynomial `G(z)`, its two-adic coefficient
content is the minimum valuation of its nonzero coefficients. This
content is additive under multiplication: after dividing out the two
contents, reduction modulo two multiplies two nonzero polynomials
over a field, and their product remains nonzero.

Fix `t in S` and substitute `T=t+alpha z`. For `q != t`, the
linear factor in (1) has content

\[
v_2(q-t)-v_2(q),
\]

since all nonzero grid differences have valuation at most `a`.
For `q=t`, its content is `a+1-v_2(t)`.

If `t=0`, or if `v_2(t)>=u`, every factor with nonzero multiplicity
has content zero. Thus `F(t+alpha z)` is integral in those cases.

Now let `i=v_2(t)<u` and put `D_i=L/2^i`. The exact content is

\[
\boxed{
A_i=m_i(a+1-i)
 +hD_i\sum_{j=1}^{a-i}\frac{j\min(j,u-i)}{2^j}.
}
\tag{3}
\]

To derive it, separate the factor at `q=t`. Points of valuation
below `i` contribute zero. For each `j>=1`, there are `D_i/2^j`
points with `v_2(q)=i+j`, and equally many points with
`v_2(t-q)=i+j`. In the latter set `v_2(q)=i`.
Pairing those contributions gives

\[
\sum_{q\ne0,t}m_{v_2(q)}
       (v_2(t-q)-v_2(q))
=D_i\sum_{j=1}^{a-i}\frac{j}{2^j}(m_i-m_{i+j}).
\]

Since `m_i-m_(i+j)=h min(j,u-i)`, equation (3) follows.
Every summand in (3) is nonnegative. Consequently

\[
F(t+\alpha z)\in\mathbb Z_{(2)}[z]
\qquad\text{for every }t\in S.
\tag{4}
\]

These are all the local integrality hypotheses in the source
Hermite lemma. The coefficient ring may be extended exactly as
in that lemma; the rational witness embeds in it.

## 3. Indistinguishability below the precision scale

The `j=1` term of (3) is present because `i<u<=a`. It gives

\[
A_i\ge\frac{hD_i}{2}\ge\frac{MD_i}{4u},
\tag{5}
\]

where `floor(M/u)>=M/(2u)` for `M/u>=1`.

Suppose that, for every retained center `t` with `v_2(t)=i<u`,
one is supplied only the residues of `F(t+alpha z)` modulo
`2^(P_t)`, with

\[
1\le P_t\le\left\lfloor\frac{MD_i}{4u}\right\rfloor.
\tag{6}
\]

Equations (3)--(5) show that **all these residue polynomials are
zero**. The same data are obtained from the identically zero
polynomial. Both inputs satisfy the degree and integrality
hypotheses, but their values at zero differ modulo two.

No recovery procedure based only on the data (6) can distinguish
them. Allowing more jet coefficients does not help, because the
entire local polynomial vanishes at that precision. In particular,
a uniform replacement of the source's order `MD_i/u` precision by
`o(MD_i/u)` across all retained centers cannot solve the same
general interpolation problem.

The result says that at least one center must escape the lower
precision envelope for universal recovery. It does not say that
every center individually needs that much precision, nor that the
source's numerical constants are optimal.

## 4. The next sufficient improvement

For the source's mode multiplicity `O(M/D_i)`, a valid precision
bound `O(MD_i/u^(1+delta))`, with fixed `delta>0`, would change
the corresponding partition estimate to

\[
\exp(O(M/u^{(1+\delta)/2})).
\tag{7}
\]

That observation is only a local sufficient condition. A complete
algorithm would also have to preserve the source's group separation,
polynomially bounded descriptions and coefficient precision, entry
query cost, and recursive wire reduction. Formula (1) shows that
such a precision improvement cannot hold for the same general
polynomial class and the same retained residue grid.

A concrete remaining problem is to construct a different deformation
with many specializations that are local under suitable balanced
groups, but whose evaluation centers avoid this two-adic loss.
Extension-field or root-of-unity centers can have unit pairwise
differences; preserving locality at enough such centers is the
additional requirement. The current source's integer-difference and
bit-flip construction does not automatically supply it.

## 5. Evidence and reproduction

```sh
python3 tools/check_tapered_hermite_precision.py --output evidence/tapered-hermite-precision.json
```

The checker expands the shifted polynomial by ordinary integer
multiplication, using its common denominator only after expansion.
It compares the actual coefficient content to (3), checks every
coefficient's integrality, and verifies `F(0)=1` and the degree
bound. A separate broader parameter sweep checks the closed
integer formulas. The finite tests support the all-size proof;
they are not an extrapolation from small examples.

See [`check_tapered_hermite_precision.py`](../tools/check_tapered_hermite_precision.py)
and [`tapered-hermite-precision.json`](../evidence/tapered-hermite-precision.json)
for the exact cases and recorded counts. No Lean development is
included. This note and checker were prepared with substantial
OpenAI assistant support; the source evaluator is credited above.
