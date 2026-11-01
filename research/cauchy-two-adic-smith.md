# Exact two-adic channels for Cauchy blocks

## Result and scope

The Cauchy matrix

\[
H^{(d)}_{ab}=\frac1{c+2(a-b)},\qquad 0\le a,b<d,
\tag{1}
\]

with odd integer `c` has the exact two-adic Smith exponents

\[
\boxed{\sigma_k=2k+2v_2(k!)=4k-2s_2(k),\quad 0\le k<d.}
\tag{2}
\]

Here `s_2(k)` is the number of ones in the binary expansion of `k`.
An explicit binomial `LDU` factorization supplies the channels, including
their units. All formulas work over the localization `Z_(2)`, over its
finite quotients, and after extension to the unramified coefficient ring
used in the source evaluator. The shift `c` may be any unit of that ring.

These blocks occur in the polynomial-separation evaluator of
[OpenAI's October 25, 2026 paper](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
at the pinned source revision `4eaee87548250b2682affda7d05d947c05554bae`.
The source uses geometric modes of linearly increasing valuation.
Equation (2) gives exact channels for a complete same-tile block and
exposes a stronger linear slope, asymptotically four. It does **not**
give a replacement factorization of the entire crossing operator with
that slope: a complete column tile also interacts with other row tiles.

Even granting the exact spectrum, the associated weighted-monomial
enumeration retains the scale

\[
\exp\!\left(\Theta\!\left(\frac{M}{\sqrt u}\right)\right)
\tag{3}
\]

when the precision is `Theta(MD/u)` and there are `Theta(M/D)`
independent blocks of dimension `Theta(D)`. This is a limitation of
enumerating all monomials permitted by the valuation budget. It is
not a lower bound for tensor contraction, for the source algorithm
under more structural changes, or for factoring.

The factoring score remains `1/2`. We do not claim that the classical
Cauchy factorization or its Smith specialization is new in the
literature. The contribution here is its exact form, implementation,
independent checks, and stated consequence for this research route.

## 1. Explicit integral factor channels

Use zero-based indices. Define

\[
\begin{split}
L_{ak}
 &=\binom ak\,
   \frac{\prod_{j=0}^k(c+2j)}
        {\prod_{j=0}^k(c+2(a-j))}, && a\ge k,\\
U_{kb}
 &=\binom bk\,
   \frac{\prod_{j=0}^k(c-2j)}
        {\prod_{j=0}^k(c+2(j-b))}, && b\ge k,\\
d_k
 &=\frac{(-4)^k(k!)^2}
        {c\prod_{j=1}^k(c^2-4j^2)}.
\end{split}
\tag{4}
\]

Set the remaining entries of `L,U` to zero. Every displayed
denominator is odd, and `L_kk=U_kk=1`. Thus `L` and `U` are
invertible integral triangular matrices over `Z_(2)`.

**Factorization theorem.**

\[
H=L\,\operatorname{diag}(d_0,\ldots,d_{d-1})\,U.
\tag{5}
\]

Here is an elimination proof that also fixes all signs. For a general
Cauchy matrix `1/(x_a+y_b)`, after eliminating indices below `k`,
the remaining entry is

\[
S^{(k)}_{ab}=
\frac{\prod_{j<k}(x_a-x_j)(y_b-y_j)}
     {(x_a+y_b)\prod_{j<k}(x_a+y_j)(x_j+y_b)}.
\tag{6}
\]

At `k=0` this is the original matrix. Subtracting
`S^(k)_ak S^(k)_kb / S^(k)_kk` from `S^(k)_ab`
inserts the next factors in (6); the needed numerator identity is

\[
(x_a+y_k)(x_k+y_b)-(x_a+y_b)(x_k+y_k)
=(x_a-x_k)(y_b-y_k).
\]

Take `x_a=c+2a` and `y_b=-2b`. Formula (6) gives the diagonal
entry `d_k` in (4). Its column and row multipliers are exactly
`L_ak=S^(k)_ak/d_k` and `U_kb=S^(k)_kb/d_k`.
This proves (5). The same identity holds for any rectangular
matrix with consecutive row and column indices: retain
`min(rows,columns)` diagonal channels.

The units in (4) show

\[
v_2(d_k)=2k+2v_2(k!).
\]

Counting the factors of two in `k!` gives

\[
v_2(k!)=\sum_{j\ge1}\lfloor k/2^j\rfloor=k-s_2(k).
\]

Moreover,

\[
\sigma_{k+1}-\sigma_k=2+2v_2(k+1)>0.
\]

Absorb the odd unit of each `d_k` into one triangular factor.
The diagonal entries are then `2^(sigma_k)` in divisibility order,
which proves the Smith assertion (2). In particular, the smallest
valuation of a nonzero `r` by `r` minor is

\[
\sum_{k=0}^{r-1}\sigma_k
=r(r-1)+2\sum_{k=0}^{r-1}v_2(k!).
\tag{7}
\]

There is also a direct check of (7). The Cauchy determinant
of any such minor has valuation

\[
r(r-1)+v_2\!\prod_{i<j}(a_j-a_i)
             +v_2\!\prod_{i<j}(b_j-b_i).
\]

For distinct nonnegative integers, each Vandermonde product is
divisible by `prod_(k=0)^(r-1) k!`: the quotient is the determinant
of the integer matrix with entries `binom(a_i,k)`. Consecutive
indices attain equality. This yields the same minimum without
using the factorization.

## 2. Finite precision and bit cost

Modulo `2^P`, a scalar multiple `2^beta H` has precisely the
nonzero diagonal channels

\[
\beta+\sigma_k<P.
\tag{8}
\]

The equality (5) still holds, since only odd denominators are
inverted. Formula (8) describes its exact Smith profile; it does
not merely count nonzero terms in a selected geometric expansion.

The matrices in (4) are constructive. Their rational numerators
and denominators have bit length polynomial in `d+log(|c|+2)`.
Alternatively, compute products and inverses directly modulo
`2^P`, after separating the known powers of two in each pivot.
This uses polynomially many bit operations in `d,P,log(|c|+2)`.
The ordinary integer implementation below uses exact rational
numbers and therefore provides both an exact identity checker
and a finite-precision channel oracle. No unproved fast Smith
routine is part of the construction.

The identity is universal in the indeterminate `c` after inverting
the odd expressions that appear in (4). Consequently it remains
valid in an unramified extension of `Z_(2)` when `c` is a unit.
Each such denominator is then a unit, while the valuation of
the explicit scalar pivot is unchanged. The code checks integer
shifts; the ring extension follows from this identity argument.

## 3. Where the block occurs in the source evaluator

Use the notation in the source's *Modes for interactions between
groups* subsection. Its packet slots have size `L=2^a`, the
deformation is evaluated at

\[
T=t+2^{a+1}\zeta,\qquad i=v_2(t),
\]

and its Cauchy matrix is

\[
W_{xy}(T)=\frac{\mathcal C(T)}
                  {\mathcal C'(0)(T-y+x)}.
\tag{9}
\]

Choose one residue `x_0` modulo `2^(i+1)` and the matching row
and column tiles

\[
x=x_0+2^{i+1}r,\qquad y=x_0+2^{i+1}s.
\]

Both tiles have size `d=L/2^(i+1)`. If the packet's two symbols
in coordinate `i` differ, equal `i`th bits put these row and
column ports in different groups. Thus this is a crossing block.
On it,

\[
T-y+x=2^i(c+2(r-s)),
\quad
c=t/2^i+2^{a+1-i}\zeta,
\tag{10}
\]

and `c` is a unit. If `zeta` is a unit, the source's complete
residue grid gives

\[
v_2(\mathcal C(T)/\mathcal C'(0))=a+1.
\]

Indeed the factor at `t` contributes `a+1`, while all other
factor valuations cancel against `mathcal C'(0)` because they
are the same as at `t`. Therefore (9) on this tile is an odd
unit times

\[
2^\beta H^{(d)},\qquad
\beta=a+1-i=\log_2d+2.
\tag{11}
\]

The result describes a complete tile with all its slots. Missing
real ports or additional tensor restrictions can lower its
effective dimensions. It gives no lower bound for an arbitrary
input with such restrictions.

## 4. Exact channel compression still leaves a square-root partition scale

Consider `g` independent copies of a dimension-`d` block with
weights

\[
w_k=\beta+\sigma_k,
\qquad 0\le k<d,\quad\beta\ge1.
\]

Let `A(P)` be the number of nonnegative integer mode tuples
whose total weighted degree is at most `P`. The next statement
concerns this explicitly defined enumeration problem.

**Mode-count theorem.** Let `u>=1`, `h=ceil(sqrt(u))`,
`d>=16h`, `1<=beta<=d/(4h)`, and

\[
P=\left\lfloor\frac{gd^2}{u}\right\rfloor.
\]

Then, with absolute positive constants,

\[
\boxed{\log A(P)=\Theta(gd/\sqrt u).}
\tag{12}
\]

For the upper bound, `sigma_k>=2k`, so each integer weight has
at most `g` variables. For every `z=exp(-t)`, `t>0`,

\[
A(P)\le e^{tP}\prod_{j\ge1}(1-e^{-tj})^{-g}.
\]

Expanding the logarithm and using
`sum_(j>=1) exp(-t a j) <=1/(ta)` gives

\[
\log A(P)\le tP+\frac{g}{t}\sum_{a\ge1}\frac1{a^2}.
\]

Take `t=sqrt(g/P)`. This is `O(sqrt(gP))`, hence
`O(gd/sqrt(u))` for the stated budget.

For the lower bound, put `K=floor(d/(8h))`. In each copy,
independently choose any subset of its first `K` modes, using
exponent zero or one. All `2^(gK)` resulting tuples fit the budget,
because `sigma_k<=4k` and

\[
g\sum_{k=0}^{K-1}w_k
\le g(\beta K+2K(K-1))
\le\frac{gd^2}{16h^2}
\le P.
\tag{13}
\]

The last inequality follows from `gd^2/u>=256`; the factor
`1/16` absorbs the floor. Also

\[
K\ge d/(16h)\ge d/(32\sqrt u).
\]

This proves the lower estimate in (12). These retained tuples
already use at most `gK<=gd` factors. A global degree cap of
`gd` cannot discard this lower-count family.

In the source regime `d` is exponential in the retained number
`u` of valuation levels, while `beta=O(log d)`. Thus the
hypotheses hold eventually. Replacing `gd` by a wire scale `M`
gives (3). Exact Smith compression changes constants in this
model but does not turn the square-root logarithmic saving into
a higher power.

This conclusion has a precise boundary. The lower count includes
every valuation-admissible monomial. An algorithm can attempt
to exploit relations among modes, special input tensors,
cancellation among monomials, a smaller Hermite precision loss,
or a different recursive aggregation. Equation (12) excludes
none of those possibilities. It is not a hardness theorem and
does not justify a claim that the half-rate factoring bound is
optimal.

## 5. Reproduction and checked evidence

Run

```sh
python3 tools/check_cauchy_two_adic_smith.py --output evidence/cauchy-two-adic-smith.json
```

The recorded run passes in under one second in the research
environment. It checks:

- 1,281 exact rational entries against the explicit rectangular
  factorization, including positive and negative odd shifts;
- 2,769 exact minors by independent rational elimination, with
  all minor sizes through six;
- 126 finite-chain-ring profiles by independent minimum-valuation
  pivoting, through dimension 24 and precision 80;
- 16,590 modular factorization entries;
- 80 weighted generating-function counts against literal tuple
  enumeration;
- 1,536 exact parameter checks of the lower-count construction.

The full implementation is
[`cauchy_two_adic_smith.py`](../tools/cauchy_two_adic_smith.py).
The independent checks are
[`check_cauchy_two_adic_smith.py`](../tools/check_cauchy_two_adic_smith.py),
with recorded output in
[`cauchy-two-adic-smith.json`](../evidence/cauchy-two-adic-smith.json).
The finite checks support the stated all-size proofs; they are
not used to extrapolate an asymptotic statement.

No Lean development is included. This note and its code were
prepared with substantial OpenAI assistant support. The source
evaluator and its polynomial-separation construction are credited
to the pinned OpenAI paper above.
