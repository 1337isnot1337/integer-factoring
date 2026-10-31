# Arithmetic Cauchy permanents and their complete squared-matrix Smith profiles

## Result and scope

The power-of-two restriction in the permanent theorem of
[the Walsh note](walsh-principal-permanents.md) can be removed. For every
integer `n >= 1`, every even integer `d`, and every odd two-adic unit `T`,

\[
 C_{ij}=\frac1{T+d(i-j)},\qquad 0\le i,j<n,
\]

satisfies

\[
 \boxed{\quad
 \mathrm{per}C\equiv n!T^{-n}
       \pmod {2^{v_2(n!)+1}},\qquad
 v_2(\mathrm{per}C)=v_2(n!)=n-s_2(n).
 \quad}                                                    \tag{1}
\]

Here `s_2(n)` is the number of ones in the binary expansion of `n`.
The argument applies over unramified two-adic coefficient rings and
coefficientwise over their formal power-series rings. Section 4 then
proves the previously conjectured individual Smith orders of the squared
Cauchy matrix, using an integral orthogonal-polynomial basis.

Consequently, every even-spacing principal packet mode of size `n` has
exact local content

\[
 n(a+1)+v_2(n!)=n(a+2)-s_2(n),\qquad L=2^a,                 \tag{2}
\]

at every odd retained center. Previously the all-size proof covered
power-of-two `n`; the non-power-of-two cases were finite observations.
This is a scalar valuation theorem. It supplies no new crossing rank
bound, no general scalar decoder, and no improvement to the factoring
score `alpha=1/2`.

## 1. A denominator bound for finite sampling moments

Fix positive integers `k_1,...,k_r`, put `K=sum_j k_j`, and let `n>=r`.
Sample an ordered tuple `(I_1,...,I_r)` uniformly from the distinct
elements of `{0,...,n-1}`. Then, for every prime `p`,

\[
 M_n(k_1,\ldots,k_r)
   =\frac1{(n)_r}\sum_{\substack{0\le i_1,\ldots,i_r<n\\
                                i_1,\ldots,i_r\ \text{distinct}}}
                          \prod_j i_j^{k_j}
 \quad\hbox{satisfies}\quad
 v_p(M_n)\ge-\left\lfloor\frac K{p-1}\right\rfloor .        \tag{3}
\]

We use `(n)_r=n(n-1)...(n-r+1)` and set `v_p(0)=+infinity`.

### Proof by equality patterns and finite orders

Introduce `r` distinguished root vertices. Attach `k_j` labeled leaves
to root `j`; there are exactly `K` leaf-to-root edges. For a fixed root
value `i_j`, the factor `i_j^(k_j)` counts the assignments of values in
`{0,...,n-1}` to these leaves for which each leaf value is strictly less
than its root value. The roots must have distinct values.

Partition all `r+K` vertices according to equality of their values.
Discard a partition if two roots occupy the same block or if its
leaf-to-root inequalities create a directed cycle, including a loop.
Every remaining partition has `v>=r` blocks and defines a partial order
`P` on those blocks. If `e(P)` is its number of linear extensions, the
number of injective assignments of values to its blocks is

\[
 e(P)\binom nv .                                           \tag{4}
\]

Choose the `v` values, order them, and use the linear extension to assign
them. Every assignment counted in the original moment has exactly one
equality partition, so these contributions do not overlap.

After division by `(n)_r`, its contribution is

\[
 \frac{e(P)\binom nv}{(n)_r}
   =\frac{e(P)}{v!}(n-r)_{v-r}.                            \tag{5}
\]

This also holds when `n<v`, since the falling factorial on the right
is then zero. Its last factor is an integer.

Let the undirected connected components of the quotient graph have
orders `v_1,...,v_c`, and let `P_j` be the induced partial orders. The
possible interleavings of their linear extensions give

\[
 \frac{e(P)}{v!}=\prod_{j=1}^c\frac{e(P_j)}{v_j!}.         \tag{6}
\]

Each component requires at least `v_j-1` distinct graph edges. The
quotient has at most the original `K` edges. Legendre's formula gives

\[
 \sum_j v_p(v_j!)
 =\sum_j\frac{v_j-s_p(v_j)}{p-1}
 \le\frac{\sum_j(v_j-1)}{p-1}
 \le\frac K{p-1}.                                         \tag{7}
\]

The integers `e(P_j)` cannot worsen this denominator bound. Equations
(5)--(7) prove (3) for each contribution, and hence for their sum.
The proof does not assume that the roots are independent samples.

## 2. A permanent congruence with arbitrary column shifts

Let `O` be an unramified complete discrete valuation ring of residue
characteristic `p`, with valuation normalized by `v_p(p)=1`. Take
`d in pO` and units `y_1,...,y_n`. Define

\[
 A_{ij}=(y_j+di)^{-1},\qquad 0\le i<n.                     \tag{8}
\]

The following statements hold.

1. If the units `y_j` all have the same residue, then

   \[
    \frac{\mathrm{per}A}{n!}
       \equiv\prod_j y_j^{-1}\pmod p .                    \tag{9}
   \]

2. When `p` is odd, the common-residue assumption is unnecessary.
   It is also unnecessary when `v_p(d)>=2`.

In particular the quotient in (9) is a unit. The theorem covers every
matrix size; column shifts need not be distinct.

### Expansion and integrality

Let `sigma` be a uniformly random bijection from the columns to the
row indices. Expanding each geometric series gives

\[
 \frac{\mathrm{per}A}{n!}
  =\left(\prod_jy_j^{-1}\right)
   \sum_{k_1,\ldots,k_n\ge0}
        (-d)^K\left(\prod_j y_j^{-k_j}\right)
        {\mathbb E}\left[\prod_j\sigma(j)^{k_j}\right],
 \qquad K=\sum_jk_j .                                     \tag{10}
\]

Zero exponents can be removed from a moment. For the remaining exponents,
(3) proves

\[
 v_p\left(d^K{\mathbb E}\prod_j\sigma(j)^{k_j}\right)
 \ge K v_p(d)-\left\lfloor\frac K{p-1}\right\rfloor
 \ge0.                                                    \tag{11}
\]

Thus each normalized term is integral. Convergence needs a separate
observation: with `n` fixed, every moment has denominator dividing `n!`.
Its term therefore has valuation at least `K v_p(d)-v_p(n!)`, which
tends to infinity. There are only finitely many exponent tuples at
each `K`. Consequently (10) converges and can be reduced term by term.
For `d=0` the same result is immediate from (8).

If `p` is odd, every term with `K>0` in (11) is divisible by `p`.
The same is true when `v_p(d)>=2`. Only the constant term survives,
proving (9) without a restriction on the unit residues.

### The characteristic-two case with one residue

It remains to allow `p=2` and `v_2(d)=1`. Some nonconstant terms in
(11) can have valuation zero, so deleting them would be incorrect.
Choose one unit `u` congruent to every `y_j` modulo two. In each integral
term of (10), replacing the unit factor `prod_j y_j^(-k_j)` by `u^(-K)`
does not change its residue. The prefactor likewise becomes `u^(-n)`.
The resulting series is the normalized permanent of the constant-column
matrix with entries `(u+di)^(-1)`. That permanent is exactly

\[
 n!\prod_{i=0}^{n-1}(u+di)^{-1} .                         \tag{12}
\]

Its normalized residue is `u^(-n)`, equal to `prod_j y_j^(-1)`.
This proves (9). Multiplying by `n!` gives the congruence modulo
`p^(v_p(n!)+1)` for the original permanent.

Taking `y_j=T-dj` proves (1). The shifts all reduce to `T`, including
over an unramified extension whose residue field has more than two
elements.

## 3. Exact local content for every even-spacing packet mode

Use the physical kernel and local root factorization from
[the Boolean Cauchy note](boolean-cauchy-smith-profile.md):

\[
 W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},\qquad
 C(T)=\prod_{q=-L}^{L-1}(T-q),\qquad L=2^a.
\]

Let `S={c,c+d,...,c+(n-1)d}` be contained in the packet, with positive
even `d`. At an odd retained center `t`, write `T=t+2Lz`. The common
numerator has the form

\[
 \frac{C(t+2Lz)}{C'(0)}=2^{a+1}z\,u(z),\qquad
 u(z)\equiv1\pmod2.                                      \tag{13}
\]

The remaining principal matrix has entries `(T+d(i-j))^(-1)`.
All these denominators are units of `O[[z]]` and reduce to one.
The proof of (9) works coefficientwise in this ring: its terms are
integral and converge in the two-adic topology. Consequently, with
`v=v_2(n!)`,

\[
 \boxed{\quad
 \lambda_S(t+2Lz)=\mathrm{per}W_{S,S}(t+2Lz)
  \equiv 2^{n(a+1)+v}z^n
              \pmod {2^{n(a+1)+v+1}} .\quad}              \tag{14}
\]

The full local polynomial therefore has exactly the content in (2).
Every unit local evaluation point attains this valuation. For a
power-of-two size, `s_2(n)=1`, and (14) recovers the earlier formula.

Examples outside the earlier size restriction are:

| Packet length `L` | Spacing `d` | Size `n` | Exact local content |
|---:|---:|---:|---:|
| 8 | 2 | 3 | 13 |
| 16 | 2 | 5 | 28 |
| 16 | 2 | 6 | 34 |
| 16 | 2 | 7 | 39 |
| 32 | 2 | 10 | 68 |
| 64 | 2 | 21 | 165 |

As in the Walsh note, the two binary columns `(1+chi_S)/2` and
`(1-chi_S)/2`, paired with `delta_0`, have original contractions one
and zero. Their twirled contractions differ by `lambda_S`. Equation
(14) gives their exact local collision depth for these arbitrary
sizes. The half-integral normalization of the individual twirled
values remains part of that construction.

## 4. Every individual Smith order of the squared matrix

Let `p` be any prime, let `d` be nonzero with `v_p(d)=e>=1`, and let
`A` be a unit of an unramified complete discrete valuation ring `O`.
For every `N>=1`, the Hankel matrix

\[
 H_N(A,d)=\left((A+d(i+j))^{-2}\right)_{0\le i,j<N}
\]

has the exact Smith orders

\[
 \boxed{\quad
 s_r=2re+2v_p(r!)+v_p(r+1),\qquad 0\le r<N.
 \quad}                                                   \tag{15}
\]

They are strictly increasing because

\[
 s_r-s_{r-1}=2e+v_p(r)+v_p(r+1)>0\qquad(r\ge1).             \tag{16}
\]

Reversing the columns of the Toeplitz squared Cauchy matrix
`((T+d(i-j))^(-2))` gives `H_N(T-d(N-1),d)`. The new parameter is
still a unit, so (15) also describes that matrix. For `p=2,d=2`,
this is exactly the profile that the Walsh note had checked only at
finite sizes.

### Leading determinants and their shifts

Write `Delta_n(A)=det H_n(A,d)`, with `Delta_0=1`. The ordinary
Cauchy determinant formula and the permanent theorem in Section 2,
combined with [Borchardt's identity, Theorem 1.1 of Singer (2004)][borchardt],
give

\[
 v_p(\Delta_n(A))
   =e n(n-1)+2\sum_{j=0}^{n-1}v_p(j!)+v_p(n!).              \tag{17}
\]

Indeed, the ordinary Cauchy determinant has numerator
`d^(n(n-1)) product_(j=0)^(n-1) (j!)^2`, up to sign, and unit
denominators. Its permanent has order `v_p(n!)`.
The matrices in question are nonsingular because both factors in
Borchardt's identity are nonzero. Formula (17) is unchanged when `A`
is replaced by `A+d`.

Leading determinant valuations alone do not usually determine the
Smith form. The rest of the proof supplies the integral change of
basis that is needed here.

### An integral orthogonal basis

Over the fraction field, define a bilinear form on polynomials by

\[
 \langle X^i,X^j\rangle=(A+d(i+j))^{-2}.
\]

All leading Gram determinants are nonzero. Gram--Schmidt therefore
gives monic orthogonal polynomials `P_r` of degree `r`, with norms

\[
 h_r=\langle P_r,P_r\rangle
       =\frac{\Delta_{r+1}(A)}{\Delta_r(A)},
 \qquad v_p(h_r)=s_r.                                    \tag{18}
\]

The determinant formula for a monic orthogonal polynomial gives

\[
 P_r(0)=(-1)^r\frac{\Delta_r(A+d)}{\Delta_r(A)}.            \tag{19}
\]

For clarity, evaluating its last row `(1,X,...,X^r)` at zero
leaves the moment minor with indices `i+j+1`, which shifts `A` to
`A+d`. Equations (17) and (19) show that every `P_r(0)` is a unit.

Multiplication by `X` is self-adjoint for this Hankel form. The monic
polynomials consequently satisfy the three-term recurrence

\[
 XP_r=P_{r+1}+a_rP_r+b_rP_{r-1},\qquad
 b_r=\frac{h_r}{h_{r-1}}\quad(r\ge1),                      \tag{20}
\]

with the last term absent for `r=0`. To justify that no lower term
occurs, pair `XP_r` with `P_j` for `j<r-1`: moving `X` to `P_j`
leaves a polynomial of degree below `r`, orthogonal to `P_r`.
Pairing with `P_(r-1)` gives the displayed value of `b_r`.

By (16), each `b_r` is integral. Evaluating (20) at zero gives

\[
 a_r=-\frac{P_{r+1}(0)+b_rP_{r-1}(0)}{P_r(0)}\in O,        \tag{21}
\]

and the same conclusion holds for `a_0` after omitting the absent
term. The denominator and all the displayed constant terms are
units. Starting with `P_0=1`, recurrence (20) now proves by induction
that every coefficient of every `P_r` lies in `O`.

For degrees below `N`, the change from monomials to these monic
polynomials is unit lower triangular over `O`. It therefore gives
the exact unimodular congruence

\[
 U H_N(A,d)U^{\mathsf T}=\mathrm{diag}(h_0,\ldots,h_{N-1}). \tag{22}
\]

Each diagonal entry is a unit times `p^(s_r)`. Their orders are
increasing, so (22) proves the Smith profile (15).

This argument is algebraic; it does not require a real positive
inner product. It also explains why the shifted determinant valuations
are needed in addition to the unshifted ones.

## 5. Necessary hypotheses and unresolved questions

Arithmetic row coordinates matter. For the nonconsecutive positions
`S={0,2,4,8}` at `T=1`, direct rational evaluation gives

\[
 \mathrm{per}\left(\frac1{1+x-y}\right)_{x,y\in S}
   =\frac{63488}{165375},\qquad v_2=11,
\]

whereas `v_2(4!)=3`. Formula (1) does not hold for every set of even
positions.

The common-residue hypothesis in characteristic two also matters
over extension fields. In the unramified quadratic extension, take
`omega^2+omega+1=0`, `d=2`, and column shifts `omega,omega^2`.
The two-by-two permanent is exactly zero:

\[
 \frac1{\omega(\omega^2+2)}+
 \frac1{\omega^2(\omega+2)}=0.
\]

Its numerator is `2 omega omega^2+2(omega+omega^2)=0`.
The two shifts have different residues, so (9) does not apply.

These theorems do not turn the complete principal-permanent span into
the smaller interval span. The high-height scalar collisions and the
grouped rank obstruction in the Walsh note still apply. The Smith
profile concerns the physical squared Cauchy matrix, not the full
Boolean operator or its grouping cut.

## 6. Reproduction

```sh
python3 tools/check_arithmetic_cauchy_factorial.py \
  --output evidence/arithmetic-cauchy-factorial.json
```

The checker independently compares the equality-partition moment
formula with direct finite sums, tests the moment denominator bound,
evaluates permanents by subset dynamic programming over prime-power
rings, and checks complete local packet polynomials at non-power-of-two
sizes. Rational Gram--Schmidt checks the integral orthogonal basis;
separate minimum-valuation elimination checks the Smith orders. The
nonconsecutive-position and mixed-residue counterexamples above are
negative controls. The checks are finite controls of the ordinary
all-size proofs, not a new factoring machine certificate.

The Cauchy packet and Walsh constructions retain the attribution in
the linked source notes. The moment lemma, its permanent application,
the integral Smith argument, and the new checker used OpenAI assistant support.

[borchardt]: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v11i1r48/pdf
