# Sharpness and limits of the Clifford compression argument

This note gives two exact controls for
[the Clifford compression theorem](clifford-compression-pressure-floor.md).
The first shows that its generic row-frame threshold is sharp.
The second gives a concrete cubic obstruction to extending the quadratic
proof by the same central-square argument. Neither changes the factoring
exponent.

Work over a field of characteristic two. On `t` binary sites, write
`n=2^t`, let `e_i=E01` and `f_i=E10` act at site `i`, and put

$$J_t=\begin{pmatrix}0&1\\1&0\end{pmatrix}^{\otimes t}.$$

The degree filtration `F_d` is the span of the ordered squarefree
Clifford monomials of degree at most `d`.

## 1. The row threshold is attained

Fix `1<=d<t` and set

$$A=e_1\cdots e_d.$$

Then

$$A J_t=E_{00}^{\otimes d}\otimes J_{t-d}. \tag{1}$$

This is an alternating bilinear form of rank `n/2^d`. Its radical is spanned
by the coordinate vectors whose first `d` bits are not all zero.
Inside the remaining coordinate subspace, keep only the vectors whose
next bit is zero. Together these vectors span a subspace `W` of dimension

$$q=n-\frac{n}{2^{d+1}}. \tag{2}$$

Equation (1) pairs a vector in the remaining block only with the vector
whose last `t-d` bits are complemented. Two kept vectors cannot be paired,
since their next bit is zero. Thus `W` is totally isotropic for (1).

Let `P` be its coordinate embedding and set `L=P^T`. These maps satisfy
`LP=I`, but

$$P^T A J_t P=L A J_t L^T=0. \tag{3}$$

Consequently the degree-`d` row frame is not injective at

$$2^{d+1}(n-q)=n.$$

The strict inequality in the generic row theorem cannot be weakened to
a non-strict inequality. In particular it cannot be replaced by the
factor-two stronger column threshold suggested by the quadratic case.

There is also an actual Smith-multiplicity loss. The lower-degree frames
are independent, since their strict generic threshold still holds.
After making those columns exact coordinate axes by regular basis
changes, the order-`d` coefficient factors through the images of the
degree-`d` constant frame vectors modulo the lower-degree span. Equation
(3) gives a nonzero kernel vector in that quotient. Therefore the row
Smith multiplicity at order `d` is strictly below `binom(2t,d)`.

## 2. A cubic with commuting square-zero first commutators

On seven Clifford pairs, define

$$\begin{aligned}
A={}&f_1f_2e_5+f_3f_4f_5
     +f_1f_3e_6+f_2f_4f_6\\
 &+f_1f_4e_7+f_2f_3f_7.
\end{aligned} \tag{4}$$

This operator belongs to `F_3`. Each summand squares to zero. In the
square of their sum, the only nonzero cross commutators come from the
three matched auxiliary pairs `(e_5,f_5)`, `(e_6,f_6)`, `(e_7,f_7)`.
Each contributes `f_1f_2f_3f_4`. Since three is one in the ground field,

$$A^2=f_1f_2f_3f_4\ne0. \tag{5}$$

For every first-layer vector `v`, the first commutator satisfies

$$[A,v]^2=0. \tag{6}$$

Here commutators are written with a plus sign in characteristic two.
To check (6), it suffices first to use the fourteen generators. The
commutators with `f_1,...,f_4` vanish. Those with the auxiliary generators
are products of two of `f_1,...,f_4`. For example,

$$[A,e_1]=f_2e_5+f_3e_6+f_4e_7,$$

whose summands commute and square to zero. The analogous formulas for
`e_2,e_3,e_4` have the same property.

Moreover all fourteen first commutators commute pairwise. The only
possible nonzero terms between two of the latter four formulas occur
twice and cancel; the other cases have no noncommuting auxiliary pair.
Linearity then proves both statements for arbitrary first-layer vectors.

Nevertheless, (5) is not central:

$$[A^2,e_1]=f_2f_3f_4\ne0. \tag{7}$$

Thus square-zero, pairwise commuting first commutators do not force a
central square in degree three. The extra Jacobi term identified in the
main note can really survive. This closes that particular proposed
induction; it does not disprove a stronger compression theorem proved
by some different argument.

The exact matrix control gives `rank(A)=38` and `rank(A^2)=8` in dimension
128. In particular this example does not contradict the still-unproved
rank assertion that `A^2!=0` might imply `rank(A)>=n/2^(d-1)`.

## Exact checks

[`tools/check_clifford_threshold_sharpness.py`](../tools/check_clifford_threshold_sharpness.py)
constructs the coordinate retractions and the cubic matrix directly.
It verifies the frame kernels, selected row Smith multiplicity losses,
the cubic square, all first-commutator squares, every pairwise
commutator, and the noncentral square. Its report is
[`evidence/clifford-threshold-sharpness.json`](../evidence/clifford-threshold-sharpness.json).
These computations use exact arithmetic over `F_2`; the displayed
algebraic arguments apply over every characteristic-two field.

This note used substantial OpenAI assistant support.
