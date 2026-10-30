# Partial-mode Weyl curves and an exact ordered endpoint construction

This construction separates adjacent mode exchanges without a multivariate
coefficient explosion.  Each partial curve has the full crossing as its
endpoint.  An ordered product using every adjacent edge once has a unique
intermediate occupation path for each entry, degree `O(h log h)`, and an
effective coefficient and factorization algorithm.  The finite profiles
below do not improve the existing factoring exponent.

## 1. Freeze the omitted modes

Work in characteristic two, with weights `w_j=3^j` and the full energy-`h`
label set `L_h`.  For a subset of modes `J`, define

\[
\begin{split}
P_J(z)(\alpha,y,u,\beta)
={}&\prod_{j\notin J}[\alpha_j=\beta_j][y_j=u_j]\,
 [\alpha_J+u_J=\beta_J+y_J]\\
&\cdot[t^{\alpha_J}]\prod_{j\in J}
(1+t_j)^{\beta_j}(1+(1+z)^{w_j}t_j)^{y_j}.
\end{split}                                                   \tag{1}
\]

At `z=1`, a nonzero entry requires `alpha_J<=beta_J`.  Frozen modes and
equal total energies give equal energies on `J`, so positive weights force
`alpha_J=beta_J`.  Conservation then gives `y_J=u_J`.  Therefore

\[
P_J(1)=E,\qquad \deg P_J\le h.                  \tag{2}
\]

No characteristic-power assumption is used.  Formula (1) can likewise be
written over any field, replacing `1+z` with `1-z`; the computations in this
note use characteristic two.

For `J={j,j+1}`, the normalized variable has degrees divisible by `3^j`.
After writing `q=z^{3^j}`, the two active weights become 1 and 3.  Thus the
commutator projector in `mismatched-commutator-projector.md` applies to
coefficients `3^j` and `3^{j+1}` of `P_J(1+z)`.  The two labels can have
different energies on these two active modes; the occupation proof permits
that difference.  This gives a clean second exchange on modes `(3,9)`.

## 2. A complete endpoint with only polynomially many coefficients

There are `r=floor(log_3 h)+1` modes.  Choose any permutation
`pi=(pi_1,...,pi_(r-1))` of the adjacent edges.  Here the listed order is
the application order, so define the row-matrix product

\[
F_\pi(z)=P_{\{\pi_{r-1},\pi_{r-1}+1\}}(z)
\cdots P_{\{\pi_1,\pi_1+1\}}(z).
\]

Equations (2) imply

\[
F_\pi(1)=E,\qquad \deg F_\pi\le(r-1)h.          \tag{3}
\]

Consequently its coefficients are an exact crossing decomposition with
at most `(r-1)h+1` pieces.  Using the same scalar variable in all factors
is essential to this count; no independent multivariate type grid is being
silently omitted.

For fixed radix, `log |L_h|=Theta((log h)^2)`.  Only the elementary two-sided
order estimate is needed here: coordinatewise occupation bounds give the
upper bound.  For the lower bound, independently vary the first half of
the modes other than the weight-one mode below `h^(1/2)/O(log h)`, and use
the weight-one mode only to fill the remaining energy.  This gives an
injective family of labels of the required size.
Thus the logarithm of the piece count is `O(log h)=o(log |L_h|)`.

Enumerating labels, evaluating all polynomial entries, and factoring the
two flattenings by exact elimination costs `|L_h|^O(1) h^O(1)`, or
`2^O(lambda)` for `lambda=ceil(log_2 |L_h|)`.  This is an effective interface;
it asserts no favorable rank profile without computing or proving one.

## 3. Unique paths from independent occupation flows

Every partial operator preserves the integer difference of the two row
labels.  A potentially nonzero entry therefore satisfies

\[
\alpha-y=\beta-u.                               \tag{4}
\]

An operation on edge `j` changes the first label by a multiple of

\[
v_j=3e_j-e_{j+1},
\]

and changes the second label by the same multiple.  If every edge is used
once, all flow amounts are fixed by the endpoints.  Set `q_-1=0` and

\[
q_j=\frac{\alpha_j-\beta_j+q_{j-1}}3,
\qquad 0\le j<r-1.                             \tag{5}
\]

Equal energies make these quantities integers: equivalently,
`q_j=(sum_(k<=j) 3^k(alpha_k-beta_k))/3^(j+1)`.
The same equality gives the final-coordinate equation
`alpha_(r-1)-beta_(r-1)=-q_(r-2)`.

For any application order, start at `(beta,u)` and add `q_j v_j` to both
labels when edge `j` is applied.  These are the only possible intermediate
states.  If any occupation becomes negative, the entry is zero.  Otherwise
the entry is the product of the two-mode polynomial entries along this
single path.  There is no sum over hidden intermediate labels.

This is stronger than endpoint verification alone: it gives a concrete
formula for every coefficient of a noncommutative product.  The edge
vectors are independent, and that independence is the reason uniqueness
survives arbitrary application order.

## 4. Two adjacent isolated exchanges have a three-state label description

The projector construction gives an isolated label exchange

\[
\tau_j:\ (a_j,a_{j+1})\leftrightarrow(a_j+3,a_{j+1}-1)
\]

between `a_j=0 mod4,a_(j+1) odd` and
`a_j=3 mod4,a_(j+1) even`.  No higher binary digit changes: this is precisely
the three-bit exchange `001 <-> 110` on the low two bits of `a_j` and the
low bit of `a_(j+1)`.

The union of the partial matchings `tau_0,tau_1` has components of size at
most three.  A vertex where both moves are active has residues

\[
(a_0\bmod4,a_1\bmod4,a_2\bmod2)
\in\{(0,3,0),(3,0,1)\}.
\]

Either move takes it to a vertex where the other move is inactive.
Therefore every degree-two vertex is the middle of a three-vertex path;
all other components are isolated vertices or edges.  On each three-path,
the two partial adjacency maps satisfy

\[
K_0K_1K_0=K_1K_0K_1=0,\qquad K_i^3=K_i.       \tag{6}
\]

This is a common three-state label description, not yet an auxiliary
one-dimensional crossing quotient.  In fact the two adjacency maps on a
three-path generate its full matrix algebra, so a one-dimensional character
cannot simply be assumed. On a path numbered `0,1,2`, the product `K_0^2 K_1^2` is the middle matrix unit. Multiplying it on either side by the two adjacency maps gives the adjacent matrix units, and their products give all remaining units.

At `h=63`, the two isolated exchanges have 69 singleton components, 18
edges, and nine three-paths.  The full partial curves act on additional
occupation states beyond these isolated exchange sectors.

## 5. Exact finite controls

The native engine `explore_partial_mode_sweep.cpp` implements (5) directly
and checks every supported entry at `z=1`.  The independent standard-library
checker builds each partial polynomial matrix separately and multiplies
the matrices with actual intermediate-state summation, then compares both
coefficient-rank arrays with the unique-path engine.

```sh
python tools/check_partial_mode_sweep.py
```

The check covers energies 9 and 27, including two different three-edge
orders at energy 27.  For application order `(0,1,2)` at energy 27, the full
profile has 77 nonzero coefficients and rank sum about 6263.451, above the
target 84.290.  The constant coefficient alone has cuts `(230,142)`.

At energy 63, a constant-coefficient screen gives:

| Application order | Constant row rank | Constant column rank | Its rank-sum contribution |
|---|---:|---:|---:|
| `(0,1,2)` | 5195 | 2400 | 3531.006 |
| `(0,2,1)` | 4879 | 2346 | 3383.214 |

Both already exceed the full target `132^sqrt(2)=997.572`, so a larger raw
profile is unnecessary for these two orders.  This does not exclude
coefficient recombinations, alternative partial curves, repeated edges,
or constructions using the shared small exchange sectors.  Those remain
the useful places to seek a better complete rank profile.

### Constant powers before a larger repeated-edge profile

Squaring the entire endpoint curve preserves its value `E` at 1 and allows
genuine interference among intermediate states.  As a targeted control,
let `A=F_(0,1,2)(0)` at energy 63.  Direct multiplication of the row blocks
and exact reshaping into the other cut give:

| Power | Row rank | Column rank | Product |
|---:|---:|---:|---:|
| 1 | 5195 | 2400 | 12468000 |
| 2 | 2568 | 2771 | 7115928 |
| 4 | 1844 | 2646 | 4879224 |
| 8 | 1782 | 3375 | 6014250 |
| 16 | 1782 | 3618 | 6447276 |

The row rank stabilizes between powers 8 and 16.  Every displayed product
is larger than `132^3`, hence also larger than `132^(2sqrt(2))`.
The constant coefficient alone therefore excludes the raw coefficient
decompositions of these five powers of the curve.  A full repeated-edge
profile was unnecessary.

This is a finite screen of the displayed powers.  Stabilization of the row
rank does not imply stabilization of the other cut or exclude every later
power or polynomial recombination.

The optional `explore_partial_sweep_powers.py` reconstructs the constant
partial matrices independently, multiplies and squares them, and reports
both ranks until the row rank stabilizes.  Its first pair of ranks matches
the native unique-path engine.  The energy-63 replay took about seven seconds.
