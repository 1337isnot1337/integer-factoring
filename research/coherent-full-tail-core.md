# A fixed physical code with a full Gaussian tail

## Result and scope

A physical code on the lowest `k` reflection modes, tensored with the
complete occupation alphabet on all remaining modes, admits a bounded
active-core calculation. Here `k` and the total number `M` of physical
modes are powers of two, with `M >= 2k`. For the compound matrix, an
arbitrary constant projector on the `k` modes uses only `4k` active
Gaussian modes after regular changes of basis. The other modes retain
their diagonal Gaussian factors.

To obtain the actual two profiles of a physical reflection restriction
from this compound calculation, the physical embedding and retraction
must also satisfy the height homogeneity and deletion invariance
conditions in [the coherent compound reduction](coherent-hot-compound-retractions.md).
These hypotheses allow the internal Gaussian factors to descend through
the physical maps. They are not consequences of the projector identity.
In particular, the result does not yet give the profiles of arbitrary
unequal-height projectors on the four lowest occupation states.

The complete implementation treats the coordinate code that forbids
joint occupation of the two lowest physical modes. Its physical
alphabet has `3*2^(M-2)` letters. Both profiles are obtained from a
rank-144 active matrix of order 256 and explicit Gaussian spectators.
Exact calculations through `M=512` do not improve the Gaussian pressure.
This is a uniform reduction and a finite computational result, with no
improvement to the factoring exponent.

## 1. The polynomial tensor and its compound matrix

Throughout this note, `M` is the number of physical reflection modes.
The original reflection interval has length `2M`, unlike the interval
notation in [the reflection construction](reflection-dual-rail-compression.md).
The deflated variable is `t=w^2`. The occupation state `A` has height

\[
h(A)=\sum_{a\in A}(2a+1).
\tag{1}
\]

Take a constant split pair `P,L` on the first `k` modes, with `LP=I`,
and tensor both maps with the identity on the other `M-k` modes. Write
`q_hot` for the rank of `P`. The resulting physical alphabet has

\[
q=q_{\rm hot}2^{M-k}.
\tag{2}
\]

Under the stated height and deletion hypotheses, the two natural
regular profiles are the two compressed compound profiles. The
one-particle matrix is

\[
B_M(t)=\begin{pmatrix}D&J\\J&D\end{pmatrix},\qquad
D_{ab}=(1+t)t^{2(a\mathbin{|}b)+1},\qquad
J_{ab}=[a=b]t^{2a+1},
\quad 0\le a,b<M.
\tag{3}
\]

The infinity cut has the map orientations in equation (5) of the
coherent compound note. These orientations cannot in general be
replaced by the zero-cut maps. The implemented coordinate code has
`P=L^T`, so its two profiles agree.

The raw reflection tensor remains polynomial. Its conservative degree
bound is `8M^2` in `w`, or `4M^2` in `t`: it is a restriction of selected
minors of the original polynomial Pascal matrix. A constant physical
restriction cannot increase this degree. Its common scalar contributes
`+M^2` to the zero orders in `t` and `-M^2` to the infinity orders. Those
scalar shifts cancel in the product of the two counting functions.
None of these physical degree or height costs is removed by the
bounded matrix calculation below.

## 2. A regular diagonalization with explicit operators

Swap the two column blocks of (3), and define

\[
F_{ab}=(1+t)t^{2(a\mathbin{\&}\neg b)},\qquad G=F^{\mathsf T}.
\]

Then

\[
B_{M,\mathrm{swapped}}=
\begin{pmatrix}I&0\\F&I\end{pmatrix}
\begin{pmatrix}J&0\\0&t^{2M}J\end{pmatrix}
\begin{pmatrix}I&G\\0&I\end{pmatrix}
=L\Delta R.
\tag{4}
\]

All triangular factors and their inverses are regular. To check the
Schur complement, first note `FJ=JG`. On the binary digits of a mode,

\[
G=(1+t)\bigotimes_{j=0}^{\log_2 M-1}
 \begin{pmatrix}1&t^{2^{j+1}}\\1&1\end{pmatrix}.
\]

Each digit factor squares to a scalar. Hence

\[
G^2=(1+t)^2\prod_j(1+t^{2^{j+1}})I
    =(1+t^{2M})I,
\]

which proves (4). The Gaussian diagonal has orders

\[
1,3,\ldots,2M-1,\qquad 2M+1,2M+3,\ldots,4M-1.
\tag{5}
\]

For an invertible one-particle map `T`, let `Gamma(T)` be its action
on the occupation algebra. A creator with coefficient column `v`
is sent to the creator with coefficients `Tv`; an annihilator with
coefficient row `a` is sent to the annihilator with coefficients
`aT^-1`. These identities hold in characteristic two without sign
choices. They transport every matrix unit, and therefore every
constant physical projector.

There is also a useful elementary sandwich identity. If `LP=I` and
`Pi=PL`, the nonzero Smith factors of `Pi A Pi` are precisely those
of `LAP`. Indeed, extend the image and kernel of the constant
projector to a constant basis; the sandwich is then the compressed
matrix with zero rows and columns added. The same argument permits
different constant row and column split pairs. This identity is the
reason to work with transported projectors on the full occupation
space rather than construct the exponentially large compressed
matrix.

## 3. Weighted regular concentration

Write a mode as `c+kB`, where `0 <= c < k` and `0 <= B < N=M/k`.
Put

\[
x=t^{2k},\qquad S'=\sum_{j=0}^{N-2}x^j,\qquad S=xS'.
\tag{6}
\]

For each residue `c`, keep the initial coordinate `B=0` fixed.
On its tail `B=1,...,N-1`, use a matrix `C` whose first column is

\[
v=\frac{(1,x,\ldots,x^{N-2})^{\mathsf T}}{S'},
\]

and whose column `B>=2` is `e_B+e_1`. Its inverse `U=C^-1`
has first row `(1,...,1)` and row `B>=2`

\[
e_B^{\mathsf T}+\frac{x^{B-1}}{S'}(1,\ldots,1).
\tag{7}
\]

The denominator is a local unit, since `S'(0)=1`. Direct
multiplication proves that these are inverse matrices. They
concentrate the two relevant tail vectors:

\[
(1,\ldots,1)C=e_1^{\mathsf T},\qquad
U(x,x^2,\ldots,x^{N-1})^{\mathsf T}=S e_1.
\tag{8}
\]

The weighted regularity is essential. For each of `C,U`, an entry
from tail column `j` to row `i` has valuation at least
`max(0,2k(i-j))`. This follows immediately from the displayed
columns and rows. Thus, after restoring the physical mode indices,
both `J^-1 C J` and `J^-1 U J` are regular. Their inverses are also
regular, by the same calculation.

Apply `U` on both mode groups on the left of `Delta`, and pair it
with `V=Delta^-1 C Delta` on the right. These changes preserve
`U Delta V=Delta` exactly. Their occupation actions are regular
unimodular matrices. They change Smith profiles by neither a
valuation shift nor an uncharged physical height.

The concentration is regular at the crossing endpoint as well:
`N` is an even power of two and `S'(1)=N-1=1` in characteristic
two. This endpoint check is useful for the generic-rank promise.
The concentration matrices are local calculation devices; they
are not additional rational maps imposed on the physical tensor.

## 4. The complete `4k` active core

The active Gaussian modes are

\[
x(c,0),x(c,1),y(c,0),y(c,1),\qquad 0\le c<k.
\tag{9}
\]

They have orders

\[
\{2j+1:0\le j<2k\}\ \cup\
\{2M+2j+1:0\le j<2k\}.
\tag{10}
\]

The spectators have the same two sets of orders for `2k <= j < M`.
To specify every transported physical operator, write `c_z,a_z`
for creation and annihilation on mode `z`, and put
`f_ca=(1+t)t^(2(c & ~a))`. For `0 <= a < k`, the row operators are

\[
\begin{aligned}
c_{x_a}&\longmapsto c_{x(a,0)}+
 \sum_{c<k}f_{ca}\bigl(c_{y(c,0)}+S c_{y(c,1)}\bigr),\\
a_{x_a}&\longmapsto a_{x(a,0)},\\
c_{y_a}&\longmapsto c_{y(a,0)},\\
a_{y_a}&\longmapsto a_{y(a,0)}+
 \sum_{c<k}f_{ac}\bigl(a_{x(c,0)}+a_{x(c,1)}\bigr).
\end{aligned}
\tag{11}
\]

After the column swap in (4), the column operators are

\[
\begin{aligned}
c_{v_a}&\longmapsto c_{x(a,0)},\\
a_{v_a}&\longmapsto a_{x(a,0)}+
 \sum_{c<k}f_{ca}\bigl(a_{y(c,0)}+x a_{y(c,1)}\bigr),\\
c_{u_a}&\longmapsto c_{y(a,0)}+
 \sum_{c<k}f_{ac}\bigl(c_{x(c,0)}+S'c_{x(c,1)}\bigr),\\
a_{u_a}&\longmapsto a_{y(a,0)}.
\end{aligned}
\tag{12}
\]

These are complete polynomial formulas. In particular, the local
unit denominators in the concentration basis cancel from these
operators. Equations (11)-(12) follow by conjugating first through
`L,R` in (4), then through `U,V`; (8) removes all tail directions
except the first one in each residue. The independent checker
performs those full `2M`-mode conjugations before comparing with
the displayed formulas.

Any operator on the first `k` physical modes is a polynomial in
their creators and annihilators. Substitution of (11)-(12)
therefore places each physical projector entirely on the active
space. The result has the form

\[
K_M\otimes\Gamma(\Delta_{\mathrm{spectator}}),
\tag{13}
\]

up to regular row and column units. Its complete nonzero Smith
orders are the orders of `K_M` plus every subset sum of the
spectator orders. The active matrix has order `2^(4k)` and generic
rank

\[
q_{\rm hot}^{\,2}2^{2k}.
\tag{14}
\]

For the split pairs under consideration, the endpoint is a
projector of this rank: `f_ca(1)=0`, the active diagonal is the
identity, and each initial mode is fixed by the concentration.
The projector rank also supplies the generic upper bound, so
this endpoint supplies the exact generic rank independently of
any Smith precision window. The spectator multiplicity gives
`q_hot^2 2^(2k) 2^(2M-4k)=q^2`, as required.

Equations (11)-(13), as identities for the compound, do not require
height or deletion invariance. Those additional hypotheses are
needed at the earlier passage from the physical reflection tensor
to (3). Removing them requires transporting the internal Gaussian
factors too; the present implementation does not do that.

## 5. The forbidden-pair code

Take `k=2` and retain the occupation states `{00,01,10}`. The
physical projector is

\[
\Pi=I+n_0n_1,\qquad n_a=c_a a_a.
\tag{15}
\]

This is a height-homogeneous coordinate downset. The same projector
is used on all four physical legs, and both cuts have the same
complete profile. Substituting (11)-(12) into the four copies of
(15), with the correct right-action order, gives the active matrix
implemented in [the constructor](../tools/coherent_full_tail_core.py).

The ambient order is 256 and the generic rank is 144. Occupation
is preserved. The ranks of its occupation blocks are the
coefficients of

\[
(1+2z)^2(1+z)^4,
\]

namely `1,8,26,44,41,20,4,0,0`. The largest ambient block has
order 70. The endpoint rank is checked in every block before
Smith elimination. No missing block or unresolved pivot is
included in a claimed complete profile.

A conservative complete coefficient-degree bound for this
active matrix is

\[
\deg_t K_M\le16M+40.
\tag{16}
\]

The sum of all eight active diagonal orders is `8M+32`. The two
creators in the first row projector have degree at most `2M-1`
each; the two annihilators in the second row projector have
degree at most 3 each. The corresponding column bounds are 7
and `2M-5`. Summing these four projector bounds adds at most
`8M+8`, proving (16). The implementation constructs all
coefficients, starts above their actual degree, and doubles the
series window if the promised number of pivots is unresolved.
Its explicit precision limit is 8192. A failed limit is an error,
not a partial profile.

For fixed `k`, the active order is constant and the displayed
coefficient degrees grow linearly in `M`. Thus exact polynomial
linear algebra gives an effective calculation polynomial in `M`
for these fixed-size cores, with constants exponential in `k`.
This profile-calculation bound does not replace the physical
alphabet size or the original tensor degree.

## 6. Physical entropy and the numerical comparison

Let `a_1,...,a_144` be the active orders. Define

\[
\begin{aligned}
Z_M(r)&=\sum_{i=1}^{144}r^{a_i},\\
S_M(r)&=\prod_{j=4}^{M-1}(1+r^{2j+1})(1+r^{2M+2j+1}),\\
A_M(t)&=(1+t+t^3)\prod_{j=2}^{M-1}(1+t^{2j+1}).
\end{aligned}
\tag{17}
\]

The complete profile count at either cut is `Z_M(r)S_M(r)`.
The complete physical height partition is `A_M(t)`, including
the full tail. The common-type pressure test is therefore

\[
2\log\bigl(Z_M(r)S_M(r)\bigr)
 <4\rho\log A_M\bigl(r^{1/\rho}\bigr),\qquad 0<r<1.
\tag{18}
\]

For comparison, the unmodified full Gaussian has profile
partition `prod_(j=0)^(2M-1) (1+r^(2j+1))` and physical
partition `prod_(j=0)^(M-1) (1+t^(2j+1))`. These are evaluated
at the same `M`; no limit or different alphabet is substituted
for the benchmark.

The following values are numerical local minima, with `r=exp(-u)`.
The search uses 101 grid points in `log(u)` for `10^-6 <= u <= 2`,
then refines each detected local minimum. It supplies neither a
global scale certificate nor a strict rate certificate.

| Physical modes `M` | Forbidden-pair pressure | Its `u` | Full Gaussian pressure |
| ---: | ---: | ---: | ---: |
| 4 | 0.7376486592 | 0.686516859 | 0.7166536597 |
| 8 | 0.7376471826 | 0.685065341 | 0.7106044045 |
| 16 | 0.7370708098 | 0.102316425 | 0.7083157559 |
| 32 | 0.7262448424 | 0.055042572 | 0.7075072113 |
| 64 | 0.7185222280 | 0.030503761 | 0.7072351452 |
| 128 | 0.7136531695 | 0.016956790 | 0.7071468790 |
| 256 | 0.7107699534 | 0.009384127 | 0.7071190462 |
| 512 | 0.7091235806 | 0.005159039 | 0.7071104679 |

The `M=128` value is not a fixed-code limit. There are two
detected local-minimum branches from `M=16` onward, and the
better branch drifts toward smaller `u`. The finite calculations
do not prove its limiting value. They show no pressure gain
over the unmodified Gaussian.

## 7. Independent controls and retained evidence

The [checker](../tools/check_coherent_full_tail_core.py) verifies
the concentration inverses, weighted regularity, and all
transported operators for several pairs `(k,M)`. It separately
compares (4) with the original one-particle matrix. For `M=4`,
it constructs the literal normalized Pfaffian tensor, selects
the twelve physical words, and uses the independent Python
Smith implementation on both parity blocks. This recovers all
144 active orders, matching both the previous compound
constructor and the new core.

The retained [controls](../evidence/coherent-full-tail-controls.json)
record 5,568 inverse-matrix entries, 5,568 weighted-regular
entries, 5,696 transported-operator entries, 1,808 canonical
pairings, and 1,008 one-particle factorization entries. The
literal control checks all 20,736 endpoint entries. The complete
receipt replay recomputes all 1,152 active orders for the eight
intervals in the [profile evidence](../evidence/coherent-full-tail-core.json),
including block ranks, complete coefficient degrees, precision
attempts, and source hashes. Numerical pressure values are
explicitly separate from these exact algebraic checks.

```sh
python3 tools/coherent_full_tail_core.py \
  --output evidence/coherent-full-tail-core.json
python3 tools/check_coherent_full_tail_core.py \
  --evidence evidence/coherent-full-tail-core.json \
  --output evidence/coherent-full-tail-controls.json
```

The next construction question is whether physical maps that
mix the low states can improve the pressure once their internal
Gaussian factors and height costs are retained. The `4k` compound
formula alone does not settle that question.
