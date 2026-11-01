# A Gaussian cube profile floor for every equally weighted Fock alphabet

The [one-mode endpoint-minor argument](cube-gaussian-profile-budget.md)
extends to any number of equally weighted modes on each physical leg.
Arbitrary interactions between those modes do not remove this obstruction.
The common-scale, unweighted Smith-profile method still has `rho > 3/5`,
which is above the cube half-score target `1/sqrt(3)`.

This statement uses the **full Fock alphabet** on each leg. Proper physical
subspaces, unequal mode heights, different decompositions, and stronger
coefficient-rank estimates are outside its scope. It is a limitation of
the stated profile estimate, not a lower bound on the ranks of arbitrary
coefficient pieces.

## 1. Setting

Let `M >= 1`, and place `M` binary modes at every vertex of a cube. Let
`K(x,y,z)` be an alternating `8M` by `8M` matrix over a field of
characteristic two. Every nonzero entry is homogeneous of the same degree
`D > 0`. At `(1,1,1)`, require

\[
K_{(v,i),(w,j)}=[w=\overline v][i=j].
\tag{1}
\]

Thus its Gaussian tensor has antipodal equality as its endpoint. The
physical alphabet consists of every subset of the `M` modes. A letter
with `r` occupied modes has energy `Dr/2`; the tensor entry degree is the
sum of the eight letter energies.

At cube axis `a`, write the `4M` one-particle crossblock Smith orders as
`s_(a,1) <= ... <= s_(a,4M)`, and set

\[
\gamma_a(k)=\sum_{i=1}^k s_{a,i}.
\]

Orders are taken at the relevant axis variable over the rational
function field in the other two variables. Polynomial within-cut
Gaussian factors are unimodular on the full Fock spaces. Consequently
the complete tensor cut profile consists of the subset sums of these
one-particle orders, with counting function

\[
Z_a(t)=\prod_{i=1}^{4M}(1+t^{s_{a,i}}).
\tag{2}
\]

No translation symmetry, Cayley form, or commutation of color matrices
is assumed.

## 2. Shared minors of sizes M and 2M

Take all `M` colors of the antipodal vertex pair `000,111`. Their square
submatrix has endpoint identity and homogeneous determinant degree `MD`.
The same minor lies in every axis crossblock, up to transposition.
Therefore

\[
\sum_a\gamma_a(M)\le MD.
\tag{3}
\]

Next fix two axes `a,b`, with remaining axis `c`. Select all colors of
the two vertices `000,e_c` as rows and their antipodes as columns. This
same `2M` by `2M` minor lies in both cuts `a,b`. Its endpoint is identity
and its determinant degree is `2MD`. Hence

\[
\gamma_a(2M)+\gamma_b(2M)\le 2MD.
\]

Adding over the three axis pairs gives

\[
\sum_a\gamma_a(2M)\le3MD.
\tag{4}
\]

Each complete crossblock determinant has degree `4MD`. Thus

\[
\sum_a\gamma_a(4M)\le12MD.
\tag{5}
\]

These bounds concern complete generic minors. Specializing the other
two variables to one can increase their valuations and does not suffice
to certify these generic profiles.

## 3. The same four-order relaxation, repeated M times

For `0<t<1`, the function `f(s)=log(1+t^s)` is decreasing and convex.
Group the orders at all three cuts into the first `M` per cut, the next
`M` per cut, and the remaining `2M` per cut. Let their respective means
be `u,v,w`. Sorted orders give `u <= v <= w`, while (3)--(5) give

\[
u\le D/3,\qquad u+v\le D,\qquad u+v+2w\le4D.
\]

The convex balancing argument from the one-mode proof now gives

\[
\frac13\sum_a\log Z_a(t)
\ge M\left[
 \log(1+t^{D/3})+\log(1+t^{2D/3})
 +2\log(1+t^{3D/2})\right].
\tag{6}
\]

In other words, the most optimistic relaxed profile consists of `M`
copies of `[D/3,2D/3,3D/2,3D/2]`. It need not be realizable by any kernel.
Put `r=t^(D/6)` and

\[
Z_*(r)=(1+r^2)(1+r^4)(1+r^9)^2.
\tag{7}
\]

Then the geometric mean of the three complete cut counting functions
is at least `Z_*(r)^M`.

## 4. Every common physical type is covered

Let `P` be any distribution on the full `2^M`-letter alphabet, and let
`Mp` be its expected occupation. Entropy subadditivity and concavity of
binary entropy give

\[
H(P)\le\sum_{i=1}^M H_2(p_i)\le M H_2(p).
\tag{8}
\]

Using this same physical type on all eight legs fixes total coefficient
degree `4DMpm` over `m` sites. The common profile method therefore has
contraction at least

\[
\frac{M\log_2 Z_*(r)-8Mp\log_2 r}{4H(P)}.
\tag{9}
\]

The exact positive-polynomial certificate in the one-mode note proves

\[
\log_2 Z_*(r)-8p\log_2r>\frac{12}{5}H_2(p)
\quad(0<p<1,\ 0<r\le1).
\]

Equations (8) and (9) imply `rho > 3/5` for every common type of positive
entropy. Types of zero entropy cannot provide the large alphabet needed
by the transfer theorem. This proves the claimed floor for every `M`.

## 5. Exact controls

The [checker](../tools/check_cube_equal_height_fock_floor.py) reconstructs
the shared generic minors of seeded kernels with one through four modes.
It checks all entry endpoints, homogeneous degrees, the common row and
column sets for both cuts, and all determinant valuation budgets. It
also replays the 355-coefficient positive-polynomial certificate used
above. The [receipt](../evidence/cube-equal-height-fock-floor.json) retains
the inputs, counts, and source hashes.

Run `python3 tools/check_cube_equal_height_fock_floor.py` to reproduce
these controls. The arbitrary-size conclusion follows from (3)--(9),
not from exhausting the sampled kernels.
