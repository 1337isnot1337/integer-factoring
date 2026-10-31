# Two-cut lattice synthesis for occupation-transfer kernels

This note gives an exact linear search method for the
[occupation-preserving operator family](coupled-identity-operator-filters.md).
It allows every operator coefficient in the selected occupation blocks,
including off-diagonal transfers. Its purpose is to search prescribed
zero and infinity Smith improvements without enumerating polynomial
coefficients individually. It does not assert that an improved factoring
bound has been found.

The implementation currently supports two or three binary identity
directions, with coefficients in \(\mathbf F_4\). It accepts every proved
geometry handled by the [Weyl reduction](identity-operator-weyl-reduction.md),
including [split-end selections](split-identity-mode-reduction.md).

## 1. Independent row lattices at the two cuts

Let \(S=(j_1,\ldots,j_r)\) be the selected directions, let \(P\) be the
ramification, and write \(q=-v\). Label the binary occupation states by
\(a,b\in\{0,1\}^r\). Only coefficients \(F_{ab}(q)\) with
\(|a|=|b|\) are present. Put

\[
 \kappa_b=P\sum_s(2j_s+1-2n)b_s.
 \tag{1}
\]

At infinity the active matrix is

\[
 F(z^{-2})\operatorname{diag}(z^{\kappa_b}).
 \tag{2}
\]

Choose any constant invertible binary matrix \(W\) that preserves
total occupation, and use \(G=WF\) as the unknown polynomial matrix.
This row change is a unit at infinity. A desired row lower bound
\(\eta_a\) in (2) is exactly the coefficient degree condition

\[
 \deg_q G_{ab}\le
 \left\lfloor\frac{\kappa_b-\eta_a}{2}\right\rfloor.
 \tag{3}
\]

Thus, after fixing a common lower Laurent power \(L\), every unknown
coefficient lies in the finite interval between \(L\) and the bound
in (3). An optional common upper power further restricts that interval.

Independently, let \(M(F)\) be the complete polynomial Weyl matrix at
zero, using the zero-unit denominator clearing in the cited reduction.
The variable there is \(q=z^2\). For each occupation block, compute
regular polynomial matrices \(U_0,V_0\) such that

\[
 U_0 M(I)V_0
\]

is diagonal. The computation uses only row and column operations with
unit determinants at zero; it does not divide by a truncated series.
Fix desired row lower bounds \(\theta_i\) for

\[
 U_0 M(W^{-1}G)V_0.
 \tag{4}
\]

The map \(F\mapsto M(F)\) is linear. Consequently, requiring every
coefficient below \(q^{\theta_i}\) in row \(i\) of (4) to vanish is
a finite homogeneous linear system in the Laurent coefficients of \(G\).
The constant infinity basis and the zero Smith basis are independent.
In particular, no zero-unit denominator or basis change is presumed
regular at infinity.

## 2. Endpoint conditions and exact \(\mathbf F_4\) feasibility

The endpoint repair needs \(E=F(1)\) to preserve the physical \(Y\)
energy filtration. The implementation imposes

\[
 E_{ab}=0\quad\text{for }a\ne b\text{ and }
 \sum_sj_sa_s\le\sum_sj_sb_s,
 \qquad E_{aa}\ne0.
 \tag{5}
\]

After ordering by energy within each occupation block, this is a
triangular matrix with nonzero diagonal. Therefore \(E\) and its
inverse give the regular constant physical endpoint repair from the
operator theorem. With the two- and three-direction geometries used
here, the energies within each occupation sector are distinct.
Condition (5) also explicitly excludes equal-energy off-diagonal
entries if a future geometry has a tie; this avoids treating nonzero
diagonal entries alone as an invertibility criterion.

The zero equations, (3), and the zero entries in (5) are all defined
over \(\mathbf F_2\). Let \(V\) be their binary solution space, and
let \(E_\Delta\subseteq\mathbf F_2^{2^r}\) be the image of the
endpoint diagonal map. Write \(\mathbf F_4=\mathbf F_2(\omega)\).
An \(\mathbf F_4\) solution has the form \(x+\omega y\), with
\(x,y\in V\). Its endpoint diagonal is nonzero in every coordinate
if and only if the supports of the two binary image vectors cover
all \(2^r\) coordinates. This is a finite exact criterion, including
the eight-coordinate case where separate nonzero projections are
insufficient over \(\mathbf F_4\).

The implementation computes the endpoint image and this covering
pair, lifts both image vectors back to the coefficient solution
space, and reconstructs the full Laurent kernel. It then independently
rechecks its filtered endpoint and every required zero row order.
No random choice or numerical rank decision occurs in feasibility.

## 3. What a feasible lattice certifies

If every row of a square matrix has order at least its specified
integer \(d_i\), the matrix factors as
\(\operatorname{diag}(q^{d_i})B\), with \(B\) regular. After a common
scalar shift makes the exponents nonnegative, the inclusion of image
lattices gives a surjection of cokernels onto that of the diagonal
matrix. Hence its sorted Smith orders are componentwise at least
the sorted \(d_i\). Equivalently, the same claim follows from the
invariant factors of a sublattice of a diagonal lattice.

Apply this separately in every occupation block of (4) and (2).
The zero targets \(\theta_i\) give \(z\)-orders \(2\theta_i\).
The Weyl reconstruction adds the explicit auxiliary shifts

\[
 -2P\sum_sj_su_s,
 \qquad u\in\{0,1\}^r,
 \tag{6}
\]

and assigns charge \(d-|u|\) when row \(i\) lies in occupation
block \(d\). The infinity targets are already \(z\)-orders.
Thus a feasible pair of lattices gives upper bounds on both active
partition sums, for every positive evaluation scale and every aggregate
charge. Restoring the unchanged spectators and the paired physical
partition supplies a sufficient pressure bound.

These are sufficient conditions for a candidate profile. Failure of
one prescribed pair of row lattices does not exclude another regular
row flag, another Laurent support, or the full operator family.
When the sufficient pressure bound is not improving, an individual
returned witness can still have better actual Smith orders; exact
profile replay is needed before discarding that witness.

## 4. Implementation and independent controls

- [weyl_operator_lattice_synthesis.py](../tools/weyl_operator_lattice_synthesis.py)
  constructs all occupation-preserving coefficient sources, the two
  independent lattices, and the endpoint image. The public entry
  points are `operator_basis`, `coefficient_system`, `synthesize`,
  and `bound_blocks`.
- [check_weyl_operator_lattice_synthesis.py](../tools/check_weyl_operator_lattice_synthesis.py)
  independently assembles literal weighted Weyl matrices and compares
  the solver with every constant two-mode \(\mathbf F_4\) kernel.
- [Recorded evidence](../evidence/weyl-operator-lattice-synthesis.json)
  enumerates all 4,096 constant coefficient tables, of which 324 have
  admissible filtered endpoints. Across three independent infinity
  row bases, 972 literal matrix calculations agree with 369 exact
  feasibility queries. Every one of the 77 feasible witnesses passes
  the separate literal-matrix test. Nine further controls cover
  positive and negative scalar Laurent powers, ramification three,
  and a split three-direction interval of length 32.

Replay with:

```sh
python3 tools/check_weyl_operator_lattice_synthesis.py \
  --output /tmp/weyl-operator-lattice-replay.json
```
