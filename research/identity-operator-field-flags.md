# Extension-field row flags for operator synthesis

This extends the [binary-coefficient lattice solver](identity-operator-lattice-synthesis.md)
to arbitrary constant \(\mathbf F_4\) row bases at infinity.
The selected polynomial kernel still contains every occupation-preserving
matrix coefficient. A fixed \(\mathbf F_{16}\) extension supplies
endpoint coefficients when an admissible full-support endpoint exists
over a sufficiently large field.

It also gives a small exact source formula at prime-power lengths.
The resulting changes improve the search procedure and establish
additional valid families. They do not assert an improved factoring
exponent.

## 1. Nonbinary infinity bases

Let \(W\) be an invertible constant \(\mathbf F_4\) matrix
preserving total occupation, and use \(G=WF\) as the unknown kernel.
The infinity degree caps remain

\[
 \deg_qG_{ab}\le
 \left\lfloor\frac{\kappa_b-\eta_a}{2}\right\rfloor.
 \tag{1}
\]

At zero, the source for the coefficient \(G_{ab}\) is the
\(\mathbf F_4\) linear combination of the original Weyl sources
given by column \(a\) of \(W^{-1}\). The prescribed zero row
orders therefore remain homogeneous linear equations, now over
\(\mathbf F_4\). The endpoint equations are obtained from
\(F(1)=W^{-1}G(1)\), not from the diagonal of \(G(1)\).

As before, the endpoint is required to be triangular in physical
\(Y\) energy with nonzero diagonal. The implementation explicitly
forbids equal-energy off-diagonal terms in this particular linear
solver. The general operator engine separately permits an arbitrary
invertible endpoint block within an equal-energy level.

There are at most eight endpoint diagonal coordinates for the current
two- and three-mode solver. Over \(\mathbf F_{16}\), feasibility of
a nonzero endpoint diagonal is equivalent to every coordinate being
a nonzero functional on the homogeneous solution space: eight proper
hyperplanes cannot cover a vector space over a sixteen-element field.
The fixed extension is independent of the interval length. For
\(\mathbf F_4\) output, the solver instead checks the full finite
endpoint image.

The packed linear solver retains both \(\mathbf F_4\) bit planes.
After reconstruction it checks the complete linear system, the
filtered endpoint, and every prescribed zero row valuation. The
[independent extension-field profiler](f16-identity-operator-profiles.md)
can then compute the full actual profiles of a returned kernel.

## 2. Small sources from the scalar Pascal contraction

For \(n=2^k\), put \(m=nP\) and \(c=1-q^m\). The
[simultaneous reflected change](prime-power-joint-identity-filters.md)
gives \(C=cI\) for every selected subset of indices. For binary
states \(a,b\), the Weyl formula gives

\[
 c^r\Psi_C(E_{ab})
 =c^{|b|}\sum_{R\subseteq a\cap b}
   |(a\setminus b)\cup R\rangle
   \langle(b\setminus a)\cup R|.
 \tag{2}
\]

Indeed, multiplication by \(Y^{1-b}\) requires the input state to
be a subset of \(b\). The derivatives \(\partial^{1-a}\)
then require every element of \(b\setminus a\) to be occupied.
The remaining input subset \(R\subseteq a\cap b\) is free,
and the output is \((a\setminus b)\cup R\). All nonzero
binary differentiation coefficients equal one.

The zero row indexed by \(s\) also receives the monomial weight
\(q^{P\sum_i j_i s_i}\). For the Gaussian kernel \(F=I\),
equation (2) is diagonal, with entry

\[
 q^{P\sum_i j_i s_i}\,c^{|s|}\,q^{m(r-|s|)}.
 \tag{3}
\]

Thus its row order is

\[
 P\sum_i j_i s_i+m(r-|s|).
 \tag{4}
\]

Sorting these diagonal entries within each occupation block supplies
a Gaussian Smith basis without elimination.

The previous polynomial Weyl compiler clears a denominator \(\Delta\)
and has \(C^{-1}=N/\Delta\). In the scalar case,
\(N=N_{00}I\) and \(\Delta=cN_{00}\). Its complete source
matrices are therefore exactly \(N_{00}^r\) times the sources
in (2). This is a common zero unit. Removing it preserves every
entry valuation and every row-lattice constraint; it makes no
assumption about regularity of that unit at infinity.

## 3. Arbitrary selected operator modes

The shared operator interfaces now accept geometries such as
`selected:3,17,29`. At a prime-power interval length, the constructor
uses

\[
 F_x=I,\qquad F_v=(1-\tau)^{n-1}I,\qquad
 \operatorname{retained}(j)=n-1-j.
 \tag{5}
\]

The selected zero forms are the independent pairs

\[
 \bigl(z^{a_j}Y_j+(1-\tau^n)V_{n-1-j}\bigr)
 \bigl(U_j-z^{-a_j}X_{n-1-j}\bigr).
 \tag{6}
\]

The proof of the reflected change is regular for both individual
linear forms. Hence it applies to \(L^aR^b\) and to off-diagonal
operator kernels, as well as to diagonal filters. The common
normalizing factors cancel because \(|a|=|b|\), exactly as in the
original operator theorem. The global polynomial kernel, its endpoint
repair, its degree bound, and its infinity profile are unchanged.

## 4. Exact implementation controls

The new components are:

- [f4_weyl_operator_lattice_synthesis.py](../tools/f4_weyl_operator_lattice_synthesis.py)
  for nonbinary infinity flags and \(\mathbf F_4\) or
  \(\mathbf F_{16}\) endpoint lifts;
- [scalar_pascal_operator_basis.py](../tools/scalar_pascal_operator_basis.py)
  for (2)--(4);
- [selected_identity_mode_filters.py](../tools/selected_identity_mode_filters.py)
  for arbitrary selected indices.

The [flag checker](../tools/check_f4_operator_lattice_flags.py) compares
all 4,096 constant two-mode \(\mathbf F_4\) coefficient tables.
There are 324 admissible filtered endpoints. Across all five
two-dimensional \(\mathbf F_4\) flags, 1,620 literal matrix
calculations agree with all 420 feasibility comparisons and all
90 returned witnesses. It also checks 213 complete source entries
against the earlier compiler, including their exact common unit.
The [record](../evidence/f4-operator-lattice-flags.json) is reproducible.

The [arbitrary-selection checker](../tools/check_selected_identity_operators.py)
compares seven full active and reduced profiles over characteristics
two, three, and five. Four original physical cuts at matrix order
256 agree for an interior pair and a ramified nonconsecutive triple;
18,981 endpoint coefficients are checked. Its
[record](../evidence/selected-identity-operators.json) also includes
invalid-geometry rejection controls.

```sh
python3 tools/check_f4_operator_lattice_flags.py \
  --output /tmp/operator-flags-replay.json
python3 tools/check_selected_identity_operators.py --extended \
  --output /tmp/selected-operators-replay.json
```
