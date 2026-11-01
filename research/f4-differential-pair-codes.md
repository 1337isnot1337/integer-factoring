# F4 differential codes and a matroid Hilbert formula

## Construction and scope

Homogeneous differential equations give deletion-closed physical codes
whose changes occur below the final height boundary. This extends the
[existing binary construction](coherent-hot-compound-retractions.md#21-homogeneous-differential-equations-give-another-construction)
to F4. It is outside the fixed-core hypothesis of the
[boundary pressure floor](hot-boundary-natural-pressure-floor.md).

The constructor and the dimension formula below are exact. Six initial
complete-profile trials did not improve the hyperoval crossing bound.
A smaller physical-code dimension is not itself a pressure improvement
or a faster factoring algorithm.

## 1. Common kernels and a coordinate retraction

Write D_a for deletion of occupation bit a. All these operators commute
and square to zero in characteristic two. Let T_1,...,T_r be nonzero
linear combinations of deletion monomials, each homogeneous in physical
height and differential degree. Inside a coordinate height downset, set

```math
V=\bigcap_j\ker T_j.
```

Every D_a preserves V, because it commutes with every T_j and decreases
height. Row reduction separately in each height/occupation cell gives
a homogeneous basis. Its leading occupation words form a deletion
downset: when a leading word contains bit a, deleting that bit gives
the leading word of the deleted vector. Integer order on the masks
preserves comparison after subtracting the same occupied bit.

Choose a reduced basis P with identity on these leading coordinates,
and let L be coordinate evaluation there. Then LP=I. The column space
of P is deletion closed, and the coordinate space of L is deletion
closed because its leading-word set is a downset. Thus both natural
coherent cut interfaces apply.

`tools/coherent_differential_kernels_f4.py` constructs this basis without
an unspecified choice of complement. It checks the differential
equations, LP=I, both deletion closures and the image/nullity counts.

## 2. Disjoint paired derivatives

Suppose the operators have the form

```math
T_j=\sum_{i=1}^m C_{ji}D_{a_i}D_{b_i},
```

where the pairs (a_i,b_i) are disjoint and all pair heights
h(a_i)+h(b_i) equal one common H. Let r(S) be the rank of the coefficient
columns indexed by S. Use x to count height and y to count occupation.
The common-kernel Hilbert polynomial is

```math
\boxed{
\sum_{S\subseteq[m]}
 (1+x^H y^2)^{|S|-r(S)}
 \prod_{i\notin S}y\bigl(x^{h(a_i)}+x^{h(b_i)}\bigr)
 \prod_{c\text{ unused}}(1+x^{h(c)}y).
}                                                        \tag{1}
```

**Proof.** The four occupation states of pair i split into the invariant
two-dimensional space spanned by its vacuum and double occupation,
and the two singleton states killed by D_ai D_bi. Choose S to be the
pairs using the former space. The inactive pairs contribute the displayed
singleton factors.

On the active pairs, write z_i for double occupation. Their space is
K[z_i:i in S]/(z_i squared), and the operators are constant-coefficient
derivations with coefficient matrix C_S. In characteristic two, any
invertible linear variable change preserves this square-zero algebra.
It reduces the independent derivations to coordinate derivatives.
Their common kernel is the square-zero polynomial algebra in
|S|-r(S) remaining variables. Each has height H and occupation two,
so its Hilbert polynomial is the first factor in (1). Unused modes
are unrestricted. The invariant direct sum over S proves the formula.

Setting x=y=1, with no unused modes, gives

```math
\dim V=2^m\sum_{S\subseteq[m]}2^{-r(S)}.                 \tag{2}
```

The weighted formula requires the common pair height. It does not infer
a crossing rank from a physical-code dimension.

## 3. A field-extension difference not visible over F2

Take pairs (3,4), (2,5), (1,6), (0,7), all of height 16, and

```math
C=\begin{pmatrix}
1&1&1&1\\
0&1&\omega&1+\omega
\end{pmatrix},\qquad \omega^2=\omega+1.
```

Every two columns are independent, so r(S)=min(2,|S|). Equation (2)
gives the full eight-mode kernel dimension

```math
16\left(1+\frac42+\frac{6+4+1}{4}\right)=92.           \tag{3}
```

Over F2 there are only three nonzero column directions in dimension
two. Four nonzero columns therefore contain a repeated pair. Relative
to the rank pattern in (3), that pair adds at least
16(1/2-1/4)=4 to (2). A zero column increases the dimension further.
Thus every binary two-equation system on these four pairs has kernel
dimension at least 96. Columns (1,0),(0,1),(1,1),(1,0) attain 96.

This is a genuine difference between the available physical modules.
It is not a better pressure bound. The 92-dimensional full eight-mode
kernel in (3) is also a different code from the 92-label low-height
core used in the boundary theorem.

## 4. Exact checks and completed profile trials

The [kernel controls](../evidence/coherent-differential-kernels-f4-controls.json)
compare three binary reduced bases, nine independent dense F4 nullities
and 1,344 literal annihilation equations. The
[Hilbert controls](../evidence/differential-pair-hilbert-controls.json)
compare every weighted cell in ten cases against the constructed
kernel bases. They also enumerate the 78 rank-two binary matrices
with four nonzero columns and reproduce the minimum dimension 96.

The [six completed profile trials](../evidence/coherent-differential-kernels-f4-screen.json)
use height caps 20, 22 and 24 on twelve modes:

- The two-equation F4 system has dimensions 48, 60 and 74, with
  numerical balanced pressures approximately 0.7115538741,
  0.7113797760 and 0.7112776576.
- The single equation with coefficients (1,1,omega,1+omega) has
  dimensions 52, 66 and 83, with numerical pressures approximately
  0.7112239692, 0.7109542579 and 0.7107709894.

Every trial retains both complete profiles and source fingerprints.
These numerical values are search diagnostics, not exact strict
certificates or global minima. Neither family beats the current
hyperoval construction. The formula helps compare coefficient choices
with the same weighted physical dimensions, but a favorable two-cut
profile remains a separate obligation.

Run `python3 tools/check_coherent_differential_kernels_f4.py` and
`python3 tools/check_differential_pair_hilbert.py` for the quick controls.
The separate search script reconstructs the longer complete-profile
trials. No new Lean work or independent acceptance is claimed.

Substantial OpenAI assistant support was used for the derivation,
implementation, computations, and text.
