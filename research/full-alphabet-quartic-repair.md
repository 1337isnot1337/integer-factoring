# Full-alphabet quartic crossings in characteristic two

## Result and scope

Changing the extreme quartic coefficient from one to a value
\(\gamma\ne1\) makes its full sixteen-dimensional residual invertible.
All four physical extreme states can then be retained on both crossing
pairs. An explicit partial correction has only five off-diagonal entries
and becomes the identity at \(\gamma=0\). Independent matching types
give polynomial crossing tensors of every admissible length.

The choices \(\gamma=\omega,\omega^2\) over \(\mathbb F_4\) give
non-Gaussian examples, certified by a failed Gaussian minor identity.
This is a new admissible construction family. No pressure witness below
\(1/\sqrt2\), or factoring exponent above \(1/2\), is asserted.

## 1. The full residual and its rank-four projector

Use the Gaussian \(G\), extreme modes \(0,T\), and
\(c=1+w\) of [the coordinate-core construction](coordinate-core-crossings.md).
For any finite mode set containing these extremes, start with

```math
K_\gamma=G\left(1+\frac{\gamma Q}{c^2}\right),
\qquad Q=x_0y_Tu_0v_T.
\tag{1}
```

The full physical extreme alphabet, on every boundary, is ordered
\((\mathrm{low},\varnothing,\mathrm{top},\mathrm{both})\).
Its height vector in units of \(T\) is \(h=(0,0,1,1)\).
Rows of the residual matrix use \(v+4y\), and columns use \(x+4u\).

After \(w=1+t^T\) and the actual occupied-index gauge, the bulk
endpoint is its crossing. The residual \(F_\gamma\) is obtained by
expanding

```math
\Gamma(\mathrm{crossing})
(1+u_Tv_0)(1+x_Ty_0)(1+\gamma x_0y_Tu_0v_T).
\tag{2}
```

Write \(F_\gamma=F_0+\gamma Q_0\) and
\(P=F_0^{-1}Q_0\). Direct binary multiplication gives

```math
F_0^2=I,\qquad P^2=P,\qquad \mathrm{rank}(P)=4.
\tag{3}
```

For an independently reproducible description, encode each binary row
of \(P\) as an integer whose bit \(j\) is column \(j\). The sixteen
row integers are

```text
1, 0, 0, 8, 0, 0, 0, 1, 0, 0, 1, 4096, 4096, 1, 8, 32768.
```

An idempotent decomposes the vector space into its image and kernel.
Thus, over every field of characteristic two,

```math
\det F_\gamma=(1+\gamma)^4,
\qquad
F_\gamma^{-1}
=\left(I+\frac{\gamma}{1+\gamma}P\right)F_0^{-1}
\quad(\gamma\ne1).
\tag{4}
```

At \(\gamma=1\), the residual has rank twelve; this is the
rank-deficient point used by the rectangular construction. For all
other coefficients it has rank sixteen.

## 2. A sparse partial correction

Put \(\delta=1+\gamma\). Define a sixteen by sixteen matrix
\(C_\gamma\) as follows. Begin with the identity, replace the four
diagonal entries numbered \(10,11,14,15\) by \(\delta^{-1}\),
and put \(\gamma/\delta\) in positions

```math
(10,0),\quad(10,7),\quad(10,13),\quad(11,12),\quad(14,3).
\tag{5}
```

All unspecified off-diagonal entries are zero. If
\(h(i)=h(i\bmod4)+h(\lfloor i/4\rfloor)\), then

```math
C_\gamma F_\gamma=U_\gamma,\qquad
(U_\gamma)_{ii}=1,\qquad
(U_\gamma)_{ij}\ne0\ (i\ne j)\Longrightarrow h(i)<h(j).
\tag{6}
```

The checker verifies (6) as a polynomial identity over
\(\mathbb F_2[\gamma]\), after multiplying by \(\delta\).
It does not infer the identity from its four field specializations.
The correction is block lower triangular in height order, with
determinant \(\delta^{-4}\). At \(\gamma=0\), it is exactly the
identity and (1) is the original Gaussian tensor.

For the complete inverse repair, replace \(C_\gamma\) by
\(F_\gamma^{-1}\); its endpoint is exactly the identity even before
matching physical types. Both repairs are available in the constructor.

## 3. Height dressing and polynomial typed powers

On output core pairs let \(H_{ii}=c^{h(i)}\), and define

```math
\widetilde K_\gamma=H^{-1}C_\gamma H K_\gamma.
\tag{7}
```

The original numerator in (1) has denominator \(c^2\). Each dressed
correction entry contributes a power between \(c^{-2}\) and \(c^2\).
Therefore \(c^4\widetilde K_\gamma\) is polynomial. This safe bound
is sufficient; the recorded bulk examples reduce to a single factor
\(c\) after exact cancellation.

The ramified, index-gauged version of (7) is simply
\(C_\gamma\widehat K_\gamma\). It is regular in \(t\), with
endpoint (6) on the core and the crossing on the bulk. Fix independent
matching full-letter types on \((X,V)\) and \((Y,U)\). Their total
occupied-index heights agree, so the ramified gauge cancels from each
restricted tensor-power entry. Every endpoint term using a strict
height decrease in (6) is impossible at equal total height. The endpoint
is consequently

```math
[x=v][y=u].
\tag{8}
```

The tensor-power entries are rational functions with no possible finite
pole other than \(w=1\). Regularity after the ramification rules out
that pole as well. Hence they are polynomials in \(w\). For each fixed
seed, their degrees grow at most linearly with the word length. The
usual effective coefficient and finite rank-factorization interfaces
therefore apply, subject to proving a strict pressure bound.

The partial correction leaves the all-low output row unchanged, and
the quartic is absent from that row. Its all-low background is still
the undeformed Gaussian. Retaining the fourth state on \(Y,U\)
changes the physical partition and the leading-pressure problem;
the rectangular leading-face certificates do not automatically apply.

## 4. Exact controls and search use

Run:

```sh
python3 tools/check_full_alphabet_quartic_repair.py --full \
  --output evidence/full-alphabet-quartic-repair.json \
  --seed-directory /tmp/full-alphabet-quartic-seeds
```

The [checker](../tools/check_full_alphabet_quartic_repair.py) constructs
entries from independent Gaussian minors, and separately expands the
full exterior product. It checks the ramified endpoint coefficient by
coefficient, polynomial division of two-copy independently typed
entries, the rank-four projector, and the symbolic partial identity.

The [recorded evidence](../evidence/full-alphabet-quartic-repair.json)
contains twenty-four cases: both repairs, coefficients \(0,\omega,\omega^2\),
and mode sets \(\{0,1\}\), \(\{0,1,2\}\), \(\{0,2,5\}\),
and \(\{0,1,2,3\}\). It checks 443,904 exterior coefficients,
443,904 ramified endpoint entries, and 1,653,600 independently typed
two-copy entries. All Gaussian controls satisfy every tested two by two
compound identity; each nonzero nonbinary case has an explicit failure
of such an identity, proving it is not a Gaussian compound tensor.

The two nonbinary coefficients are exchanged by the field automorphism
\(a\mapsto a^2\). This automorphism preserves every weighted minor
valuation. Thus their pressure functions agree at the same physical
weights, and only one needs numerical optimization. Floating searches
remain exploratory and do not certify a factoring exponent.

## 5. Fixed-size profile calculation for long intervals

For a full interval, the quartic and either correction are polynomials
in the same extreme creation and annihilation operators used in the
[bounded-defect proof](rectangular-defect-factorization.md). Removing
the three-state projections does not introduce any new operator.
Consequently the seven active modes at infinity and nine at zero
still suffice. All core coordinate states are now retained, giving
full-rank matrices of orders \(128\) and \(512\) for \(M\ge8\),
with the same explicit scalar spectator factors.

The full-rank promise follows from the all-length typed endpoint.
If either local square flattening of order \(N\) had rank less than
\(N\), its tensor-power ranks would be exponentially too small to
contain the crossing on two uniform physical type classes, whose
dimension is \(N^m\) up to a polynomial factor. Regular invertible
core-coordinate maps preserve this rank. This argument supplies the
rank promise independently of truncated series calculations.

Allowing weights on all four extreme physical cores uses the
[ten-mode extension](rectangular-full-core-factorization.md), now with
all \(1024\) coordinate states retained. The smaller calculation has
23 continuous parameters: four bulk charges, three weights on each of
\(V,Y,U\) at each place, and \(\log s\). The full-core version adds
six independent \(X\) weights, for 29 parameters.

The [bounded program](../tools/full_alphabet_quartic_profiles.py) and
[profile checker](../tools/check_full_alphabet_quartic_profiles.py)
were compared with the independent original four-mode constructor.
The [profile evidence](../evidence/full-alphabet-quartic-profiles.json)
records 589,824 exact entries over \(\mathbb F_{16}\) and eighteen
complete quarter-weight profiles, including a Gaussian control.
Both independent full-tensor Gibbs roots agree with their reduced
calculations to less than \(3\cdot10^{-16}\), and all 52 derivative
coordinates agree with finite differences to less than
\(4\cdot10^{-10}\). The latter comparisons are numerical
implementation controls, separate from the exact matrix identities.

```sh
python3 tools/check_full_alphabet_quartic_profiles.py \
  --output evidence/full-alphabet-quartic-profiles.json
python3 tools/search_rectangular_core_pressure.py --M 8 \
  --full-alphabet-gamma 2 --seconds 3600 \
  --output /tmp/full-alphabet-M8-search.json
```
