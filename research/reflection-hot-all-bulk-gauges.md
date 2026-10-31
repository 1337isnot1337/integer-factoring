# All physical-leg bulk orderings with a bounded reflection core

## Result and scope

Two explicit regular cross pivots remove the ordering restriction from
the [separate-leg bulk construction](reflection-hot-leg-bulk-gauges.md).
Every nonnegative choice of the eight physical-leg bulk heights is
supported, independently at the two cuts. The active support remains
\(\min(M,4k)\) for a power-of-two hot-core size \(k\).

Relative hot-state heights must still satisfy the exact hot-cold
regularity inequalities from that construction. Internal factors
joining two hot variables are retained, and the tensor modification
must be an even exterior polynomial. The theorem does not cover an
arbitrary coefficient mask or every possible physical state gauge.
It supplies a broader exact evaluator, not an improved exponent.

## 1. A second regular cross pivot

At either cut, reverse the logical mode indices and write the natural
cross matrix in the boundary order
\((\ell_1,\ell_2)\mid(r_1,r_2)\) as

\[
K=\begin{pmatrix}DT&D\\D&DT\end{pmatrix},\qquad
D_i=z^{P(2+4i)},\quad 0\le i<M/2.
\tag{1}
\]

The compression theorem gives

\[
T^2=B I,\qquad B=1+F,\qquad F=z^{2MP},\qquad
S=DTD^{-1},\qquad S^2=B I.
\tag{2}
\]

Both \(T\) and \(S\) are regular matrices, and \(B\) is a unit
at the local place. The previous identity-pivot factors were

\[
L_I=\begin{pmatrix}0&I\\I&S\end{pmatrix},\qquad
R_I=\begin{pmatrix}I&T\\0&I\end{pmatrix},\qquad
L_IKR_I=\operatorname{diag}(D,FD).
\tag{3}
\]

The second factorization is

\[
\boxed{
L_T=\begin{pmatrix}I&0\\S/B&I\end{pmatrix},\qquad
R_T=\begin{pmatrix}T/B&I\\0&T\end{pmatrix},\qquad
L_TKR_T=\operatorname{diag}(D,FD).
}
\tag{4}
\]

For example, the lower-left block of \(L_TK\) is
\(DT^2/B+D=0\); its lower-right block is \(FDT/B\).
Multiplying by \(R_T\) gives (4). Every displayed factor and
its inverse is regular, because the only denominator is a power
of the local unit \(B\).

With uniform per-leg bulks, (4) remains regular whenever

\[
b_{\ell_2}\ge b_{\ell_1},\qquad
b_{r_2}\ge b_{r_1}.
\tag{5}
\]

Its two weighted off-diagonal blocks acquire exactly these
nonnegative differences. The lower diagonal half receives left
bulk \(b_{\ell_1}\) and right bulk \(b_{r_1}\); the upper half
receives \(b_{\ell_2}\) and \(b_{r_2}\).

The cross matrix (1) is unchanged by simultaneously interchanging
its two row groups and its two column groups. Applying that coordinate
permutation when \(b_{r_1}>b_{r_2}\), then choosing (3) or (4)
according to the left ordering, covers all four orderings. These are
permutations of the locally gauged boundary variables. The physical
labels, defect coefficients, and their original natural heights are
transported explicitly; no global tensor relabeling is assumed.

In every branch the lower-half spectator shift is

\[
\min(b_{\ell_1},b_{\ell_2})+
\min(b_{r_1},b_{r_2}),
\tag{6}
\]

and the upper-half shift is the corresponding sum of maxima.

## 2. The hot CAR support stays bounded

The change from the identity-pivot diagonal basis to (4) is

\[
P=L_TL_I^{-1}
 =\begin{pmatrix}S&I\\0&S/B\end{pmatrix},\qquad
P^{-\mathsf T}
 =\begin{pmatrix}S^{\mathsf T}/B&0\\I&S^{\mathsf T}\end{pmatrix}.
\tag{7}
\]

A creation vector is multiplied by \(P\), and its annihilation
vector by \(P^{-\mathsf T}\). This preserves all canonical
anticommutation pairings. It also preserves the diagonal Gaussian
when the companion column change is included.

Put \(s=z^{4P}\), \(c=1+z^{2P}\), and let
\(0\le j<k\) index a selected hottest logical bit. The only
spreading vectors in (7) are columns of \(S\) and
\(S^{\mathsf T}\):

\[
S_{rj}=c\,s^{r-(j\mathbin{\&}r)},\qquad
(S^{\mathsf T})_{rj}=c\,s^{j-(j\mathbin{\&}r)}.
\tag{8}
\]

Here \(\&\) is bitwise intersection. On a residue class
\(r=a+k\ell\), with \(k\) a power of two and \(j<k\), the
first sequence is a scalar times \((s^k)^\ell\), and the
second is constant. The additional factors \(B^{-1}\), \(F\),
and the selected diagonal order are independent of \(\ell\).

Thus the exact same two regular residue sweeps used for the cross-only
hot-state core concentrate every creation and annihilation vector
onto two positions per residue in each diagonal half. Their Gaussian
companions remain regular, and their bulk weights cancel because
the sweeps stay within a half. The resulting support has at most
\(4k\) modes, or all \(M\) modes when the interval is smaller.

The exact eight physical generators per selected bit are implemented
in
[`reflection_hot_alternate_coordinates.py`](../tools/reflection_hot_alternate_coordinates.py).
That file applies the existing support concentration and coefficient
exporter through explicit data records. It introduces no new native
elimination routine.

## 3. Weighted vacuum coordinates and complete coefficients

The normalized hot creation and annihilation operators in the new
basis remain regular under (5). The left spreading tails gain
\(b_{\ell_2}-b_{\ell_1}\); on the column side they gain
\(b_{r_2}-b_{r_1}\). After a simultaneous group swap the same
statement uses the transported physical labels.

A principal unit vacuum minor always works, by the projector argument
in the preceding note. The implementation may instead reuse a
nonprincipal unit minor whose row and column labels have equal total
occupations in each diagonal half. This is sufficient: diagonal bulk
conjugation multiplies its determinant by

\[
z^{b_A(\sum |I_A|-\sum |J_A|)
    +b_B(\sum |I_B|-\sum |J_B|)}=1.
\tag{9}
\]

The two exact integer equalities in (9) are checked before a cached
basis is accepted. Its coordinate height is still the physical hot
state height plus the diagonal-half bulk occupation of the selected
vacuum row or column. The old parity-coordinate exporter therefore
computes complete coefficient matrices for this basis as well.

Powers of \(B\) are cleared only as local unit denominators in the
profile matrix. They are not discarded from coefficient-width bounds.
The source builder doubles its requested precision until both the
ordinary program degree bound and the actual retained entry width
certify a complete polynomial matrix. Only then is that matrix passed
to the adaptive Smith engine. If the complete source exceeds the
supported precision, the construction rejects it; widening an
elimination workspace is never used to replace missing coefficients.

The global tensor, its common Laurent clearing, its ordinary degree,
and its repaired endpoint are exactly those of the original exterior
defect construction. No local denominator in (4) is inserted into
that global tensor.

## 4. Interface and verification

[`reflection_hot_all_bulk_gauges.py`](../tools/reflection_hot_all_bulk_gauges.py)
exposes `AllBulkHotStatePressure`, with the same offset and bulk
arguments as the separate-leg evaluator. It chooses the appropriate
branch separately at each cut and caches that branch's complete
coefficient matrices. Physical paired-type energies and all spectators
are accounted for by the preceding separate-leg formulas.

The [independent checker](../tools/check_reflection_hot_all_bulk_gauges.py)
and [saved extended replay](../evidence/reflection-hot-all-bulk-gauges.json)
verify fifty complete original weighted cuts and 398,080 endpoint
coefficients. These include all four orderings, independently chosen
branches at the two cuts, nonlinear two-bit states, and a ramification
of three. Forty-eight exact generator identities independently check
the change (7) against the previous physical operators.

Twenty-four weighted-basis cases verify 256 normalized creation and
annihilation operators and 48 unit vacuum minors. The support checks
extend through \(M=64\), and all 68 analytic derivatives agree with
centered finite differences with maximum error below
\(1.55\times10^{-10}\).

The source-expansion control starts with a requested precision of 128
at \(M=128\). It increases the retained source window to 2,048 and
the builder precision to 4,096. The program degree bound is 2,310,
the cleared local-unit denominator has degree 1,524, and the largest
actual entry coefficient has degree 1,782. These are deliberately
recorded as separate quantities. The final source is certified complete
before any weighted Smith calculation.

The shared coordinate and support helpers were also replayed through
the full [original hot-state regression](../evidence/reflection-hot-state-gauges.json):
56 complete weighted cuts and 328,960 endpoint coefficient checks
passed. Its separate adaptive-elimination control still widens the
Smith workspace from 512 to 4,096 while retaining a complete source
whose maximum coefficient width is fourteen bits.

To rebuild the full record, run:

```sh
python tools/check_reflection_hot_all_bulk_gauges.py --extended \
  --output evidence/reflection-hot-all-bulk-gauges.json
```
