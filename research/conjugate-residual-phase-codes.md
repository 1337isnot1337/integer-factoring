# Conjugate residual codes and descent to a binary crossing

## Status

This note constructs three explicit eight-letter polynomial crossings. They
start with two conjugate copies of a rational tensor over \(\mathbb F_4\),
repair a residual endpoint phase, and then descend to \(\mathbb F_2\).
Every polynomial division and every endpoint entry is checked exactly.

These are new candidates for coefficient-rank and physical-type optimization.
No improvement to the accepted factoring exponent is claimed. Numerical
optimization is deliberately separate from the construction certificate.

## 1. The residual phase

Use the sparse one-species Gaussian from
[the interval construction](schur-one-species-interval.md). Let the mode set
\(S\) contain \(0\) and \(T=\max S>0\), and let the ambient power of two
\(M\) satisfy \(M>T\). Put \(c=1+w\). Consider

\[
 \mathcal T(w)=G(w)
 \left(1+\frac{x_0y_Tu_0v_T}{c^2}\right).
 \tag{1}
\]

The general sparse kernel is polynomial for this index set. The argument
below uses inclusion support and index bounds; it does not require closure
of \(S\) under subsets.

To describe the endpoint, introduce \(t\) by \(w=1+t^T\), and multiply
each physical coefficient by
\(t^{\text{output index sum}-\text{input index sum}}\). All Gaussian
entries become regular in \(t\). The quartic has index increase \(2T\),
so its coefficient becomes one. The only additional Gaussian edges that
survive at \(t=0\) are \(u_Tv_0\) and \(x_Ty_0\). The extreme-mode
endpoint is the tensor

\[
 R=\Gamma(\text{crossing})
       (1+u_Tv_0)(1+x_Ty_0)(1+x_0y_Tu_0v_T).
 \tag{2}
\]

The bulk modes already have the crossing endpoint. The polynomiality and
the local form (2) are also explained in
[the localized quartic note](localized-quartic-profiles.md).

Order the four extreme-mode states as
\((\varnothing,\{0\},\{T\},\{0,T\})\). For scalars \(p,r\), use
the input columns

\[
 B_pe_0=|\{0\}\rangle,\qquad
 B_pe_1=|\varnothing\rangle+p|\{0,T\}\rangle
\]

and output rows

\[
 C_r^{(0)}=\langle\{0\}|,\qquad
 C_r^{(1)}=\langle\varnothing|+r\langle\{0,T\}|.
 \tag{3}
\]

Here core bit zero is the low-only state; core bit one is the vacuum/both
superposition. Direct expansion of (2) gives

\[
 (C_r\otimes C_r)R(B_p\otimes B_p)
 =\mathrm{diag}(1,a,a,1),\qquad a=1+pr,
 \tag{4}
\]

when outputs are ordered \((v,y)\) and inputs \((x,u)\). In tensor
notation, the nonzero entries still require \(x=v\) and \(y=u\), and
their scalar factor is \(a^{x\mathbin\oplus u}\).

## 2. Two conjugate copies repair the phase

Let \(\omega\in\mathbb F_4\) satisfy \(\omega^2=\omega+1\), and write
\(\overline z=z^2\) for coefficient conjugation. Choose \(p=1,r=\omega\),
so \(a=\omega^2\). Use the conjugate code in a second copy.

Keep both copies' bulk occupations independent, but require their logical
core bits to agree on each physical leg. The two endpoint phases multiply
to

\[
 a^{x\oplus u}\overline a^{x\oplus u}=1,
\]

because \(a\overline a=a^3=1\). The grouped endpoint is therefore an
ordinary crossing. With \(k=|S|-2\) bulk modes, its logical alphabet has
size \(2^{2k+1}\): two independent bulk blocks and one repeated core bit.

This step uses two copies, not three identical copies. It preserves all
bulk freedom and reduces only the independent core freedom.

## 3. The original variable suffices after a common type restriction

Applying (3) directly to (1) gives a rational tensor in the original
variable \(w\). Its only possible pole is at \(w=1\). A uniform height
identity controls that pole.

For one logical label \(\ell\), let \(b(\ell)\) be its core bit and
\(n_i(\ell)\) its bulk occupations. Define the possibly half-integral height

\[
 h(\ell)=\sum_{i\in S\setminus\{0,T\}}
              \left(i-\frac T2\right)n_i(\ell)
              +\frac T2 b(\ell).
 \tag{5}
\]

To verify the identity, consider any nonzero physical term contributing to
a coded entry. Write \(B_{\rm in},B_{\rm out}\) for the total numbers of
both-occupied extreme states on its input and output legs. Physical charge
conservation holds for both \(G\) and its quartic perturbation. Since a
logical core state has particle count
\(1-b+2[\text{both}]\), charge conservation gives

\[
 2(B_{\rm out}-B_{\rm in})
   =\Delta b-\Delta N_{\rm bulk},
\]

where \(\Delta\) means output minus input. The physical index increment
is consequently

\[
 E=\Delta\!\sum_i i n_i+T(B_{\rm out}-B_{\rm in})
   =h(y)+h(v)-h(x)-h(u).
 \tag{6}
\]

It is independent of which components of the superpositions (3) were used.
On physical support it is an integer, even when (5) uses half-integers.

Thus the regular, energy-gauged coded tensor is exactly \(t^E\) times
the rational coded tensor evaluated at \(w=1+t^T\). For the grouped
construction, add the two heights in (5). In a tensor power restricted to
one common logical word type, the total height difference is zero. The
entire product therefore equals a polynomial in \(t\), with the crossing
endpoint furnished by Section 2.

The original product has no possible poles except at \(w=1\). Regularity
after \(w=1+t^T\) excludes that pole, so the original product is itself
a polynomial in \(w\). Its degree is bounded linearly in the number of
copies, since the original rational seed has fixed numerator and denominator
degrees. This gives a polynomial typed crossing without replacing the
original variable by a degree-\(T\) parameter map in the final construction.

The common denominator must be handled consistently when evaluating local
Smith profiles. A two-copy numerator has denominator \(c^4\). Its
zero-place orders are unchanged by this denominator, whereas every
infinity order increases by four. This restoration is justified by the
global polynomiality argument, rather than by discarding a denominator.

## 4. An eight-letter local polynomial crossing

Take \(S=\{0,1,2\}\) and \(M=4\). Its only bulk mode has
\(i-T/2=0\), so a single-copy height is simply its core bit. Encode a
grouped physical label as

\[
 \ell=\beta_1+2\beta_2+4b,
 \qquad \beta_1,\beta_2,b\in\{0,1\}.
 \tag{7}
\]

The grouped energy increment is \(2\Delta b\). Since \(c=t^2\), the
entrywise normalization

\[
 K_{xyuv}(w)=
 c^{b(y)+b(v)-b(x)-b(u)}\mathcal T^{\rm grouped}_{xyuv}(w)
 \tag{8}
\]

is already a local polynomial tensor. Its endpoint at \(w=1\) is the
strict eight-letter crossing, before any type restriction.

For an explicit construction, let \(N\) be the one-copy coded numerator
of \(c^2G+x_0y_Tu_0v_TG\). Each grouped entry is a product of an entry
of \(N\) and a conjugated entry of \(N\). Equation (8) divides this
product by

\[
 c^{4-\Delta b},\qquad
 \Delta b=b(y)+b(v)-b(x)-b(u).
 \tag{9}
\]

All 4,096 divisions are exact. All 4,096 endpoint entries equal
\([x=v][y=u]\). The polynomial degree is 48.

The same finite construction, with the same normalization (8), also gives
strict polynomial eight-letter crossings for \(S=\{0,1,3\}\) and
\(S=\{0,2,3\}\), again with ambient \(M=4\) and degree 48. Those two
additional local statements are independently certified by exact polynomial
division and exhaustive endpoint checks. The midpoint-height argument above
is the symbolic local proof for \(S=\{0,1,2\}\); a general local-polynomial
claim for arbitrary \(S\) is not needed here.

The sign in (8) is material. Using input minus output instead fails 29
polynomial divisions for \(S=\{0,1,2\}\). The checker retains this negative
control.

## 5. Galois descent gives a binary tensor of the same degree

Let \(\pi\) swap the two bulk-copy blocks in every physical label and fix
the repeated core bit. The grouped construction satisfies

\[
 \overline K=\pi^{\otimes4}K.
 \tag{10}
\]

This semilinear symmetry yields a concrete change of physical basis. For
each fixed label of \(\pi\), keep its basis vector. For each two-cycle
\((a,b)\), use the two input basis columns

\[
 e_a+e_b,\qquad \omega e_a+\omega^2e_b.
 \tag{11}
\]

Let the resulting matrix be \(B\), and use \(C=B^{-1}\) on each output
leg. On the two-cycle, the matrices are

\[
 B_{\{a,b\}}=\begin{pmatrix}1&\omega\\1&\omega^2\end{pmatrix},
 \qquad
 C_{\{a,b\}}=\begin{pmatrix}\omega^2&\omega\\1&1\end{pmatrix}.
 \tag{12}
\]

Their determinant is one. The columns in (11) are fixed by
\(\pi\circ\overline{(\cdot)}\), so

\[
 \overline B=\pi B,\qquad \overline C=C\pi.
\]

For the transfer matrix with input legs \((x,u)\) and output legs
\((v,y)\), equation (10) therefore gives

\[
 \overline{(C\otimes C)K(B\otimes B)}
   =(C\otimes C)K(B\otimes B).
 \tag{13}
\]

Every coefficient of the transformed tensor lies in \(\mathbb F_2\).
The input and output maps are inverse, so the crossing endpoint is
preserved. The constant basis changes preserve every unrestricted
coefficient-flattening rank. They do not raise the polynomial degree.
For the three explicit seeds above, the descended degree remains 48.

Although unrestricted ranks are unchanged, the new binary physical basis
can expose different diagonal gauges and different coordinate type
restrictions. This is why optimization in both bases is useful.

## 6. Reproduction and independent checks

Run:

```sh
python tools/check_conjugate_residual_phase_codes.py --output evidence/conjugate-residual-phase-codes.json
```

To export the complete \(\mathbb F_4\) and \(\mathbb F_2\) tensors, add
`--seed-directory /tmp/conjugate-phase-seeds`. The four-index order is
`x+8*y+64*u+512*v`. The original tensors use two bits per polynomial
coefficient; the descended tensors use ordinary binary polynomial bits.

The checker verifies all 256 residual-code entries for all choices
\(p,r\in\mathbb F_4\), every polynomial division, both strict endpoints,
the conjugation/swap identity, the inverse basis matrices, and selected
coefficient ranks in both flattenings before and after descent. The evidence
records canonical coefficient hashes, degree, support size, and the
wrong-sign negative control.

An independent bit-plane implementation produced exactly the same 4,096
binary coefficient masks for each of the three descended seeds. Separate
native \(\mathbb F_4\) validation checked rational-weight Smith queries
against an independent fraction-free oracle before numerical optimization.
These construction checks do not turn a floating-point optimization result
into an exponent certificate.
