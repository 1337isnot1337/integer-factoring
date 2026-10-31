# Signed residual codes over arbitrary fields

## Scope

This note studies a finite, signed exterior residual tensor. It proves that
the three-letter constant-code obstruction extends to this tensor over
every field, and gives an exact three-conjugate phase repair over
\(\mathbb F_{27}\).

These are endpoint statements. This note does not claim a complete
odd-characteristic polynomial crossing construction, a rank-pressure
certificate, or an improved factoring exponent. In particular, any future
transfer to a factoring construction must retain the stated local Fock
ordering and verify its polynomial endpoint and tensor-power conventions.

## 1. A sign-safe definition

Order the one-particle input modes as
\((X_0,X_T,U_0,U_T)\), and output modes as
\((V_0,V_T,Y_0,Y_T)\). Let

\[
 D_0=I_4+E_{0,3}+E_{2,1}.
\]

For a matrix \(D\), write \(\Gamma(D)\) for its action on the full
exterior algebra. Its entry between two ordered subsets of the same size
is the corresponding determinant; entries of different sizes are zero.
Creation and annihilation operators use the usual sign
\((-1)^{\text{number of occupied earlier modes}}\).

Define the integer matrix

\[
 R=\Gamma(D_0)+
   a^\dagger_{V_T}a^\dagger_{Y_T}\,
   \Gamma(D_0)\,b_{U_0}b_{X_0}.
 \tag{1}
\]

Thus the rightmost annihilation operator removes \(X_0\) first. The
definition fixes the signs before reduction modulo a field characteristic.
All entries of \(R\) belong to \(\{-1,0,1\}\).

View each pair of modes on one physical leg as a four-state space ordered
\((\varnothing,\{0\},\{T\},\{0,T\})\). The tensor entry is
\(R_{v+4y,x+4u}\). The unperturbed block identity is therefore the ordinary
crossing \([x=v][y=u]\); no additional swap sign is introduced by this
physical-leg convention.

The inverse one-particle matrix is
\(D_0^{-1}=I_4-E_{0,3}-E_{2,1}\). Its contraction minor with input indices
\((X_0,U_0)\) and output indices \((V_T,Y_T)\) is

\[
 \det\begin{pmatrix}0&-1\\-1&0\end{pmatrix}=-1.
\]

If \(N=(R-\Gamma(D_0))\Gamma(D_0)^{-1}\), the signed
creation-annihilation relations give

\[
 N^2=-N,\qquad \operatorname{rank}N=4.
 \tag{2}
\]

Hence \(I+N\) is a rank-twelve projection, and \(R\) has rank twelve
over every field. The checker verifies (2) as an integer matrix identity.
In characteristic two, (1) agrees with the residual endpoint used in
[the binary quartic code construction](residual-quartic-codes.md).

Ordinary tensor powers of the four-index tensor (1) are well defined over
every field. No identification with a globally reordered exterior algebra
is needed for the finite code statements below. This distinction avoids
silently introducing signs when moving fermion modes between sites.

## 2. The three-letter obstruction extends to every field

Flatten \(R\) with rows \((x,v)\) and columns \((y,u)\). Its rank is eight.
Its column space, viewed as a space of four-by-four matrices, is

\[
 W_-=
 \left\{
 \begin{pmatrix}
 t&b&0&0\\
 0&t&e&f\\
 c&d&t&-b\\
 0&c&g&s
 \end{pmatrix}:t,b,c,d,e,f,g,s\in F
 \right\}.
 \tag{3}
\]

Put \(J=\operatorname{diag}(1,1,1,-1)\), and let \(W\) be the matrix
space with the same pattern but \(+b\) in position \((2,3)\). Then
\(W_-=WJ\), and the row space of the flattening is \(W^TJ\).

These are field-independent spanning statements. The checker verifies all
integer pattern equations and exhibits coordinate determinants \(\pm1\)
for eight selected columns and eight selected rows. Thus no exceptional
prime is hidden in the identification.

### Compression lemma over an arbitrary field

Every compression \(A\mapsto LAK\) of \(W\), where
\(L:F^4\to F^3\) is surjective and \(K:F^3\to F^4\) is injective,
has image dimension at least five.

The proof is the same matrix-space argument as in
[the characteristic-two obstruction](residual-quartic-three-letter-obstruction.md),
with signs retained. Write \(\ker L=Fa\) and
\(\operatorname{im}K=\ker b^T\). The killed matrices are exactly
\(ar^T+sb^T\). The parameterization has one-dimensional redundancy

\[
 (r,s)=(\lambda b,-\lambda a).
\]

The four zero-entry constraints, two diagonal equalities, and two repeated
entry equalities defining \(W\) impose at least four independent linear
conditions on \((r,s)\). For completeness, the rank argument separates
the following cases.

- If \(a_0b_0\ne0\), the four zero-entry equations give four distinct
  pivots. If only one of \(a_0,b_0\) is nonzero, two zero-entry pivots and
  the diagonal/repeated-entry equations give four pivots. Transposition
  and exchange of coordinates one and two identify the two cases.
- If \(a_0=b_0=0\) and both zero-entry pivots \(r_0,s_0\) occur, the
  remaining equations contain two nonzero blocks on disjoint variable
  sets. Each contributes at least one more rank.
- If exactly one of those pivots occurs, the nonzero remaining coordinate
  of \(a\) or \(b\) gives three additional independent equations. If
  neither occurs, \(a_2b_1\ne0\), and four independent equations remain.

Replacing an equality by a subtraction changes none of these nonzero
pivots. The solution space therefore has dimension at most four, and
quotienting its one-dimensional redundancy bounds the compression kernel
on \(W\) by three. Since \(\dim W=8\), the image has dimension at least
five. The statement also holds for \(W^T\), and multiplication by the
invertible diagonal matrix \(J\) does not change it.

### Theorem

No four constant physical linear maps restrict (1) to a crossing on three
letters, over any field.

Indeed, a successful physical map has rank three. Factor the rank-eight
crossing-cut matrix through an eight-dimensional space. The two compressed
factors each retain rank at least five by the lemma, so their product has
rank at least \(5+5-8=2\). A three-letter crossing has rank one at this
cut. This is a contradiction.

The theorem concerns a single copy and constant physical maps. It does
not exclude phase gates, temporal composition, polynomial maps, or codes
using additional modes.

## 3. The signed two-letter phase family

Use the input columns

\[
 B_pe_0=|\{0\}\rangle,\qquad
 B_pe_1=|\varnothing\rangle+p|\{0,T\}\rangle
\]

and output rows

\[
 C_r^{(0)}=\langle\{0\}|,\qquad
 C_r^{(1)}=\langle\varnothing|+r\langle\{0,T\}|.
\]

An integer symbolic expansion of (1) gives

\[
 (C_r\otimes C_r)R(B_p\otimes B_p)
  =\operatorname{diag}(1,1+pr,1+pr,1+2pr).
 \tag{4}
\]

The last diagonal entry explains an important difference from
characteristic two. Two nontrivial members of this family cannot be inverse
phase gates in odd characteristic. If their parameters are \(a,b\), the
two required equations are

\[
 (1+a)(1+b)=1,\qquad(1+2a)(1+2b)=1.
\]

Subtracting twice the first expanded equation from the second gives
\(2ab=0\). Thus \(ab=0\) and then \(a+b=0\), so both parameters are zero.

## 4. Three conjugates over \(\mathbb F_{27}\) repair the phase

Let \(\alpha\) satisfy

\[
 \alpha^3+\alpha^2-1=0
 \tag{5}
\]

over \(\mathbb F_3\). The cubic has no root in \(\mathbb F_3\), so it
defines \(\mathbb F_{27}\). Its elementary symmetric functions are
\(e_1=2,e_2=0,e_3=1\). Consequently,

\[
 \operatorname{Norm}(1+\lambda\alpha)
  =1+2\lambda+\lambda^3=1
  \qquad (\lambda=1,2).
 \tag{6}
\]

Take the three phase gates (4) with \(p=1\) and
\(r=\alpha,\alpha^3,\alpha^9\). Repeat each logical core bit in all
three copies. Every diagonal phase multiplies to one by (6), so the
grouped two-letter endpoint is an exact crossing.

This is an exact finite residual-code identity. A complete construction
with independent bulk occupations would still need its polynomial
regularity, physical-type transfer, and effective coefficient-rank bounds
verified with the odd-characteristic sign conventions. Those additional
claims are not part of this note.

## Exact reproduction

```sh
python tools/check_signed_residual_codes.py --output evidence/signed-residual-codes.json
```

The checker constructs (1) with literal signed creation and annihilation
matrices, verifies (2) over the integers, identifies both crossing-cut
spaces with unimodular coordinate witnesses, expands (4) symbolically,
and implements \(\mathbb F_{27}\) directly to check all grouped endpoint
entries. It also records the three flattening ranks over
\(\mathbb F_2,\mathbb F_3,\mathbb F_5,\mathbb F_7\).

The first of those ranks is 15 over \(\mathbb F_3\), rather than 16 as
over the other three checked fields. The column rank remains 12 and the
crossing-cut rank remains eight; this extra finite rank drop does not
evade the three-letter obstruction.
