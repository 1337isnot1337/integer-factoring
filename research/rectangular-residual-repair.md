# Repairing an invertible rectangular residual core

## Result and scope

An invertible twelve-dimensional corner of the residual quartic endpoint
gives a rational crossing family with four extreme core states on one
crossing pair and three on the other. A joint output correction repairs
the local ramified endpoint exactly. Independent matching word types on
the two crossing pairs then give polynomial tensors in the original
variable, with the exact rectangular crossing endpoint.

The construction works for arbitrary finite mode sets containing their
minimum index \(0\) and maximum index \(T>0\). It retains all four
extreme states on one pair. The correction may affect the two flattening
profiles, which must be evaluated for the resulting tensor. This note
does not assert a pressure or factoring-exponent improvement.

## 1. The twelve-dimensional residual corner

Use the Gaussian \(G\), quartic \(Q=x_0y_Tu_0v_T\), and variables
\(c=1+w\), \(w=1+t^T\) from the
[coordinate-core construction](coordinate-core-crossings.md).
Start from \(G(1+Q/c^2)\). On the extreme modes, use the following
coordinate states:

\[
\begin{aligned}
(x,v)&:\quad (\text{low},\ \varnothing,\ \text{top},\ \text{both}),\\
(u,y)&:\quad (\text{low},\ \varnothing,\ \text{top}).
\end{aligned}
\tag{1}
\]

The corresponding core height lists, in units of \(T\), are

\[
h_X=(0,0,1,1),\qquad h_U=(0,0,1).
\tag{2}
\]

All bulk modes remain free. If \(b=2^{|S|-2}\), the four alphabet
dimensions are \((4b,3b,3b,4b)\), in the order \((x,y,u,v)\).
Both flattenings therefore have dimension \(12b^2\).

After the regular index gauge, the bulk endpoint is its crossing. The
extreme residual is the twelve-dimensional matrix \(F\) obtained by
restricting

\[
\Gamma(\text{crossing})
(1+u_Tv_0)(1+x_Ty_0)(1+x_0y_Tu_0v_T)
\tag{3}
\]

to (1). Rows use \(v+4y\), and columns use \(x+4u\). The following
table gives every nonzero column in each row of \(F\) and \(F^{-1}\).
All displayed entries are one; omitted entries are zero.

| Row pair | Nonzero input pairs in \(F\) | Nonzero source pairs in \(F^{-1}\) |
|---|---|---|
| \((0,0)\) | \((0,0),(3,1),(2,2)\) | \((0,0),(3,1)\) |
| \((1,0)\) | \((1,0),(2,1)\) | \((1,0),(2,1)\) |
| \((2,0)\) | \((2,0)\) | \((2,0)\) |
| \((3,0)\) | \((3,0)\) | \((3,0)\) |
| \((0,1)\) | \((0,1),(1,2)\) | \((0,1),(1,2)\) |
| \((1,1)\) | \((1,1)\) | \((1,1)\) |
| \((2,1)\) | \((2,1)\) | \((2,1)\) |
| \((3,1)\) | \((3,1),(2,2)\) | \((0,0),(2,2)\) |
| \((0,2)\) | \((0,2)\) | \((0,2)\) |
| \((1,2)\) | \((1,2)\) | \((1,2)\) |
| \((2,2)\) | \((0,0),(2,2)\) | \((0,0),(3,1),(2,2)\) |
| \((3,2)\) | \((3,2)\) | \((3,2)\) |

Literal multiplication verifies \(F^{-1}F=FF^{-1}=I\). In particular,
the inverse is not obtained by merely adding the off-diagonal part of
\(F\) to the identity.

## 2. A correction in the original variable

Write \(K\) for the rational tensor after restriction (1). Its numerator
is the restricted \(c^2G+QG\). On output core pairs define

\[
H_{(v,y),(v,y)}=c^{h_X(v)+h_U(y)}.
\]

Apply the joint output map

\[
J(c)=H^{-1}F^{-1}H,
\qquad \widetilde K=J(c)K.
\tag{4}
\]

Thus an entry of \(F^{-1}\) from source \((v',y')\) to target
\((v,y)\) is multiplied by

\[
c^{h_X(v')+h_U(y')-h_X(v)-h_U(y)}.
\tag{5}
\]

Every exponent is an integer between \(-2\) and \(2\). A universally
valid polynomial numerator for \(\widetilde K\) therefore has common
denominator \(c^4\). No claim that this denominator is minimal is
needed for the construction.

The direction of the conjugation in (4) follows from the energy gauge.
If \(\widehat K\) is the regular ramified tensor, then

\[
\widehat K=H_{\rm out}K H_{\rm in}^{-1},
\]

including the bulk index factors. Since \(F^{-1}\) acts only on the
core, the bulk factors commute with it. Consequently the energy-gauged
version of \(\widetilde K\) is precisely

\[
F^{-1}\widehat K.
\tag{6}
\]

It is regular in \(t\), and its endpoint is the identity on the core
pairs and the crossing on the bulk modes.

## 3. Polynomiality and independent paired types

Fix a type class \(\mathcal A\) of words in the \(4b\)-letter alphabet
and an independent type class \(\mathcal B\) in the \(3b\)-letter
alphabet. Restrict \(x,v\) to \(\mathcal A\) and \(u,y\) to
\(\mathcal B\). Then every entry of
\(\widetilde K^{\otimes m}\) is a polynomial in \(w\), and its value
at \(w=1\) is

\[
[x=v][u=y].
\tag{7}
\]

The index heights cancel separately on each crossing pair. Thus the
product of the ramified coefficient gauges is one on the restricted
entries. Equation (6) proves regularity after \(w=1+t^T\). The only
possible finite pole in the original rational tensor is at \(w=1\),
so this regularity proves polynomiality in \(w\) itself. Its endpoint
is (7), because the corrected ramified endpoint is already the exact
crossing before any type restriction.

This proof differs from an argument that simply discards residual
transitions of one energy sign: \(F\) has transitions in both
directions. The explicit inverse in (4) removes them first.

One can obtain equal physical dimensions by tensoring with the copy
having the two crossing pairs exchanged. The resulting type-class
dimension on every leg is \(|\mathcal A||\mathcal B|\). The two
flattening ranks multiply, so this pairing is compatible with the
ordinary rank-cost transfer. The common denominator must still be
included when computing the local Smith profiles: division by \(c^d\)
adds \(d\) to every infinity order and leaves all zero orders unchanged.
The [paired-type pressure theorem](paired-type-crossing-pressure.md)
gives the corresponding entropy normalization and independent Gibbs
optimization.

## 4. Exact finite seeds and replay

The [checker](../tools/check_rectangular_residual_repair.py) constructs
the numerator by minors of the direct one-particle matrix, applies the
full correction (4), cancels its actual common power of \(c\), and
records the remaining denominator. The following finite cancellations
are exact; they are not asserted as a general formula for all mode sets.

| Modes \(S\) | Ambient \(M\) | Alphabet dimensions \((d_X,d_U)\) | Numerator degree after cancellation | Remaining denominator |
|---|---:|---:|---:|---:|
| \(0,1\) | 2 | \((4,3)\) | 9 | \(1\) |
| \(0,1,2\) | 4 | \((8,6)\) | 26 | \(c\) |
| \(0,1,2,3\) | 4 | \((16,12)\) | 34 | \(c\) |

Run from the repository root:

```sh
python tools/check_rectangular_residual_repair.py \
  --output evidence/rectangular-residual-repair.json
```

The [evidence](../evidence/rectangular-residual-repair.json) includes a
gapped-set control \(S=\{0,2,5\}\), \(M=8\). It verifies the literal
residual corner, both inverse identities, every local ramified endpoint,
and every entry of every independent two-copy type pair. For the largest
listed seed this means 36,864 local coefficient checks and 136,896 typed
polynomial divisions and endpoint checks.

The `--seed-directory` option exports complete rational numerators with
their denominator degrees, logical heights, and coefficient hashes.
Their local `strict_endpoint` flag is false: the established assertion
is the polynomial crossing after the paired-type restriction. All
arithmetic in the checker is exact over \(\mathbf F_2\).

## 5. A polynomial two-parameter inverse core

There is an unrestricted parameter extension over every characteristic-two
field. Replace the third input \(u\)-core vector by

\[
|\text{top}\rangle+\alpha|\text{both}\rangle
\]

and the third output \(y\)-core covector by

\[
\langle\text{top}|+\beta\langle\text{both}|.
\tag{8}
\]

Keep the other states in (1). Both terms in each superposition have the
same index height \(T\). The energy gauge and the regularity proof
therefore apply unchanged.

Let \(F_0\) be the matrix in Section 1. With the same numerical indices
\(v+4y\), define

\[
\begin{aligned}
A={}&E_{0,9}+E_{3,10}+E_{7,9}+E_{10,9}+E_{11,8},\\
B={}&E_{4,10}+E_{7,3}+E_{8,11}+E_{9,10}+E_{10,3}.
\end{aligned}
\tag{9}
\]

The residual with (8) satisfies the exact polynomial identity

\[
F(\alpha,\beta)=F_0(I+\beta B)(I+\alpha A).
\tag{10}
\]

This can be checked by comparing the four coefficient matrices in
\(1,\alpha,\beta,\alpha\beta\). Relative to \(F_0\), those
three nonconstant matrices are respectively

\[
\begin{aligned}
F_\alpha&=E_{0,9}+E_{3,10}+E_{11,8},\\
F_\beta&=E_{8,11}+E_{9,10}+E_{10,3},\\
F_{\alpha\beta}&=E_{8,8}+E_{9,9}+E_{10,10}.
\end{aligned}
\tag{11}
\]

The matrices in (9) have

\[
A^2=E_{3,9},\qquad B^2=E_{4,3}+E_{9,3},
\qquad A^3=B^3=0.
\]

Consequently

\[
\det F(\alpha,\beta)=1,
\quad
F(\alpha,\beta)^{-1}
=(I+\alpha A+\alpha^2A^2)
 (I+\beta B+\beta^2B^2)F_0^{-1}.
\tag{12}
\]

Thus no exceptional parameter values or endpoint phase conditions are
needed. Use the height-dressed version of this inverse in (4). The
paired-type polynomial crossing theorem applies to all \(\alpha,\beta\).
It also applies after polynomial substitutions
\(\alpha=\alpha(w)\), \(\beta=\beta(w)\): the inverse remains
polynomial in these parameters and its value at the endpoint is the
required inverse residual. Laurent-polynomial substitutions can be used
after a common power of \(w\) clears their poles at zero. That clearing
power equals one at \(w=1\).

The checker verifies all four binary corners directly from the exterior
endpoint, the four coefficient identities in (10), and the nilpotence
identities. These finite matrix identities prove (10)–(12) as polynomial
identities over every characteristic-two field.

For the doubly mixed corner \(\alpha=\beta=1\), the following finite
seeds have an additional exact cancellation: their entire \(c^4\)
denominator cancels entrywise.

| Modes \(S\) | Ambient \(M\) | Alphabet dimensions \((d_X,d_U)\) | Polynomial degree | Positive-height endpoint extras |
|---|---:|---:|---:|---:|
| \(0,1\) | 2 | \((4,3)\) | 9 | 10 |
| \(0,1,2\) | 4 | \((8,6)\) | 25 | 73 |
| \(0,1,2,3\) | 4 | \((16,12)\) | 33 | 512 |

These are local polynomial tensors with a triangular endpoint. Every
extra endpoint entry strictly increases the actual total index height,
so independent paired types eliminate it. Their denominators, endpoint
extras, full ramified endpoints, and all two-copy paired-type products
are checked exactly. The entrywise cancellation in this table is a
finite result; the general theorem uses the safe rational denominator
and does not assume the cancellation for larger mode sets.

## 6. A smaller repair that preserves the low-mode background

For the original corner \(\alpha=\beta=0\), full inversion is not
necessary. Put

\[
C_0=I+E_{10,0}+E_{10,7}.
\tag{13}
\]

Then \(C_0^2=I\), and \(U=C_0F_0\) has identity diagonal and only
strictly height-decreasing off-diagonal entries. Explicitly,

\[
U=I+E_{0,7}+E_{0,10}+E_{1,6}+E_{4,9}+E_{7,10}.
\tag{14}
\]

The dressed correction is just

\[
J_{\rm partial}(c)=H^{-1}C_0H
 =I+c^{-2}E_{10,0}+c^{-1}E_{10,7}.
\tag{15}
\]

Its output row for the all-low core state \((v,y)=(0,0)\) is unchanged.
Thus it preserves that Gaussian background exactly. This can be useful
when optimizing rare core excitations.

The regular ramified tensor now has endpoint \(U\), instead of the
identity. Every extra transition in (14) strictly lowers total height.
Independent paired types have equal input and output total height, so
they eliminate every term using an extra transition. The same argument
as in Section 3 proves polynomiality and the exact typed crossing.
There is no need to multiply by \(U^{-1}\).

The checker verifies (13)–(15), their height inequalities, and all local
ramified and two-copy independent-type checks for the three finite mode
sets in the first table. These tensors are exported separately with
`repair_kind` equal to `partial`; their numerator degrees and remaining
denominators are recorded exactly in the evidence.

For completeness, the mixed family has principal determinants
\(1,(1+\alpha\beta)^2,1\) on the successive height subspaces of
dimensions \(4,10,12\). This follows directly from (11): the middle
principal block is upper triangular, with its last two diagonal entries
equal to \(1+\alpha\beta\). Thus the same block-elimination principle
works whenever \(\alpha\beta\ne1\). The doubly mixed binary point
\(\alpha=\beta=1\) has a singular middle principal block, so (13)
does not supply its repair; the full polynomial inverse (12) remains
valid there.

### An accessible non-Gaussian control

The partial repair retains a non-Gaussian interaction even when there are
no bulk modes. For \(S=\{0,1\}\), \(M=2\), its local tensor is
polynomial, has vacuum coefficient one, and has the exact Wick residual

\[
\begin{aligned}
\mathcal W_{x_0y_Tu_0v_T}
={}&K_\varnothing K_{x_0y_Tu_0v_T}
 +K_{x_0y_T}K_{u_0v_T}
 +K_{x_0u_0}K_{y_Tv_T}
 +K_{x_0v_T}K_{y_Tu_0}\\
={}&(1+w)^4\ne0.
\end{aligned}
\tag{16}
\]

Every coefficient in (16) belongs to the retained rectangular alphabet:
neither the \(u\) nor the \(y\) leg uses its excluded both-occupied
state. Therefore no choice of the missing physical coefficients can
complete this tensor to a fermionic Gaussian. This is an obstruction on
the actually available coordinates, not an arbitrary zero completion.

The checker verifies (16) in the physical variable order
\((x_0,x_T,y_0,y_T,u_0,u_T,v_0,v_T)\), where the quartic has mask
153. It also checks all 41 fully accessible quartic Wick relations;
five fail. This supplies an exact non-Gaussian control for the partial
repair, independently of its numerical rank pressure.
