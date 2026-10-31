# Signed global Pascal crossings over every field

## Result and scope

The one-species Pascal construction extends to every field. The extension
uses a global Gaussian tensor before imposing the common physical type.
This convention retains the fermionic signs in odd characteristic and
gives an ordinary crossing endpoint.

Its normalized infinity Smith orders are still the union of the original
mode set and its greedy Pascal complement. The complement is now computed
over the chosen field. Consequently the same explicit pressure formula
can be searched in any fixed prime characteristic. Initial intervals
retain their previous profiles; noninterval sets can have different
profiles in different characteristics.

This is a constructive extension, with complete signed endpoints and both
Gaussian rank-profile interfaces. It does not establish a pressure witness
below \(1/\sqrt2\) or a stronger factoring exponent. In particular,
ordinary tensor powers of unsigned local determinants are not the
construction in this note.

## 1. A polynomial kernel without a characteristic restriction

Let \(K\) be a field and let \(S\subseteq\mathbb Z_{\ge0}\) be finite
and nonempty, with \(n=|S|\). Choose an even integer
\(M>2\max S\). All matrices below are indexed by \(S\). Put

~~~math
P_{ab}=\binom ab,\qquad
R_{ab}=(-1)^{a-b}\binom ab.                            \tag{1}
~~~

Binomial coefficients are reduced in \(K\), and vanish when \(b>a\).
Both matrices are lower unitriangular. For a general noninterval set,
no claim that \(R=P^{-1}\) is needed.

Define the signed Delannoy polynomial

~~~math
\begin{aligned}
O_{ab}(t)
&=(-1)^{a+b}[X^aY^b]\frac1{1-t(X+Y+XY)}\\
&=(-1)^{a+b}\sum_{k=\max(a,b)}^{a+b}
  \binom ka\binom a{k-b}\,t^k.                        \tag{2}
\end{aligned}
~~~

In particular, \(t^{\max(a,b)}\) divides \(O_{ab}\), and its degree
is at most \(a+b\). Its coefficient at degree \(a\) is
\((-1)^{a+b}\binom ab\).

With row boundaries \((V,Y)\) and column boundaries \((B,X)\), set

~~~math
\boxed{\quad
\mathcal D_S(w)=
\begin{pmatrix}
\mathrm{diag}((1+w)^{a+1})P^{\mathsf T}&w^M R\\
I&w^{M-1}O(w^{-1})
\end{pmatrix}.
\quad}                                                \tag{3}
~~~

Every entry is polynomial: the exponent \(M-1\) clears all degrees
in (2). At \(w=-1\), the upper-left block vanishes and the upper-right
block is \(R\), since \(M\) is even. The lower-right block need not
vanish. The common physical type will remove its endpoint contribution.

## 2. Global mode order makes the endpoint an ordinary crossing

Fix a copy count \(m\) and mode counts \(0\le k_a\le m\). A physical
label is a word of \(m\) subsets of \(S\), with mode \(a\) occupied
exactly \(k_a\) times. There are

~~~math
D_m=\prod_{a\in S}\binom m{k_a}                       \tag{4}
~~~

such labels. Every label has the same total charge
\(k=\sum_a k_a\).

Replace each scalar entry of (3) by that entry times \(I_m\). The result
is a global one-particle matrix \(\mathcal D_S^{[m]}\), with mode/copy
indices ordered in the same way at each physical boundary. For labels
\((x,y,u,v)\), select rows from \((v,y)\) and columns from \((u,x)\).
Each side has exactly \(2k\) selected indices. Define

~~~math
\Gamma_m(w;x,y,u,v)
=(-1)^k\det\mathcal D_S^{[m]}(w)[(v,y),(u,x)].          \tag{5}
~~~

**Endpoint theorem.**

~~~math
\boxed{\quad
\Gamma_m(-1;x,y,u,v)=[x=v][y=u].
\quad}                                                \tag{6}
~~~

**Proof.** At the endpoint, each of the \(k\) selected \(V\) rows can
use only an \(X\) column. There are exactly \(k\) such columns, so
they are all consumed. The remaining \(Y\) rows must use the \(B\)
identity block. Thus the lower-right block contributes no matching, and
\(y=u\) is necessary.

The other endpoint factor is a minor of \(R\otimes I_m\). Every
nonzero entry matches an input mode to a mode at least as large. The
common mode counts equalize the sum of mode indices on \(x\) and \(v\).
Every matched increase must therefore be zero. Equal modes are connected
by the identity on the copy indices, so \(x=v\). Conversely, when the
labels agree, the surviving principal minor is one.

Interchanging the selected \(B\) and \(X\) column blocks contributes
\((-1)^{k^2}\). The factor \((-1)^k\) in (5) cancels it, because
\(k^2\equiv k\pmod2\). This proves (6) over every field.

### Independent matching types on the two crossing pairs

The same argument permits mode counts \(k^X_a\) on \((X,V)\) and
different counts \(k^Y_a\) on \((Y,B)\). Put
\(K_X=\sum_a k^X_a\) and \(K_Y=\sum_a k^Y_a\). There are
\(K_X\) selected \(V\) rows and exactly that many \(X\) columns,
so the upper-right block consumes all of them at the endpoint.
The remaining identity block gives \(y=u\), and the pairwise matching
mode counts again force \(x=v\).

The endpoint determinant sign is \((-1)^{K_XK_Y}\). For the Gaussian
exterior coefficient, the Pfaffian sign contributes
\((-1)^{\binom{K_X+K_Y}{2}}\), giving the constant endpoint scalar

~~~math
\beta=(-1)^{\binom{K_X}{2}+\binom{K_Y}{2}}.             \tag{6a}
~~~

Multiplying the whole typed tensor by \(\beta\) normalizes the endpoint
to an ordinary crossing. All physical boundary charges are fixed, so
every boundary-block reordering used for either cut still has one constant
sign. Both rank-profile arguments below remain valid.

### Why a global construction is necessary

In odd characteristic, the local selected determinants generally have
endpoint signs depending on two local occupation counts. Multiplying
those local tensors does not turn that sign into a function of the common
global type.

For one mode and two copies, take the common type with one occupied
position. The diagonal endpoint with both labels occupied in the same
copy has local determinant product \(-1\); placing their occupations in
different copies gives \(+1\). No common scalar corrects both entries.
The checker includes this negative control over \(\mathbb F_3\).

Equation (5) avoids that problem by defining the global matrix first.
It is not identified with an ordinary tensor power of local determinants.

## 3. Signed Gaussian profiles at both cuts

For an alternating matrix \(\Omega\), define its Gaussian exterior
element by its principal Pfaffians. Equivalently, in ordered Grassmann
variables \(e_i\), it is the product

~~~math
\prod_{i<j}(1+\Omega_{ij}e_i e_j).                     \tag{7}
~~~

Its degree-two factors commute. This definition uses no division by
factorials, so it is valid in every characteristic.

Apply it to

~~~math
\Omega=
\begin{pmatrix}0&\mathcal D_S^{[m]}\\
-(\mathcal D_S^{[m]})^{\mathsf T}&0\end{pmatrix},
~~~

in boundary order \((V,Y,B,X)\). A coefficient with \(r\) selected
plus variables and \(r\) selected minus variables is
\((-1)^{r(r-1)/2}\) times the selected determinant. Here \(r=2k\),
so this sign is \((-1)^k\). The coefficient is exactly (5).

For either crossing cut, reorder the four boundary blocks so each side
is contiguous. Each physical boundary has charge \(k\). Every exchanged
pair of boundary blocks therefore contributes the fixed sign
\((-1)^{k^2}\). After restriction, the reordered Gaussian coefficient
matrix differs from (5) by a single constant sign. Thus either cut can
be analyzed in its own contiguous Grassmann order.

At zero, use the cut \((X,Y)\mid(V,B)\). The internal blocks are the
two diagonal blocks of (3); the cross blocks are \(w^M R\) and \(I\),
up to transposition and signs. If their internal gauged entries are
regular, the corresponding internal Gaussian multiplication maps are
invertible over the local power-series ring: they have scalar fermion
degree zero equal to one. Removing them leaves the direct sum of exterior
powers of the cross matrix.

At infinity, use \((X,B)\mid(V,Y)\). Every minus variable is on one
side and every plus variable on the other, so there are no internal
blocks. The flattening is the exterior compound of the whole kernel,
up to harmless diagonal signs on exterior degrees.

These are ordinary matrix rank statements after the physical restriction.
They do not identify an ungraded local tensor power with a globally
reordered exterior algebra. For \(m\) copies, each one-particle Smith
order is repeated \(m\) times, since the global cross matrix is the
corresponding one-particle matrix tensored with \(I_m\).

## 4. The greedy Pascal profile over an arbitrary field

Let \(D=\mathrm{diag}(t^a)\) and
\(Q=\mathrm{diag}((1+t)^a)\). Consider

~~~math
F_S(t)=
\begin{pmatrix}
QP^{\mathsf T}D&DR\\
D&(1+t)O
\end{pmatrix}.                                        \tag{8}
~~~

For \(c\ge0\), let
\(u_c=(\binom ca)_{a\in S}\). Let \(c_1<\cdots<c_n\) be the greedy
independent columns among \(c\notin S\), computed over \(K\).
Their existence and the bound \(c_n\le2\max S+1\) follow from
[translated consecutive Pascal bases](pascal-complement-convolution.md).

**Profile theorem.** The Smith orders of (8) are

~~~math
\boxed{\quad \mathrm{sort}(S\cup\{c_1,\ldots,c_n\}).
\quad}                                                \tag{9}
~~~

**Proof.** Regular block elimination gives

~~~math
F_S\sim\mathrm{diag}(D,K_S),\qquad
K_S=DR-(1+t)QP^{\mathsf T}O.                           \tag{10}
~~~

The last column elimination is regular because row \(a\) of \(O\)
is divisible by \(t^a\).

The generating function (2) gives the formal identity

~~~math
(1+t)^{a+1}\sum_{c\ge0}\binom ca O_{cb}(t)
=t^aR_{ab}.                                          \tag{11}
~~~

To verify it, substitute \(X=-1-Z\) into
\((1-t(X+Y+XY))^{-1}\), obtaining
\((1+t+tZ(1+Y))^{-1}\), and extract \([Z^aY^b]\).
The substitution is valid \(t\)-adically: the coefficient at index
\(c\) has valuation at least \(c\).

Removing the indices in \(S\) from (11) gives

~~~math
K_S=U\,\mathrm{diag}(t^c:c\notin S)\,V,
\quad U_{a,c}=(1+t)^{a+1}\binom ca,
\quad V_{c,b}=t^{-c}O_{cb}(t).                        \tag{12}
~~~

The constant column of \(U\) is \(u_c\). The constant row of \(V\)
is \((-1)^c u_c^{\mathsf T}\mathrm{diag}((-1)^b)\).

For each size \(j\), the first \(j\) greedy columns form the unique
independent set of minimum sum of indices. Cauchy--Binet therefore bounds
every \(j\)-minor of (12) below by valuation \(c_1+\cdots+c_j\).
Choose the same row and column index set on which these greedy columns
have a nonzero minor. Its leading contribution is a nonzero signed
square of that minor. The minimum-weight index set is unique, so no
other term can cancel it. The determinantal minimum is exactly
\(c_1+\cdots+c_j\). Thus \(K_S\) has orders \(c_1,\ldots,c_n\),
proving (9).

For \(S=\{0,\ldots,n-1\}\), every field gives
\(C=\{n,\ldots,2n-1\}\). For a noninterval set the field matters.
For example,

| Modes \(S\) | Field | Greedy complement |
|---|---|---|
| \(\{1,2\}\) | \(\mathbb F_2\) | \(\{3,5\}\) |
| \(\{1,2\}\) | \(\mathbb F_3\) | \(\{4,5\}\) |
| \(\{0,1,2,4\}\) | \(\mathbb F_2\) | \(\{3,5,6,7\}\) |
| \(\{0,1,2,4\}\) | \(\mathbb F_3\) | \(\{3,5,7,8\}\) |
| \(\{0,1,2,4\}\) | \(\mathbb F_5\) | \(\{3,5,6,9\}\) |

Larger complement orders lower the displayed infinity profile sum. This
comparison alone does not provide a below-threshold pressure witness.

## 5. Integer gauges and the same pressure criterion

Set \(w=-u^2\), so the polynomial endpoint occurs at \(u=1\).
The substitution doubles all valuations and changes coefficient signs;
it does not change coefficient ranks. List heights in boundary order
\((B,X,V,Y)\). Use

~~~math
\begin{aligned}
h_0(a)&=(0,\ 2a-2M+1,\ 0,\ 2a+1),\\
h_\infty(a)&=(2a-n+1,\ 2M-n-1,\ 2a-n+2,\ -n).
\end{aligned}                                         \tag{13}
~~~

At zero, (2) makes the lower-right internal entry regular, since its
order is at least \(2(M-1-a-b)\). The upper-left internal block is
regular as well. Both cross blocks have orders
\(\{2a+1:a\in S\}\).

At infinity, remove the common scalar of order \(1-2n\). Regular unit
row and column factors reduce the remaining matrix to (8), with
\(t=-u^{-2}\). Thus, writing \(\nu\) for the union list (9), the
infinity orders are \(2\nu+1-2n\). The restored physical heights are

~~~math
G_a=8a-4(n-1).                                        \tag{14}
~~~

Here the reduction uses only units in the original field. To see its
signs explicitly, put \(E=\mathrm{diag}((-1)^a)\). After removing the
common scalar, multiply the \(V\) rows by \(-E\), the \(B\) columns
by \(E\), and the \(X\) columns by \(-I\). The block matrix becomes
\(\bigl(\begin{smallmatrix}(1+t)QP^{\mathsf T}D&DR\\D&O\end{smallmatrix}\bigr)\).
Multiplying the \(V\) rows by \((1+t)^{-1}\) and the \(X\) columns
by \(1+t\) gives (8). No square root of \(-1\) is adjoined.

Adding a common height \(h\) to the two plus boundaries at infinity
adds \(h\) to every infinity order and \(2h\) to each physical height.
Therefore the exact criterion for an independent product type is

~~~math
\begin{aligned}
\mathcal B_{S,h}(s)
&=\prod_{a\in S}(1+s^{-(2a+1)})^2
  \prod_{\nu\in S\cup C_K}(1+s^{-(2\nu+1-2n+h)}),\\
\mathcal Z_{S,h,R}(s)
&=\prod_{a\in S}
 (1+s^{-[8a-4(n-1)+2h]/R}),\\
\mathcal B_{S,h}(s)&<\mathcal Z_{S,h,R}(s)^R,
\qquad R=4\rho.                                      \tag{15}
\end{aligned}
~~~

This is the same formula as in the binary Pascal note, with \(C_K\)
computed in the selected field. General independent mode heights are
also permitted, provided the internal blocks remain regular. Their
restoration is constant on the common global type, just as in
[fixed-mode crossings](fixed-mode-crossings.md).

For the independent matching types above, split the restored heights as

~~~math
G_X(a)=\sum_{e\in\{0,\infty\}}(h_{e,X}(a)+h_{e,V}(a)),
\qquad
G_Y(a)=\sum_{e\in\{0,\infty\}}(h_{e,Y}(a)+h_{e,B}(a)).   \tag{15a}
~~~

Let \(\mathcal B(s)\) be the product of the two profile sums for these
gauges. The [paired-type transfer](paired-type-crossing-pressure.md),
applied to this signed global construction, gives the criterion

~~~math
\boxed{\quad
\mathcal B(s)<[Z_X(s)Z_Y(s)]^R,\qquad R=2\rho,
\quad Z_D(s)=\prod_{a\in S}(1+s^{-G_D(a)/R}).
\quad}                                                  \tag{15b}
~~~

The two optimal mode-probability vectors can differ. Tensoring the
rectangular crossing with the exchanged copy restores equal physical
label dimensions, keeping the two coefficient indices separate. The
rank products multiply and the physical logarithmic dimensions add,
which accounts for the normalization \(R=2\rho\).

For the displayed gauges (13), including the common shift \(h\),
both pair heights equal \(4a-2(n-1)+h\). Their optimizing distributions
therefore agree, and (15b) becomes (15). Independent height choices
enlarge the paired-type search beyond that symmetric specialization.

## 6. Effectivity and factoring transfer

Fix a finite mode set, polynomial kernel, field \(\mathbb F_p\), and
heights before the copy count grows. The global selected determinants
have matrix size \(O(m)\) and polynomial degree \(O(m)\). Each is
computable by polynomial arithmetic over the fixed field. Enumerating
all typed entries and computing coefficient rank factorizations costs
\(\exp(O(m))\) bit operations. The number of coefficients is polynomial
in \(m\), and the label count satisfies

~~~math
\log D_m=m\sum_a H(k_a/m)+O(n\log(m+1)).               \tag{16}
~~~

The two-cut coefficient transfer uses the repeated one-particle orders
from Section 3. Its geometric coefficient weights cancel between zero
and infinity, leaving the standard restored pressure. A strict rational
witness to (15), or to the corresponding independent-gauge criterion,
therefore supplies effective crossing pieces with the same rate parameter
\(\rho\). Rounding the type counts changes their logarithmic dimension
by \(O(\log m)\). For an arbitrary binary label width \(\lambda\),
choose a type with at least \(2^\lambda\) labels and restrict to the first
\(2^\lambda\) labels in a
fixed enumeration. The strict pressure margin absorbs the rounding and
padding costs.

The factoring implication uses two identified dependencies. The local
[scale-one transfer, Sections 1--3](scale-one-modular-transfer.md) constructs
the finite-field clique evaluator from effective crossing pieces. Scalar
recovery over each fixed prime field is proved in
[PR #31 at revision b1ba452](https://github.com/1337isnot1337/integer-factoring/blob/b1ba452fa4a3d19fe6fea3822c53613f3b582104/research/fixed-prime-factor-transfer.md).
That source is external to this checkout; it must not be confused with
an existing local file. The [transfer source audit](fixed-prime-transfer-source-audit.md)
records the precise theorem, its multiplier and grouping dependencies,
and the replayed controls. The construction theorem in Sections 1--5
does not depend on that recovery theorem.

This completes the construction and sign interfaces in the stated
one-species family. It does not extend the signed residual-code claims
to arbitrary non-Gaussian local tensor powers. No witness with
\(\rho<1/\sqrt2\) is supplied here.

## 7. Exact reproduction

```sh
python3 tools/check_all_field_pascal_crossings.py \
  --output evidence/all-field-pascal-crossings.json
```

The standard-library checker works over \(\mathbb F_2,\mathbb F_3,
\mathbb F_5,\mathbb F_7\). It verifies 1,348 complete global typed
entries, comparing literal signed Pfaffians with determinants and checking
the ordinary crossing endpoint. It also checks unequal matching types,
including the nontrivial endpoint scalar (6a), over all four fields.
It computes both full Gaussian cut
profiles directly for eight small cases and compares them with the
predicted exterior subset sums. The normalized Pascal orders are checked
in 28 field/mode cases, with independent block-Toeplitz ranks for the
smallest matrices. The negative control in Section 2 verifies why the
global convention is necessary.
