# An explicit 220-dimensional constant-profile crossing family

The two-site tensor power of the four-state Schur Gaussian crossing has an
explicit affine family of $2^{220}$ polynomial crossings over
$\mathbb F_2$. Every member has exactly the same two Smith profiles at the
prescribed tensor-product gauges.

The family contains 40 independent directions beyond products of one-site
endpoint-vanishing sections. These directions are valid new crossings.
Their normalized residue matrices share a strict triangular form, which
prevents a pressure gain at these gauges.

This is a statement about an explicit family and two fixed gauges. It does
not classify every tensorized section, exclude a gain after changing
physical gauges, or bound all polynomial crossing paths. It gives no
improvement to the factoring exponent.

## 1. The one-site lattices

Let $G(z)$ be the four-state polynomial crossing in
[the Schur Gaussian reduction](schur-gaussian-digit-reduction.md). Its physical
indices are $(x,y,u,v)\in\{0,1,2,3\}^4$ and

$$
G(1)_{xyuv}=[x=v][y=u].
$$

Use the gauges from
[polynomial crossing interpolation](polynomial-crossing-interpolation.md).
The four rows of each displayed height array correspond to $x,y,u,v$:

$$
h_0=
\begin{pmatrix}
0&0&0&0\\
0&0&0&0\\
0&1&1&2\\
0&-1&-1&-2
\end{pmatrix},
\qquad
h_\infty=
\begin{pmatrix}
0&1&1&2\\
0&1&1&2\\
0&0&0&0\\
0&2&2&4
\end{pmatrix}.
$$

At zero flatten as $(x,y)\mid(u,v)$ and substitute $z=\tau$.
At infinity flatten as $(x,u)\mid(y,v)$ and substitute $z=\tau^{-1}$.
The gauge at either endpoint multiplies an entry by
$\tau^{h_{\mathrm{end}}(x,y,u,v)}$, where the four leg heights are added.

Both matrices have Smith orders

$$
\sigma=(0,1,1,2,3,3,4,4,4,4,5,5,6,7,7,8).
\tag{1}
$$

Write $D=\operatorname{diag}(\tau^{\sigma_i})$. Fix regular invertible
left and right changes $U_e,V_e$ at each endpoint. Their finite jets
suffice to verify the row divisibility statements below. The checker
constructs 32 coefficients of these changes and checks
$U_eG_eV_e\equiv D\pmod{\tau^{32}}$, together with invertibility of their
constant terms.

For an integer $t$, define $S_t$ to be the Laurent polynomial tensors $T(z)$ such
that both matrices

$$
D^{-1}U_0T_0V_0,
\qquad
D^{-1}U_\infty\bigl(\tau^tT_\infty\bigr)V_\infty
\tag{2}
$$

are regular at $\tau=0$. Here $T_e$ includes the stated physical gauge and
flattening. The right changes preserve the row ideals, so only the left
changes are needed to solve for these spaces.

Because all $\sigma_i$ are nonnegative, membership forces every entry of
$T_0$ and $\tau^tT_\infty$ to be regular. Consequently an entry with
physical tuple $\lambda$ has precisely the finite exponent window

$$
-h_0(\lambda)\le k\le h_\infty(\lambda)+t.
\tag{3}
$$

Imposing the remaining row divisibility conditions is therefore an exact
finite binary linear system. The checker finds:

| Space | Monomial variables | Constraint rank | Dimension |
| --- | ---: | ---: | ---: |
| $S_{-2}$ | 769 | 768 | 1 |
| $S_{-1}$ | 1024 | 1018 | 6 |
| $S_0$ | 1280 | 1259 | 21 |
| $S_1$ | 1536 | 1480 | 56 |

The subspace $S_0+zS_0$ has dimension 36 in $S_1$. Choose the 20
complement representatives returned by the checker's deterministic binary
elimination.

The unique nonzero element of $S_{-2}$ is

$$
a(z)=z^2 E^{(v)}_{3,0}G(z),
\tag{4}
$$

where the matrix unit acts on the physical $v$ leg. Formula (4) describes
the section; the construction also obtains it independently as the
one-dimensional kernel of the system.

## 2. The two-site perturbations

Combine each pair of four-state labels as $a+4b$, obtaining physical
alphabet size 16. Let

$$
G^{(2)}(z)=G(z)\otimes G(z),
\qquad q(z)=1+z.
$$

First form the 252 generators

$$
q\,(S_{-1}\otimes S_0)
\quad\text{and}\quad
q\,(S_0\otimes S_{-1}),
\tag{5}
$$

using every pair of basis elements. Their span has dimension 180.
For the 20 selected representatives $b_j\in S_1$, adjoin

$$
q\,a\otimes b_j,
\qquad q\,b_j\otimes a.
\tag{6}
$$

The 292 resulting generators have rank 220. All of their Laurent
expressions simplify to polynomials. Each vanishes at $z=1$.
Thus the affine space

$$
G^{(2)}+\operatorname{span}_{\mathbb F_2}\{\Delta_1,\ldots,\Delta_{292}\}
\tag{7}
$$

contains exactly $2^{220}$ distinct polynomial crossings.

Every product in (5) or (6) uses section twists $t+u=-1$. At zero, $q$
is a unit and both section matrices satisfy their row conditions. At
infinity,

$$
q(\tau^{-1})\tau^{-t}\tau^{-u}
=(1+\tau)\tau^{-1-t-u}=1+\tau.
\tag{8}
$$

It follows that every generator preserves the two tensor-product row
lattices. After the corresponding regular changes, its matrix is divisible
row by row by $D\otimes D$.

## 3. A common triangular residue certificate

For a section $T\in S_t$, let $R_e(T)$ be the constant term of the
corresponding normalized matrix in (2). For a generator
$q\,T\otimes T'$ with twists summing to $-1$, (8) shows that the
normalized constant term is

$$
R_e(T)\otimes R_e(T').
\tag{9}
$$

The checker constructs all 292 such $256\times256$ binary matrices at
each endpoint. They admit a common strict triangular permutation. The
successive support layers have sizes

$$
1,\ 8,\ 28,\ 56,\ 70,\ 56,\ 28,\ 8,\ 1.
\tag{10}
$$

The evidence file gives every index in these layers. The checker verifies
that every nonzero matrix entry goes from a later layer to an earlier
layer. This alone makes every linear combination nilpotent.

A stronger independent calculation iterates the common image operation

$$
W_0=\mathbb F_2^{256},
\qquad
W_{j+1}=\sum_i W_jR_e(\Delta_i),
\tag{11}
$$

using row vectors. At both endpoints the dimensions are

$$
256,\ 240,\ 176,\ 80,\ 16,\ 0.
\tag{12}
$$

Every product of five residue matrices is therefore zero. In particular,
if $N$ is the residue of any perturbation in (7), then $N^5=0$ and

$$
(I+N)^{-1}=I+N+N^2+N^3+N^4.
\tag{13}
$$

After division by the prescribed row factors, the baseline has constant
matrix $I$. Every perturbed matrix consequently has invertible constant
term $I+N$. Its normalized matrix is an invertible power-series matrix.
Both Smith profiles are exactly

$$
\operatorname{sort}\{\sigma_i+\sigma_j:1\le i,j\le16\}.
\tag{14}
$$

The multiplicities are:

| Order | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Multiplicity | 1 | 4 | 6 | 8 | 17 | 24 | 22 | 28 | 36 |

The multiplicities for orders $9,\ldots,16$ are the reverse of those for
orders $7,\ldots,0$.

The 40 new directions alone have a smaller nilpotence certificate.
At each endpoint their residues have the form
$A\otimes B+C\otimes A$, where $A^2=0$ and
$AB_jA=0$ for each of the 20 selected sections. Every such combination
has cube zero. The checker verifies these local $16\times16$ identities
as an additional control.

## 4. Reproduction and scope

Run:

~~~sh
python tools/check_tensorized_smith_sections.py \
  --output evidence/tensorized-smith-sections.json
~~~

The [checker](../tools/check_tensorized_smith_sections.py) reconstructs
the seed, its Smith changes, all four section spaces, the quotient
representatives, and all tensor generators. It checks every endpoint
evaluation, the exact polynomial span ranks, all row-lattice conditions,
the common triangular permutation, and the common image dimensions.
The [evidence](../evidence/tensorized-smith-sections.json) records the
construction hash and source hashes.

No random search or floating-point computation is used in this
certificate. The finite linear calculations prove the claims for every
member of the displayed affine family.

The rank calculation shows that the 40 directions in (6) enlarge the
inherited span. It does not assert that the resulting 220-dimensional
space is the entire endpoint-vanishing section space of the tensorized
lattices. The profile obstruction also leaves changes of physical gauges,
changes of row lattices, and other polynomial crossing families open.
