# Complete Gaussian Smith profiles from a Clifford matrix

## Result and scope

An invertible Gaussian operator on `M` fermionic modes has a complete
`2^M`-dimensional Smith profile that can be recovered from a `2M`-dimensional
matrix. The reduction remains valid when its creation and annihilation
factors have poles. It therefore removes the internal regularity cone in
the earlier [independent Gaussian gauge search](../tools/search_reflection_gaussian_gauges.py).

For the reflection-compressed family, this gives an explicit evaluator for
all `4M` independent physical mode heights at the two places. Every entry
of the smaller matrix has at most two Laurent monomials. The construction
keeps both complete cut profiles and both physical type partition functions.

The reduction is a uniform theorem over characteristic two, proved below
by elementary integral operations. It is a constructive search interface,
not a new factoring exponent. A finite unsuccessful search cannot exclude
the whole gauge family. General occupation-dependent state heights and
non-Gaussian operators do not satisfy this Clifford reduction.

## 1. The Clifford Smith theorem

Let `k` be any field of characteristic two, `O=k[[z]]`, and `K=k((z))`.
Write

\[
\mathcal F=\bigwedge K^M.
\]

For each mode let `C_i` be creation and `A_i` be annihilation in the usual
subset basis. They obey

\[
C_i^2=A_i^2=0,\quad C_iA_j+A_jC_i=\delta_{ij},
\quad C_iC_j=C_jC_i,\quad A_iA_j=A_jA_i.
\tag{1}
\]

Identify the span of these `2M` operators with pairs `(x,y)` through
`sum x_i C_i + sum y_i A_i`. Its square is

\[
q(x,y)\,I,\qquad q(x,y)=\sum_i x_i y_i.
\tag{2}
\]

Let `F` be an invertible operator whose conjugation preserves this span,
and let `B` be the matrix of that conjugation in creation-then-annihilation
order. Equation (2) shows that `B` preserves the quadratic form `q`.
Preservation of only the alternating polar form would not suffice here.

**Theorem.** The Smith orders of `B` have the form

\[
-a_1,\ldots,-a_M,a_M,\ldots,a_1,
\qquad a_1\ge\cdots\ge a_M\ge0.
\tag{3}
\]

If

\[
\bar\nu=2^{-M}v_z(\det F),\qquad
\mu=\bar\nu-\frac12\sum_{i=1}^M a_i,
\tag{4}
\]

then the complete Smith multiset of `F` is

\[
\boxed{\left\{\mu+\sum_{i\in S}a_i:S\subseteq\{1,\ldots,M\}\right\}.}
\tag{5}
\]

In particular, `mu` is an integer. The average `bar nu` can be a
half-integer; no rounding is allowed in (4).

### 1.1 Every integral operation used below has an integral Fock lift

An integral mode change `L in GL_M(O)` acts as `Lambda(L)` on Fock space
and as `diag(L,L^{-T})` on the Clifford span. Its Fock matrix and inverse
both have entries in `O`.

For `i != j`, each of

\[
I+b C_iC_j,\qquad I+b A_iA_j,\qquad I+b C_iA_j\quad(b\in O)
\tag{6}
\]

is its own inverse and is an integral unimodular Fock matrix. Their
conjugations are the corresponding elementary orthogonal shears.
The operator `C_i+A_i` is also its own inverse. Its conjugation exchanges
the `i`th creation and annihilation coordinates and fixes the others.
Thus a particle-hole exchange has an integral lift even in characteristic
two. We work on the full Fock space; preservation of each parity summand
is not required for the reducing operations.

Every primitive isotropic vector `v=(x,y)` over `O` can be sent to the first
creation vector `e_1` by these operations. First exchange a creation and
annihilation coordinate if necessary so that some entry of `x` is a unit.
An integral mode change then sends `x` to `e_1`. Isotropy forces the first
entry of the new `y` to be zero. An alternating lower shear sends the
remaining entries of `y` to zero. Its matrix has entries
`Q_{i1}=Q_{1i}=y_i`, diagonal zero, and is lifted by the corresponding
product of annihilation factors in (6).

### 1.2 Extract one hyperbolic pair

The minimum entry valuation of `B` is `-a` for some `a>=0`: its determinant
is a unit, so all entries cannot have positive valuation. Choose a column
containing an entry of valuation `-a`. This column is isotropic. Permute
input hyperbolic pairs, and exchange its two coordinates if needed, to
make it the first column. Its rescaling by `z^a` is primitive isotropic.
The preceding normalization, applied on the left, makes that column
`z^{-a}e_1`.

Let `f_1` be its hyperbolic partner and `W` the orthogonal complement of
their span. Polar-form preservation forces the `f_1` row to have just one
nonzero entry, namely `z^a` in the `f_1` column. The `e_1` row on `W` is a
linear functional `beta`, all of whose coefficients have valuation at
least `-a`.

Choose `v in W(O)` such that

\[
\langle v,w\rangle=z^a\beta(w)\qquad(w\in W).
\]

The required integral right operation is the Eichler map

\[
e_1\mapsto e_1,\qquad
w\mapsto w+\langle v,w\rangle e_1,\qquad
f_1\mapsto f_1+v+q(v)e_1.
\tag{7}
\]

The last quadratic term is essential. Direct substitution in
`q(alpha e_1+w+gamma f_1)=alpha gamma+q(w)` proves that (7) preserves
`q`, including in characteristic two. It kills the `e_1` row on `W`.
Polar orthogonality now forces the `W` component of the `f_1` column to
vanish, and the isotropy of that column forces its `e_1` component to
vanish. The resulting matrix splits off
`diag(z^{-a},z^a)`.

There is no unproved integral-lifting step in (7). These maps satisfy
`H_v H_w=H_{v+w}`. For a basis vector of `W`, `H_{b e_j}` is lifted by
`I+b C_1C_j`, and `H_{b f_j}` by `I+b C_1A_j`. A product gives the desired
lift for every `v`.

All reducing matrices and their inverses are integral. The remaining
orthogonal block still has entries of valuation at least `-a`. Induction
therefore extracts (3), with nonincreasing `a_i`. The reducing operations
lift to unimodular Fock matrices.

### 1.3 Recover the scalar and the full spectrum

After these lifts, the conjugation matrix is
`diag(z^{-a_i},z^{a_i})`. The exterior operator
`Lambda(diag(z^{-a_i}))` has the same conjugation matrix. Their quotient
commutes with every `C_i` and `A_i`, and hence is a scalar.

To see this last assertion without an irreducibility assumption, the
operator `product_i A_i C_i` projects to the vacuum. Multiplying it on
the left and right by products of creation and annihilation operators
produces all matrix units. Their commutant consists of scalars.

The reduced Fock operator consequently has diagonal orders
`c-sum_{i in S}a_i` for an integer `c`. Replacing `S` by its complement
gives (5). Averaging the orders identifies `mu` with (4), and also proves
its integrality. This completes the theorem.

The argument is an explicit split quadratic Cartan decomposition. General
Cartan decompositions over complete DVRs are proved, for example, in
[Alper--Halpern-Leistner--Heinloth, Theorem 1.3 and Remark 3.5](https://sites.math.washington.edu/~jarod/papers/iwahori.pdf).
The proof above supplies the characteristic-two quadratic and Fock-lift
details needed for this application directly.

## 2. A Gaussian kernel gives its Clifford matrix explicitly

Partition an alternating kernel into left and right variables:

\[
\begin{pmatrix}P&K\\K^T&Q\end{pmatrix},
\]

where `K` is invertible. The Fock matrix whose entries are its selected
Pfaffians factors as

\[
F=s\,U_P\,\Lambda(K)\,V_Q,
\quad
U_P=\prod_{i<j}(I+P_{ij}C_iC_j),\quad
V_Q=\prod_{i<j}(I+Q_{ij}A_iA_j).
\tag{8}
\]

This identity is valid over `K`; it does not require `P` or `Q` to be
integral. It follows by expanding the products: the remaining cross
matchings are the minors of `K`, with determinant and permanent equal
in characteristic two.

Conjugating the generators in (1) gives

\[
B=
\begin{pmatrix}I&P\\0&I\end{pmatrix}
\begin{pmatrix}K&0\\0&K^{-T}\end{pmatrix}
\begin{pmatrix}I&0\\Q&I\end{pmatrix}
=
\begin{pmatrix}
K+PK^{-T}Q&PK^{-T}\\K^{-T}Q&K^{-T}
\end{pmatrix}.
\tag{9}
\]

Both triangular factors preserve the quadratic form because `P` and `Q`
are alternating. Moreover `det U_P=det V_Q=1`, and
`det Lambda(K)=(det K)^{2^{M-1}}`. Thus

\[
\bar\nu=v_z(s)+\frac12v_z(\det K).
\tag{10}
\]

Physical mode gauges multiply `F` on the left and right by exterior
diagonal operators. If their height vectors are `l` and `r`, they change
`B` to

\[
\operatorname{diag}(z^l,z^{-l})\,
B\,\operatorname{diag}(z^r,z^{-r})
\tag{11}
\]

and add `(sum l+sum r)/2` to (10). No sign or regularity restrictions on
these heights are necessary.

Rational heights are covered after ramification. Real weighted profiles
follow by continuity: every determinantal valuation is the minimum of
finitely many affine functions of the row and column heights, so the
finite Smith orders and both sides of (5) vary continuously.

## 3. Closed matrices for the reflection-compressed family

Let `M>=4` be a power of two, `m=M/2`, `F=w^{2M}`, `u=1+F`, and let `J`
be the `m`-dimensional reversal matrix. Define

\[
A_{ij}=(1+w^2)w^{4(i\mathbin{\&}\neg j)},\quad
C_{ij}=(1+w^2)w^{4(i\mathbin{\&}j)},\quad
D_{ij}=(1+w^2)w^{4(m-1-(i\mathbin{|}j))}.
\tag{12}
\]

The bitwise operations use the `log_2 m` low bits. These matrices obey

\[
C=AJ,\quad D=JA,\quad A^T=JAJ,\quad A^2=uI,\quad CD=uI.
\tag{13}
\]

For the last identity for `A`, write its bit factors as
`[[1,1],[w^{4*2^r},1]]`. Their squares are scalar. The product, including
the factor `(1+w^2)^2`, is `(1+w^4)^m=1+w^{2M}`.

For the zero cut `(x,y)|(u,v)`, the cross matrix is
`K_0=[[A,I],[I,A]]`, with inverse `F^{-1}K_0`. Substitution of the full
internal blocks from the original Laurent kernel into (9) gives

\[
B_0=
\begin{pmatrix}
F^{-1}A&F^{-1}I&F^{-1}C&0\\
F^{-1}I&F^{-1}A&0&F^{-1}C\\
D&0&F^{-1}A^T&F^{-1}I\\
0&D&F^{-1}I&F^{-1}A^T
\end{pmatrix}.
\tag{14}
\]

For the infinity cut `(x,u)|(y,v)`, the cross matrix is
`K_infty=[[F^{-1}C,I],[I,D]]`, with inverse
`[[FD,FI],[FI,C]]`. In the variable `w`, its Clifford matrix is

\[
B_\infty=
\begin{pmatrix}
C&FI&A&0\\
FI&FD&0&A^T\\
A^T&0&FD&FI\\
0&A&FI&C
\end{pmatrix}.
\tag{15}
\]

Use `w=z^{-1}` when taking its infinity valuations. The cross determinants
are `det K_0=w^{M^2}` and `det K_infty=w^{-M^2}`. The original common
scalar is `w^{M^2/2}`. If the physical heights at the corresponding place
are `h_i`, the two complete Fock means are therefore

\[
\bar\nu_0=M^2+\frac12\sum_i h_i,
\qquad
\bar\nu_\infty=\frac12\sum_i h_i.
\tag{16}
\]

These means include the common scalar, whose zero and infinity orders
cancel in the combined pressure. Omitting the cross determinant in (16)
would give an incorrect gain.

## 4. Pressure and optimization

Let `beta_1,...,beta_{2M}` be all signed Clifford Smith orders after
multiplying ordinary valuations by the profile scale. Equation (5) gives
the complete cut partition function as

\[
\log Z=-\bar\nu+
\frac12\sum_{j=1}^{2M}
\log\left(e^{\beta_j/2}+e^{-\beta_j/2}\right).
\tag{17}
\]

The implementation evaluates (17) stably. It then uses the same two
physical type partitions and positive Gibbs root as the earlier gauge
engine. Optimizing the actual root avoids mistaking a zero-entropy limit
of absolute pressure for an improved rate.

The derivative of the Clifford term with respect to each signed order is
`tanh(beta_j/2)/4`. Native weighted Smith pivot records therefore give
piecewise exact parameter derivatives. At profile walls these are
one-sided choices; a finite difference control must use a step smaller
than the distance to a wall. The search rescales physical heights by
`M*scale` and optimizes the logarithm of scale to reduce stiffness.

## 5. Verification and limits

The independent checker constructs the original Gaussian kernel and
compares its rational block formula (9) with (14)--(15). It also compares
the subset-sum prediction to complete original polynomial Fock matrices
under integer gauges, including negative internal regularity margins.
Small cases use a separate exact rational DVR elimination, so the central
claim is not checked only by the native series oracle. Toy cases include
odd determinant orders, creation poles, and particle-hole exchanges.

Run

```sh
python3 tools/check_reflection_unrestricted_gaussian.py \
  --output evidence/reflection-unrestricted-gaussian.json
python3 tools/search_reflection_unrestricted_gaussian.py \
  --interval 32 --restarts 20 --seconds 3600 \
  --output evidence/reflection-unrestricted-gaussian-search-M32.json
```

The full theorem is the algebraic argument in Sections 1--3. Finite
controls check the implementation and normalization. The native series
window must resolve the promised full rank; unresolved pivots raise an
error. Floating searches remain exploratory, and a proposed strict
improvement requires a separate exact pressure certificate.

The construction uses the reflection kernel at revision
`c2e951fa1f391edd97ba587d15ec13d3d6c8892d`. The proof, implementation,
and exact controls used substantial OpenAI assistant support. No Lean
source was changed.
