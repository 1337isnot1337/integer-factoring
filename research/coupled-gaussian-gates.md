# Growing coupled Gaussian gates with complete valuation costs

## Result and scope

This note extends the complete Clifford profile interface to quadratic gates
that couple different physical letters. Such a gate can change the determinant
of a flattening. The formula below tracks that change and every rational
denominator in the reduced matrix.

The construction is effective and preserves the strict crossing endpoint.
It is a new seed interface, not an improved factoring bound. It supports a
number of coupled modes that grows with M. Numerical searches using this
interface remain separate from the exact identities.

The coefficient field in
[`reflection_gaussian_coupled_gates.py`](../tools/reflection_gaussian_coupled_gates.py)
is F2. The complete Smith theorem used here is proved in
[`clifford-gaussian-smith-profiles.md`](clifford-gaussian-smith-profiles.md).

## 1. Physical gates and the endpoint

There are m=M/2 occupation modes on each physical letter x,y,u,v. Let U and V
each be a creation or annihilation operator on distinct physical modes. Put

\[
 G=1+g(w)UV,\qquad g(w)=w^a+w^b,\quad a\ne b.         \tag{1}
\]

The operators U,V commute and square to zero. Hence G has inverse itself on
the complete tensor space. It is the identity at w=1. Any finite ordered word
of these gates applied to the reflection tensor still satisfies

\[
 T_{xyuv}(1)=\mathbf1_{x=v}\mathbf1_{y=u}.             \tag{2}
\]

A block in the implementation chooses two different physical letters, one
operator kind on each letter, and identity pairs, reversal pairs, all pairs,
or one specified pair. The gates within one such block commute. Different
blocks may use different kinds and need not commute.

Although each gate is invertible on the complete tensor space, it need not
act as independent row and column changes on a particular flattening. Its
cut determinant cost must therefore be computed.

## 2. A rank-two Clifford update

Fix one of the two cuts, with M row modes and M column modes. Let F be its
invertible Fock matrix, and let B be the conjugation matrix on the vector
space spanned by C_1,...,C_M,A_1,...,A_M. Write

\[
 \langle r,s\rangle=rs+sr
\]

for the polar pairing, and write flip(i)=i+M modulo 2M.

Suppose the physical gate has one row operator u and one column operator v.
Column action is transposed, so

\[
 F'=F+g uFv^T=(1+g u w)F,\qquad w=Fv^TF^{-1}.       \tag{3}
\]

Both u and w are isotropic Clifford vectors. Put

\[
 \tau=\langle u,w\rangle,\qquad d=1+g\tau.
\]

For X=uw, one has X^2=tau X and therefore

\[
 (1+gX)^{-1}=1+\frac gdX.                            \tag{4}
\]

Here d is nonzero. Indeed, F is regular and invertible at w=1 by (2), while
g(1)=0, so d(1)=1. For any Clifford vector a, direct multiplication in (4)
gives

\[
 (1+gX)a(1+gX)^{-1}
 =a+g\langle w,a\rangle u
       +\frac gd\langle u,a\rangle w.               \tag{5}
\]

If u is the basis vector e_i and w is column j of B, orthogonality of B turns
the first pairing in (5) into the coordinate with index flip(j). Thus

\[
 \boxed{B'=B+gE_{i,\operatorname{flip}(j)}
       +\frac gd B_{*,j}B_{\operatorname{flip}(i),*},
 \quad d=1+gB_{\operatorname{flip}(i),j}.}            \tag{6}
\]

The transpose convention matters: physical column creation uses an
annihilation index j, and physical column annihilation uses a creation index.

If both physical operators lie on the same side of the cut, they are distinct
orthogonal canonical vectors. Their gate has tau=0 and is the usual left or
right Clifford shear. On the column side the transpose flips both operator
kinds. These same-side gates have determinant one on the cut.

## 3. The determinant mean changes

The size of F is 2^M. If tau=0, then X^2=0 and det(1+gX)=1. If tau is nonzero,
X acts as tau on the image of u, since uwu=tau u. Its image is contained in
that subspace. A canonical creation or annihilation operator has rank
2^(M-1). Consequently X has eigenvalues tau and 0, each with multiplicity
2^(M-1). In both cases,

\[
 \boxed{\det(1+gX)=d^{2^{M-1}}.}                    \tag{7}
\]

At either valuation place s, the mean Fock Smith order changes by

\[
 \boxed{\operatorname{mean}_s(F')-
         \operatorname{mean}_s(F)=\frac12v_s(d).}     \tag{8}
\]

Start from means M^2 at zero and 0 at infinity for the original reflection
tensor. Apply (8) for every cross-cut elementary gate, separately at the two
cuts. Same-side gates add zero. At the end, arbitrary physical mode heights
h add (sum h)/2 at their place.

Let a_1,...,a_M be the positive half of the final Bogoliubov Smith spectrum.
The complete Fock orders are

\[
 \left\{\operatorname{mean}_s(F)
        -\frac12\sum_i a_i+\sum_{i\in S}a_i:
          S\subseteq\{1,\ldots,M\}\right\}.          \tag{9}
\]

All scalar normalizations are included in the mean. In particular, reusing
the original reflection means after a cross-cut gate would generally give
an incorrect pressure.

## 4. Rational reduced matrices do not remove pole costs

Equation (6) can create rational entries even though the original tensor
coefficients remain Laurent polynomials. Write B=N/D with a common nonzero
polynomial denominator D. At a fixed place, the kth determinantal valuation
of B is that of N minus k v_s(D). Hence

\[
 v_i(B)=v_i(N)-v_s(D),\qquad
 v_0(D)=\operatorname{ord}_0D,\quad
 v_\infty(D)=-\deg D.                               \tag{10}
\]

The implementation clears D exactly over F2[w], applies the finite matrix
oracle to N, and subtracts (10) from every signed order and every active
minor valuation. It repeats the oracle at a larger series window when a
promised pivot is unresolved. The maximum supported window is a computational
limit, not a criterion for favorable cancellation.

The original coefficient tensor, its common Laurent clearing, the rational
denominator of B, and the determinant factors d are separate quantities.
Equations (8)--(10) account for each of them.

## 5. Degree and effective construction bounds

In a block with fixed kinds on two different physical letters, each physical
mode can occur at most once in a nonzero product of quadratic factors. Thus
at most m factors occur, including for an all-pairs block. If the block has
powers a_j<b_j, a word of blocks has Laurent clearing at most

\[
 C=m\sum_j\max(0,-a_j).
\]

The original reflection degree is at most 3M^2/2. The cleared coefficient
tensor therefore has degree at most

\[
 \frac{3M^2}{2}
   +m\sum_j\bigl(\max(0,-a_j)+\max(0,b_j)\bigr).     \tag{11}
\]

Multiplication by w^C preserves the endpoint and adds opposite scalar shifts
at zero and infinity. A bounded number of blocks with powers O(M) gives
degree O(M^2). Direct coefficient construction still uses 2^O(M) work.

## 6. Independent controls

The checker applies every gate directly to the original tensor coefficient
array, without using the rank-two formula. It checks every crossing endpoint,
the clearing and degree bounds, both complete Fock profiles under independent
integer mode heights, and the changed determinant means. Small controls also
compare full original Fock conjugation matrices and use exact rational DVR
elimination without a series truncation.

For creation-only gates, a separate calculation adds g to the appropriate
entries of the original alternating kernel and recomputes its Gaussian block
factorization. This independently checks (6) and (8). The controls include
all four creation/annihilation combinations and noncommuting words.

```sh
python tools/check_reflection_gaussian_coupled_gates.py \
  --output evidence/reflection-gaussian-coupled-gates.json
```

The exact identities provide a scalable search interface. A finite numerical
screen, a favorable local optimization status, or a pressure approaching zero
does not by itself prove an improved asymptotic bound.
