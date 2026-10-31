# A first-jet theorem for arbitrary square-zero physical frames

## Result

For every power-of-two reflection interval `M>=4`, every square-zero
matrix on its full physical alphabet of size `q=2^(M/2)` produces a
polynomial coupled frame with a simple pole before cancellation.
The matrix can have any possible rank over a coefficient field of
characteristic two. The tensor retains the strict crossing endpoint
and the original degree. The four-letter and sixteen-letter examples
retain degrees 24 and 96, respectively.

The general frame argument depends only on four low coefficients.
A matching-parity argument proves those hypotheses for the entire
reflection family. It does not extrapolate from a finite catalogue.
The result supplies a larger constructive search family;
it does not assert an improved pressure or factoring exponent.

## 1. The general local hypothesis

Use the Bell normalization from the
[coupled rational-frame construction](coupled-rational-physical-frames.md).
View the normalized tensor in the third flattening, with rows `(x,v)`
and columns `(y,u)`. It is an operator `M(w)` on `End(F^q)`.

Suppose that, over a field of characteristic two,

\[
 M(w)=I+w^2(\operatorname{ad}_H+cI)+w^4R(w),
 \qquad \operatorname{ad}_H(X)=HX+XH,                 \tag{1}
\]

where `R(w)` is polynomial. This hypothesis says exactly that the
constant operator is identity, the coefficients at degrees one and
three are zero, and the coefficient at degree two has the displayed
form. No restriction is imposed on higher coefficients.

Let `G` be any square-zero matrix: `G^2=0`. Set

\[
 A(w)=I+w^{-1}G.
 \tag{2}
\]

Act with `A` on `x,u` and with `A^(-T)` on `y,v`, as before. Since
the characteristic is two, `A^(-1)=A`. In the third flattening, the
action is

\[
 M\longmapsto LML^{-1},\qquad
 L=A\otimes A^{-T}=\operatorname{Ad}_A.              \tag{3}
\]

Here `Ad_A(X)=AXA^(-1)`. Both `L` and `L^(-1)` have degree at most
two in the scalar `f=w^(-1)`.

## 2. Why every pole cancels

The identity term in (1) is fixed by (3). Conjugation of a commutator
operator is another commutator operator:

\[
 L\operatorname{ad}_H L^{-1}
       =\operatorname{ad}_{AHA^{-1}}.
 \tag{4}
\]

The matrix inside (4) is

\[
 AHA^{-1}=H+f(GH+HG)+f^2GHG.                        \tag{5}
\]

Multiplication by `w^2` therefore removes all possible poles from the
degree-two term in (1). The scalar operator `cI` is fixed. The
remaining expression `w^4 L R(w) L^(-1)` has at most four powers of
`f`, so its factor `w^4` removes all of them. The complete transformed
tensor is polynomial.

This proof works for any rank of `G`, including rank `q/2`. It also
works for every coefficient field of characteristic two satisfying
(1); no enumeration of square-zero matrices is needed.

All introduced powers are nonpositive. Since the physical frame is
identity at infinity, the leading degree coefficient is unchanged.
The strict crossing endpoint is preserved by the paired dual maps at
`w=1`. Every physical matrix has determinant one, so both complete
cost-cut determinants remain unchanged.

The first-jet hypothesis is a condition on the source tensor. It is
not asserted to persist after an arbitrary transformation. Composing
additional frames on a modified tensor still requires either a fresh
first-jet certificate or the complete Laurent polynomiality test.

## 3. The first-jet hypothesis holds at every reflection interval

The [uniform reflection construction](reflection-dual-rail-compression.md)
gives both a selected-minor representation and an explicit Laurent
Gaussian formula. They supply two uniform arguments.

### 3.1. Uniform Bell factorization at zero

At zero, the one-particle matrix in that construction is

\[
 D(0)=\begin{pmatrix}G_0&0\\I&H_0\end{pmatrix},
 \qquad G_0H_0=I.
\]

Every retained row and column subset contains `M/2` indices in each
physical block. The zero upper-right block forces its selected
`M`-minor to be the product of the corresponding `M/2`-minor of
`G_0` and `M/2`-minor of `H_0`.

The matrices `G_0,H_0` preserve the reflection symplectic form. The
selected-minor homology representation proved in Section 2 of the
uniform construction is multiplicative. Its two selected-minor
matrices are therefore inverses. In the literal implementation order,
the selected rows are `(v,y)` and columns are `(u,x)`. Since the
one-particle factors are symmetric, the factor on `(x,y)` is
`Q=rho(H_0)`, while the factor on `(u,v)` is `rho(G_0)=Q^(-1)`.
For the polynomial diagonal blocks, the corresponding factors are
`rho((1+w)H_M)` and `rho((1+w)G_M)`. Thus the constant factorization is

\[
 T(0;x,y,u,v)=Q_{x,y}(Q^{-1})_{u,v}.
\]

This proves that the Bell normalization exists at every interval in
the stated family, independently of the finite coefficient checks.

### 3.2. A coefficientwise identity modulo four

In the explicit Laurent Gaussian kernel, every base exponent is a
multiple of four. This includes the scalar power `M^2/2`, the powers
`+2M,-2M`, and all bitwise exponents. Exactly the four cross-group
blocks `xy,xu,yv,uv` have an additional factor `1+w^2`.

Partition the physical variables into `{x,v}` and `{y,u}`. The four
marked blocks are exactly the edges crossing that partition. In any
perfect matching of a selected variable set, the number of crossing
edges has parity

\[
 r=\bigl(\operatorname{popcount}(x)+
          \operatorname{popcount}(v)\bigr)\pmod2.
\]

Pairs wholly within the first group remove two vertices and do not
change this parity. Since `(1+w^2)^2=1+w^4`, every Pfaffian summand
has the same factor `(1+w^2)^r` times a Laurent polynomial in `w^4`.
Consequently every coefficient tensor entry has the form

\[
 T(w;x,y,u,v)=(1+w^2)^r U_{x,y,u,v}(w^4).
\]

The complete source tensor is polynomial by its selected-minor
construction. The factor `1+w^2` is a unit at zero, so the Laurent
polynomial `U` cannot have a negative exponent. Thus this is an
identity of polynomials.

Let `P` be diagonal physical-letter parity. The identity gives

\[
 T_1=T_3=0,\qquad T_2=(P_x+P_v)T_0.
\]

After applying `Q^(-1)` on `x` and `Q^T` on `v`, the constant tensor
is the Bell product. The two parity operators become `H` and `H^T`,
where `H=Q^(-1) P Q`. In the third flattening this is exactly
`ad_H`. Hence (1), with `c=0`, holds at every power-of-two reflection
interval `M>=4`. Section 2 proves the all-size square-zero pole-one
theorem.

### 3.3. Finite exact controls

In the four-letter and sixteen-letter examples, the normalized
operator additionally has the simple checked form

\[
 H_{s,s}=\operatorname{popcount}(s)\pmod2,
 \qquad H_{s,t}=0\quad(s\ne t).                     \tag{6}
\]

The verifier checks every entry of the coefficients at degrees zero,
one, two, and three. There are 1,024 such entries in the four-letter
case and 262,144 in the sixteen-letter case. The matrices in (6) have
two and eight diagonal ones, respectively.

The general extraction routine also permits a scalar `c`: its value
is the `(0,0;0,0)` coefficient at degree two. It then extracts `H`
with `H_(0,0)=0` and checks the entire proposed adjoint operator.
Changing `H` by a scalar matrix does not change `ad_H`.

The verifier also checks the full polynomial identity above for every
entry of both finite tensors. These controls check the conventions
and implementation. The all-size assertion follows from Sections
3.1--3.2; it does not require the additional simplification `H=P`
at larger intervals.

## 4. Involutions give inexpensive full-rank generators

For any involutive permutation matrix `P`,

\[
 G=I+P\qquad\text{satisfies}\qquad G^2=0.           \tag{7}
\]

Each two-cycle contributes rank one. A perfect matching of the
sixteen letters therefore gives rank eight. Fixed points can be
retained to obtain smaller ranks. At the crossing endpoint,
`A(1)=P`, so the paired maps remain invertible.

The construction also accepts arbitrary square-zero `F4` matrices,
including conjugates of these generators by full physical bases.
If a generator is supplied in the original `x` coordinates, the
implementation converts it to Bell coordinates by `Q^(-1) G Q`.
Restoring the original coordinates afterwards recovers the requested
generator on the `x` leg and its paired dual on `v`.

These full-letter generators are not restricted to one-particle
exterior-power changes. That scope alone is not a pressure claim;
every candidate still requires its complete weighted cost.

## 5. Implementation and independent controls

The API constructs all five auxiliary-scalar coefficient tensors,
retains one common `w^4` shift, and tests the final polynomial before
dividing. Thus the implementation also directly verifies the
polynomiality guaranteed by the theorem. A matrix with `G^2!=0` is
rejected, because its paired inverse would require a different
formula.

Independent dense physical-leg contractions check the resulting full
polynomial coefficients. The verifier also checks the involution,
original/Bell coordinate conversion, recovery of the earlier
rank-one implementation, nonbinary conjugated generators, and
negative controls for both a nonnilpotent generator and a nonzero
linear jet. Selected four-letter complete weighted profiles are
compared with the independent exact rational local oracle.

```sh
python3 tools/check_square_zero_coupled_frames.py \
  --output evidence/square-zero-first-jet-frames.json
```

Files:

- [`square_zero_coupled_reflection_frames.py`](../tools/square_zero_coupled_reflection_frames.py)
- [`check_square_zero_coupled_frames.py`](../tools/check_square_zero_coupled_frames.py)
- [`square-zero-first-jet-frames.json`](../evidence/square-zero-first-jet-frames.json)
