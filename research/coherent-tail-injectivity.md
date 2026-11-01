# All-field injectivity of the surviving full-tail band

## Result and scope

The leading-band rank result for constant `F4` maps on two hot modes
has a structural proof over every field of characteristic two. It
holds for every fixed power-of-two number `k` of hot modes.

More precisely, in the `T=t^(2M)` leading band of the full physical
tail construction, **every nonzero joint functional on the two hot
physical groups is injective on the `2^(2k)` leading input channels**.
This holds independently on the row side and the transposed column
side. A single pair of nonzero physical output functionals is already
enough; product structure is not needed for the linear-algebra lemma.

It follows that every fixed, nonzero physical split code on finitely
many modes, tensored with a complete Gaussian tail, retains the full
surviving Gaussian band. For a code of dimension `q_hot < 2^k`, the
near-critical pressure correction from the preceding note therefore
has the positive entropy loss

\[
L=\log(2^k/q_{\rm hot}).
\tag{1}
\]

The result removes the field-size and two-mode restrictions from
[the finite leading-band certificate](coherent-physical-tail-low-band.md).
It does not give an all-scale pressure floor, complete high-band
profiles, or an improved factoring exponent. The limit keeps the
physical support and label gauges fixed as the Gaussian tail grows.

## 1. A general coefficient-extraction lemma

Work in the characteristic-two occupation algebra: the variables
commute and each variable squares to zero. Let `z_1,...,z_d` be
physical variables, `b_1,...,b_d` their partner variables, and
`f_1,...,f_d` independent free variables. Let `G_0(z,f)` be any
polynomial with

\[
G_0(0,f)=1.
\]

Take linear forms

\[
\ell_j=(Bz)_j+(Cf)_j,\qquad C\in\operatorname{GL}_d(K).
\]

Define the encoding

\[
\mathcal E(g)=G_0(z,f)\prod_{i=1}^{d}(1+z_i b_i)\,g(\ell_1,\ldots,\ell_d)
\tag{2}
\]

for any occupation polynomial `g` in `d` input variables. A physical
functional has coefficients `alpha_A` on the monomials `z_A`.

**Lemma.** Every nonzero physical functional `alpha` makes
`(alpha tensor I) E` injective.

Choose a support set `A` of largest cardinality with `alpha_A != 0`.
Extract the coefficient of the partner monomial `b_A` from the
conditional output. The corresponding partner factor in (2) is
exactly `z_A`. Any further physical variable from `G_0` or `g(ell)`
either repeats a variable of `z_A`, making the product zero, or raises
the physical degree above the largest degree seen by `alpha`.
Consequently the extracted coefficient is exactly

\[
[b_A]\bigl(\alpha\otimes I\bigr)\mathcal E(g)=\alpha_A\,g(Cf).
\tag{3}
\]

Substitution by `C` is an automorphism of the occupation algebra.
Equation (3) is therefore an explicit left inverse: extract `b_A`,
divide by `alpha_A`, then substitute by `C^-1`. This proves the
lemma over the whole coefficient field `K`. The functional may have
arbitrary coefficients in that field.

The dimension of the partner space is relevant. It lets a highest
physical-degree coefficient survive without being mixed with the
logical input variables. Merely observing generic rank at a few
field specializations would not establish this statement.

## 2. The hot matrices and their exact identities

Let `K=F(t)`, where `F` is any field of characteristic two. For a
power of two `k`, put

\[
x=t^{2k},\qquad \sigma=(1+x)^{-1},\qquad \lambda=x\sigma,
\]

and define `k`-by-`k` matrices

\[
F_{ab}=(1+t)t^{2(a\mathbin{\&}\neg b)},\qquad
H_{ab}=(1+t)t^{2(a\mathbin{\&}b)},\qquad
J_{ab}=[a+b=k-1].
\tag{4}
\]

The digit factorizations give

\[
FJ=H,\qquad FH=(1+x)J,\qquad HF^{\mathsf T}=(1+x)J.
\tag{5}
\]

In particular, `F,H,J` are invertible over `K`; their determinants
are local units at `t=0`. The scalars `x` and `lambda` are nonzero
in `K`. The injectivity proof uses their invertibility in the
rational-function field; it does not call them local units.

Use the complete physical operator formulas from
[sections 2-3 of the preceding note](coherent-physical-tail-low-band.md).
Set `T=0`. The terminal two `y` blocks disappear from every selected
operator, leaving six blocks of `k` Gaussian variables. Denote them
by

\[
(h,p,q,r,v,s)=(x0,x1,xp,xe,y0,y1).
\tag{6}
\]

The leading Gaussian input is the full occupation algebra on `h,p`;
all other Gaussian modes begin in the vacuum. There are `2^(2k)`
input channels.

## 3. The row-side normal form

The physical row creators, written as vectors of multiplication
operators, are

\[
C_x=h+F^{\mathsf T}(v+xs),\qquad C_y=v.
\]

Choose new creation coordinates

\[
u=h+F^{\mathsf T}(v+xs),\quad
P=p+\sigma F^{\mathsf T}v,\quad Q=q,\quad
R=r+F^{\mathsf T}v,\quad S=s,
\tag{7}
\]

and keep `v` itself. This is an invertible triangular change over
`K`, regular with regular inverse at zero. Direct substitution in
the physical annihilators, using (5), gives

\[
\begin{aligned}
A_x&=\partial_u+Hv+JR,\\
A_y&=\partial_v+Hu+HQ+HR+x(1+x)JS.
\end{aligned}
\tag{8}
\]

For example, the derivative terms
`partial_v+F partial_h+sigma F partial_p+F partial_r` become
`partial_v` in the new coordinates. The two copies of
`HF^T v` in the creation part cancel.

Conjugating by the occupation Gaussian whose quadratic kernel is

\[
u^{\mathsf T}Hv+u^{\mathsf T}JR+
v^{\mathsf T}H(Q+R)+x(1+x)v^{\mathsf T}JS
\tag{9}
\]

turns (8) into the standard derivatives in `u,v`, while leaving
their creators unchanged. Thus, in these physical coordinates,
the leading input is a Gaussian times an occupation polynomial in
the original logical forms

\[
h=u+F^{\mathsf T}v+xF^{\mathsf T}S,\qquad
p=P+\sigma F^{\mathsf T}v.
\tag{10}
\]

Now change only the partner variables:

\[
b_u=JR,\qquad b_v=H(Q+R).
\tag{11}
\]

This is invertible and does not change `P,S`. The Gaussian from
(9) factors as

\[
G_0(u,v,S)\prod_{a<k}(1+u_a b_{u,a})(1+v_a b_{v,a}),
\tag{12}
\]

where `G_0(0,0,S)=1`. The logical forms (10) contain no partner
variables. Setting the physical variables `u,v` to zero leaves

\[
(h,p)=(xF^{\mathsf T}S,P),
\tag{13}
\]

an invertible linear map on the `2k` free variables `P,S`.
The coefficient-extraction lemma with `d=2k` applies. Every
nonzero joint functional on the physical row groups preserves all
`2^(2k)` input channels.

## 4. The transposed column-side normal form

For a right action, transpose the coefficient matrix and exchange
the transported creator and annihilator roles. The resulting
physical creators are again pure multiplication operators:

\[
C_v=h+F^{\mathsf T}(v+\lambda s),\qquad C_u=v.
\]

Use the new coordinates

\[
u=h+F^{\mathsf T}(v+\lambda s),\quad
P=p+F^{\mathsf T}v,\quad Q=q,\quad
R=p+r+F^{\mathsf T}v,\quad S=s,
\tag{14}
\]

and retain `v`. The physical annihilators become

\[
\begin{aligned}
A_v&=\partial_u+Hv+JR,\\
A_u&=\partial_v+Hu+HR+\sigma HQ+xJS.
\end{aligned}
\tag{15}
\]

The Gaussian kernel is now

\[
u^{\mathsf T}Hv+u^{\mathsf T}JR+
v^{\mathsf T}(HR+\sigma HQ)+xv^{\mathsf T}JS.
\tag{16}
\]

Choose partners

\[
b_u=JR,\qquad b_v=HR+\sigma HQ.
\tag{17}
\]

They are independent of the free variables `P,S`. The original
logical forms are

\[
h=u+F^{\mathsf T}v+\lambda F^{\mathsf T}S,\qquad
p=P+F^{\mathsf T}v.
\tag{18}
\]

Their physical-zero part `(lambda F^T S,P)` is invertible.
The same coefficient-extraction lemma therefore proves injectivity
for every nonzero joint functional on the transposed column side.
This proves the opposite side directly; it does not assume equality
of the two physical profiles.

## 5. Consequences for physical restrictions and pressure

For a constant physical split pair `LP=I`, every retraction row
and every embedding column is nonzero. The ambient height
corrections multiply their coordinates by nonzero monomials in
`t`. The new physical label gauges multiply entire rows or columns
by nonzero monomials. They remain nonzero functionals over `K`.

Apply sections 3 and 4 to one pair of physical outputs on each
side. The left leading factor has full column rank `2^(2k)`,
the right leading factor has full row rank `2^(2k)`, and the
central leading Gaussian diagonal is invertible. Hence the
leading active core has rank exactly

\[
2^{2k}
\tag{19}
\]

at both cuts, over every characteristic-two field. This statement
does not require the physical maps to preserve height, occupation,
parity, or deletion invariance.

An arbitrary fixed finite support is contained in an initial
power-of-two block. Extending the code by identities on unused
modes in that block preserves the ratio `q_hot/2^k`. Thus the
result is not restricted to contiguous original supports.

For fixed maps and bounded label gauges, the complete bivariate
core has only finitely many coefficient functions in `t,T`.
The same determinantal argument as in the preceding note makes
its Smith orders eventually affine in `2M`. Exactly `2^(2k)`
orders have slope zero. The Gaussian active block has the same
number of slope-zero orders. The ratio of their counting
functions therefore has logarithm
`O(u)+O(exp(-2Mu))` as `u -> 0` and `Mu -> infinity`.

The physical hot alphabet, however, has `q_hot` letters in place
of `2^k`. Repeating the pressure calculation in the preceding
note gives, for `rho_0=1/sqrt(2)` and bounded `s`,

\[
\begin{aligned}
u&=\frac{\rho_0}{2M}(\log M+s),\\
\rho_M(u)
&=\rho_0+\frac{6\rho_0}{\pi^2M}
 \left[\log(2^k/q_{\rm hot})(\log M+s)+e^{-s}\right]
 +o(M^{-1}).
\end{aligned}
\tag{20}
\]

Thus a proper fixed code has a positive leading correction on
this scale branch. The loss cannot be removed by a larger
constant coefficient field or an unequal-height mixture on the
same fixed physical support.

The proof is not uniform in a growing `k`. It does not exclude
physical spaces whose support grows with the reflection interval,
constraints that couple the code to the tail, or a different
finite scale. Such changes alter hypotheses used in the fixed
partner-space reduction and the asymptotic error bounds.

## 6. Exact controls

The proof is algebraic; the finite controls do not replace its
all-field argument. The [checker](../tools/check_coherent_tail_injectivity.py)
verifies (5) as exact binary polynomial identities for
`k=1,2,4,8`. It independently conjugates both sides of the
existing physical operator interface through (7) and (14),
checking the annihilator normal forms and the logical forms.

The checker also constructs the complete encoded matrices in
the partner coordinates for `k=1,2`. It uses joint functionals
of every possible largest occupation degree, including mixtures
with coefficients outside `F4`. For each one it checks (3)
entry by entry, rather than checking only the resulting rank.

The [retained controls](../evidence/coherent-tail-injectivity-controls.json)
record 255 exact polynomial identity entries, 20,400 coordinate
identity entries, and 2,656 exact partner-coefficient entries in
sixteen joint-functional cases. These complement the original
full-kernel and literal physical tensor controls in the preceding
note.

```sh
python3 tools/check_coherent_tail_injectivity.py \
  --output evidence/coherent-tail-injectivity-controls.json
```
