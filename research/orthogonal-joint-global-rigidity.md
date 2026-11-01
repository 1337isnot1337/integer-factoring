# Global orthogonal rigidity for the five-letter joint gauge

## Result and scope

Let `B(z)` be the retained five-letter joint-energy matrix from
[joint energy conjugations](joint-energy-conjugations.md), with physical
half-heights `(0,5,3,1,4)`. Let `C0=B[-8]`.
Over **any field of characteristic two**, suppose a constant matrix `U`:

1. preserves the exact integer joint-energy blocks;
2. satisfies `U^T U=I`;
3. makes the opposite realignment of `U C0 U^T` have rank one.

Then **`U C0 U^T=C0`**. Furthermore, its action on the complete polynomial
crossing has a two-parameter normal form

```math
U(a,b)=I+aP+bQ,\qquad a,b\in K.                 \tag{1}
```

In the original joint index `x+5*y`, `P` has all-one two-by-two blocks
on `(1,5)`, `(2,10)`, `(3,15)` and is zero elsewhere. `Q` has one such
block on `(4,20)`. Thus `P^2=Q^2=PQ=QP=0`.
The parameters in (1) determine distinct tensors.

This is a **global** statement, not an identity-based formal statement.
It does not cover nonorthogonal gauges, variable-dependent gauges,
or gauges that increase the leading realignment rank.
It does not prove a stronger factoring bound. The score remains `1/2`.

## 1. Reduce the leading coefficient to five by five

Reorder the physical labels by heights `(0,1,3,4,5)`. The leading operator
acts on five-by-five matrices as `Z -> X Z^T X`, where

```text
X = [1 0 0 0 0]
    [0 1 1 0 1]
    [0 1 0 0 1]
    [0 0 0 1 0]
    [0 1 1 0 0].
```

It is symmetric, invertible, and its fourth power is the identity.
The energy-zero block of an orthogonal `U` is the scalar one.
Consequently the vacuum anchor of `C=U C0 U^T` is one.
Symmetry of `C` and rank one of its realignment imply uniquely

```math
C_{(x,y),(v,u)}=A_{x,u}A_{v,y},\qquad A_{0,0}=1.
```

To see this directly, factor the realignment into two vectors and use
its anchored row and column. Transposing the joint operator interchanges
these vectors while transposing their physical indices. The anchor fixes
the scalar ambiguity. Thus `C(Z)=A Z^T A`.

Diagonal entries of this congruence are squares. Injectivity of squaring
in a field gives `vec(A)=U vec(X)`; this does not require a perfect field.
Every characteristic-two orthogonal matrix fixes the all-one vector:
its row sums have square one and therefore equal one. Apply this within
each energy block. All components of `vec(X)` are zero or all-one vectors,
except at energies four and six. The general possible matrix is therefore

```text
A = [1 0 0 r 0]
    [0 1 u 0 p]
    [0 v w 0 1]
    [s 0 0 1 0]
    [0 q 1 0 0],
```

with `r+s+u+v=0` and `p+q+w=0`.

## 2. Fourth powers and energy-block invariants force A=X

Since `C^4=I`, write `L=A A^T` and `R=A^T A` to obtain
`L^2 Z R^2=Z` for every matrix `Z`. Thus `L^2=lambda I` and
`R^2=lambda^-1 I`. Also `R` and `L` are similar, so `lambda^2=1`.
In characteristic two this forces `lambda=1`.

The even-height block of `A` is `[[1,r],[s,1]]`.
The off-diagonal entry of `(A A^T)^2` on this block is `(r+s)^3`.
It follows that `r=s`; its diagonal equation then gives `r^4=0`.
Hence `r=s=0`, `v=u`, and `w=p+q`.
The remaining block is

```text
B = [1 u p]
    [u w 1]
    [q 1 0].
```

Because `U` preserves energy, each diagonal energy block of `C` is similar
to the corresponding block of `C0`, even though `C` need not itself preserve
energy. Comparing characteristic polynomials at energy four gives
`u^2+w=1`. At energy six it gives

```math
pq+w^2=1,\qquad w(pq+u)=0.                    \tag{2}
```

The exact characteristic polynomials are

```math
(\lambda^2+1)(\lambda^2+(u^2+w)^2),
\qquad
\lambda^3+(pq+w^2)^2\lambda+w^2(pq+u)^2.
```

Their source values are `(lambda+1)^4` and `lambda^3+lambda`.

If `w=0`, equations (2) and `p+q=w` force `p=q=1` and then `u=1`.
If `w` is nonzero, they imply `pq=u=1+w^2` and `w^3=1`.
For `w=1`, the two choices `(p,q)=(0,1),(1,0)` with `u=0` both fail
`(B B^T)^2=I` by direct multiplication.

For `w!=1`, we have `w^2+w+1=0`, `u=w`, `p=w+q`, and
`q^2+wq+w=0`. The bottom-right entry of `(B B^T)^2-I` is, after
using `w=p+q`,

```math
p^2+q^4+q^2u^2+u^2.
```

Reducing by these equations gives `1+wq`. Its vanishing forces `q=w^2`.
But substituting that value into `q^2+wq+w` gives one, a contradiction.
Only `A=X` remains. This proves the global leading-coefficient claim.

## 3. Quotient the remaining orthogonal centralizer

The exact binary linear equations for energy preservation and commutation
with `C0` define the ten-dimensional algebra `L` in the earlier note.
Every basis element is symmetric, and all basis elements commute.
These facts hold after extension to every characteristic-two field.
Thus orthogonality of `U` in `L` is equivalent to `U^2=I`.

The square-zero subspace `N` of `L` has dimension four. In the published
binary basis of `L`, a basis of `N` has coordinate bitmasks
`6,24,64,384`. All sixteen products between these basis matrices vanish.
The square map has binary coefficients; its kernel extends to every field
because its coordinate equations are squares of binary linear equations.
Therefore every remaining orthogonal unit is `I+R` with `R` in `N`.

Write the basis as `N0,N1,N2,N3`. The commutators with `B[-7]` have rank two,
and the kernel is spanned by `N0+N3` and `N1+N2`.
Both kernel directions commute with **every** coefficient of `B(z)`.
Since `N*N=0`, this two-dimensional subgroup acts trivially, leaving the
unique representatives `P=N0`, `Q=N1` in (1).

## 4. Verification and what remains open

Run:

```sh
python3 tools/check_orthogonal_joint_global_rigidity.py
```

The standard-library-only checker reconstructs the source, verifies all
625 leading-factor entries, the fourth power, the centralizer equations,
all sixteen zero radical products, the two-dimensional stabilizer, and
the ordinary polynomial identities used above. Polynomial powers are
**not** Booleanized. Its receipt is
[`orthogonal-joint-global-rigidity.json`](../evidence/orthogonal-joint-global-rigidity.json).

Run `python3 tools/check_orthogonal_joint_finite_profiles.py` for the
separate complete experiment over F4, F8 and F16. It finds five profile
classes, distinguished by `(a,b)=(0,0)`, exactly one zero, equal nonzero
parameters, or distinct nonzero parameters. All 336 tested pairs are
prefix-dominated by the reference profile. These finite experiments are
**not** an all-field Smith-profile theorem. Proving that dominance over
arbitrary extensions, or escaping the stated orthogonal restrictions,
remains open.
