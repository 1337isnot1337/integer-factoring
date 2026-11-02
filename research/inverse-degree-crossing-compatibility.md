# Inverse degree is an additional crossing compatibility constraint

## Result and scope

The abstract seven-letter optimizer of the
[shared-minor relaxation](shared-minor-smith-floor.md) is **not realizable**
by a scalar polynomial crossing, even after arbitrary integer ramification.
Its degree and profiles satisfy the determinant and shared-minor budgets,
but fail an inverse-degree compatibility condition.

Here is the general ordinary theorem, over every field. Let the scalar
polynomial crossing K(z) have actual degree M and endpoint E. Write mu
for its row/infinity orders and nu for its column/zero orders. Assume

```math
\det K_{\rm col}(z)=c z^{nM},\quad c\ne0,\qquad
\mu_1=0<\mu_2,\qquad n=q^2.
```

Then

```math
\boxed{\mu_2\le2(\nu_n-M).}                              \tag{1}
```

Thus the formal degree 12r pattern

```math
\mu=\nu=(0,\underbrace{7r,\ldots,7r}_{6},
                  \underbrace{13r,\ldots,13r}_{42})
```

cannot occur for any positive integer r: (1) would say 7r<=2r.
This rules out a tempting apparent post-half witness before any search.
It does not exclude every seven-letter path or a nearby, asymmetric
profile. No stronger factoring exponent is established.

The nilpotent-resolvent constructions in
[PR #31](https://github.com/1337isnot1337/integer-factoring/pull/31) motivated
the inverse question. The argument here recognizes a constraint directly
from the two profiles, rather than assuming the path was built as a
resolvent. Substantial OpenAI Codex assistance was used. No Lean files change.

## 1. Reverse the column matrix and identify the inverse degree

Put H(w)=w^M K(w^-1), and let P be the permutation matrix that converts
the column-cut convention into the ordinary operator on Mat_q. Thus

```math
H_{\rm col}(w)=T(w)P,\qquad
\Phi(T)(x,y,u,v)=T(E_{vy})_{xu},\qquad T(1)=I.
```

The determinant hypothesis makes H_col, and hence T, unimodular:
their determinants are nonzero constants. In particular T(0) is invertible
and L(w)=T(w)^-1 is a polynomial matrix. Let d be its actual degree.

The largest pole order of K_col(z)^-1 at zero is nu_n. But

```math
K_{\rm col}(z)^{-1}=z^{-M}P^{-1}L(z^{-1}),
```

whose largest entry pole is M+d. Therefore d=nu_n-M.

## 2. A long rank-one jet survives inversion

The hypothesis mu_1=0<mu_2=k means the reshuffled matrix Phi(T(w))
has rank exactly one modulo w^k. More explicitly, a unit pivot and the
vanishing of all two-by-two minors modulo w^k give

```math
T(w)(X)=A(w)X^T B(w)\pmod{w^k}.
```

Because T(0) is invertible as an operator, both A(0) and B(0) are
invertible q-by-q matrices. Hence A and B can be inverted modulo w^k.
The inverse map is again a transpose channel:

```math
L(w)(Y)=B(w)^{-T}Y^T A(w)^{-T}\pmod{w^k}.
```

Thus every two-by-two minor of Phi(L(w)) is divisible by w^k. But those
minors are polynomials of degree at most 2d. If k>2d, each is the zero
polynomial. This contradicts Phi(L(1))=Phi(I)=E, whose row flattening has
rank q^2>1. Therefore k<=2d, proving (1).

This proof uses a scalar local alphabet and a unimodular reversed column
matrix. It does not assert the same inverse argument for rectangular
auxiliary cuts or nonmonomial column determinants.

## 3. Affine inverses recover the single-operator resolvent

If d=1, write L(w)=A+wB. Since L(1)=I, A+B=I. Set S=-A^-1 B.
Then A(I-S)=I and

```math
T(w)=(I-wS)^{-1}(I-S).
```

Unimodularity forces det(I-wS)=1, so S is nilpotent by
Cayley--Hamilton. This is exactly the existing single-operator resolvent
case, obtained as a consequence of the inverse degree.

There is an even stronger first-jet control in this case. If Phi(I-S)
has row rank one, I-S is an invertible transpose channel U. Its square
is an invertible left-right multiplication, so Phi(U^2) has full row
rank n. The first coefficient of Phi(T) is Phi(U-U^2). Quotienting by
the constant row and column spaces removes at most two dimensions from
Phi(U^2). Hence at least n-2 Smith orders of Phi(T) equal one.
This explains why a correct nilpotent column spectrum alone can be
misleading. It agrees with PR #31's partial-transpose power-rank analysis.

## 4. A concrete noncommuting quadratic-inverse constructor

For nilpotent matrices S,T that are strictly upper triangular in one
fixed basis, put N=S+T and define

```math
R(w)=(I-wS-w^2T)^{-1}(I-N).
```

The inverse is a finite polynomial sum because every product of enough
strictly upper triangular matrices vanishes. The endpoint is R(1)=I,
so Phi(R) gives a complete crossing. Moreover

```math
R(w)^{-1}=(I-N)^{-1}(I-wS-w^2T)
```

has degree at most two. No commutativity is required. Thus this constructor
escapes the affine-inverse recognition while retaining a computable exact
endpoint and a bounded inverse degree. It still obeys (1): when its
constant reshuffling has rank one, its second Smith order cannot exceed four.

The [exact implementation](../tools/check_quadratic_inverse_paths.py)
starts on Mat_7(F2) with U(X)=J_7(1)X^T and N=I-U. A deterministic
kernel flag has levels 1 through 16, with dimensions 4,3,...,3. In this
basis split N by input level into S and T. Choose six of levels 2 through
16 to have weight one, and nine to have weight two. There are 5,005 choices.

Every path has actual degree 24. A longest nonzero product has 15 edges,
visits all sixteen kernel levels, and has total weight 6+2*9=24. Every
shorter path skips a positive weight, so cannot cancel that leading
coefficient. The row constant has rank one, and the entire polynomial
resolvent identity and physical endpoint are checked explicitly for the
best retained witness.

This first concrete family does not improve the bound. Quotienting its
first coefficient by the constant row and column spaces gives at least
20 order-one channels, across every one of the 5,005 paths. An independent
two-step Toeplitz rank checks the best witness. The
[record](../evidence/quadratic-inverse-seven-letter-paths.json) includes the
complete multiplicity histogram and deterministic best weights; a
[three-letter control](../evidence/quadratic-inverse-three-letter-control.json)
separately checks the constructor without asserting the same rejection.

## 5. Why two jets suffice to reject this family

Let c>=20 be the number of order-one infinity channels. The polynomial
inverse has degree at most two, so every zero order is at most M+2.
For 1<=s<=3/2,

```math
B(s)\ge \frac{49}{s^2}(1+c/s)\ge\frac{8428}{27}.
```

For s>=3/2, at least one zero channel and one infinity channel have
order zero, so B(s)>=s^24>=(3/2)^24. The checker verifies that both
rational bounds exceed 7^(20/7) by raising to the seventh power. Thus
the complete 5,005-path family has rho>5/7; no full Smith-profile search
is needed for the exclusion.

Run `python3 tools/check_quadratic_inverse_paths.py`. Its construction,
finite screen, endpoint and inverse controls use only exact F2 arithmetic.
Remaining constructive routes include more general noncommuting pencils,
nonmonomial determinant paths, additional physical gauges, and direct
coefficient recombinations. A new candidate must satisfy inverse as well
as determinant and shared-minor constraints. The
[nearby asymmetric profile](seven-letter-asymmetric-profile.md) is one
such feasible profile: it satisfies the inverse-degree bound and crosses
the profile threshold, but it is not constructed algebraically.
