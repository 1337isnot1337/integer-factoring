# An all-degree small-alphabet Smith pressure floor

## Result and scope

Let a polynomial local crossing have physical alphabet size q>=2 over any
field. Arbitrary fixed auxiliary dimension is allowed, with both auxiliary
indices placed on one physical side of each cut. Suppose both cut matrices
have full generic rank n=q^2, as follows from a crossing endpoint with a
valid scalar boundary contraction. Then its **bare, ungauged** Smith
pressure, at every scale s>=1 and every polynomial degree, satisfies

```math
B(s)\ge F_n:=\frac{4n^2}{(n+1)^{(n+1)/n}}.
```

For q=2,3,4 this implies

```math
\rho_{\rm bare}=\frac{\log B(s)}{4\log q}
>\frac{577}{816}>\frac1{\sqrt2}.
```

Consequently no polynomial local path on these alphabets can improve the
factoring exponent through this particular bare-profile pressure criterion,
regardless of degree or fixed auxiliary dimension. There is no fixed
constant-term or parity restriction. Removing a common monomial factor
does not change the pressure.

This is a lower bound on a sufficient rank certificate, not a lower bound
on actual coefficient ranks or factoring time. It does not cover independent
physical or auxiliary gauges in the two cuts, cancellation-module
refinements, arbitrary recombination of coefficients, or larger alphabets.
A block of three binary positions already has q=8 and is outside the
small-alphabet exclusion. The formula F_n remains valid for larger q,
but then need not exceed the half-rate threshold.

The argument generalizes the preceding [finite polynomial catalogues](extended-binary-polynomial-catalogues.md).
The polynomial/Smith interface and degree-one obstructions were developed
in [PR #30](https://github.com/1337isnot1337/integer-factoring/pull/30) and
[PR #39](https://github.com/1337isnot1337/integer-factoring/pull/39).
This all-degree proof used substantial OpenAI Codex assistance. It has no
new Lean formalization and does not change the proved factoring score 1/2.
PR #39 also contains a [stronger scalar binary theorem with independent
physical heights](https://github.com/1337isnot1337/integer-factoring/blob/research/coupled-frontier-c480f802/research/universal-binary-smith-floor.md).
The simple ungauged argument here covers q=3,4 as well, and fixed auxiliary
matrices through rectangular minors. It does not replace that gauged theorem.

## 1. Normalize the polynomial and its two cuts

Divide every entry by the largest common power z^h. Let M be the actual
maximum degree of the resulting polynomial matrix tensor. Its constant
tensor is nonzero. In either cut, transpose if necessary so the rectangular
matrix has n columns, while both auxiliary indices stay on the row side.
The rows may number n times the square of the auxiliary dimension.

Write the zero Smith orders as nu_i. At infinity use the regularization
t^M K(t^-1), with orders mu_i, and put d_i=M-mu_i. The pressure is

```math
U(s)=\sum_{i=1}^n s^{d_i},\qquad
V(s)=\sum_{i=1}^n s^{-\nu_i},\qquad B(s)=U(s)V(s).
```

The cuts may be interchanged, provided the convention is fixed. For an
original common factor z^h, the zero orders increase by h and the actual
degree increases by h, while the regularized infinity matrix is unchanged.
The two pressure factors acquire opposite scalar powers, which cancel.

If the endpoint is K(1)=E times the auxiliary identity, contraction against
vectors with pairing one gives the n-by-n permutation flattening of E.
More generally it is enough for that scalar contraction to equal E.
It is a linear map on the auxiliary side and cannot increase rank. Hence
both rectangular cut matrices have full generic column rank n.

## 2. Four degree and valuation facts

The normalized matrices have the same polynomial entries in different
arrangements. In particular, they have the same maximum degree M and the
same nonzero constant tensor. Determinantal-divisor identities give

```math
\max_i d_i=M,\qquad \sum_i d_i\ge0,\qquad
\min_i\nu_i=0,\qquad \sum_i\nu_i\le nM.                 \tag{1}
```

For the first identity, some infinity entry has nonzero constant
coefficient because some original entry has degree M. Thus the smallest
mu is zero. For the second, reversing a full n-by-n minor replaces its
infinity valuation by nM minus its original degree. Therefore
sum_i d_i is the maximum degree of a nonzero full minor, which is
nonnegative. This argument works for rectangular matrices: the Smith
determinantal divisor uses the minimum valuation among all full minors.

The third identity follows from a nonzero constant entry. Finally, a
nonzero full minor has degree at most nM. The minimum valuation among
such minors is at most that degree; this proves the last inequality.
No relation between the two cuts beyond their common entry degrees is used.

## 3. Arithmetic-geometric mean bound

Remove one maximal term M from the list d_i. The remaining n-1 numbers
sum to at least -M. Applying the arithmetic-geometric mean inequality
to their positive s-powers gives

```math
U(s)\ge s^M+(n-1)s^{-M/(n-1)}.
```

Remove one zero from the list nu_i. Their remaining sum is at most nM,
so the same inequality, with negative exponents, gives

```math
V(s)\ge1+(n-1)s^{-nM/(n-1)}.
```

Set x=s^(nM/(n-1))>=1. Multiplication gives the scalar bound

```math
B(s)\ge\frac{(x+n-1)^2}{x^{(n+1)/n}}.
```

The logarithmic derivative of this last expression is

```math
\frac{(n-1)(x-(n+1))}{nx(x+n-1)}.
```

Thus its minimum on x>=1 occurs at x=n+1 and equals F_n. This proves
the theorem, including M=0, where x=1. The argument has no degree cutoff
and no assumption about the characteristic or auxiliary dimension.

## 4. Exact comparison with the half-rate threshold

Approximate values are provided only for readability:

| q | F_(q^2) | log(F_(q^2))/(4 log q) |
| --- | --- | --- |
| 2 | 8.5598759037 | 0.7743974703 |
| 3 | 25.0861433189 | 0.7332695228 |
| 4 | 50.4600849912 | 0.7071338414 |
| 5 | 84.4049130427 | 0.6890022872 |

The q=4 exclusion is close, so numerical rounding is insufficient.
Use a=577, b=816. The integer inequality 2a^2>b^2 gives
a/b>1/sqrt(2). For each q=2,3,4, the checker verifies exactly

```math
(4n^2)^{nb}>(n+1)^{(n+1)b}\,q^{4an},\qquad n=q^2.
```

Taking positive roots proves F_n>q^(4a/b), so the pressure certificate
has rho>a/b. All arithmetic is on integers; no logarithm approximation
enters the exclusion. The q=5 row illustrates the limit of this proof:
the lower bound there lies below the critical threshold and supplies
neither an obstruction nor an improved construction.

## 5. A stronger all-degree scalar binary parity theorem

When the physical labels are bits, the tensor is scalar, K(0)=I,
K(1)=E, and every odd-parity entry vanishes, the sharper bound

```math
B(s)\ge256/27
```

holds at every degree over every field. In characteristic two the
elementary pencil attains it. This closes all higher-degree parity paths,
not just the finite cubic and quartic catalogues.

Use entry order a,b,c,d,e,f,g,h from the extended-catalogue note. The
row cut has two two-by-two blocks, with maximum entry degrees u and v.
Their determinants are nonzero polynomials, because the endpoint blocks
are invertible. If their degrees are delta_1,delta_2>=0, the restored
infinity degrees are (u,delta_1-u) and (v,delta_2-v). Thus, with
M=max(u,v),

```math
U(s)\ge s^u+s^{-u}+s^v+s^{-v}\ge s^M+s^{-M}+2.
```

At zero the first column block has a unit entry, while the second has
minimum entry valuation m. If the determinant valuations are r and t,
the four zero orders are (0,r,m,t-m). Each determinant has degree at
most 2M, so r,t<=2M. Arithmetic-geometric mean gives

```math
V(s)\ge1+s^{-2M}+2s^{-M}=(1+s^{-M})^2.
```

For y=s^-M, multiplying yields B(s)>=(1+y)^4/y, whose minimum on
0<y<=1 is 256/27 at y=1/3. Neither a Gaussian identity nor a bound on
the number of intermediate polynomial coefficients is needed.

## 6. Checks and the next constructive target

Run `python3 tools/check_all_degree_smith_floor.py`. It verifies the exact
small-alphabet comparisons and the constraints (1) on every recorded
profile in the preceding catalogues, including all 3,999 quartic profiles.
The ordinary proof above supplies the unbounded-degree conclusion;
the finite profiles are implementation controls.

Larger local alphabets, independent gauges, and actual coefficient-rank
or cancellation improvements remain meaningful constructive targets.
This theorem explains why simply lengthening a small local polynomial
path while retaining its bare ungauged profile cannot improve the score.
