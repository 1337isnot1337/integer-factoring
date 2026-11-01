# Finite-cover degree compression pays physical state cost

A finite algebra can compress a polynomial tensor's degree while retaining
its coefficients. The new state must appear in the crossing interface.
For the Artin--Schreier construction below, the scalar result is a linear
pencil on unequal paired alphabets. The
[unequal-alphabet degree-one theorem](linear-smith-barrier.md) already
excludes its bare balanced target. No large profile computation is needed.

This is an exact construction and cost invoice, not a factoring improvement
or a universal prohibition on finite-cover methods. Weighted physical types,
asymmetric supporting angles, different actual coefficient bounds and
nonlinear outputs are outside the conclusion. The complete score is **1/2**.

## 1. A finite free algebra compresses the entry degree

Fix an odd prime `p` and let

```
B = F_p[w],
A = B[u]/(u^p-u-w).
```

This algebra is free of rank `p` over `B`, with basis `1,u,...,u^(p-1)`.
Its multiplication matrix `J=rho(u)` has constant subdiagonal entries,
and its last column records `u^p=u+w`:

```
J e_b = e_(b+1)                       (0<=b<p-1),
J e_(p-1) = w e_0 + e_1,
J^p-J = w I.                                             (1)
```

For `0<=i<p`, every entry of `J^i` has degree at most one in `w`.
Indeed a product `u^i u^b` has degree at most `2p-2`, so one use of the
relation reduces it to `u^(i+b-p+1)+w u^(i+b-p)` when needed. No second
reduction is necessary. The coefficient of `w` in `J^i e_0` is zero.

Thus every matrix polynomial `F(u)` of degree at most `p-1` becomes a
matrix polynomial `rho(F)(w)` of degree at most one. Entry-degree compression
alone does not supply scalar tensor coefficients at the old alphabet size.
The matrices `rho(a)` act on a real `p`-dimensional state space.

## 2. Make the state and strict endpoint explicit

Let `F` be a `q^2`-by-`q^2` matrix polynomial. Its regular representation
acts on `A tensor V tensor V`, of dimension `p q^2`. Assign the algebra
state to one matched physical pair. The paired alphabets now have sizes

```
q_X=p q,             q_Y=q,             n=q_X q_Y=p q^2.    (2)
```

Suppose `E=rho(F)(1)` is invertible. The fully scalar polynomial matrix

```
C(w)=rho(F)(w) E^(-1)                                    (3)
```

has degree at most one and `C(1)=I_n`. Relabel its columns from `X*Y`
to `Y*X`; equivalently set

```
K(w;x,y,u',v)=C(w)[(x,y),(v,u')].                         (4)
```

Then `K(1;x,y,u',v)=[x=v][y=u']` is a literal rectangular crossing.
The joint inverse in (3) may change its opposite cut. No favorable opposite
profile is assumed. The degree-one theorem applies to the resulting tensor
regardless of that change.

There is a simple sufficient regularity test. If `F(u) G(u)=(u^p-u)^2 I`,
then regular representation gives

```
rho(F)(w) rho(G)(w)=w^2 I.
```

Hence `E` is invertible, with explicitly acquired inverse `rho(G)(1)`.
No root is supplied as advice and no element of a function field is a free
scalar. Fixed `p,q` and the full coefficient tables are seed data; all
`p` states in (2) participate in tensor powers and physical entropy costs.

The leading coefficient in `w` need not be invertible: for source degree
less than `p`, it kills the state `e_0` on every matrix-input coordinate.
One cannot identify that leading fibre with a strict endpoint for free.
Equations (3)--(4) provide the charged regular endpoint instead.

## 3. The bare linear criterion already excludes this proposal

For any scalar degree-one crossing with paired dimensions `(q_X,q_Y)`,
both profiles have `n=q_X q_Y` integral orders and total order at most `n`.
The endpoint-slice lemma gives `r+c>=max(q_X,q_Y)` for their unit ranks.
Discrete convexity therefore proves

```
B(s) >= n max(q_X,q_Y),
rho_bare >= log(n max(q_X,q_Y))/(2 log n) >= 3/4.           (5)
```

For (2), this is

```
rho_bare >= (2 log p+3 log q)/(2 log p+4 log q) > 3/4
```

when `p>1,q>=2`. It exceeds the strict balanced target `1/sqrt(2)`.
The apparent spectral gain cannot be evaluated at the old physical entropy.
No bare thermal or field-size sweep on this representation is warranted.

The same reasoning is a useful gate for any scalar linearization: if the
fully accounted output is a degree-one crossing on paired alphabets of
sizes at least two, apply (5) before a profile search. This statement is
about the displayed Smith bound, not actual coefficient ranks.

## 4. Controls and credit

Run `make verify-linear-smith`. In addition to the existing complete binary
catalogue, the [checker](../tools/check_linear_smith_barrier.py) verifies
72 actual unequal pencils with independent complete block-Toeplitz ranks,
and (1) plus the degree bound for `p=3,5,7,11`. The
[receipt](../evidence/linear-smith-barrier.json) retains the profiles, field
sizes, seed and source hashes. The ordinary proof supplies the universal
claim. The earlier endpoint-slice and convexity arguments retain credit;
this contribution supplies their unequal-pair extension and the finite-cover
invoice. Substantial OpenAI assistant support; no new Lean or independent
mathematical review is claimed.

The [complete replay receipt](../evidence/finite-cover-full-verification.json)
pins a clean successful default gate plus twenty-six focused targets at
`b9ffbafe981eafe9d15cb585e824c184ed42888a`, in 834.33 seconds. This includes
the unequal-pair/companion controls and preceding main/PR63 integration;
later additions are outside the source pin.
