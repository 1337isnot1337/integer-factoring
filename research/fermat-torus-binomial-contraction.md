# A division-free contraction for Fermat-torus dyadic boxes

## Result and scope

One frequency star of the Fermat transform can be contracted without
enumerating its exponentially many frequency values. A product of
Frobenius binomials on the norm-one torus reduces to a polynomial number
of binary-carry coefficient queries. The resulting arithmetic circuit
uses only addition and multiplication and works over every commutative
coefficient algebra. It permits zero coefficients and zero divisors.

This directly covers products of dyadic geometric input boxes. It does
not prove that the result of eliminating one frequency has the same
product representation in the next frequency. That remaining obligation
is stated in Section 5.

## 1. The torus-binomial theorem

Fix `m>=1`, put `M=2^m+1`, and let an explicit field `E` of characteristic
two contain a primitive `M`th root of unity `zeta`. Define
`T={zeta^j:0<=j<M}`, the `M` roots in this coefficient field. For example,
`F_(2^(2m))` contains such a root. Let `A` be any commutative algebra
over that field. Consider

```math
F(V)=cV^h\prod_{i=1}^r(a_i+b_iV^{2^{e_i}}),
\qquad a_i,b_i,c\in A,\quad e_i\ge0,\quad h\in\mathbb Z.    \tag{1}
```

There is an explicitly constructible division-free arithmetic circuit
of size

```math
O\bigl((r+1)^3(m+\log(r+2))\bigr)                         \tag{2}
```

which computes

```math
\sum_{v\in T}F(v).                                      \tag{3}
```

The sum uses these `M` scalar roots, even if the coefficient algebra
has additional roots of the same equation. Its topology depends only
on `m,h,e_i`, not on the coefficients.
Exponent reduction and circuit generation have polynomial bit cost
in those integer input lengths and the displayed circuit size.
The caller separately charges the arithmetic and representation cost
of the coefficient algebra.

An explicit root field is part of the setup, rather than a free input
oracle. For the example field, exhaustive irreducible-polynomial and
root searches have the coarse deterministic bound `2^O(m)`, including
trial division of the known integer `M`. Once the field is represented,
the contraction itself uses the polynomial number of operations in
(2). It does not repeat that field construction for each assignment.

## 2. Reduce the star to a short list of coefficients

On the torus,

```math
v^{2^m}=v^{-1},\qquad v^{2^{2m}}=v.                      \tag{4}
```

Reduce each `e_i` modulo `2m`. Every nonconstant factor then has
exponent `2^j` or `-2^j` for some `0<=j<m`. Rewrite a negative
factor as

```math
a+bV^{-2^j}=V^{-2^j}(b+aV^{2^j}).                       \tag{5}
```

This introduces no coefficient inversion. In particular, (5) remains
valid if either coefficient is zero or is a zero divisor.

After collecting these shifts, (1) becomes

```math
cV^H P(V),\qquad
P(V)=\prod_i(A_i+B_iV^{2^{j_i}}),\qquad
L=\sum_i2^{j_i}\le r2^{m-1},                             \tag{6}
```

with `H` reduced modulo `M`. The formal degree bound is enough even
when leading coefficients vanish.

Character orthogonality gives

```math
\sum_{v\in T}v^k=
\begin{cases}1,&M\mid k,\\0,&M\nmid k.\end{cases}          \tag{7}
```

The first case equals one because the odd integer `M` is one in
characteristic two. In the second case, the geometric sum vanishes
in the root field, hence in every coefficient algebra over it.
Therefore (3) equals

```math
c\sum_{\substack{0\le K\le L\\K+H\equiv0\pmod M}}
                [X^K]P(X).                              \tag{8}
```

There are at most `r+1` queried coefficients. The bit length of every
query is `O(m+log(r+2))`. Neither the range of `K` nor the torus is
enumerated.

## 3. A fixed binary-carry circuit

At bit position `j`, form the polynomial

```math
E_j(Y)=\prod_{i:j_i=j}(A_i+B_iY).                        \tag{9}
```

Its coefficient indexed by `s` is the sum of the weights of choices
contributing `s` copies of `2^j`. The empty product is one.

For a target `K`, a state is an incoming carry `q` in `[0,r]`.
For each coefficient `[Y^s]E_j`, keep exactly the transitions

```math
q+s\equiv\operatorname{bit}_j(K)\pmod2,\qquad
q'=\frac{q+s-\operatorname{bit}_j(K)}2.                   \tag{10}
```

Start with state zero and weight one. Process every exponent position
and enough additional positions to clear the formal degree bound `L`.
At the end, keep only carry zero.

Every choice of constant or nonconstant term in the original product
has one ordinary binary carry sequence. It is retained precisely when
its total degree is `K`, with its original product of coefficients.
This proves the coefficient formula over any commutative ring.

There are at most `r+1` states, and every local polynomial has degree
at most `r`. A dense, fixed transition schedule takes
`O((r+1)^2(m+log(r+2)))` additions and multiplications per query.
Together with (8), this proves (2). The schedule does not branch on
whether a coefficient vanishes. It can therefore be unrolled as a
division-free arithmetic circuit.

## 4. Exact dyadic input boxes

In any commutative algebra of characteristic two,

```math
S_{2^w}(Z)=\sum_{d=0}^{2^w-1}Z^d
          =\prod_{j=0}^{w-1}(1+Z^{2^j}).                 \tag{11}
```

Consequently

```math
S_{2^w}(\gamma V^{2^e})
 =\prod_{j=0}^{w-1}
       (1+\gamma^{2^j}V^{2^{e+j}}).                     \tag{12}
```

For a product of dyadic boxes, take `r` to be the sum of their binary
widths. Repeated squaring constructs all coefficient powers in
polynomially many algebra operations. The theorem then contracts the
whole shared-frequency star jointly.

There is no exceptional ratio value in this construction. For example,
`S_(2^w)(1)=0` when `w>=1`, and (11) produces that value directly.
The formula also handles `gamma=0`, width zero, and coefficient
algebras with nilpotents.

An interval `[a,a+2^w)` with `a>=0` contributes an extra scalar `gamma^a` and
monomial phase `V^(a2^e)`. Those fit (1), with repeated squaring and
integer modular exponent reduction charged explicitly. Arbitrary
products of many interval decompositions are not replaced by a free
single term.

The implementation is
[fermat_torus_binomial_contraction.py](../tools/fermat_torus_binomial_contraction.py).
It exposes both the general binomial circuit and the dyadic-box
specialization.

## 5. What a second elimination still needs

For fixed values of every other network variable, the Fermat transform's
power-of-two twists have exactly the Frobenius form used here. The
theorem removes the exponential enumeration of that one frequency.

The coefficient circuit in (8) is generally a sum of many products in
the remaining variables. It need not remain a single product of
Frobenius binomials. A repeated contraction must bound the size and
structure of that new representation or supply another charged
algorithm for it.

Likewise, if a remaining multiplication-gate label controls an exponent
in (1), fixing that label for one query does not remove its full physical
domain. The fixed-exponent arithmetic circuit theorem does not supply
a factorization of the jointly labeled gate.

Thus this is an exact local contraction primitive, complementary to the
[rational-star](frobenius-rational-trace-contraction.md) and
[relative Frobenius](relative-frobenius-trace-contraction.md) results.
It does not by itself improve the full factoring exponent.

## 6. Exact controls

The [checker](../tools/check_fermat_torus_binomials.py) compares complete
torus sums with direct finite-field evaluation, including composite
torus orders. It separately evaluates the original geometric sums in
the dyadic-box controls. It tests the circuit over
`E[epsilon]/(epsilon^2)` with genuine zero divisors.

Large-exponent controls use `m=20,37,64,127` and an independent complete
monomial expansion with a small number of factors. This verifies the
signed exponent reduction and modular coefficient selection without
enumerating those large tori. A final control confirms that the
arithmetic gate counts are independent of coefficient values.

The [record](../evidence/fermat-torus-binomial-contraction.json) contains
the exact counts and source hashes.

```sh
python3 tools/check_fermat_torus_binomials.py \
  --output evidence/fermat-torus-binomial-contraction.json
```
