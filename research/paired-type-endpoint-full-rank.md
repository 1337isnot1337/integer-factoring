# Paired-type crossing endpoints force full generic rank

## Proposition

Let `F` be any field. Let the two physical alphabets have fixed sizes
`A,B >= 2`, and let

```math
K(w;x,y,u,v)\in F(w),
\qquad x,v\in[A],\quad y,u\in[B].
```

Fix an endpoint `a` of the projective line over `F`. Local entries of
`K` may have poles at `a`. Suppose the following statement holds for
every positive integer `m` and every pair of physical word types.

Let `X` be the words in `[A]^m` of the first chosen type, and let `Y`
be the words in `[B]^m` of the second chosen type. Restrict
`K(w)^(tensor m)` to `x,v in X` and `y,u in Y`. Assume every entry
of this restricted tensor is regular at `a`, with specialization

```math
\left.\left(K(w)^{\otimes m}\big|_{X\times Y\times Y\times X}\right)
\right|_{w=a}=[x=v][y=u].                                \tag{1}
```

The notation in (1) means specialization **after the restriction and
removal of any removable singularities**. It does not assume that the
individual local tensor `K(a)` is defined.

Then each of the two local cut matrices

```math
(x,y)\mid(u,v),\qquad (x,u)\mid(y,v)
```

has generic rank exactly `AB` over `F(w)`.

### Proof

Fix either cut, and let its local generic rank be `r`. Up to fixed
coordinate permutations, the same cut of the tensor power is the
Kronecker power of the local matrix. Its rank is exactly `r^m`.
Restricting the physical words selects a submatrix, so the generic rank
of the retained matrix is at most `r^m`.

The retained matrix is square, of order `|X||Y|`. Its regular
specialization in (1) is a permutation matrix. Therefore its generic
rank is at least `|X||Y|`: a minor that vanishes identically as a
rational function still vanishes after any regular specialization.
Consequently

```math
r^m\ge |X||Y|.                                           \tag{2}
```

There are at most `(m+1)^A` types of length-`m` words over `[A]`.
Their classes partition all `A^m` words. Choose a largest class, and
do the same for `[B]`. Then

```math
|X|\ge\frac{A^m}{(m+1)^A},\qquad
|Y|\ge\frac{B^m}{(m+1)^B}.
```

Equation (2) now gives

```math
r\ge\frac{AB}{(m+1)^{(A+B)/m}}.
```

Letting `m` tend to infinity gives `r>=AB`. The local matrix is square
of order `AB`, so equality follows. The argument applies to each cut
separately. It uses no assumption on the characteristic of `F`. □

## A weaker hypothesis suffices

The proposition does not require every type if a suitable sequence of
type classes is already known. Suppose the regular crossing endpoint
holds along an unbounded sequence of lengths, for paired word sets
`X_m,Y_m` such that

```math
\frac1m\log|X_m|\longrightarrow H_X,\qquad
\frac1m\log|Y_m|\longrightarrow H_Y.
```

The same proof gives

```math
\operatorname{rank}_{F(w)}K_{\mathrm{cut}}
\ge \exp(H_X+H_Y).                                      \tag{3}
```

Balanced uniform types have `H_X=log A` and `H_Y=log B`, so they
alone force full rank. A proof that only supplies one nonuniform type
with smaller entropy yields the weaker bound (3); it must not claim
full rank from that bound alone.

## Application to removable-pole occupation gates

The [charged quartic construction](reflection-charged-quartic-tiles.md)
starts from the reflection Gaussian and applies an ordered word of
occupation gates. Write `m_0=M/2`, `L=m_0-1`, and `chi=1+w`.
A supported gate has coefficient

```math
g_j(w)=\frac{\gamma_j w^{a_j}}{\chi^{s_j}},\qquad
\Delta_j>L s_j,
```

where `Delta_j` is its paired occupation-charge increase. The charge
proof establishes Laurent polynomiality and the exact crossing endpoint
for every retained pair of full-letter word types. Laurent powers of
`w` are units at the endpoint `w=1`.

The proposition therefore proves full local generic rank for these
rational tensors, including words with genuine local endpoint poles.
A separate nonzero determinant specialization is not a mathematical
hypothesis of their Smith-profile construction. A finite-field
specialization remains useful as an independent implementation check.

The conclusion also survives restriction to any fixed nonempty physical
alphabet on each matched pair, provided each alphabet has at least two
letters. Paired word types on the retained symbols still have total
paired charge zero. The same polynomiality and endpoint proof applies,
and largest type classes now have entropy rates `log A,log B` for the
retained alphabet sizes. Intermediate occupation tuples used to compute
gate coefficients may lie outside those alphabets; this does not change
the argument.

## Scope

This proposition concerns scalar tensors and their complete local cut
matrices. It proves neither favorable Smith orders nor an improved
crossing parameter. It also does not turn an incomplete finite-window
Smith computation into a certificate: every claimed finite profile still
needs all its pivots resolved or another exact verification.

If additional auxiliary indices are placed on both sides of a cut, its
ambient order can exceed `AB`. The same rank argument can then give
`rank >= AB`, but it does not imply full rank in that larger ambient
matrix. An auxiliary word construction must justify its contraction and
rank comparison separately.

The proof is an ordinary algebraic argument, independently reviewed in
this research session. No Lean formalization or new factoring exponent
is asserted.
