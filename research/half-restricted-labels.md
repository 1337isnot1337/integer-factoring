# Restricted labels and the intermediate span

This note tests one proposed extension of the Appendix B crossing pieces.
Keep the original, unprojected operators, but encode labels only by
partitions of `h` whose parts are at most `L`. The crossing identity
continues to hold on any such injectively encoded subset. The issue is
the size of the intermediate spans, not the identity.

The source is Claude (Anthropic), *A 2^o(n) Time Algorithm for 3-SAT*,
in `anthropics/formal-math` at
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.
The operator and factor definitions are in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean):
`op`, `angle`, `aFam`, and `PieceFacts.row`.
Their proofs are assembled in
[GeneralPieces/Main.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Main.lean).

## Exact counterexample to the literal product-span bound

Write `p_{<=L}(e)` for the number of partitions of `e` whose parts are
at most `L`. A tempting bound is that the intermediate image has at most
`p_{<=L}(2h)` dimensions, as if it consisted only of products of two
low-mode inputs. This bound is false.

Take `h=L=2`. The two labels are `(1,1)` and `(2)`. At central piece
index `s=0`, the row factor's term `i=4` applies the operator
`a^<4>[2]` to the input state. The source definitions give

```math
a^{(1,1)}=t_1^2,\qquad a^{(2)}=-2t_1-t_2.
```

In `a^<4>`, the weight-one term is multiplied by `binom(3,4)=0`.
The weight-two terms are multiplied by `binom(4,4)=1`. Thus the two
operators used here are `(t_1^2)[2]` and `(-t_2)[2]`.

Direct use of the source's normal-ordered definition gives the four images

```math
\begin{aligned}
(t_1^2)[2]t_1^2 &= t_1^4+4t_1t_3,\\
(t_1^2)[2]t_2 &= t_1^2t_2+4t_4,\\
(-t_2)[2]t_1^2 &= -t_1^2t_2,\\
(-t_2)[2]t_2 &= -t_2^2.
\end{aligned}
```

For example, in the second line the two mode-index choices `(-4,2)`
and `(2,-4)` each contribute `2t_4`, because `nabla_2 t_2=2`.
The choice `(-1,-1)` contributes `t_1^2t_2`.

On the four output columns `(t_1^4,t_4,t_1^2t_2,t_2^2)`, their
coefficient matrix is

```math
\begin{pmatrix}
1&0&0&0\\
0&4&1&0\\
0&0&-1&0\\
0&0&0&-1
\end{pmatrix}.
```

Its determinant is four. The intermediate span therefore has dimension
four. However, the weight-four part of `Q[t_1,t_2]` has basis
`t_1^4,t_1^2t_2,t_2^2` and dimension `p_{<=2}(4)=3`.
Mixed creation and annihilation create an additional independent
high-mode output; low-mode input labels do not confine the image to
the low-mode product space.

## Scope and reproduction

This is a counterexample to the literal `p_{<=L}(2h)` span bound.
It does not exclude a larger compressed intermediate space, a bound
`p_{<=L}(h)^(1+o(1))` for a particular growing-parameter regime, or a
different factorization of the complete crossing piece. In particular,
the four images above are a factor block. Their independence alone is
not a lower bound on the rank after composition with the other factor
and summation over `i`.

The [exact checker](../tools/check_half_restricted_labels.py) constructs
the input labels and the five weight-four output monomials, evaluates
every coefficient using the previously audited `afam` and `opcoef`, and
computes the four-by-four determinant over the integers. It uses no
floating-point rank calculation.

```sh
python3 tools/check_half_restricted_labels.py
```

The [recorded counterexample](../evidence/half-restricted-labels.json)
contains the four output vectors, the minor, its determinant, and hashes
of both the checker and the shared operator implementation. Python's
standard library suffices. No old rank suite or Lean build is needed
for this check.
