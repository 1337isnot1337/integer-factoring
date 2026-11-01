# Nonlinear label maps give a genuinely new layer family, but no gain yet

## Result and scope

This note records an effective crossing family that escapes the common
linear-code convolution algebra. An explicit Toffoli relabeling produces
noncommuting polynomial layers and a non-Gaussian tensor. Its first
complete finite test does **not** improve the crossing rank bound.

At binary width three, the entire coefficient symmetry span generated
by the original pencil, its Toffoli image, and their three-layer product
has dimension thirteen over `F_2`. Exhausting all `8191` nonzero vectors
proves that its optimal maximum rank product is **224**, exactly the
original pencil's value. Arbitrary piece counts inside that finite
`F_2` span are covered. Extension-field recombinations are not covered.

There is no new factoring exponent, no successful finite rank-sum seed,
and no asymptotic rank theorem for the new nonlinear family.

## 1. A uniform polynomial family

Let `E_2` be the two-label swap and `I_4` the row identity. In
characteristic two, set

```math
p(z)=I_4+z(E_2+I_4),\qquad P_\lambda(z)=p(z)^{\otimes\lambda}.
                                                               \tag{1}
```

The elementary nilpotent satisfies `(E_2+I_4)^2=0`, so

```math
P_\lambda(z)^2=I,\qquad P_\lambda(1)=E_{2^\lambda}. \tag{2}
```

Choose a uniformly computable permutation `pi_lambda` of the binary
labels, for example a reversible circuit with polynomially many
CNOT and Toffoli gates. Apply that same permutation to all four tensor
legs. In the row arrangement this is conjugation by the permutation
matrix `Pi_lambda tensor Pi_lambda`; write the resulting pencil as
`Q_lambda(z)`.

Both cut ranks of each constant tensor are preserved by this label
change, and the crossing itself is fixed. Hence

```math
Q_\lambda(z)^2=I,\quad Q_\lambda(1)=E_{2^\lambda},\quad
T_\lambda(z)=P_\lambda(z)Q_\lambda(z)P_\lambda(z)
                                                               \tag{3}
```

satisfies `T_lambda(1)=E_(2^lambda)` and has degree at most `3lambda`.
Its coefficients therefore form an exact scale-one crossing
decomposition with polynomially many pieces.

All matrices have dimension `2^(2lambda)`. Direct coefficient
construction, multiplication, and elimination take `2^O(lambda)`
operations. The joint coefficient span of any fixed number of these
families and their eight crossing symmetries has dimension `O(lambda)`.
Over `F_2`, exhaustive span optimization and retained rank
factorizations also take `2^O(lambda)` operations. Thus the family has
the effectiveness required for the existing transfer framework.
If the crossing lies in the span of the candidates below a chosen
rank cutoff, retain an independent subset of those candidates that
spans the same space, together with the crossing expression.
This uses at most `O(lambda)` pieces despite the exhaustive search.

What is missing is a uniform rank-product bound below the half
threshold. A fixed gate on three coordinates, extended by the identity,
does not establish such a bound. An asymptotic construction must specify
and analyze the full sequence of label maps.

## 2. Why linear relabelings collapse

For a common linear label change, additive conservation remains
`x+y=u+v`. Put `S=x+y` and `w=x+u`. Each fixed-`S` row block is a
convolution operator on the abelian group `F_2^lambda`.

All common-linear relabelings therefore lie in the same direct sum
of commutative group algebras. Their row pencils commute. In view of
(2), the palindromic product then collapses:

```math
P_\lambda Q_\lambda P_\lambda=Q_\lambda.             \tag{4}
```

Products of distinct linear bases can also be rewritten as linear
code-section frames; they are already represented by the general
code-section family. Linear maps are not a way around that family.
The nonlinear label map below avoids this commutativity argument.

## 3. Explicit non-Gaussian and noncommuting controls

Use external-bit order
`x_0,...,x_(lambda-1),y_0,...,u_0,...,v_0,...`.
A normalized even Gaussian tensor with vacuum coordinate one must
have each four-bit coordinate equal to the Pfaffian of its two-bit
coordinates. In characteristic two this is the sum of the three
pairing products.

For the two-bit CNOT map
`(a_0,a_1) -> (a_0,a_1+a_0)`, the common relabeling of (1) has

```math
Q_{\{0,1,3,4\}}=z+z^2,
\qquad \operatorname{Pf}(Q_{\{i,j\}})_{i,j\in\{0,1,3,4\}}=0.
                                                               \tag{5}
```

Thus it is non-Gaussian over `F_2(z)`. This is a scope control only:
that small CNOT relabeling is already among the earlier finite label
permutations, and its layers still commute as in (4).

For three bits, use the Toffoli permutation

```math
\pi(a_0,a_1,a_2)=(a_0,a_1,a_2+a_0a_1),
```

whose label table is `(0,1,2,7,4,5,6,3)`. At the four-bit coordinate
`{0,1,6,10}`, its relabeled pencil has coordinate zero and Pfaffian
`z+z^2`. It too is non-Gaussian over the rational function field.

Unlike the CNOT case, this Toffoli image does not commute with the
original row pencil. The `z^2` coefficient of `PQ+QP` is one at
row labels `(x,y)=(0,3)` and column labels `(u,v)=(1,2)`.

Direct polynomial multiplication gives the exact identity

```math
PQP=Q+(z^2+z^4)D,\qquad (r(D),c(D))=(2,8).          \tag{6}
```

The correction vanishes at `z=1`, as it must. It is outside the
eleven-dimensional coefficient symmetry span of `P,Q` alone; adding
the product and its symmetries enlarges the span to dimension thirteen.
This is a new low-rank direction in this finite construction, but that
fact alone does not decompose the crossing more efficiently.

## 4. Complete finite outcome

The checker generates all coefficients of `P`, `Q`, and `PQP`, applies
every one of the eight crossing-preserving leg permutations, and
constructs their exact thirteen-dimensional span. It then computes
both matrix ranks of every nonzero vector.

The vectors with product below `224` span a ten-dimensional subspace
that does not contain the crossing. The vectors with product at most
`224` span the whole space and do contain it. A retained binary
functional annihilates every strictly cheaper vector and takes value
one on the crossing. This proves the finite optimum without any
restriction on the number of pieces inside the span.

The original coefficient ranks

```text
(64,1), (28,8), (8,28), (1,64)
```

already attain that optimum. The nonlinear construction has therefore
passed the scope and effectiveness checks but has not passed a useful
rank-improvement test.

```sh
python tools/check_nonlinear_layer_probe.py \
  --output evidence/nonlinear-layer-probe.json
```

The replay uses only the Python standard library. Its polynomial,
Pfaffian, commutator, coefficient-span, rank, and separator checks are
exact. None of them proves an asymptotic contraction bound.

The elementary coefficient construction is developed in the pinned
[PR #23 source](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/binary-coefficient-entropy.md).
The [code-section normal forms](code-section-normal-forms.md) give the
linear frame correspondence used above. Substantial OpenAI assistant
support was used in this exploration and exact replay.
