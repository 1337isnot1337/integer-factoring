# Exterior monodromy and occupation-graph rank certificates

## A finite exact crossing family

Work in characteristic two. Fix an external alphabet of size `d`,
put `D=d^2`, and let `S` be an `r`-element subset of its ordered-pair
alphabet. The auxiliary space is `Lambda^r(F^D)`, of dimension
`binom(D,r)`. Write `c_q` and `a_p` for exterior creation and
annihilation operators. Their product `c_q a_p` preserves this space.
For `p=(x,y)` and `q=(v,u)` set

```math
K(z;p,q)=[p=q]I+(1+z)c_q a_p,
\qquad
P_{t,S}(z)=\langle S|\prod_{j=1}^t K(z;p_j,q_j)|S\rangle.
\tag{1}
```

All products keep their original position order. Equation (1) has
degree at most `t`, and `K(1;p,q)=[p=q]I`. Thus

```math
P_{t,S}(1)=E_{d^t}.                                \tag{2}
```

For fixed `d,S` this is a uniform finite-auxiliary crossing family
with at most `t+1` coefficient pieces and entry-construction cost
polynomial in `t`. Dense exact scalar ranks and factors cost `d^O(t)`.
The construction uses neither a cross-half commutation assertion nor
free amortization of an arbitrary interleaved matrix product.

In the row-vector convention of (1), an off-diagonal local letter
`p!=q` can act only when `q` is occupied and `p` is empty. It then
replaces `q` by `p` and contributes `1+z`. A diagonal letter leaves
the occupation unchanged and contributes `z` if that letter is
occupied, and one otherwise. The final occupation must equal `S`.
Consequently every nonzero entry is a monomial `z^a(1+z)^b`.
This exact exclusion walk is convenient for construction and checks.

The one-particle case is the matrix-unit monodromy from the companion
note. Higher exterior powers change the actual crossing coefficients;
they are not a division by a vanishing exterior dimension.

## 1. A triangular row restriction for every occupation

Order the pair alphabet with every letter outside `S` before every
letter inside it, and order words lexicographically. At the first
position where row word `p` differs from column word `q`, a nonzero
entry requires `p` outside `S` and `q` inside `S`. All preceding
diagonal letters have left the initial occupation unchanged. Therefore
the entire row polynomial matrix is triangular. Its diagonal at a
word containing `k` letters in `S` is exactly `z^k`.

In particular, for any coefficient functional `L` with `L(z^k)!=0`,

```math
\operatorname{rank}_{row} L(P_{t,S})
 \ge {t\choose k}r^k(D-r)^{t-k}.                    \tag{3}
```

To obtain the displayed minor, retain just the words having exactly
that occupation count. A triangular matrix with nonzero constant
diagonal is invertible over every field. The binomial expression
counts labels as an ordinary integer; it is not reduced modulo two.

## 2. A directed graph for the column restriction

Define a directed graph `G_S` on the pair alphabet. For column-cut
row letter `a=(x,u)` and column letter `b=(v,y)`, put an edge

```math
b\longrightarrow a
 \quad\Longleftrightarrow\quad
 a\ne b,\ (v,u)\in S,\ (x,y)\notin S.              \tag{4}
```

Let `T` be a set on which this graph is acyclic. Put `n=|T|` and
`s=|T intersect S|`, and use a topological ordering of `G_S[T]`.
Restrict all column row and column words to `T^t`.

At a first differing position, a nonzero entry requires exactly the
edge in (4). The row letter therefore occurs after the column letter
in the topological order. The restricted polynomial matrix is lower
triangular. Its diagonal is again `z^k` on words having `k` letters
in `S`. Thus

```math
\operatorname{rank}_{col} L(P_{t,S})
 \ge {t\choose k}s^k(n-s)^{t-k}
 \quad\hbox{if }L(z^k)\ne0.                        \tag{5}
```

Suppose any number of tensors from the full coefficient span sum to
the crossing. Write them as `L_j(P_(t,S))`. For every `k` represented
by a row diagonal, (2) forces `sum_j L_j(z^k)=1`. At least one piece
therefore meets both (3) and (5). The rank sum of the decomposition
satisfies the exact finite bound

```math
\boxed{\mathcal S_t\ge
 \max_{0\le k\le t}{t\choose k}
       (\sqrt{rs})^k
       (\sqrt{(D-r)(n-s)})^{t-k}.}                 \tag{6}
```

Terms with a zero label count are interpreted by the usual binomial
polynomial convention. Summing the `t+1` terms proves the simpler
asymptotic certificate

```math
\mathcal S_t\ge
 \frac{(\sqrt{rs}+\sqrt{(D-r)(n-s)})^t}{t+1}.
\tag{7}
```

These statements allow arbitrary coefficient recombinations and
arbitrary characteristic-two extension fields. They are restrictions
on this family with its fixed occupation boundary, not on all
exterior auxiliary constructions.

## 3. Rectangles and all binary occupations

If `S=A times B` is a rectangle, the entire graph `G_S` is acyclic.
Indeed for a letter `a=(x,u)` put
`f(a)=[x in A]-[u in B]`. An edge `b=(v,y)->a` in (4) has
`v in A`, `u in B`, and either `x notin A` or `y notin B`.
Consequently `f(a)<f(b)`. A strict potential forbids directed cycles.
Taking `T` to be the full pair alphabet in (7) gives

```math
\mathcal S_t\ge D^t/(t+1),\qquad \rho\ge1.         \tag{8}
```

Complementing `S` reverses every edge of (4). Thus the same full-graph
conclusion holds whenever the complementary occupation is a rectangle.
The empty and full occupations already give a single crossing piece.

For `d=2`, the only remaining occupations are the two diagonal
two-element subsets. Their graph is a directed four-cycle. Removing
one vertex gives `n=3`, with either `s=1` or `s=2`. Equation (7) yields

```math
\mathcal S_t\ge(2+\sqrt2)^t/(t+1).                 \tag{9}
```

Thus every binary fixed-occupation family is excluded from improving
the target exponent. The weakest leading restriction supplied here is
`rho>=log_2(2+sqrt(2))/2`, approximately `0.88578`.
For the diagonal occupation at `t=7`, the single `k=3` term of (6),
with the favorable three-letter restriction, has squared lower bound
`2,508,800`. This already exceeds `2^20`, whereas
`128^(2sqrt(2))<2^20`. It rejects arbitrary recombinations at that
finite width, without relying on asymptotics or decimal ranks.

## 4. The smallest ternary half-occupation candidates

There are 126 four-element subsets of the nine-letter pair alphabet
when `d=3`. The exact checker supplies, for every such `S`, a
seven-letter acyclic set `T` containing exactly three occupied letters.
Every witness includes an explicitly verified topological order.
Complementation covers the five-particle occupations.

These finite certificates and the all-length theorem (7) prove

```math
\mathcal S_t\ge
 (\sqrt{12}+\sqrt{20})^t/(t+1)                     \tag{10}
```

for every ternary half-occupation family. The base is approximately
`7.93624`, well above `3^(sqrt(2))`. No large rank scan is required.
For example, occupation mask 267 has a witness consisting of the
seven letters in mask 127. Letter `(x,y)` is encoded by `x+3y`.

The exact witnesses establish this bounded classification; no claim
about arbitrary larger alphabets is inferred from it.

## 5. Reproducible profiles and controls

The native matrix-unit checker also evaluates the exterior walk.
At seven binary positions, for a rectangular two-particle occupation,
the exact row and column ranks are respectively

`[8192,10880,12160,13440,13440,12160,10880,8192]`

and

`[8192,10880,12400,13600,13600,12400,10880,8192]`.

For the diagonal occupation the row profile is unchanged, while the
column profile is

`[8630,13656,15208,15620,15620,15208,13656,8630]`.

The scalar endpoint is independently checked by explicit exterior
operators. Full small cut matrices verify both native rank profiles.
The graph checker verifies the triangular restrictions and every
ternary witness using exact integers and binary field arithmetic.
The evidence file is `exterior-monodromy-verification.json`.

A coherent boundary, such as an initial diagonal occupation and a
final sum of the two perfect-matching occupations, still satisfies
the crossing endpoint because the initial occupation has coefficient
one in the final vector. It is outside the fixed-boundary hypothesis
above. Its seven-site tested column ranks were larger, so that example
does not supply an improvement. More general normalized exterior
boundary functionals remain a separate constructive interface.
