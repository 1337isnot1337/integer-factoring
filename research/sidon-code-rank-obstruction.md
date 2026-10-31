# Sparse additive graph codes do not rescue the elementary crossing

## Result and scope

A subsequent [active-difference theorem](difference-support-rank-floor.md)
extends this argument to arbitrary codes and uniformly excludes the
multiplication graph family left open by this note.

Let C be a binary word code of size q such that each nonzero xor difference
comes from at most one unordered pair of distinct codewords. For the
elementary crossing restricted to C, every decomposition of its endpoint
into tensors from the full multivariate coefficient span satisfies

```
sum_i sqrt(rank_row(A_i) rank_col(A_i)) >= q(q-1)/sqrt(2).
```

This holds over every field of characteristic two and for any number of
pieces. In particular, these codes fail the repeated-block rank-sum target
q^sqrt(2) for q >= 4. With only s pieces, some piece has rank product at
least q^2(q-1)^2/(2s^2). For q exponential in the label width and s
polynomial in that width, this also excludes every fixed rho < 1 in the
ordinary maximum-rank-product interface.

This is a restriction on one crossing's coefficient span, not a lower bound
for factoring, general physical codes, or order-sensitive auxiliary seeds.
It is an ordinary proof with finite checks, not independent mathematical
review or new Lean formalization. The complete factoring score stays 1/2.

## Proof by literal blocks

Use independent variables z_j at the word positions. The local tensor is
I + z_j(E+I). A nonzero word entry requires x xor y = u xor v and
x intersection y = u intersection v. The difference hypothesis implies
that, when x != y, the unordered pairs {x,y} and {u,v} coincide. Thus the
row matrix is a direct sum of singleton diagonal blocks and 2-by-2 blocks,
one for each unordered pair {x,y}.

For any coefficient functional L, write D = support(x xor y) and set

```
a_xy = L(product_(j in D) (1+z_j)),
b_xy = L(product_(j in D) z_j).
```

The corresponding row block is [[a_xy,b_xy],[b_xy,a_xy]]. If b_xy is
nonzero, this block has rank at least one. Let t(L) count those active
unordered pairs. Consequently r(L) >= t(L).

In the column cut, the off-diagonal labels (x,u), x != u, have exactly one
possible nonzero entry: column (u,x), of value b_xu. These labels form a
permuted diagonal matrix, independent of the block on diagonal labels
(x,x). Each active unordered pair supplies two nonzero diagonal entries.
Therefore c(L) >= 2t(L), and sqrt(r(L)c(L)) >= sqrt(2)t(L).

For a decomposition sum_i A_i = kappa E_C with kappa != 0, the sum of the
b_xy values is kappa for every unordered pair. Each pair must therefore
be active in at least one summand. Summing gives

```
sum_i t(L_i) >= binom(q,2),
sum_i sqrt(r_i c_i) >= sqrt(2) binom(q,2).
```

The maximum-product statement follows by pigeonhole and r_i c_i >= 2t_i^2.
For q >= 4 the sum bound is already greater than q^(3/2): this reduces
to (q-1)^2 > 2q. Since 3/2 > sqrt(2), the required threshold is missed.
No positivity of coefficients, semisimplicity, or characteristic-zero
lifting is used in the argument.

## Cubic field graphs meet the hypothesis at every size

In any fixed binary basis of F_(2^k), consider C = {(a,a^3): a in F_(2^k)}.
For a nonzero difference h, the derivative is

```
(a+h)^3+a^3 = h*a^2+h^2*a+h^3.
```

The linearized part has kernel {0,h}, since h*a*(a+h)=0. Every derivative
fiber therefore has at most two elements, precisely the two orientations
of one unordered pair. The code meets the difference hypothesis. Its
restricted elementary crossing is excluded uniformly, not only in the
small tested fields. This explanation is self-contained and does not
assume a general theorem about APN functions.

## Exact diagnostic controls

Run from the repository root:

```
python3 tools/check_sidon_code_crossings.py > /tmp/sidon-code-crossings.json
cmp /tmp/sidon-code-crossings.json evidence/sidon-code-crossings.json
```

The script checks 16 finite codes and all 752 binary coefficient functionals.
Fifteen codes meet the difference hypothesis. Every polynomial entry is
constructed twice (subset expansion and dense binomial expansion); endpoints
are literal crossings, ranks are checked with reversed row order, and every
active-pair inequality is checked. Complete finite rank-sum optimization
uses integer square-root lower weights at precision 2^-24 and verifies all
shortest-path dual inequalities. The decimal upper values are descriptive;
the lower-bound comparisons use exact integer squares.

The two cubic graphs (k=2,3) have optimal rank-sum intervals near 13.74597
and 50.33202. The distinct multiplication graph {(a,b,ab)} at k=2 is NOT
a difference-unique code, so the theorem does not cover it. Its tested
six-bit, sixteen-symbol coefficient span nevertheless has lower rank sum
greater than 191.27121, versus target approximately 50.45251. This single
basis and finite case does not exclude its general family.

## Consequence for construction choices

Low collision counts are not automatically helpful: they isolate swap
entries into independent column blocks while retaining separate row blocks.
A new physical-code construction must permit enough nontrivial mixing to
avoid that proof, or change the local polynomial/auxiliary mechanism.
Do not enlarge the cubic-graph catalogue in search of an escape already
excluded by the all-size identity above. Higher-additive-energy codes and
language-dependent auxiliary endpoints remain outside this result.

Predecessors: [the arbitrary-code entry formula](downset-code-crossings.md)
and [the fixed-block rank-sum contract](auxiliary-coefficient-span.md).
The code-generation experiment was suggested during the director checkpoint
in [Discussion 51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19323155).
This note and checker were developed with OpenAI assistant support.
