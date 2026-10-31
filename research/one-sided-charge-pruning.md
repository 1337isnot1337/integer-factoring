# One-sided charge changes cannot improve common-type crossings

## Result and scope

A perturbation that moves a conserved charge strictly in one direction
vanishes after a common word-type restriction. Before restriction, it
can only increase each convex Smith-profile sum, for either crossing cut
and for arbitrary diagonal physical gauges.

This gives a search rule: test charge-neutral interactions or collections
of charge changes with a nontrivial nonnegative zero-sum relation. A
single nonzero charge direction cannot improve a typed crossing.

The statement applies to arbitrary fields and to several simultaneous
integer charges. The exact finite checks below use characteristic two.
It is a restriction on the stated construction, not a lower bound on
factoring. The known factoring exponent remains `1/2`.

The derivation and checks were developed with substantial OpenAI
assistant support.

## 1. A charge compatible with both cuts

Let `A` be a finite physical alphabet and let `c:A -> R` be a charge.
Choose real constants `a,b`. For a four-leg tensor entry, define

~~~math
Q(x,y,u,v)=a c(x)+b c(y)-b c(u)-a c(v).                 \tag{1}
~~~

The crossing endpoint `E(x,y,u,v)=[x=v][y=u]` has charge zero.
Suppose a scalar tensor over a field or a Laurent-series field has the
decomposition

~~~math
T=T_0+R,\qquad
\mathrm{supp}T_0\subseteq\{Q=0\},\qquad
\mathrm{supp}R\subseteq\{Q>0\}.                  \tag{2}
~~~

Replacing `Q` by `-Q` includes strictly negative changes. Polynomial
crossing paths usually also impose `R(1)=0`, so the endpoint stays fixed.
The pruning argument itself does not need that endpoint condition.

For the row cut `(x,y)|(u,v)`, the row and column grades are

~~~math
a c(x)+b c(y),\qquad b c(u)+a c(v).                    \tag{3}
~~~

For the column cut `(x,u)|(y,v)`, they are

~~~math
a c(x)-b c(u),\qquad a c(v)-b c(y).                    \tag{4}
~~~

In each case, the row and column grade multisets agree. Sorting by these
grades makes `T_0` block diagonal with square blocks and makes `R`
strictly block triangular. Diagonal row or column weights leave this
support unchanged. The same holds for the physical gauges used in
[type-optimized pressure](type-optimized-crossing-pressure.md).

## 2. Smith-profile monotonicity

Work over a discretely valued field. Suppose the chosen square
flattening of `T_0` is generically invertible. Let its ascending Smith
orders be `nu_1,...,nu_N`, and let those of `T` be
`mu_1,...,mu_N`, after the same arbitrary diagonal gauges.
Then

~~~math
\sum_{j=1}^k\mu_j\le\sum_{j=1}^k\nu_j
\quad(1\le k<N),\qquad
\sum_{j=1}^N\mu_j=\sum_{j=1}^N\nu_j.                 \tag{5}
~~~

Indeed, the minimum valuation of a nonzero `k`-minor of a block diagonal
matrix is attained by a minor that chooses equally many rows and columns
inside each grade block. Its determinant is a product of within-block
minors. The identical row and column sets in `T` give a block triangular
minor with the same diagonal blocks, so its determinant is unchanged.
The minimum over all minors of `T` can only be smaller. At full size,
the triangular determinant is exactly the original determinant. The
determinantal characterization of Smith orders gives (5).

Thus `mu` majorizes `nu` when the lists are ordered in descending order.
For every convex real function `f`,

~~~math
\sum_j f(\mu_j)\ge\sum_j f(\nu_j).                    \tag{6}
~~~

In particular, for every `s>0`,

~~~math
\sum_j s^{-\mu_j}\ge\sum_j s^{-\nu_j}.               \tag{7}
~~~

Apply this separately at zero and infinity, with either orientation of
the two cuts. The restored physical-height term is the same, so the
product Smith pressure cannot improve. This includes arbitrary common
physical types and arbitrary diagonal physical gauges. Singular
flattenings require a separate statement about their nonzero invariant
factors; (5)--(7) are asserted here only in the invertible square case.

There is also a direct rank version. Apply any common coefficient-linear
functional to every entry of `T_0` and `T`. The resulting matrices retain
the same block diagonal and strictly triangular pattern. Every nonzero
minor witnessing the rank of the former remains unchanged in the
latter. Thus the rank cannot decrease in either cut. This assertion
does not require generic invertibility.

## 3. Exact equality after type restriction

Take `m` local factors, allowing position-dependent tensors satisfying
(2) with the same charge. Restrict all four external words to the same
letter-count type. Then

~~~math
\left.\bigotimes_{j=1}^m(T_{0,j}+R_j)\right|_{\rm type}
=\left.\bigotimes_{j=1}^mT_{0,j}\right|_{\rm type}.     \tag{8}
~~~

To see this, fix four words of the common type. Their total charge is
zero, because the four sums of `c` over the words agree. Expand the
product on the left. A nonzero term using only base entries has total
charge zero. A nonzero term using at least one perturbing entry has
strictly positive total charge, since all other local charges are
nonnegative. The latter is impossible for the fixed boundary words.
Every such term is therefore zero before summation.

The proof still holds when each local entry is an auxiliary matrix and
the tensor word uses their ordered product, with arbitrary fixed or free
auxiliary boundaries. It does not commute those matrices or assume
their products are nonzero. It only excludes every local charge pattern
that contains a perturbation.

The conclusion also holds when the four types differ but have equal
total charge in the combination (1). A common type is a convenient
sufficient condition, not a necessary one.

## 4. Several charges and a practical pruning criterion

Let `Q` take values in `Z^r`, with each coordinate of the additive form
(1), and with the base supported at zero. If a
linear functional `ell` is strictly positive on every perturbation
charge vector, apply the preceding proof to `ell(Q)`.

More generally, a perturbative tensor word can survive a zero-charge
restriction only if its multiset of nonzero local charge vectors obeys

~~~math
\sum_{v\in V} n_v v=0,\qquad n_v\in\mathbb Z_{\ge0},
\qquad \sum_v n_v>0.                                  \tag{9}
~~~

For a finite set of integer vectors, absence of such a relation is
equivalent to the existence of a strictly positive separating linear
functional. One direction follows by applying the functional to (9).
For the other, separate zero from the compact convex hull of `V`.
If zero belongs to that hull, the rational linear system defining its
convex representations has a rational solution; clearing denominators
gives (9). This is a necessary survival criterion, not a sufficient
criterion for a rank or pressure improvement.

For two copies of the four-letter Gaussian seed, retain each copy's
integer charge separately. Each copy has four positive and four negative
fermions among its four physical legs. There are exactly

~~~math
2\binom42^2+(4\cdot4)^2=328                            \tag{10}
~~~

quartic exterior monomials neutral in both copy charges. A pure positive
quartic has copy-occupation vector `(k,4-k)`, with counts
`1,16,36,16,1` for `k=0,...,4`. A positive and negative quartic can cancel
both charges only when these copy occupations agree, leaving

~~~math
1^2+16^2+36^2+16^2+1^2=1810                           \tag{11}
~~~

opposite-charge pairs. These counts reduce the search space; they do
not assert an improving construction.

## 5. Why opposite directions remain open

The strictly one-sided hypothesis is essential. Over `F_2((z))`, take

~~~math
A_0=I_2,\qquad
A=\begin{pmatrix}1&1+z\\1+z&1\end{pmatrix}.             \tag{12}
~~~

The off-diagonal changes have opposite grade directions, and `A(1)=I_2`.
At zero the profiles are `(0,0)` and `(0,2)`, respectively. Hence
`1+s^{-2}<2` for `s>1`: the zero-end profile improves. This is only a
matrix control. It does not claim an improvement of the two-cut,
two-end crossing pressure.

## Verification

~~~sh
python3 tools/check_one_sided_charge_pruning.py \
  --output evidence/one-sided-charge-pruning.json
~~~

The checker verifies both cuts and both places, independently chosen
physical gauges, every tested coefficient functional, common-type
polynomial identities, ordered noncommuting auxiliary products, the
quartic counts, and the opposite-direction negative control. The
ordinary proof establishes the unrestricted statements above; finite
checks only guard the formulas and implementation conventions.
