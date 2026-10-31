# Fixed smooth-square search reduces to variable-cap square decisions

The [universal smooth-square theorem](universal-smooth-square.md) supplies
one factored square that every class in the square subgroup represents.
This note reduces finding its representation to a precise decision
problem. It needs at most O(m log(T+1)) adaptive queries. Every queried
norm is still a square, uses the original split-prime support, and has
prime exponents no larger than the original ones. The target remains in
the original square subgroup throughout.

The decision oracle is an explicit unresolved assumption. In particular,
an algorithm that answers YES only for the original guaranteed instances
does not supply this oracle. Negative instances and unequal prime
exponents are essential to the reduction. The finite checker below uses
complete small-group enumeration to implement the oracle; it does not
provide an efficient decision algorithm for large discriminants.

## 1. The exact decision contract

Fix an imaginary quadratic order discriminant Delta, distinct unramified
split rational primes p_1,...,p_m, and chosen roots modulo 4p_i. These
define canonical prime-form classes g_1,...,g_m. Let H be their subgroup.
Fix an integer T>=1 and put

```math
M=\prod_{i=1}^m p_i^T.
```

For a nonnegative integer cap vector t=(t_1,...,t_m), write

```math
M(t)=\prod_{i=1}^m p_i^{t_i},
\qquad
\mathcal B(t)=\left\{\prod_{i=1}^m g_i^{2e_i}:
                         -t_i\le e_i\le t_i\right\}.
```

The norm-to-word argument in the universal-square note applies one prime
at a time. It gives the exact equivalence

```math
[Q]\in\mathcal B(t)
\quad\Longleftrightarrow\quad
Q(x,y)=M(t)^2\text{ for some integers }x,y.                 \tag{1}
```

Representations may be nonprimitive. Indeed, their coordinate gcd removes
a square factor of M(t)^2. The remaining primitive norm has an even
exponent at each split prime. Its two conjugate orientations recover an
exponent 2e_i with |e_i|<=t_i. Conversely, a bounded word constructs a
primitive represented norm with exponents 2|e_i|. Scaling both coordinates
by the product of p_i^(t_i-|e_i|) supplies M(t)^2.

Assume an exact oracle D(Q,t) for (1), for every cap vector satisfying
0<=t_i<=T. The inputs include the split-prime list and its roots, so no
factorization of M(t) is requested. A zero cap omits that prime. The
oracle must return the correct boolean answer even when no representation
exists.

The desired input is a target h=[Q] promised to belong to B(T,...,T).
For the existing application, the universal-square theorem supplies this
promise for every h in H^2. The reduction does not query or assume a
uniform decision algorithm for classes outside H^2.

## 2. Recovering the word while keeping square norms

Maintain an offset vector a, a residual target h', and caps t, with

```math
h=\left(\prod_i g_i^{2a_i}\right)h',
\qquad h'\in\mathcal B(t).                                \tag{2}
```

Initially a=0, h'=h, and every cap is T. Process the coordinates in order.
For a current cap t_i=t>1, set

```math
s=\lfloor t/2\rfloor,
\qquad r=\lceil t/2\rceil.
```

The two integer intervals centered at -s and s with radius r cover the
old interval:

```math
[-t,t]=[-s-r,-s+r]\ \cup\ [s-r,s+r].                      \tag{3}
```

For even t they overlap at zero. For odd t they overlap at -1,0,1.
Both child intervals have odd cardinality. This choice preserves the
symmetric exponent-box form and hence keeps every queried norm a square.

Replace only the i-th cap by r to obtain t'. Query

```math
D(h'g_i^{2s},t').                                         \tag{4}
```

If the answer is YES, replace h' by h'g_i^(2s) and a_i by a_i-s.
If the answer is NO, replace h' by h'g_i^(-2s) and a_i by a_i+s.
In either case use caps t'. Equation (3), the old positive promise, and
the exact oracle answer show that the chosen child still has a solution.
No query for the right child is necessary. Equation (2) is preserved.

When t_i=1, test the three possible remaining exponents as follows.
Set the i-th cap to zero. First query the shifted target h'g_i^2,
which tests exponent -1. If that fails, query h', which tests exponent
zero. If both fail, the old promise forces exponent +1, so use h'g_i^-2.
Update a_i by -1, zero, or +1 respectively. This takes at most two queries.

Each child interval lies inside its parent interval. Thus every recovered
offset satisfies |a_i|<=T. At the end all caps are zero, so B(0,...,0)
contains only the identity. Equation (2) therefore gives

```math
h=\prod_i g_i^{2a_i},\qquad |a_i|\le T.                   \tag{5}
```

The explicit CRT and quadratic-form reduction construction in
[universal smooth squares](universal-smooth-square.md) then recovers
integer coordinates representing M^2. It also verifies the final
identity and the represented norm exactly.

## 3. Query and bit costs

The cap update t -> ceil(t/2) reaches one in ceil(log_2 T) steps. Hence
the number of queries is at most

```math
m\bigl(\lceil\log_2T\rceil+2\bigr).                     \tag{6}
```

The bound also holds when m=0: the promised target is then the identity
and the represented norm is one. The implementation does not make an
initial decision query, because positivity of the initial instance is a
promise supplied by the preceding theorem.

Let b be the discriminant bit length and let L be the bit length of M^2.
Every query has a reduced target form with O(b) bits per coefficient.
Its explicit norm has at most L bits; alternatively, its factored input
uses the same prime list and smaller binary caps. Raising a reduced form
to a shift uses binary exponentiation with O(log(T+1)) compositions.
The polynomial composition and reduction routine is described in
[factored norm search](factored-norm-search.md). All non-oracle work,
including the final coordinates, has bit cost polynomial in b, L, and m.

Consequently, if the variable-cap exact decision problem can be solved
within a uniform bit bound C(b,L,m), the original fixed-square search has
cost

```math
O\!\left(m\log(T+1)\,C(b,L,m)\right)+\operatorname{poly}(b,L,m),
```

with the cap input size absorbed by L. A randomized, approximate, or
one-sided promise oracle would need additional guarantees; none is
silently substituted here.

This reduction sharpens the remaining arithmetic obligation. An exact
decision algorithm on all the smaller smooth squares would supply the
required bounded words, dominant relation rows, and class-number
multiple. The conditional generator result in
[zero-free class generators](zero-free-class-generators.md) addresses a
different obligation. It makes H equal to the full class group under its
explicit analytic inputs; it does not implement this decision oracle.

## 4. Implementation and finite evidence

The implementation is
[smooth_square_decision_reduction.py](../tools/smooth_square_decision_reduction.py).
Its callback is deliberately named `decide`, and its contract requires
unequal and zero caps as well as negative instances. It checks the final
group identity before constructing and verifying the coordinates.

The independent finite oracle in
[check_smooth_square_decision_reduction.py](../tools/check_smooth_square_decision_reduction.py)
enumerates each complete small subgroup, determines the exact orders of
the squared generators, and forms the reachable class set for each cap
vector. This full finite information is charged to the oracle. The
checker then recovers and normalizes every target witness in the square
subgroup, using the cap fixed by the existing universal-square evidence.

Run the replay with

```console
python3 tools/check_smooth_square_decision_reduction.py > evidence/smooth-square-decision-reduction.json
```

The checked [evidence](../evidence/smooth-square-decision-reduction.json)
contains 255 discriminants and 1,258 recovered target representations.
It verifies 15,046 decision queries, including 2,818 negative answers.
All three query types occur. The largest observed number of queries for
one target is 25, below its bound of 28. The tested exponent caps reach
228. An always-YES callback is rejected on a promised positive input,
confirming that original-instance existence cannot replace the required
decision information.

These finite checks validate the executable reduction and its witnesses.
They do not establish a fast decision algorithm or an improved factoring
exponent.
