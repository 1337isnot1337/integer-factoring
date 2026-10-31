# The complete Boolean Cauchy Smith profile

## Status

An elementary row division gives every invariant factor of the full
Boolean substitution operator at every nonzero grid center. The result
holds at all packet lengths that are powers of two. At either nearest
center, and at an odd local evaluation point, the inverse loses exactly
one bit of precision, independent of the packet length.

There is also a uniform inverse-entry algorithm with time
`exp(O(L^2))` and polynomial space. This cost fits the entry-query budget
in the quantitative contraction proof. The remaining obstruction to
inserting this inverse directly is its interaction between groups: its
normalized reduction modulo two already has group-cut rank
`2^(L-|t|)`. In particular, the one-bit inverse at `t=1` has rank
`2^(L-1)` across the source's two-group cut.

These statements concern the full operator, whose inputs and outputs
are all Boolean packet states. They do not provide a scalar decoder
from one deformed contraction. No improved factoring exponent is
claimed. The canonical scalar collisions in
[the decoder note](canonical-packet-decoder-frontier.md) remain valid.

The separate [arithmetic Cauchy theorem](arithmetic-cauchy-factorial-law.md)
proves every invariant factor of the physical squared Cauchy matrix at
arbitrary sizes. That smaller matrix has a different profile from the
full Boolean substitution operator studied here; its theorem does not
remove this operator's grouped-rank cost.

## 1. Definitions and coefficient rings

Let `L=2^a`, `a>=1`, and use the Cauchy kernel from
[the quantitative contraction source][quantitative]:

\[
 \mathcal S=\{-L,-L+1,\ldots,L-1\},\qquad \alpha=2L,
 \qquad C(T)=\prod_{q\in\mathcal S}(T-q),
\]

\[
 W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},\qquad 0\le x,y<L.
 \tag{1}
\]

The quotient in (1) is a polynomial. Work first over
`O=Z_(2)`, the rationals with odd denominator. The argument also works
over any unramified extension with residue field of characteristic two,
including the source's `O=Z_(2)[U]/(f)`. An additional truncated variable
`V` is handled by extending scalars after the factorization below.

For an `L` by `L` matrix `A`, let `S(A)` be its Boolean monomial
substitution matrix. Its rows and columns are subsets of `[L]`, and

\[
 S(A)_{I,J}=[X_J]\left(\prod_{x\in I}(AX)_x
                     \bmod(X_y^2-X_y)_{y<L}\right).
 \tag{2}
\]

Here `X_J` is the squarefree monomial indexed by `J`. Boolean reduction
is essential: repeated occurrences of the same column variable reduce
to one occurrence. In general `S(A)S(B)` is not `S(AB)`.

The corresponding matrix in the Boolean value basis is

\[
 K(A)_{Y,X}=\prod_{x=0}^{L-1}l_{Y_x}((AX)_x),\qquad
 l_0(v)=1-v,\quad l_1(v)=v.
 \tag{3}
\]

Let `Z` be the Boolean zeta matrix, `Z_(X,I)=1` when `I` is a subset
of `X`, and zero otherwise. Its inverse is the integral Boolean
Möbius matrix, so

\[
 K(A)=(Z^{\mathsf T})^{-1}S(A)Z^{\mathsf T}.
 \tag{4}
\]

Thus the value and monomial operators have the same Smith orders.
Both basis changes are tensor products of two by two integral unit
matrices. They also preserve ranks across cuts that assign each Boolean
bit, on either side of the operator, to a specified group.

## 2. Physical row normalization

Fix a nonzero `t` in `S`, and put

\[
 b=|t|,\qquad
 B_t=\{x: x+t\notin[0,L)\},\qquad
 \sigma_t(x)=(x+t)\bmod L.
 \tag{5}
\]

There are exactly `b` boundary rows. On the disk `T=t+2Lz`, define

\[
 W'_{xy}(z)=
 \begin{cases}
  W_{xy}(t+2Lz)/(2z),&x\in B_t,\\
  W_{xy}(t+2Lz),&x\notin B_t.
 \end{cases}
 \tag{6}
\]

Every quotient in (6) is a polynomial in `O[z]`, and

\[
 W'(z)\equiv P_{\sigma_t}\pmod 2.
 \tag{7}
\]

The congruence is a polynomial identity: the right side is independent
of `z`.

**Proof.** Write `w_t=C'(t)/C'(0)`. The factorial formula for this ratio
and the oddness of every binomial coefficient in row `2L-1` show that
`w_t` is a two-adic unit. The local root product is

\[
 \frac{C(t+2Lz)}{C'(0)}
   =2Lz\,w_t\prod_{r\in\mathcal S\setminus\{t\}}
                    \left(1+\frac{2Lz}{t-r}\right).
 \tag{8}
\]

Every nonzero difference `t-r` has valuation at most `a`, so every
nonconstant factor coefficient in (8) is even. If `q=y-x=t`, the
corresponding kernel entry is the product in (8) without `2Lz` and
reduces to one. If `q!=t`, cancellation of its root factor gives

\[
 \frac{W_{xy}(t+2Lz)}{2z}
  =\frac{Lw_t}{t-q}
       \prod_{r\in\mathcal S\setminus\{t,q\}}
                    \left(1+\frac{2Lz}{t-r}\right).
 \tag{9}
\]

This polynomial is integral. Its reduction is nonzero exactly when
`v_2(t-q)=a`. The only possibilities are `t-q=+L` or `-L`.

For a nonboundary row the unreduced entry with `q=t` is its unique
odd entry. All other entries are divisible by `2z`. For a boundary
row `q=t` does not occur, and the unique column satisfying
`t-q=+L` or `-L` is `y=sigma_t(x)`. Formula (9) makes that normalized
entry odd and all other normalized entries even. This proves (6)
and (7). The endpoint `t=-L` is included; then all rows are boundary
rows and the normalized permutation is the identity. \(\square\)

## 3. Exact Boolean Smith factors

Define diagonal entries and a unit operator by

\[
 d_I(z)=(2z)^{|I\cap B_t|},\qquad D_t(z)=\mathrm{diag}_I d_I(z),
 \qquad U_t(z)=S(W'(z)).
 \tag{10}
\]

The physical row factors pull out of every input monomial, even after
Boolean reduction. Therefore

\[
 S(W(t+2Lz))=D_t(z)U_t(z).
 \tag{11}
\]

Equation (7) gives

\[
 U_t(z)\equiv P_{\sigma_t,\mathrm{subsets}}\pmod 2.
 \tag{12}
\]

In particular, at every integral evaluation point `z=zeta`,
`det U_t(zeta)` is a unit. No genericity assumption on `zeta` is used.

**Theorem.** Suppose `zeta!=0` and `h=v_2(zeta)>=0`. The Smith orders
of either full Boolean operator `S(W(t+2Lzeta))` or
`K(W(t+2Lzeta))` are

\[
 \boxed{\quad j(1+h)\quad\text{with multiplicity}
          \binom bj\,2^{L-b},\qquad 0\le j\le b.\quad}
 \tag{13}
\]

At `zeta=0` the rank is exactly `2^(L-b)`. More precisely, (11) is
an equivalence to a diagonal matrix with that many unit entries and
all remaining entries zero.

**Proof.** Right multiplication by `U_t(zeta)^(-1)` is unimodular.
The valuation of `d_I(zeta)` is `(1+h)|I intersect B_t|`. There are
`binom(b,j) 2^(L-b)` subsets with intersection size `j`. Equation (4)
transfers the factorization to the value basis. At zero, precisely
the entries with empty boundary intersection survive. \(\square\)

For `t=+1` or `-1`, and `zeta` a unit, exactly half the Smith orders
are zero and half are one. The inverse needs exactly one additional
bit when acting on the entire packet vector. This is compatible with
much larger scalar reconstruction losses: a single bilinear scalar
does not reveal the packet vector on which this inverse acts.

The factorization is over `O` before reduction modulo a precision.
After extending scalars to `O[V]/(V^ell)` it remains an exact
invertible-equivalence factorization, although Smith terminology over
that non-domain is unnecessary.

## 4. Permanent content and parity-character contractions

The top squarefree coefficient in (2) can have no repeated column.
Consequently it is the permanent, and row normalization gives

\[
 \mathrm{per}\,W(t+2Lz)=(2z)^b\,\mathrm{per}\,W'(z)
                 \equiv2^bz^b\pmod{2^{b+1}}.
 \tag{14}
\]

The complete local permanent polynomial has content exactly `b`.
The same row argument gives the corresponding determinant congruence,
up to its irrelevant odd sign.

Let `A` be an arbitrary row tensor and let `A_tilde` be its multilinear
extension. Its top coefficient is

\[
 a_{[L]}=\sum_Y(-1)^{L-|Y|}A(Y).
\]

The parity character `chi(X)=(-1)^|X|` annihilates every Boolean
monomial of degree below `L`. Substitution cannot increase degree.
It follows exactly that

\[
 \sum_X(-1)^{|X|}\widetilde A(WX)
       =(-1)^L a_{[L]}\,\mathrm{per}\,W.
 \tag{15}
\]

For the even lengths under consideration the sign is positive. Thus
the coefficient of `z` at either nearest center, divided by two and
reduced modulo two, extracts `sum_Y A(Y)` modulo two when the column
tensor is the parity character. This is a controlled scalar functional.
It does not extract `sum_Y A(Y)B(Y)` for an arbitrary column tensor.
Lower-degree components of that general contraction still mix.

## 5. A normalized inverse without a precision-dependent matrix search

For nonzero `z`, define in the value basis

\[
 R_t(z)=(2z)^bK(W(t+2Lz))^{-1}.
\]

An integral expression, also defined at zero, is

\[
 R_t(z)=(Z^{\mathsf T})^{-1}
 U_t(z)^{-1}\,\mathrm{diag}_I(2z)^{b-|I\cap B_t|}
 Z^{\mathsf T}.
 \tag{16}
\]

The entries are rational functions with denominator `det U_t(z)`,
which evaluates to an odd unit at every integral point. They can also
be viewed in the two-adic completion of `O[z]`: after removing the
permutation in (12), the inverse is a convergent geometric series in
an even polynomial matrix. Formula (16) satisfies

\[
 K(W)R_t=R_tK(W)=(2z)^bI.
 \tag{17}
\]

This inverse is not obtained by substituting the ordinary linear
inverse `W^(-1)` into (2). Boolean substitution does not preserve
composition of arbitrary linear maps.

### 5.1 Uniform entry computation

For an integral `zeta` and target precision `p`, an entry of (16) can
be computed using

\[
 \exp(O(L^2))\,\mathrm{poly}(p,e,\ell,L)
       \quad\text{time and}\quad
 \mathrm{poly}(p,e,\ell,L)\quad\text{space}.
 \tag{18}
\]

Here `e` is the degree of the unramified extension; `ell` may be
omitted if the truncated variable is absent. Input coordinate lengths
are included in the polynomial factor. The algorithm has these steps.

1. **Evaluate the normalized physical kernel.** Use (8) and (9),
   canceling known powers of two in each rational scalar before
   reduction. All remaining scalar denominators are odd. The
   integers involved in the factorial ratio and root products have
   `O(L log L)` bits. Each entry takes a polynomial number of ring
   operations; the value at `zeta=0` uses the removable quotient in
   (9), so no division by `zeta` occurs. More generally this evaluation
   may form the rational constants exactly and only then reduce
   modulo `2^p`. That avoids cancellation of a two-power denominator
   after finite precision has already been discarded.

2. **Query a Boolean monomial entry.** For any `I,J`, Boolean Möbius
   inversion writes `U_(I,J)` as a signed sum over subsets of `J` of
   the product of the row forms indexed by `I`. Enumerating those
   at most `2^L` states uses `exp(O(L))` time and polynomial space.
   This is an implicit entry oracle for the `N=2^L` dimensional
   matrix `U`, not a stored `N` by `N` array.

3. **Query a determinant and an adjoint entry.** The division-free,
   bounded-fan-in determinant and adjoint circuits of
   [Berkowitz (1984)][berkowitz] work over any commutative ring.
   They have polynomial size in `N`, depth `O(log^2 N)`, and
   logspace-uniform gate construction. A depth-first evaluation
   recomputes shared subcircuits instead of storing them. A
   constant-fan-in circuit of depth `O(log^2 N)` then unfolds to
   at most `exp(O(log^2 N))=exp(O(L^2))` operations. Gate identifiers,
   the recursion stack, and ring values use polynomial space in
   `L,p,e,ell`. The logspace gate generator can be simulated within
   polynomial time in `N` per queried gate, which is absorbed by
   the same time bound. Cramer's formula uses the transposed adjoint
   entry with the usual row/column convention.

4. **Divide only by a unit.** By (12), `det U` reduces to one in
   characteristic two. Its inverse modulo `2^p` can be found by
   repeated precision doubling, starting from one modulo two.
   The update `x <- x(2-(det U)x)` applies over the extension ring
   and after adjoining `V`. It costs polynomially many bit
   operations in the ring representation. No extra `p`-dependent
   loss is needed to compute an entry of `U^(-1)`.

5. **Apply the diagonal and zeta matrices.** At most `N^2` terms
   suffice for an arbitrary requested entry of (16). Each term
   contains a nonnegative power of `2zeta`, so it is computed
   without division. These loops preserve (18).

The dense inversion routines in the accompanying finite checker are
independent controls of (16). They do not implement the polynomial-space
determinant circuit algorithm. Formula (18) is a complexity consequence
of the cited circuit theorem and the explicit entry oracle above.

### 5.2 Comparison with the source budget

The quantitative proof sets

\[
 u_d=\left\lfloor\varepsilon\frac{\lambda}{2^d}\right\rfloor,
 \qquad a_d=2\sum_{j=0}^d u_j,\qquad L=2^{a_d},
 \qquad\lambda=\lceil\log_2s\rceil,
 \qquad\varepsilon=\frac1{128\cdot1024}.
\]

Its parameter inequalities imply

\[
 L\le s^{1/1000},\qquad 2L\le\Delta=\lceil s^{9/10}\rceil.
\]

The source permits implicit tensor entries costing `exp(O(s^(9/10)))`
time and polynomial space. Its actual packet tensors use
`exp(O(L log s))` time. Since `L^2<=s^(2/1000)`, (18) fits the permitted
entry budget whenever the precision and coefficient representation
have the source's polynomial size. A new inverse tensor has `2L`
Boolean ports, so its degree also fits the source degree cap.

The description can record the original packet data, `t`, `zeta`, and
its creation precision. Each later query recomputes its creation-time
residue and takes the same canonical representative. It need not call
an inherited entry routine, so this description does not create a
recursive stack of composed tensor-entry computations.

These entry and degree facts are necessary, but they do not establish
the edge-count or grouping invariant needed by the recursion.

## 6. Where direct inverse insertion reconnects the groups

At a unit `zeta`, the identity packet can be written as

\[
 \delta_{Y,X}=(2\zeta)^{-b}
      \sum_Z K(W(t+2L\zeta))_{Y,Z}R_t(\zeta)_{Z,X}.
 \tag{19}
\]

This replaces one packet by the deformed packet plus a new `2L`-port
inverse tensor. It restores the original contraction algebraically.
It does not preserve the source's factorization into smaller groups.

To see the problem in the exact source geometry, suppose the two
symbols of a packet differ in coordinate `i=v_2(t)<a`. Original row
and column orientations in that coordinate are opposite. The source
groups vertices by their symbol in coordinate `i`. An original row
slot `x` and a column slot `y` are in the same group exactly when
their `i`th bits differ. Since adding `t` flips that bit, the leading
shift of `K(W)` is local between its row and intermediate-column ports.

Both the intermediate ports `Z` and the original ports `X` in `R_t`
have the column orientation. Thus their groups now agree when their
`i`th bits agree. The inverse cyclic shift flips that bit, so it
connects these two groups. Giving the intermediate ports the row
orientation instead would move the same problem back into `K(W)`.
An intermediate Boolean wire has one group assignment, not separate
assignments at its two ends.

This is visible already modulo two. Equation (16) becomes

\[
 R_t(z)\bmod2=(Z^{\mathsf T})^{-1}
 P_{\sigma_t,\mathrm{subsets}}^{-1}
 \mathrm{diag}_I\mathbf1[B_t\subseteq I]\,Z^{\mathsf T}.
 \tag{20}
\]

It is independent of `z`. In the monomial basis the diagonal factor
is a tensor product: a rank-one projection on each of the `b`
boundary bit factors, and an identity on each of the other `L-b`
bit factors. The inverse permutation sends every matched bit wire
across the two groups. A projected wire contributes cut rank one;
an identity wire contributes cut rank two. The zeta changes act
separately on each bit and preserve the cut rank. Therefore

\[
 \boxed{\quad \mathrm{rank}_{\mathbb F_2}
     (R_t\text{ flattened across the two source groups})
             =2^{L-b}.\quad}
 \tag{21}
\]

For `t=+1` or `-1`, any direct sum of products of one tensor per
group must use at least `2^(L-1)` terms even for the reduction of the
normalized inverse modulo two. A procedure that expands each packet
this way independently therefore spends an exponential-in-`L` number
of choices per packet. The source may have `Theta(M/L)` packets,
so this direct replacement can cost `exp(Theta(M))` across a call.
The small per-entry inverse cost does not remove that cost.

At an outer retained center, `b=L-1`, the leading cut rank is only
two. Its normalization instead costs `L-1` bits per packet. This
does not alone rule out a useful structured higher-precision
decomposition. The exact local physical matrix, its higher-order
crossings, and the tensor product precision loss all have to be
accounted for; a small rank of (20) is not enough.

Equation (21) is an obstruction to this specified direct inverse
insertion and product expansion. It is not a lower bound for all
contraction algorithms, all inverse factorizations, or scalar
reconstruction methods that exploit particular input tensors.

## 7. Reproducible finite controls

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tools/check_boolean_cauchy_smith.py \
  --output /tmp/boolean-cauchy-smith-profile.json
```

The checker builds the physical kernel from exact rational global
polynomials, independently constructs Boolean-value and monomial
operators, and compares their integral change of basis. Small matrices
also undergo minimum-valuation elimination and dense unit inversion.
The tested zero points use exact derivatives to evaluate removable
quotients. Nonzero even points test the valuation multiplier in (13).

It checks the complete local permanent polynomial, not only its value
at one point. Arbitrary signed row tensors test (15). Bit-packed
flattenings independently compute the cut ranks in (21). The recorded
counts, cases, and source hashes are in
[boolean-cauchy-smith-profile.json](../evidence/boolean-cauchy-smith-profile.json).

[quantitative]: https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex
[berkowitz]: https://doi.org/10.1016/0020-0190(84)90018-8
