# A degree-linear spectral trace interface

## Constructive reduction

For the implicit pole Cauchy matrix `M`, form

\[
 S=MM^{\mathsf T}.
 \tag{1}
\]

Modulo two, `S` is the diagonal indicator of the true factor rows.
The number `D` of factor pairs can be recovered from only `r` supplied
traces

\[
 \mathrm{Tr}(S^j)\pmod{2^r},\qquad r\le j\le2r-1,
 \qquad r=\mathrm{bits}(m),
 \tag{2}
\]

where `m` is the matrix dimension. The largest numerical matrix power
is `2r-1`, not an exponential in `r`.

This is a complete reduction with a precise missing primitive. The
global traces in (2) have not been evaluated at a qualifying compressed
cost. The note also identifies explicit rank-two displacement generators
and a degree-three leading inverse model; neither fact is used to assume
that the missing traces are inexpensive.

The strongest complete factoring score remains unchanged.

Increasing the displacement also gives a simpler single-trace contract:
with sufficient padding, `Tr(M M^T)` alone is the factor count modulo
the required power of two. The low-degree Hermite interface applies
when the original displacement is retained.

## 1. An integral matrix with the desired unit subspace

Let

\[
 t=2^kw,\quad z=N+t,\quad w^2+w+1=0,\qquad
 M_{a,b}=\frac{t}{z-ab}.
 \tag{3}
\]

Use the existing square-domain and displacement contract. Thus every
nonzero integer difference `N-ab` has valuation below `k`; the sufficient
explicit condition is `2^k > max_{a,b}|N-ab|`. The matrix
`M` is integral, and modulo two its entries are one at `ab=N` and zero
elsewhere. Distinct positive nodes ensure at most one such entry in
each row or column. Therefore

\[
 S\bmod2=\mathrm{diag}(1_{\text{factor row}}),
 \qquad S^2-S\equiv0\pmod2.
 \tag{4}
\]

The transpose in (1) is the ordinary transpose, not quadratic
conjugation. The argument takes place over the unramified quadratic
integer ring and its quotients modulo powers of two. Merely assuming
distinct node residues would not suffice: without the stated valuation
bound, additional quadratic unit entries could survive modulo two, and
their ordinary transpose product need not be idempotent.

## 2. Padding reduces the contract to one trace

Suppose `B` is an upper bound for `v_2(N-ab)` over the nonzero
differences in the two domains. Such a bound is available from their
encoded endpoints: if the largest absolute difference is `U>0`, take
`B=bits(U)-1`. If every difference is zero, the count is immediate.
For a requested precision `r>=1`, choose

\[
 k\ge B+\lceil r/2\rceil.
\]

Every exact-match entry of `M` equals one. Every other entry has
valuation `k-v_2(N-ab)>=ceil(r/2)`. Therefore

\[
 \boxed{\mathrm{Tr}(MM^{\mathsf T})
       =\sum_{a,b}M_{a,b}^2\equiv D\pmod{2^r}.}
\]

The padding also ensures `k>B`, so no additional unit entries occur.
With the stronger choice `k>=B+r`, the entry sum itself satisfies

\[
 \sum_{a,b}M_{a,b}\equiv D\pmod{2^r}.
\]

These are positive signs for the convention `M=t/(z-ab)` in (3).
The required `k` remains polynomial in the endpoint bit lengths and
`r`. With `r=bits(m)`, either supplied residue recovers the exact count.
The unresolved operation is evaluating that singular joint-kernel
aggregation without enumerating its nodes. Displacement padding does
not supply an evaluator by itself.

## 3. Hermite projection uses degree `2r-1`

Use the integer polynomial

\[
 P_r(X)=X^r\sum_{i=0}^{r-1}{r+i-1\choose i}(1-X)^i.
 \tag{5}
\]

It is zero modulo `X^r` and one modulo `(1-X)^r`. Consequently

\[
 P_r(X)^2-P_r(X)
      \text{ is divisible by }X^r(X-1)^r.
\]

Substituting `S` and using (4) shows that

\[
 E_r=P_r(S)\pmod{2^r}
 \tag{6}
\]

is an idempotent matrix. Its reduction modulo two is the same diagonal
factor-row projector as in (4).

### Its trace is the ordinary integer rank

Here is an elementary local-ring argument. For the proof, order the
factor rows first. Form a change-of-basis matrix whose first `D` columns
are `E_r e_i` and whose remaining columns are `(I-E_r)e_j`. Its reduction
modulo two is the identity, so it is invertible over the finite local
ring. Idempotence gives

\[
 E_r C=C\,\mathrm{diag}(I_D,0).
\]

Thus `E_r` is similar to that diagonal idempotent, and

\[
 \mathrm{Tr}(E_r)=D\pmod{2^r}.
 \tag{7}
\]

The proof does not require knowing the factor rows computationally.
They are used only to establish the trace identity.

Since `0<=D<=m<2^r`, the residue in (7) determines the integer `D`
without ambiguity. Expanding (5) gives

\[
 P_r(X)=\sum_{j=r}^{2r-1}\gamma_jX^j,
 \qquad
 \gamma_{r+i}=(-1)^i{2r-1\choose r+i}{r+i-1\choose i}.
 \tag{8}
\]

The coefficients have `O(r)` bits. Taking the indicated integer linear
combination of the traces in (2) therefore recovers (7) with polynomial
additional bit cost.

The [supplied-trace interface](../tools/pole_cauchy_spectral_trace.py)
constructs the required degrees and coefficients and decodes valid
trace data. For dimension `2^200`, it requests 201 traces, with largest
numerical degree 401. This small query description does not evaluate
those traces.

### Comparison with repeated squaring

The power `S^(2^(r-1))` also agrees with `E_r` modulo `2^r`.
To see this, use the decomposition into image and kernel of `E_r`.
The restriction of `S` to the image is `I+2A`, and its restriction to
the kernel is `2B`. Repeated squaring makes the first block the identity
and the second block zero to the required precision. For `r=1`, (4)
already gives the assertion.

Repeated squaring requires only `r-1` explicit matrix multiplications,
but a symbolic power representation could otherwise hide exponential
numerical degree. The Hermite polynomial removes that degree issue:
only the linear range (2) is needed for a proposed trace evaluator.

## 4. Explicit rank-two displacement

Let `D_A` and `D_B` be the diagonal node matrices, let `1` be the all-one
vector, and let `beta_b=1/b`. Directly from (3),

\[
 D_A M-zMD_B^{-1}=-t\,1\beta^{\mathsf T}.
\]

Transposing and subtracting the two products gives, with `u=M beta`,

\[
 \boxed{D_A S-SD_A=t(u1^{\mathsf T}-1u^{\mathsf T}).}
 \tag{9}
\]

This is an explicit displacement of rank at most two. For powers,

\[
 [D_A,S^j]=\sum_{i=0}^{j-1}S^i[D_A,S]S^{j-1-i},
 \tag{10}
\]

so the displacement rank is at most `2j`.

Formula (10) is not a construction-time bound for its generators. It
requires vectors such as `S^i u` and `S^i 1`. Constructing those vectors
implicitly and evaluating their necessary scalar products remain part
of the computation. Explicit dense multiplication still depends on the
numerical matrix dimension.

## 5. The diagonal retains the global information

For the scalar function

\[
 u(a)=\sum_{b\in\mathcal B}\frac{t}{b(z-ab)},
\]

equation (9) determines off-diagonal entries as divided differences:

\[
 S_{a,a'}=t\frac{u(a)-u(a')}{a-a'}\quad(a\ne a'),
 \qquad S_{a,a}=t u'(a).
 \tag{11}
\]

Equivalently, if `L_B(Y)=P_B'(Y)/P_B(Y)` for
`P_B(Y)=prod_b(Y-b)`, and `H_B=sum_b 1/b`, then

\[
 u(a)=\frac tz\left(H_B+L_B(z/a)\right).
 \tag{12}
\]

A single value or derivative in (12) is an affine-product or
one-variable rational-sum calculation. The written degree of `P_B`,
however, is `m`. The general explicit-rational-weight compiler cannot
receive that polynomial at the cost of its short product description.

The displacement relation alone leaves the diagonal entirely free:
adding any diagonal matrix to `S` does not change its commutator with
`D_A`. The special formula (11) fixes that diagonal, but its trace is
the unresolved global sum

\[
 \mathrm{Tr}(S)=\sum_{a,b}\frac{t^2}{(z-ab)^2}.
 \tag{13}
\]

At sufficiently large displacement, (13) itself is the factor count
modulo the requested power of two. Thus the first matrix-power trace
already contains a joint product-pole aggregation problem.

There is also a direct warning inside the rank-two generator:

\[
 \sum_a\frac{u(a)}a
   =\frac tz\left(
       \left(\sum_a\frac1a\right)\left(\sum_b\frac1b\right)
       +\sum_{a,b}\frac1{z-ab}\right).
 \tag{14}
\]

The separated term is inexpensive, but the remaining term is the
original global pole trace. Fast pointwise generator evaluation does
not supply this scalar product for free.

## 6. A concrete leading inverse partial involution

There is additional exact structure on the symmetric power-of-two
domains

\[
 \mathcal A=\mathcal B=\{1,3,\ldots,2m-1\},\qquad
 m=2^s,
\]

with `k>s+1`. Put `Q=2m`, `L=k-s-1`, and let `f(a)` be the unique
representative of `N/a mod Q` in this domain. Then

\[
 K:=w\,2^L M^{-1}\pmod2
\]

is exactly the zero-one matrix with

\[
 K_{f(a),a}=1
 \quad\Longleftrightarrow\quad
 v_2(N-a f(a))=s+1,
 \tag{15}
\]

and every other entry zero. The zero difference is excluded from (15).

To prove this, every Vandermonde derivative on the domain has valuation
`2m-2-s`. This follows from the factorial expression for an arithmetic
progression and the complementary binary digits of `i` and `m-1-i`.
For each row, its unique reciprocal partner modulo `Q` has capped
valuation

\[
 d_a=\min(k,v_2(N-a f(a)))\ge s+1,
\]

using `d_a=k` for an exact factor pair. Every other cross difference has
valuation at most `s`. The explicit Cauchy inverse formula therefore
gives

\[
 v_2((M^{-1})_{f(a),a})=d_a-k,
\]

while an unmatched entry has valuation at least `s+2-k`. After scaling
by `2^L`, only the entries in (15) remain. Their normalized real unit
factors all reduce to one, and the sole quadratic unit is `w^(-1)`.
The extra factor `w` produces the zero-one values claimed.

The map `f` is an involution, and the retention condition is symmetric
in its two nodes. Hence

\[
 K^2=\mathrm{diag}(K1),\qquad K^3=K.
 \tag{16}
\]

The [node interface](../tools/pole_cauchy_spectral_trace.py) computes
one partner and its retention bit in polynomial time. Enumerating all
partners still costs `m` queries.

### Low-degree resolvents are accessible, but multiplicities remain

The minimal polynomial of `K` divides `T(T-1)^2`. Its resolvent has the
exact expression

\[
 (I-xK)^{-1}
    =I+\frac{x}{1-x^2}K+\frac{x^2}{1-x^2}K^2.
 \tag{17}
\]

Moreover,

\[
 u^{\mathsf T}K^2v=(u\mathbin{\circ}v)^{\mathsf T}K1.
\]

Thus a weighted bilinear resolvent reduces to two weighted `K`
contractions and an ordinary one-variable scalar product. For explicitly
written integral rational weights, those contractions are supplied by
the preceding inverse engines. The scalar product is a one-variable
rational sum: the valuation partition reduces it to polynomial power
sums on its leaf intervals, with all degrees charged.

If `rho` is the number of retained nodes, the characteristic polynomial
over the residue field is nevertheless

\[
 T^{m-\rho}(T-1)^\rho.
 \tag{18}
\]

The trace of the resolvent only sees multiplicities in characteristic
two; complete two-cycles disappear from that trace. A small minimal
polynomial and fast bilinear resolvents do not determine the integer
`rho`. Writing (18) with an unknown exponent does not evaluate it.

Furthermore, (15) selects a first nonzero inverse layer, not the original
unit Smith directions: true factor pairs have `d_a=k` and disappear
from (15) when `L>0`. Continuing through all layers must preserve the
restricted reciprocal sets and count their multiplicities. No bound
removing that dimension dependence is asserted here.

## 7. Exact controls

The [checker](../tools/check_pole_cauchy_spectral_trace.py) verifies 48
local-ring projections after nontrivial basis changes, including
agreement with the repeated-square exponent `2^(r-1)`. Thirty-six pole
matrices recover their literal factor-pair counts from the supplied
trace interface. Exact quadratic arithmetic checks 386 displacement
entries. A further 36 controls verify the minimally padded one-trace
identity and the more strongly padded entry-sum identity.

Twenty leading inverse matrices through dimension sixteen agree with
the predicted partial involution in all 1,364 entries. Their cubic
identities and resolvent formula are checked independently. The
[evidence](../evidence/pole-cauchy-spectral-trace.json) also records the
dimension-`2^200` trace contract, without pretending to evaluate its
201 required traces.

```sh
python3 tools/check_pole_cauchy_spectral_trace.py \
  --output evidence/pole-cauchy-spectral-trace.json
```
