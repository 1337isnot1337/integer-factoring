# Every invariant factor of an arithmetic Cauchy power

## Result and scope

Let `p` be a prime, `T` a `p`-adic unit, and `d` nonzero with `e=v_p(d)>=1`.
For every positive integer `m`, define

\[
H_n^{(m)}(T,d)=\bigl((T+d(i+j))^{-m}\bigr)_{0\le i,j<n}.
\]

Write `(m)_r=m(m+1)...(m+r-1)`, with `(m)_0=1`. The complete Smith order list is

\[
\boxed{s_r=2re+v_p(r!)+v_p((m)_r),\qquad 0\le r<n.}
\tag{1}
\]

The same orders hold for `((T+d(i-j))^{-m})`: reverse its columns and replace
`T` by the still-unit value `T-d(n-1)`. The proof also works over an unramified
finite extension of `Z_p`. No matrix-size or power-of-two restriction is imposed.

For `m=1`, this is the arithmetic Cauchy profile. For `m=2`, it recovers the
[squared-Cauchy theorem](arithmetic-cauchy-factorial-law.md), since
`v_p((2)_r)=v_p(r!)+v_p(r+1)`. The proof uses Cauchy--Binet and a
factorial-normalized alternant bound, without a permanent identity for higher powers.

This supplies exact spectral data for additional kernels. It does not by itself
give a scalar factoring decoder, a cheaper grouped operator, or an exponent above
`1/2`. Those interfaces remain necessary. The proof and finite controls used
OpenAI assistant support.

## 1. A factorial-normalized finite-difference bound

Let `delta` have valuation at least one, and put

\[
z_i=\frac{i}{1+\delta i},\qquad f_r(x)=\left(\frac{x}{1+\delta x}\right)^r.
\]

For the forward difference operator `Delta`, all nonnegative integers `i,r` satisfy

\[
v_p\left(\frac{\Delta^i f_r(0)}{r!}\right)\ge-\frac{r-i}{p-1}.
\tag{2}
\]

Zero has valuation infinity; nonintegral bounds mean the appropriate integer ceiling.
The Stirling generating function gives

\[
\frac{i!S(k,i)}{k!}=\operatorname{coeff}_{x^k}(e^x-1)^i.
\]

Each coefficient summand is `1/(j_1!...j_i!)`, for positive `j_l` summing to `k`.
The factorial bound `v_p(j!)<=(j-1)/(p-1)` therefore gives

\[
v_p\left(\frac{i!S(k,i)}{k!}\right)\ge-\frac{k-i}{p-1}.
\tag{3}
\]

The cases `i=0` and `r=0` are direct. For `r>=1`, expand

\[
f_r(x)=\sum_{k\ge r}(-\delta)^{k-r}\binom{k-1}{r-1}x^k.
\tag{4}
\]

Using `Delta^i x^k|_(x=0)=i!S(k,i)`, a nonzero summand of
`Delta^i f_r(0)/r!` has valuation at least

\[
(k-r)e+v_p(k!/r!)-\frac{k-i}{p-1}\ge-\frac{r-i}{p-1},
\]

because `k>=r` and `e=v_p(delta)>=1`. This proves (2). Convergence is separate:
before applying (3), each summand has valuation at least `(k-r)e-v_p(r!)`,
which tends to infinity.

The generating function is [DLMF 26.8.12](https://dlmf.nist.gov/26.8.E12).
Its displayed composition expansion supplies the entire valuation argument.

## 2. The arithmetic alternant bound

Take increasing nonnegative integers `r_0<...<r_(n-1)`, and define

\[
K=\sum_j(r_j-j),\quad D=\sum_j\bigl(v_p(r_j!)-v_p(j!)\bigr),\quad
B=\left\lfloor\frac K{p-1}\right\rfloor.
\tag{5}
\]

Both `K,D` are nonnegative because `r_j>=j`. Write

\[
S_R(z)=\frac{\det(z_i^{r_j})_{i,j}}{\prod_{i<j}(z_j-z_i)}.
\tag{6}
\]

This Schur polynomial has integer coefficients and degree `K`; only these
properties are used. Polynomiality also follows by dividing the alternating
integer polynomial by the distinct linear factors `z_j-z_i`.
At `z_i=i/(1+delta i)`,

\[
\prod_{i<j}(z_j-z_i)=\left(\prod_{j=0}^{n-1}j!\right)
\prod_{i=0}^{n-1}(1+\delta i)^{-(n-1)}.
\tag{7}
\]

The last product is a unit. Taking successive forward differences of rows
is a unit-triangular operation, so

\[
S_R(z)\frac{\prod_j j!}{\prod_j r_j!}
=u\det\left(\frac{\Delta^i f_{r_j}(0)}{r_j!}\right)_{i,j},
\tag{8}
\]

where `u` is a unit. In each determinant term, the bounds in (2) sum to
`-K/(p-1)`, independently of the permutation. Thus `v_p(S_R(z))>=D-B`.
The arguments are integral and `S_R` has integer coefficients, giving also
nonnegative valuation. Together,

\[
\boxed{v_p(S_R(z))\ge\max(0,D-B).}
\tag{9}
\]

The arithmetic arrangement of nodes matters; (9) is not asserted for arbitrary
nodes in `Z_p`.

## 3. A unique leading Cauchy--Binet term

Normalize `delta=d/T`; the scalar `T^(-m)` changes no Smith orders. Put
`a_i=delta*i/(1+delta*i)=delta*z_i`. The identity

\[
1+\delta(i+j)=(1+\delta i)(1+\delta j)(1-a_i a_j)
\tag{10}
\]

makes the normalized matrix unit-diagonally equivalent to

\[
K_{ij}=(1-a_i a_j)^{-m}=\sum_{r\ge0}c_r a_i^r a_j^r,
\qquad c_r=\binom{m+r-1}{r}=\frac{(m)_r}{r!}.
\tag{11}
\]

The coefficients `c_r` are integers and `a_i` are divisible by `p`. Finite
Cauchy--Binet on partial sums, followed by the convergent limit, gives

\[
\det K=\prod_{i<j}(a_j-a_i)^2
\sum_R\left(\prod_jc_{r_j}\right)S_R(a)^2.
\tag{12}
\]

The leading set `R_0=(0,1,...,n-1)` has `S_(R_0)=1`. Every other set has
`K>=1`. Homogeneity and (9) give

\[
v_p(S_R(a))\ge eK+\max(0,D-B).
\tag{13}
\]

Also,

\[
\frac{\prod_jc_{r_j}}{\prod_jc_j}
=\frac{\prod_j\prod_{h=j}^{r_j-1}(m+h)}{\prod_j(r_j!/j!)}
\]

has valuation at least `-D`. Each nonleading term in (12), relative to its
leading term, therefore has valuation at least

\[
2eK+2\max(0,D-B)-D\ge2eK-B\ge K\ge1.
\tag{14}
\]

The first inequality holds both for `D<=B` and `D>=B`. The bound tends to
infinity with `K`, and only finitely many exponent sets have any given `K`.
This proves convergence of the normalized sum and rules out leading cancellation.
Consequently

\[
v_p(\det H_n^{(m)})=n(n-1)e+
\sum_{r=0}^{n-1}\left[v_p(r!)+v_p((m)_r)\right].
\tag{15}
\]

The additional diagonal and Vandermonde unit factors are all `1` modulo `p`.
Thus one also obtains

\[
\boxed{\frac{\det H_n^{(m)}(T,d)}
{d^{n(n-1)}T^{-mn-n(n-1)}\prod_{r=0}^{n-1}r!(m)_r}\equiv1\pmod p.}
\tag{16}
\]

## 4. Integral orthogonal polynomials give every Smith order

Repeat the integral-basis argument of the squared-Cauchy theorem for the moment
functional `L(x^k)=(T+dk)^(-m)`. Its leading Hankel determinants are nonzero by
(15). Gram--Schmidt over the fraction field gives monic orthogonal polynomials `P_r` and norms

\[
h_r=L(P_r^2)=\frac{\det H_{r+1}^{(m)}(T,d)}{\det H_r^{(m)}(T,d)}.
\tag{17}
\]

The constant-coefficient determinant formula is

\[
P_r(0)=(-1)^r\frac{\det H_r^{(m)}(T+d,d)}{\det H_r^{(m)}(T,d)}.
\tag{18}
\]

Since `T+d` remains a unit, (15) makes every `P_r(0)` a unit. Equation (17)
gives `v_p(h_r)=s_r`, and

\[
s_r-s_{r-1}=2e+v_p(r)+v_p(m+r-1)>0\qquad(r\ge1).
\tag{19}
\]

Orthogonality gives the monic three-term recurrence

\[
P_{r+1}(x)=(x-b_r)P_r(x)-c_r'P_{r-1}(x),\qquad c_r'=h_r/h_{r-1}.
\tag{20}
\]

Pairing with earlier orthogonal polynomials removes all lower-degree terms;
pairing with `P_(r-1)` gives the stated coefficient. By (19), `c_r'` is integral.
Evaluating at zero and using the units in (18) proves `b_r` integral as well.
The initial `P_0=1` and `P_1=x-L(x)/L(1)` are integral. Induction proves that
every `P_r` has integral coefficients.

This monic change of basis is unit triangular over `Z_p` and diagonalizes the
moment Gram matrix to `diag(h_0,...,h_(n-1))`. Their valuations increase by (19),
proving the individual Smith orders (1), rather than only their sum.

## 5. A tunable first gap, with its limits

For `m=p^a`,

\[
v_p((p^a)_r)=a+v_p((r-1)!),\qquad 1\le r\le p^a.
\]

Hence

\[
s_0=0,\qquad s_r=a+2re+v_p(r!)+v_p((r-1)!)\quad(1\le r\le p^a).
\tag{21}
\]

This gives an arbitrarily large gap after the first unit invariant factor,
followed by an explicit full profile. A factoring-network substitution still
needs a valid realization, scalar recovery, normalization, and grouped cost.
No improved evaluator is asserted here.

## 6. Exact controls

Run

```sh
python3 tools/check_arithmetic_cauchy_powers.py \
  --output evidence/arithmetic-cauchy-powers.json
```

The checker verifies finite-difference and alternant bounds with exact rational
arithmetic, checks (16) independently, compares complete modular Smith profiles
with (1), and checks integral orthogonal bases. The
[recorded evidence](../evidence/arithmetic-cauchy-powers.json) separates these
finite controls from the all-size proof above.
