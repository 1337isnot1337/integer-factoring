# Globally normalized Bezout pairs on circuit-isometric domains

## Constructive result and its boundary

The [balanced Bezout compiler](balanced-bezout-modular.md) constructs
complete normalized pairs on two separated arithmetic progressions.
Here both index sequences may be perturbed by arbitrary-degree integer
arithmetic circuits. The complete pair, exact optimal dyadic scale, and
short coefficient supports remain computable in polynomial bit time.
The numerical interval count and scale may be exponential in their
binary encodings.

Let `g_A,g_B` be supplied division-free integer polynomial circuits in
one variable. Define

$$
\phi(i)=i+2g_A(i),\qquad \psi(i)=i+2g_B(i),
$$

$$
A_i=a+D_A\phi(i),\qquad B_i=b+D_B\psi(i),
\qquad 0\le i<m. \tag{1}
$$

Assume `m,h,N,a,b,D_A,D_B` are positive integers, `h>=1`, and

$$
a\equiv b\equiv1\pmod{2^{h+1}},\qquad
v_2(D_A)=v_2(D_B)=h+1,\qquad
N\equiv1+2^h\pmod{2^{h+1}}. \tag{2}
$$

Write

$$
P(X)=\prod_{i<m}(X-A_i),\qquad
Q(X)=\prod_{i<m}(N-B_iX). \tag{3}
$$

There is a unique rational pair `U_0,V_0`, each of degree less than
`m`, with `U_0P+V_0Q=1`. Its exact least common dyadic clearing
exponent, even allowing higher-degree Bezout pairs, is

$$
\boxed{r=h(2m-1)-s_2(m-1).} \tag{4}
$$

Put `U=2^r U_0`, `V=2^r V_0`. Given numerical precision `p>=1`, the
algorithm returns **every coefficient** of these two polynomials modulo
`2^p`, represented as

$$
U(X)\equiv\sum_{j=0}^{J_U}u_j(X-1)^{m-1-j},\qquad
V(X)\equiv\sum_{j=0}^{J_V}v_j(X-1)^{m-1-j},
\qquad J_U,J_V<4p-1. \tag{5}
$$

All omitted coefficients vanish modulo `2^p`. Both leading coefficients
are odd. Runtime is polynomial in the written circuits, the binary
scalar parameters, and `p`. The formal circuit degrees, numerical
progression length, actual node heights, and numerical clearing exponent
are not expanded.

The root classes in (2) are already separated: every `A_i B_j` is one
modulo `2^(h+1)`, whereas `N` is `1+2^h`. Thus this algorithm evaluates
globally normalized Bezout data on a proved coprime class. It does not
identify an unknown factor pair, evaluate the unrestricted
[factor-count scalar](bezout-count-values.md), or improve the complete
factoring exponent `1/2`.

The new ingredient is the [degree-free circuit product and moment
engine](circuit-box-products.md). It supplies the short product jets
for these nonlinear nodes, including reciprocal odd units. The
normalization and terminal-window argument follows the earlier
[balanced denominator](balanced-bezout-denominator.md) and
[complete modular](balanced-bezout-modular.md) proofs; those sources
already establish the exact denominator for general ordered dyadic
isometries. This note gives the missing compiler for the stated circuit
class. The work is AI-assisted ordinary mathematics and code, without
new Lean formalization.

## 1. The circuit maps certify the isometry

For any integer polynomial `g` and distinct integers `i,j`,
`g(i)-g(j)` is divisible by `i-j`. Consequently

$$
\frac{(i+2g(i))-(j+2g(j))}{i-j}
=1+2\frac{g(i)-g(j)}{i-j}
$$

is odd. Both index maps in (1) therefore satisfy

$$
v_2(\phi(i)-\phi(j))=v_2(i-j),\qquad
v_2(\psi(i)-\psi(j))=v_2(i-j). \tag{6}
$$

This proof does not require the coefficients or values of `g` to be
small, and does not expand its circuit. It proves that all `A_i` and
all `B_i` are distinct. Each `B_i` is odd, hence nonzero. Signed circuit
coefficients are allowed; positivity of every generated node is not
needed for the Bezout statement.

Interleave the two root sets by

$$
x_{2i}=A_i,\qquad x_{2i+1}=N/B_i,\qquad
z_i=(x_i-1)/2^h. \tag{7}
$$

All `z_i` belong to `Z_2`. Their parities agree with their indices.
For two even indices, (2) and (6) give

$$
v_2(z_{2i}-z_{2j})=1+v_2(i-j)=v_2(2i-2j).
$$

For two odd indices, use that `N,B_i,B_j` are units:

$$
v_2\left(\frac{N/B_i-N/B_j}{2^h}\right)
=1+v_2(i-j).
$$

Opposite parities have difference valuation zero. Thus

$$
\boxed{v_2(z_i-z_j)=v_2(i-j)\quad(i\ne j).} \tag{8}
$$

The rational roots are distinct, including across the two classes,
so the pair (3) is coprime over `Q`.

## 2. The exact scale and terminal Newton window

Let `f` be the parity function on `Z_2`, and put `K=2m-1`. The
interpolating polynomial `F` with `F(A_i)=0` and `F(N/B_i)=1` has
Newton coefficients

$$
F(X)=\sum_{k=1}^{K}c_k\prod_{i<k}(X-x_i),\qquad
c_k=2^{-hk}f[z_0,\ldots,z_k]. \tag{9}
$$

The ordered-isometry divided-difference identity is

$$
v_2(f[z_0,\ldots,z_k])=k-1-v_2(k!). \tag{10}
$$

For completeness, the short unit formula in section 3 proves (10)
directly: its leading term is one and all later terms are even after
the stated normalization. Therefore the most negative coefficient
valuation in (9) is

$$
-r=-(h-1)K-1-v_2(K!).
$$

Since `s_2(2m-1)=s_2(m-1)+1`, this is exactly (4). All Newton basis
polynomials are integral, and the leading monomial coefficient has
valuation `-r`. Hence `F` has coefficient content `-r`.

The interpolation identities are `F=U_0P` and `1-F=V_0Q`.
Multiplication by `P` or `Q` preserves minimum coefficient valuation,
since their leading coefficients are units. Thus both `U_0` and `V_0`
have exact content `-r` and degree `m-1`. A hypothetical higher-degree
Bezout pair with a smaller clearing exponent can be reduced modulo
the monic polynomial `P`; the reduction preserves integrality and
would contradict the content of its unique remainder `V_0`. This
proves optimality among all polynomial pairs.

For the normalized Newton coefficient, define

$$
\delta_k=v_2(2^r c_k)
=(h-1)(K-k)+v_2(K!)-v_2(k!). \tag{11}
$$

A product of `K-k` consecutive positive integers contains at least
`floor((K-k)/2)` even factors. Therefore

$$
K-k\ge2p\quad\Longrightarrow\quad\delta_k\ge p. \tag{12}
$$

Only the last `2p` possible Newton coefficients survive modulo `2^p`.
The implementation computes their exact gaps by binary digit counts
and stops when the gaps reach `p`.

## 3. Huge-order coefficient units need only short jets

Define the odd unit

$$
D_k=(-1)^{k+1}\frac{k!}{2^{k-1}}f[z_0,\ldots,z_k]. \tag{13}
$$

The parity Mahler series is

$$
f(z)=\sum_{n\ge1}(-1)^{n+1}2^{n-1}\binom zn.
$$

It follows from the binomial expansion at nonnegative integers and
extends uniformly to `Z_2`, where every binomial polynomial takes
integral values. For an order-`k` divided difference, multiplying each
Lagrange term by `k!` makes it integral: (8) gives denominator valuation
`v_2(i!)+v_2((k-i)!)<=v_2(k!)`. Terms of degree below `k` vanish.
After the normalization (13), terms indexed by `n=k+j` have valuation
at least `j`.

Set

$$
E_n(T)=\prod_{0\le i<n}(1-iT),\qquad
G_k(T)=\prod_{i=0}^{k}(1-z_iT),
$$

$$
H_{k,j}=[T^j]\frac{E_{k+j}(T)}{G_k(T)},\qquad
B_{k,j}=\prod_{s=1}^j(k+s),\quad B_{k,0}=1.
$$

Expanding the falling factorial and using the complete-homogeneous
formula for monomial divided differences gives

$$
\boxed{D_k\equiv\sum_{j=0}^{p-1}(-1)^j2^j
\frac{H_{k,j}}{B_{k,j}}\pmod{2^p}.} \tag{14}
$$

The same Lagrange argument proves that each displayed quotient is
integral. The term `j=0` is one, proving that `D_k` is odd and hence
proving (10).

Retain the additional precision

$$
g_k=\max_{0\le j<p}\bigl(v_2(B_{k,j})-j\bigr)
=s_2(k)-\min_{0\le j<p}s_2(k+j). \tag{15}
$$

The `j=0` term makes this guard nonnegative. The equality uses
`v_2(B_(k,j))=j+s_2(k)-s_2(k+j)`. It charges the weight `2^j` before
choosing precision and is sharper than retaining the entire maximum
denominator valuation. At precision `W=p+g_k`, the numerator divided
by `2^v_2(B_(k,j))` is known to at least `p-j>=1` bits. Invert the
remaining odd denominator and multiply by `2^j`; all `p` final bits
are correct. The retained numerical precision is polynomial in `p`
and the binary length of `k`, even when `k` is enormous.

The odd part `u_k` of `k!` is computable without forming a factorial.
Split `1,...,k` into its exact-valuation classes; each odd part is an
odd arithmetic progression. There are `O(log(k+1))` such classes,
and the circuit product engine computes each unit product. Finally,

$$
2^r c_k=(-1)^{k+1}2^{\delta_k}D_k u_k^{-1}. \tag{16}
$$

Moving from `k` to `k-1` updates `u_k` by division by the odd part of
`k`. Only one complete factorial-unit calculation is needed.

## 4. Circuit product jets close the nonlinear-node cost

The remaining data in (14) are short jets of `G_k`. Its even-index
nodes are the integer circuit values

$$
z_{2i}=\frac{a-1}{2^h}+\frac{D_A}{2^h}\phi(i). \tag{17}
$$

Its odd-index nodes are

$$
z_{2i+1}=
\frac{(N-b)/2^h-(D_B/2^h)\psi(i)}
{b+D_B\psi(i)}. \tag{18}
$$

All displayed divisions by `2^h` are exact divisions of written
integer constants. The remaining denominator in (18) is always odd.
At requested coefficient precision `W>=1`, replace its inverse by

$$
y^{-1}\equiv y^{2^{W-1}-1}\pmod{2^W}
\quad\text{for odd }y. \tag{19}
$$

For `W=1`, the exponent is zero, which correctly gives one. Ordinary
binary powering writes (19) with `O(W)` multiplication nodes. Its
formal degree can be exponential in `W`, but the degree-free product
compiler charges only the written circuit. Substitution into (18)
therefore gives the correct rational root value modulo `2^W` at every
index, without materializing either the exact root or its denominator.

For a root circuit `q` and an interval of `M` indices, reverse its two
global packets:

$$
\prod_i(1-q(i)T)
\equiv\prod_{c=0}^1
\sum_{j=0}^{J_c}(-2)^j e_{c,j}T^j(1-cT)^{m_c-j}
\pmod{2^W}. \tag{20}
$$

To obtain a jet through degree `d`, take
`J_c=min(m_c,W-1,d)`. Terms at larger `j` cannot contribute. Each
coefficient of `(1-T)^(m_c-j)` through degree `d` is a small-index
binomial coefficient with an encoded large upper argument. It has
polynomial construction cost. The underlying normalized moments use
only `max(1,W-1)` bits by the weighted Newton precision theorem.

Apply (20) to (17) and (18) on their respective prefix counts, then
multiply their two jets. The same routine computes `E_n` using `q(i)=i`.
This proves the complete polynomial construction cost of (14),
including the nonlinear and reciprocal circuits. No numerical index
interval is scanned.

The leading coefficient of `Q` is also needed. It is
`prod_i(-B_i)`, obtained by evaluating the global product of the
integer circuit `B_i` at zero. Every `B_i` is a unit, so this product
has no singular normalization loss.

## 5. The complete output has short support

Put `Y=X-1` and `H=2^r F`. Every `x_i-1` is divisible by `2^h`.
The coefficient at distance `j` from the top of each Newton basis
product has valuation at least `hj`. Combining this with (12), the
reverse polynomial

$$
H^*(T)=T^K H(1+1/T)
$$

has no nonzero coefficient modulo `2^p` beyond

$$
J_H=2p-1+\left\lfloor\frac{p-1}{h}\right\rfloor. \tag{21}
$$

The implementation uses the last actually retained Newton order,
which can make this support smaller. It constructs one initial
root-product jet, then removes the at most `2p-1` boundary linear
factors as the order decreases. These are explicit boundary queries,
not an enumeration of either node set.

For

$$
P^*(T)=T^mP(1+1/T),\qquad Q^*(T)=T^mQ(1+1/T),
$$

the roots in (2) imply

$$
v_2([T^j](P^*)^{\pm1})\ge(h+1)j,\qquad
v_2([T^j](Q^*)^{\pm1})\ge hj. \tag{22}
$$

Both constant coefficients are units; their inverse-series recurrences
preserve these bounds. The exact interpolation identities give

$$
U^*=H^*/P^*\pmod{T^m},\qquad
V^*=-H^*/Q^*\pmod{T^m}. \tag{23}
$$

The reversed term `2^r T^K` from `2^r-H=VQ` lies outside the retained
range `T^m`. Equations (21)--(23) give

$$
J_U\le J_H+\left\lfloor\frac{p-1}{h+1}\right\rfloor,
\qquad
J_V\le J_H+\left\lfloor\frac{p-1}{h}\right\rfloor.
$$

Each resulting list has fewer than `4p` coefficients, including any
retained zeros. No other coefficient can survive. This proves the
global output statement (5), not merely a local value or Taylor jet.

## 6. Complete cost and remaining obligation

All scalar encodings, count arithmetic, isometry checks, and the
formula (4) have polynomial bit cost. There are fewer than `2p`
Newton coefficients. For each one, (15) allocates polynomially many
precision bits and requests only a degree-`p-1` product jet. Unit
inversion compilation adds `O(W)` nodes, and the univariate quotient
used by the product engine has degree at most `2W`. Its intermediate
coefficient arrays therefore stay polynomial in the written inputs
and numerical precision.

The exact integer values `g_A(i),g_B(i)` may have exponential bit
length. The algorithm evaluates their circuits only in finite
coefficient quotients or at a bounded number of modular boundary
indices. It never writes those exact values. The numerical exponent
`r` is returned in binary and is never used as an allocation size.

This constructs a larger class of complete normalized Bezout pairs
than the arithmetic progression compiler. The separation hypothesis
in (2) remains part of the theorem. For the unrestricted factoring
domains, a pair of roots may meet or remain close through unknown
levels; the uniform ordered-isometry certificate and its explicit
parity labels are then absent. No step here decides those labels or
the resulting factor support.

## 7. Exact implementation controls

Files:

- [Compiler](../tools/circuit_box_bezout.py).
- [Independent checker](../tools/check_circuit_box_bezout.py).
- [Receipt](../evidence/circuit-box-bezout.json).

The checker builds 60 complete pairs independently by rational
Lagrange interpolation and exact polynomial division. It checks all
1,140 interleaved distances, every optimal scale, 180 complete
normalized pairs at three precisions, all 1,080 coefficients, and
462 normalized divided-difference units. Signed nonlinear examples
are included. Separate controls compare 45 reversed product jets and
48 normalized factorial units with literal finite calculations.

An actual large-domain computation uses `m=2^80+17`, precision six,
and two 204-node perturbation circuits representing `i^(2^200+1)`.
It returns nonzero complete normalized pairs, with 11 and 14 retained
coefficient positions, at an encoded optimal exponent exceeding
`2^81`. Its maximum numerical working precision is six bits and it
uses three explicit boundary queries.

The large control compares two output precisions, checks its factorial
unit by an independent residue-bucket calculation, and evaluates the
complete Bezout identity at eight quadratic unit points. The latter
reference enumerates only 64 residue classes of the index, raises the
corresponding factors to their exact encoded multiplicities, and
obtains nonzero values for both output polynomials. It does not expand
the original domains or independently reconstruct their rational
coefficient arrays; that all-size claim rests on the proof above.

A second large control uses `m=2^80` with the same nonlinear circuits.
It exercises the large factorial-ratio loss in (15): the stated
algorithm uses up to 86 working bits for six output bits. Charging the
coefficient weights reduces the earlier unweighted guard from 90 bits.
Eight Newton terms return the complete
pair. Independent quadratic evaluations and agreement with a separate
four-bit compilation check this case as well.

```sh
python3 tools/check_circuit_box_bezout.py --output /tmp/circuit-box-bezout.json
cmp /tmp/circuit-box-bezout.json evidence/circuit-box-bezout.json
```

The complete factoring score remains `1/2`.
