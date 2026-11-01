# Global products of unrestricted-degree circuit values

## Constructive result and scope

Fix a positive number of variables `d`. Let an integer polynomial
`q(X_1,...,X_d)` be supplied by a division-free arithmetic circuit, using
integer constants, addition, and multiplication. Its formal degree need
not be bounded or expanded. Let

$$
\mathcal B=\prod_{i=1}^d [L_i,U_i)\cap\mathbb Z^d,
\qquad 0\le L_i\le U_i,
\qquad M=|\mathcal B|.
$$

Given a numerical precision `p>=1`, a uniform deterministic algorithm
constructs an ordinary addition/multiplication circuit for

$$
R_q(Z)=\prod_{a\in\mathcal B}(Z-q(a))
\quad\text{in }(\mathbb Z/2^p)[Z]. \tag{1}
$$

The output has `O(p+1+log(M+1))` nodes. Every output constant has at most
`p` bits, and every large power is expanded into binary-powering nodes.
For fixed `d`, generation has bit cost polynomial in the written input
circuit, the endpoint bit lengths, and `p`. Neither the numerical box
cardinality nor the potentially enormous integer values of `q` are
enumerated or materialized.

This extends the [bilinear packet construction](bilinear-pole-packets.md),
which obtains its moments by expanding bounded-degree polynomial powers.
The present construction removes that degree restriction. Its universal
integer-function quotient is different from the functional quotient in
the [Cauchy circuit inverse moments](pole-cauchy-circuit-moments.md): here
the quotient relations vanish pointwise at every integer input.

The result does not normalize the product at a singular integer center,
compute a rational gcd, or improve the complete factoring exponent.
In particular, `q(u,v)=uv` already has small degree; this extension does
not remove that case's existing singular-normalization problem. The
complete score remains `1/2`. This is an ordinary proof and implementation,
developed with AI assistance from repository base
`941cc1d8f56e6e80ec45938f1d660532de2f5046`; no new Lean proof is claimed.

## 1. A small quotient preserves every integer value

For an integer precision `P>=1`, let

$$
D=\min\{j\ge0:v_2(j!)\ge P\},
\qquad F_D(X)=\prod_{i=0}^{D-1}(X-i). \tag{2}
$$

The identity `v_2(j!)=j-s_2(j)` computes `D` with `O(P)` small integer
steps. In particular `D<=2P`: the `P` even factors of `(2P)!` already
contribute valuation at least `P`.

For every integer `x`,

$$
F_D(x)=D!\binom{x}{D}\equiv0\pmod{2^P}. \tag{3}
$$

The binomial coefficient is an integer also for negative `x`. Thus every
polynomial reduction in

$$
\mathcal A_P=
(\mathbb Z/2^P)[X_1,\ldots,X_d]/
(F_D(X_1),\ldots,F_D(X_d)) \tag{4}
$$

preserves its value modulo `2^P` at every integer tuple. Since every
defining polynomial is monic, the quotient has the explicit basis

$$
X_1^{e_1}\cdots X_d^{e_d},\qquad 0\le e_i<D. \tag{5}
$$

There are at most `(2P)^d` coefficients, each of `P` bits. Ordinary
polynomial multiplication followed by monic reduction computes every
circuit operation in polynomial time for fixed `d`. No factorization or
division by a nonunit is involved.

This is a statement about integer **values**. Reducing `q` in (4) does
not assert an identity for `q` in the unrestricted polynomial ring. The
global polynomial identity (1) will follow separately from exact
elementary-symmetric identities and controlled modular precision.

## 2. Normalize parity without dividing in the quotient

Split the box into at most `2^d` input parity boxes. On the box indexed
by `e=(e_1,...,e_d)` in `{0,1}^d`, write

$$
a_i=e_i+2x_i.
$$

The `x_i` again range over binary-described nonnegative integer
intervals. The output parity `c=q(e) mod 2` is constant on this box.
There is an integer polynomial `r_e` with

$$
q(e+2X)=c+2r_e(X). \tag{6}
$$

A circuit for `r_e` has only a constant-factor increase in node count.
At every gate retain a pair `(c,r)` denoting `c+2r`, with `c` equal to
zero or one. A variable gives `(e_i,X_i)` and a constant `z` gives
`(z mod 2,(z-(z mod 2))/2)`. If two inputs have pairs `(c,u),(f,v)`,
their sum gives

$$
\left((c+f)\bmod2,\quad
u+v+\left\lfloor\frac{c+f}{2}\right\rfloor\right), \tag{7}
$$

and their product gives

$$
(cf,\;cv+fu+2uv). \tag{8}
$$

These are exact integer polynomial identities. Evaluate the resulting
normalized circuit directly in (4). In particular, the algorithm never
tries to divide an arbitrary quotient element by two.

## 3. Complete power sums have polynomial cost

For `c=0,1`, put

$$
S_c=\{a\in\mathcal B:q(a)\equiv c\pmod2\},
\quad m_c=|S_c|,
\quad s_{c,j}=\sum_{a\in S_c}\left(\frac{q(a)-c}{2}\right)^j.
\tag{9}
$$

The counts `m_c` are exact sums of the parity-box cardinalities. Compute
successive powers of the reduced `r_e` in (4). Their coefficients have
the form (5), so summation separates into one-variable power sums:

$$
\sum_{x\in I_1\times\cdots\times I_d}
X_1^{e_1}\cdots X_d^{e_d}
=\prod_{i=1}^d\sum_{x_i\in I_i}x_i^{e_i}. \tag{10}
$$

All exponents here are less than `D`, even if the original circuit has
exponential formal degree. Exact interval power sums through degree
`D-1` use the triangular identity

$$
t^{j+1}=\sum_{i=0}^{j}\binom{j+1}{i}
\sum_{x=0}^{t-1}x^i. \tag{11}
$$

Solve for the degree-`j` sum by exact integer division by `j+1`, then
subtract two prefix sums. With `b`-bit endpoints, every integer in this
calculation has `O(bD+D log D)` bits. Reducing these sums modulo `2^P`
and applying (10) gives every `s_{c,j} mod 2^P` through the requested
degree. The number of terms in the original box causes no precision
loss: a sum of errors divisible by `2^P` is still divisible by `2^P`.

## 4. The coefficient weights remove the Newton guard

Set

$$
J=\min(M,p-1),\qquad P=\max(1,p-1). \tag{12}
$$

Use the preceding algorithm to compute the moments through degree `J`
modulo `2^P`. Let `e_{c,j}` be the elementary symmetric
functions of the normalized roots in (9). Newton's identities say

$$
e_{c,0}=1,\qquad
j e_{c,j}=\sum_{i=1}^j(-1)^{i-1}e_{c,j-i}s_{c,i}. \tag{13}
$$

At stage `j`, retain `e_{c,j}` modulo `2^(P-v_2(j!))`. Every earlier
coefficient on the right is known to at least
`P-v_2((j-1)!)` bits. Hence the numerator is known to that precision.
It is divisible by `2^v_2(j)` by the exact identity (13). Divide its
canonical residue by that power of two and invert the odd part of `j`.
The result has precisely the stated retained precision. For `1<=j<=J`
it remains positive, since

$$
P-v_2(j!)=P-j+s_2(j)\ge s_2(j)\ge1.
$$

The product does not need `e_(c,j)` itself to `p` bits. It needs the
weighted coefficient `(-2)^j e_(c,j)`. Multiplication by that weight
raises the error valuation to

$$
P+j-v_2(j!)=P+s_2(j)\ge p. \tag{13a}
$$

Thus `p-1` normalized moment bits suffice for every positive packet
coefficient when `p>=2`. The constant coefficient is exactly one. At
`p=1` there are no positive packet coefficients and only the exact
root-parity counts matter; the implementation uses its harmless minimum
working precision of one.

This uniform moment precision cannot be reduced: for the singleton
even root `2r`, the coefficient of `Z^0` in `Z-2r mod 2^p` determines
`r mod 2^(p-1)`. The construction therefore has the sharp worst-case
moment precision, rather than a Newton factorial guard. The generic
guarded elementary-function helper is retained separately for checks.

No step assumes that an even integer has an inverse modulo a power of
two. The implementation keeps each changing precision explicitly and
performs the final weighted multiplication only afterward.

## 5. A global polynomial circuit

For `J_c=min(m_c,p-1)`, ordinary expansion gives the polynomial identity

$$
R_q(Z)\equiv
\prod_{c=0}^1
(Z-c)^{m_c-J_c}
\sum_{j=0}^{J_c}(-2)^j e_{c,j}(Z-c)^{J_c-j}
\pmod{2^p}. \tag{14}
$$

Every omitted coefficient has a factor `2^j` with `j>=p`. Empty packets
contribute one. Repeated values of `q` are allowed and are counted with
their multiplicities.

Each short sum in (14) uses `O(p)` ordinary nodes by Horner evaluation.
The remaining powers use `O(1+log(M+1))` nodes by binary powering.
Together with the constants and two packet multiplications, this proves
the claimed output size. Every computed constant is reduced modulo
`2^p` before writing the output.

Equation (14) is equality of every coefficient in
`(Z/2^p)[Z]`, so its circuit remains valid under substitution into any
commutative algebra over `Z/2^p`. The intermediate integer-function
quotient was used to compute the moments, not to weaken this final
polynomial identity to a pointwise statement.

For a written circuit of `S` nodes, the generator performs `2^d` circuit
evaluations and at most `p` quotient multiplications per parity box.
Naive quotient arithmetic uses a polynomial number of operations in
`D^d`, with `D<=2p`. Prefix sums, count arithmetic, and (13) also have
polynomial bit cost in `S,b,p` for fixed `d`. The constants of the input
circuit are read at their actual written lengths and reduced during
evaluation. Large hidden integer values at intermediate gates are not
constructed.

The dependence on `d` is not claimed polynomial: the explicit basis and
parity splitting cost up to `(2p)^d` and `2^d`, respectively. Encoding a
general many-variable Boolean problem does not evade these costs.

## 6. Unit traces retain the same extension

Suppose `Z-c` is a unit for each nonempty root-parity packet in a
commutative `Z/2^p`-algebra. Then

$$
\sum_{a\in\mathcal B}\frac1{Z-q(a)}
\equiv
\sum_{c=0}^1\sum_{j=0}^{p-1}
\frac{2^j s_{c,j}}{(Z-c)^{j+1}}
\pmod{2^p}. \tag{15}
$$

This is a finite geometric expansion with an error divisible by `2^p`
at every node. Use the exact count `m_c` for the term `j=0`. For every
`j>=1`, the factor `2^j` makes moments modulo `2^(p-1)` sufficient;
the same working precision (12) therefore applies. For example, both
denominators are units at `Z=w`, where
`w^2+w+1=0` in the unramified quadratic extension. Truncated jets at
`w` are also allowed, charging their numerical jet order.

At the factor-count center `N+2^k w`, one root-parity denominator is
nonunit. Neither (15) nor raw division by (14) is justified there. The
checker requires unit inverses and rejects a singular center. The
existing exponentially large content examples therefore still apply.

## 7. Implementation and exact controls

The [constructor](../tools/circuit_box_products.py) accepts `var`,
`const`, `add`, and `mul` input nodes. It writes the same ordinary
`x`, `const`, `add`, and `mul` product-node format as the earlier packet
constructors. The input and output indices are explicit, and invalid
forward references are rejected.

Run the [independent checker](../tools/check_circuit_box_products.py):

```sh
python3 tools/check_circuit_box_products.py --output /tmp/circuit-box-products.json
cmp /tmp/circuit-box-products.json evidence/circuit-box-products.json
```

The [retained receipt](../evidence/circuit-box-products.json) records
2,000 literal null-polynomial tests, 216 quotient-product point checks,
216 generic guarded Newton checks, 216 weighted packet checks at the
reduced precision against ordinary elementary expansion, and
76 complete polynomial-coefficient comparisons in one through three
variables. The same small controls compare normalized moments and exact
quadratic unit traces with literal sums. Their role is to check the
implementation; the all-size statement is proved above.

A separate large control uses a 161-node circuit with formal degree
greater than `2^81` and a box with more than `2^204` integer points. At
four-bit output precision the constructor uses three moment bits and a
quotient of dimension 16, and writes 448 ordinary output nodes. An
independent reference counts only 256 pairs of residue classes modulo 16,
then expands binomial
multiplicities with exact integer coefficients. It agrees with both
complete packet coefficient lists and exact packet multiplicities.
The original box is not enumerated. The control establishes no
singular-normalization or factoring claim.
