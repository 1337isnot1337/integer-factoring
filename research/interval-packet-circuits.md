# Small global modular circuits for binary-described interval products

## Constructive theorem

Given positive interval endpoints `L<=U`, an integer `N`, and precision
`p>=1`, there is a deterministic polynomial-bit-time algorithm that
writes ordinary arithmetic circuits for

```math
P_{[L,U)}(X)=\prod_{a=L}^{U-1}(X-a),\qquad
Q_{N,[L,U)}(X)=\prod_{a=L}^{U-1}(N-aX)
\quad\text{in }(\mathbb Z/2^p)[X].                \tag{1}
```

If the integer inputs have at most `n` bits, each circuit has `O(n+p)`
nodes. Its operations are addition and multiplication; each constant
has at most `p` bits. Generation, including every constant, has cost
polynomial in `n,p`. Exponentiation is expanded into ordinary circuit
nodes by binary powering. No unit or parity condition on `N` is needed.

This is **global polynomial equality** in the displayed ring, so the
same circuit is valid under substitution into every commutative
`Z/2^p`-algebra. In particular it extends the earlier
[unit-point jet evaluator](interval-unit-jets.md) to arbitrary points
and jets, including nonunit points.

The inputs to the [orbit-gcd question](interval-orbit-gcd.md) therefore
have inexpensive global modular circuits. Recovering their common
factor over `Q` still requires an additional elimination algorithm.
The strongest complete factoring score remains `1/2`.

## 1. Opposite-parity pairs share a quadratic base

First suppose `U-L=2h`. Pair `a=L+t` with
`a'=L+2h-1-t`, for `0<=t<h`. Their sum

```math
S=2L+2h-1
```

is odd, so their product is even. Put
`c_t=(L+t)(L+2h-1-t)/2`, an integer. The paired factors are

```math
(X-a)(X-a')=Y_P+2c_t,\qquad Y_P=X^2-SX,
```

and

```math
(N-aX)(N-a'X)=Y_Q+2c_tX^2,\qquad Y_Q=N^2-SNX.
```

Let `e_j` be the elementary symmetric functions of `c_0,...,c_(h-1)`.
Expanding in these pair contributions and reducing modulo `2^p` gives

```math
P=\sum_{j=0}^{J}2^j e_jY_P^{h-j},\qquad
Q=\sum_{j=0}^{J}2^j e_jX^{2j}Y_Q^{h-j},
\qquad J=\min(h,p-1).                            \tag{2}
```

Every omitted coefficient is divisible by `2^p`. This statement
concerns polynomial coefficients; it uses no equality test on a
finite set of values. If the interval length is odd, apply (2) to
all but its last element and multiply the last linear factor.
The empty interval gives one.

## 2. All needed symmetric functions have polynomial cost

The integer-valued quadratic sequence is

```math
c_t=\frac{L(L+2h-1)+(2h-1)t-t^2}{2}.
```

For each `j<=J`, expand the numerator's `j`th power as a polynomial
of degree `2j`, sum its monomials using the exact prefix power sums,
and divide the resulting integer by `2^j`. This computes
`s_j=sum_(t<h)c_t^j` without visiting any `t`. The division is exact
because every `c_t` is an integer. Newton's identities then give

```math
e_0=1,\qquad
e_k=\frac1k\sum_{j=1}^k(-1)^{j-1}e_{k-j}s_j.
```

These too are exact integer divisions, before modular reduction.
Since `h` and each endpoint have `O(n)` bits, `s_j` and `e_j` have
`O(np)` bits for `j<p`. Computing the power sums, quadratic powers,
and Newton recurrences takes polynomially many operations on integers
of polynomial bit length. Reduce each `2^j e_j` modulo `2^p` only
after these exact computations.

For the circuit, factor `Y^(h-J)` out of (2). Horner evaluation of
the remaining expression costs `O(J)` nodes, including the successive
powers of `X^2` for `Q`. Binary powering of `Y^(h-J)` takes `O(n)`
nodes. The quadratic base and final unpaired factor use a constant
number of nodes. This proves the node and complete generation bounds.

## 3. A large artificial common factor can survive modular compression

The circuits remove one cost obligation. They do not make a modular
gcd an exact rational gcd. In fact, large artificial common factors
can persist at every precision below an exponential input degree.

Choose positive `h,k`, put `A=2h+1`, `B=2k+1`, and set

```math
N=AB,\qquad I=[1,A),\qquad J=[1,B).
```

Every positive product in this rectangle is less than `N`, so its
rational gcd (1) from the companion note is one. However the bases
in (2) are

```math
Y_P=X(X-A),\qquad
Y_Q=N^2-BNX=-NB(X-A).
```

Equation (2) implies both modular inputs are divisible by

```math
(X-A)^{\max(0,\min(h,k)-p+1)}\pmod{2^p}.           \tag{3}
```

Thus merely asking for a common divisor of the compressed polynomials
modulo `2^p` can return a polynomial of enormous degree, although the
true rational gcd is one. This is an obstruction to that particular
modular test, not a factoring lower bound. It identifies the next
obligation: remove artificial multiplicity or certify a rational
common factor using the arithmetic boundary relations, with all costs
charged. Increasing precision to the numerical interval length would
lose the desired complexity bound.

## 4. Exact evidence and circuit generation

```sh
python3 tools/check_interval_packet_circuits.py --output /tmp/interval-packets.json
cmp /tmp/interval-packets.json evidence/interval-packet-circuits.json
python3 tools/interval_packet_circuits.py --lower 1 --upper 1025 --precision 16 --output /tmp/packet.json
make verify-arithmetic
```

The [generator](../tools/interval_packet_circuits.py) writes explicit
addition/multiplication node lists. The
[checker](../tools/check_interval_packet_circuits.py) passed 256 paired
symmetric-function controls and 1,024 **complete coefficient**
comparisons with dense products at precisions 1 through 16. It also
checks (3) directly by translating dense polynomials in 204 small
cases whose queried rectangles have no factors.

The [receipt](../evidence/interval-packet-circuits.json) records four
large circuit hashes and sizes. For degree `2^200+12` at precision
48, the `P` and `Q` circuits have 546 and 641 nodes respectively.
Four evaluations agree with the independent parity/unit-jet method;
two evaluations at the zero jet agree with a separate direct
Newton formula for the first coefficients of `Q`. The long cases
are exact compressed evaluations, not dense coefficient expansions.
