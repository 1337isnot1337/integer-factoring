# A nonlinear level-five engine using four binary quadratic forms

The [Tate-coordinate construction](tate-coordinate-coefficients.md) leaves one
binary coordinate after its cheap trace at phase order five. This note
expresses the full original-index coefficient by four positive binary
quadratic-form parity counts and square tests. It is an exact nonlinear
representation; its present evaluator takes O(sqrt(N) poly(log(N+1))) time.
This remains exponential in the input bit length, and one fixed phase order
does not recover arbitrary factors. The complete factoring exponent stays 1/2.

The order-five group identities, Artin--Schreier equations, divisor pairing
and binary quadratic counts are standard tools. This is a specialization of
the preceding coefficient construction with exact source checks, developed
with OpenAI Codex assistance. No new Lean theorem or independent mathematical
acceptance is claimed.

## 1. Normalize the two conjugate torsion coordinates

The intrinsic coefficient algebra is

```math
A_5=\mathbb F_2[c]/(c^2+c+1),\qquad
c\longmapsto Z^2+Z^3\text{ in }R_5.
```

Write X(t)=X_5(Delta(t)), Y(t)=sigma(X(t)), with sigma(c)=c^2, and let
B=X+Y. Then B is the known scalar trace and B(0)=1. Frobenius on the
constant phase takes a Tate point P to 2P; since 5P=0, its x-coordinates
cycle between X and Y.

On E_v:y^2+xy=x^3+v, the doubling formula is

```math
x(2P)=x(P)^2+v/x(P)^2.
```

For completeness, its tangent slope is y/x+x. Squaring and adding that
slope, and then using the curve equation, gives the displayed formula.
Consequently

```math
X^2Y=X^4+v,\qquad Y^2X=Y^4+v.
```

Adding gives XY B=B^4. Since B is a unit, XY=B^3, and the two conjugates
satisfy

```math
X^2+BX+B^3=0.                                          \tag{1}
```

Set U=X/B and V=U+c. Then U(0)=c, V(0)=0, and

```math
V^2+V=B+1.                                             \tag{2}
```

These identities hold for the full formal series, not only for sampled
coefficients.

## 2. The Artin--Schreier solution is a unary-square series

Put

```math
V_*(t)=\sum_{r\ge1}\bigl(t^{r^2}+t^{5r^2}\bigr).
```

Every positive integer is uniquely a power of two times an odd integer.
Thus, with H(t)=Delta(t)+Delta(t^5),

```math
V_*(t)=\sum_{j\ge0}H(t^{4^j}),\qquad
V_*(t)^2+V_*(t)=\sum_{j\ge0}H(t^{2^j})=B(t)+1.           \tag{3}
```

The last equality is the positive-degree trace expansion in the preceding
note. Equation (2) has a unique zero-constant formal solution: its even
coefficient is determined by the coefficient at half its index and its odd
coefficient directly by the right side. Therefore V=V_*. In particular,

```math
B=1+V+V^2,\qquad X=cB+BV,
\qquad BV=V+V^2+V^3.                                   \tag{4}
```

The constant of X is c. At every positive index, its trace coordinate is
the coefficient of V+V^2 and its remaining scalar coordinate is the
coefficient of V+V^2+V^3.

## 3. Four quadratic parity counts give the original coefficient

For positive a,b,N let R_(a,b)(N) be the parity of the positive integer
solutions of a x^2+b y^2=N. Frobenius gives V^3=V(t)V(t^2), so

```math
[t^N]V^3=
R_{1,2}(N)+R_{1,10}(N)+R_{5,2}(N)+R_{5,10}(N)\pmod2.     \tag{5}
```

Each coefficient of V is just the XOR of the tests that N is a positive
square and that N/5 is a positive square integer. Coefficients of V^2 are
zero at odd indices and are the coefficient of V at N/2 at even indices.
Equations (4)--(5) therefore give both coordinates of X at every index.
At odd N this is the desired signed divisor value; at even N it is the
signed odd-divisor value in the Tate series. Factoring returns factor two
before calling an odd-index engine.

Enumerate y up to floor(sqrt((N-a)/b)), and test whether (N-b y^2)/a is a
positive square integer. This constructs R_(a,b)(N) in O(sqrt(N)) integer
steps. The four fixed forms give complete time O(sqrt(N) poly(log(N+1)))
and O(log(N+1)) bits of working storage. There is no divisor-list oracle
in this evaluator. Its input-length cost is still exponential and does
not improve the accepted asymptotic factoring algorithm.

## 4. The same remaining bit isolates a divisor-count residue

There is an alternative description that identifies the arithmetic hidden
by the parity counts. Assume N is odd and put M=N/5^(v_5(N)). In the
embedding above, 1 is the sum of the two opposite residue pairs, while c
is the pair {2,3}. Write the coefficient as A_N+c B_N. Then

```math
A_N=\#\{d\mid M:d=1\text{ or }4\pmod5\}\pmod2,
\qquad B_N=[M\text{ is a square}].                       \tag{6}
```

Pair d with M/d. If M is a square, the only fixed divisor is sqrt(M), so
A_N is one exactly when sqrt(M) is 1 or 4 modulo five. If M is a nonsquare
quadratic residue modulo five, the involution preserves both residue
classes and has no fixed divisor, giving A_N=0. If M is a quadratic
nonresidue, the involution exchanges the two classes, giving

```math
A_N=\tau(M)/2\pmod2=[\tau(M)=2\pmod4].                  \tag{7}
```

Thus a full value is obtainable with at most one explicitly charged query
for tau(M) modulo four. The query occurs only for nonsquare quadratic
nonresidues. This description is not a fast unrestricted tau oracle.

If M=product p^e, its divisor count is two modulo four exactly when it has
one odd exponent, that exponent is one modulo four, and all other exponents
are even. This follows by taking the 2-adic valuation of product(e+1).
It explains why prime, square and general composite cases cannot be merged.

For example M=33 and M=363 are both nonsquare composites, both three modulo
five. Their divisor counts are respectively zero and two modulo four, and
their remaining coordinate is respectively zero and one. Merely combining
primality, the square test and M modulo five does not implement (7).

## 5. A squarefree promise changes the problem

If M is promised squarefree and M>=2, then tau(M)=2^omega(M). Equations
(6)--(7) simplify the entire phase-five value to

```math
[M\text{ prime}]\,[M\bmod5\in\{2,3\}]\cdot1.
```

For M=1 the value is 1+c. Under this promise, primality testing suffices;
deterministic polynomial-time primality testing is supplied by
[Agrawal--Kayal--Saxena](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf).
The checker uses literal small primality tests, not a claimed implementation
of a fast unrestricted squarefree decomposition.

In particular, a product of two distinct primes other than five gives zero
at phase five. This is a limitation of that one phase, not evidence about
the full collection of residue levels. The squarefree promise is not
established by testing that the input is a nonsquare or is not a perfect
power, and is not assumed for the unrestricted factoring endpoint.

## 6. Reproduction and next target

Run `make verify-q5-quadratic`, or

```sh
python3 tools/check_q5_trace_zero.py
```

The [checker and explicit evaluators](../tools/check_q5_trace_zero.py) compare
the nonlinear formulas with literal divisor residues on 8,192 odd indices.
The [receipt](../evidence/q5-quadratic-coefficients.json) separately labels
the quadratic enumeration cost, 4,095 supplied small tau queries, 6,919
squarefree-promise controls, 64 Tate/Delta prefix controls, 946 distinct-prime
semiprimes, twelve large cheap controls and eight large supplied-factor
controls. The large supplied-factor controls do not implement the
unrestricted divisor-count query. Their known prime exponents supply it.

The nonlinear construction succeeds at this fixed phase but has the wrong
time bound and cannot replace the other residue levels. The next case is
phase seven, whose trace-zero part has two binary coordinates. Extending
the level-five polynomial-in-unary-squares identity requires a new argument;
the scalar trace identity alone does not supply those coordinates.
