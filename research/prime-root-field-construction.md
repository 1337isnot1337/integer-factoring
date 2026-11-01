# A prescribed binary subfield and a primitive prime-order root

## Result and role

Let `Q` be an odd prime and `w>=1`. One can deterministically construct
an explicit field containing both `F_(2^w)` and a primitive `Q`th root
of unity in

\[
 2^{O(w)}+\operatorname{poly}(Q,w)                         \tag{1}
\]

bit operations. If the binary subfield is already represented, only the
polynomial term is needed. The resulting field has binary degree

\[
 d=wf\le w(Q-1),\qquad f=\operatorname{ord}_Q(2^w).        \tag{2}
\]

Its arithmetic has polynomial bit cost in `Q,w`. The algorithm does not
search through all binary polynomials of degree `d`, and does not factor
the integer `Q-1`.

For the [adaptive prime constructor](width-sensitive-prime-constructor.md),
`w=Theta_epsilon(log n)` and `Q=n^O_epsilon(1)`. Thus this root-field
preparation and each field operation cost polynomial time in `n`.
For the earlier [narrow-prime schedule](prime-transform-boundary-constructor.md),
the same bound is `exp(O((log n)^2))`, below its written-table budget.
These are preparation bounds, not a new network contraction exponent.

The proof is a specialized deterministic Berlekamp construction. The
fixed-algebra idea is described in
[Berlekamp's own account](https://math.berkeley.edu/~berlek/poly.html).
The argument below gives the complete specialization and its cost;
it replaces a scan over all subfield constants by binary trace tests.

## 1. Explicit Frobenius fixed vectors

Put `b=2^w` and work in

\[
 R=F_b[X]/(X^Q-1).                                        \tag{3}
\]

This algebra is reduced because `Q` is odd and the derivative of
`X^Q-1` is `X^(Q-1)`. The `b`th-power Frobenius permutes its monomial
basis by `X^j -> X^(bj mod Q)` and fixes every coefficient in `F_b`.

Compute the orbits `O` of multiplication by `b` on `Z/QZ`. They provide
the explicit basis

\[
 B_O(X)=\sum_{j\in O}X^j                                 \tag{4}
\]

of the `F_b`-linear fixed algebra `{z in R:z^b=z}`. Their supports are
disjoint, and being fixed is exactly the condition that coefficients
are constant along these orbits. No nullspace computation is necessary.

Factorization over an algebraic closure, or the usual finite-field CRT,
identifies this fixed algebra with one copy of `F_b` for each irreducible
factor of `X^Q-1`. In particular, for any two different irreducible
factors, some basis vector `B_O` has different values on their components.

## 2. Binary traces split all components

Let `gamma_1,...,gamma_w` be the given binary basis of `F_b`. For every
orbit and every basis element compute in (3)

\[
 T_{O,j}(X)=\sum_{h=0}^{w-1}
                    (\gamma_j B_O(X))^{2^h}.             \tag{5}
\]

These satisfy `T_(O,j)^2=T_(O,j)`, so their values on each field
component are either zero or one. The calculations require only
coefficient squaring and the permutation `X^i -> X^(2i mod Q)`.

The tests (5) separate any two different irreducible components. To
prove this, choose an orbit vector whose component values differ by
`z!=0` in `F_b`. The absolute trace pairing
`(gamma,z) -> Tr_(F_b/F_2)(gamma z)` is nondegenerate. An elementary
proof suffices: the nonzero polynomial
`Y+Y^2+...+Y^(2^(w-1))` has degree below `b`, so it cannot vanish at
all `b` field elements; its values lie in `F_2`, so one value is one.
Multiplication by `z` is a bijection, and expanding the corresponding
`gamma` in the binary basis shows that some `gamma_j` has trace pairing
one with `z`. Its test in (5) distinguishes the two components.

Consequently gcds with the polynomials (5) split any reducible product
of irreducible factors. Only these zero/one tests are needed; there is
no enumeration of the `b` possible coefficient values.

## 3. Extract one root field

Compute `f=ord_Q(b)` by repeated multiplication modulo `Q`, stopping
at the first return to one. This uses at most `Q-1` steps. Every
nonidentity `Q`th root has Frobenius orbit length exactly `f`, so every
irreducible factor of

\[
 \Phi_Q(X)=1+X+\cdots+X^{Q-1}                             \tag{6}
\]

over `F_b` has degree `f`.

Start with `G=Phi_Q`. Process the tests (5) in a fixed order. Whenever
`gcd(G,T_(O,j))` is a proper nonconstant factor, replace `G` by that
gcd. Stop as soon as `deg G=f`. At most one pass through all tests is
needed. Indeed, after a test is processed its value is constant on
every irreducible component retained in `G`; later replacements retain
that property. If the final `G` had two components, the separating test
from Section 2 would contradict this conclusion.

The resulting `G` is monic and irreducible. Set

\[
 E=F_b[X]/(G),\qquad \zeta=X\bmod G.                      \tag{7}
\]

Then `zeta^Q=1` and `zeta!=1`. Since `Q` is prime, its order is exactly
`Q`. The prescribed subfield is the explicit constant subfield in (7).
The tower basis consisting of the binary basis of `F_b` times
`1,zeta,...,zeta^(f-1)` is an explicit binary representation of `E`.

This construction also covers `f=1`: it extracts a linear factor and
returns a primitive root already lying in the subfield.

## 4. Complete bit cost

An irreducible binary polynomial of degree `w` can be found by exhaustive
search, testing each candidate with the standard Frobenius criterion.
That costs `2^O(w)` bit operations and gives the ground-field basis.
This is the only exponential search in (1).

There are at most `Q` orbit vectors and at most `Qw` trace tests.
Each trace is computed by `w` squarings in a cyclic algebra of dimension
`Q`; those squarings are coefficient squarings and index permutations.
Polynomial gcds over `F_b`, with degrees at most `Q`, have polynomial
bit cost by the Euclidean algorithm. Even a coarse implementation
therefore uses `poly(Q,w)` bit operations for all remaining preparation.

Arithmetic in (7) uses polynomial multiplication and reduction modulo
`G`, with ground-field operations of polynomial bit cost in `w`.
Inversion uses the extended Euclidean algorithm. Alternatively,
exponentiation to `2^(wf)-2` gives a polynomial-operation inverse of
every nonzero element, since that exponent has only `wf` bits.
There is no uncharged search in the potentially much larger degree
`wf`.

The integer `Q` is assumed prime. Trial division can certify this in
`poly(Q)` time if it was not already certified by the multiplication
constructor. This remains within (1).

## 5. Exact replay and scope

The companion checker
[`check_prime_root_field_construction.py`](../tools/check_prime_root_field_construction.py)
constructs the orbit vectors and binary trace splitters, extracts a
factor, and independently checks its irreducibility by Frobenius gcd
tests. It also checks the exact primitive-root order and the degree-one
cases. Run:

```sh
python tools/check_prime_root_field_construction.py --output evidence/prime-root-field-construction.json
```

The committed [evidence](../evidence/prime-root-field-construction.json)
records 40 constructed root fields, 40 independent Frobenius
irreducibility checks and 43 additional irreducibility gcd checks,
165 exact trace idempotents, 41 proper gcd splits, and 1,080 products
in the prescribed subfield. Five models have their primitive root
already in that subfield. All comparisons passed.

When used inside a coefficient algebra, the caller must still specify
how additional selector fields or nilpotent factors are adjoined and
charge their dimensions. This result supplies a concrete root field;
it does not treat a compositum or tensor-product algebra as unit cost.
