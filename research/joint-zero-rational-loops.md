# Zero-loop reachability and a coupled two-place pencil

An all-field obstruction prevents a single square-zero simple-pole gauge
from taking the first two zero Smith orders from `(8,9)` to `(8,10)` or
better. The leading rank-one coefficient is allowed to change. The proof
uses exact affine and Frobenius consequences, with no finite-field
exhaustion or solver assumption.

A separate identity supplies an eleven-parameter strict polynomial family
coupling the zero and one places. Its first finite controls remain
uncompetitive. Neither result improves the factoring score `1/2`.

## 1. Necessary equations for a simple pole

Use `A(t)=sum_j A_j t^j`, `R(t)=D_H A(t)=sum_j R_j t^j` from
[the refined Newton analysis](joint-graded-newton-ceiling.md). The reference
has `R_j=0` for `j<8`. Let `P` be any constant energy-block matrix with
`P^2=0`. The gauge `U=I+P/t` has inverse itself, and

```math
A_U=A+t^{-1}[P,A]+t^{-2}PAP.
```

Polynomiality requires `P A_0 P=0` and `[P,A_0]+P A_1 P=0`. Because
`P` commutes with `D_H`, the same formula applies to `R_U`. Keeping its
smallest order at eight requires

```math
P R_8 P=0,\qquad [P,R_8]+P R_9 P=0.
```

Together with `P^2=0`, these are quadratic equations in 69 energy
parameters. Pure linear row-space consequences and pure-square
Frobenius roots reduce the necessary parameter space

```text
69 -> 41 -> 36.
```

All monomials are formal: linear terms and squares occupy different
columns. The argument never uses `x^2=x` or finite-field reduction.
For binary coefficients, a pure-square equation is the square of the
corresponding linear equation over any characteristic-two field, and
Frobenius is injective.

## 2. A rank-one tangent system reaches a contradiction

The first two coefficients of the normalized zero matrix are

```math
C_0=R_8+[P,R_9]+P R_{10}P,\qquad
C_1=R_9+[P,R_{10}]+P R_{11}P.
```

The `(0,0)` entry of `R_U` equals that of `R`: index zero is a singleton
energy block and `P_00^2=0` forces `P_00=0`. Thus `(C_0)_00=1`.

If the second Smith order were at least ten, `C_0` would have rank one
and `C_1` would be tangent to the rank-one locus. Equivalently there
would be vectors `u,v,a,b` satisfying

```math
C_0=u v^T,\qquad C_1=u b^T+a v^T,\qquad
u_0=v_0=1,\quad a_0=0.
```

This equivalence follows by eliminating the unit `(0,0)` pivot: every
remaining Schur-complement entry must vanish through degree one. It
permits every rank-one leading coefficient with that pivot, rather
than fixing the reference coefficient.

The 36 generator parameters and 97 free vector coordinates give 133
parameters. Exact affine and Frobenius elimination gives

```text
133 parameters, 1,317 equations -> 60 parameters
 60 parameters,   472 equations -> contradiction.
```

The first stage has 73 linear consequences. The second has 37 linear
and seven pure-square consequences; their affine system contains the
constant equation `1=0`. Substitution preserves every quadratic
monomial, including repeated variables. The retained evidence records
the consequence rows and the complete terminal generator basis.

Therefore, whenever this gauge is polynomial and its first zero order
stays eight, its second zero order is at most nine. Higher poles,
products of moves, other starting tensors and trades that worsen the
first order remain outside this obstruction.

## 3. A family coupling the zero and one places

Let `N` be the square-zero endpoint generator with all-ones blocks on
joint index pairs `(4,20)`, `(7,11)`, `(8,16)`. The checker verifies
`N^2=0`, every `N A_j N=0`, and the exact factorization

```math
[N,A(t)]=t^8(t+1)^4(K+t^4L),
```

where `K=[N,A_8]` and `L=[N,A_16]`. Consequently every polynomial `q`
with `degree(q)<=11` and `q(1)=0` gives

```math
U_q=I+\frac{q(t)}{t^8(t+1)^4}N,\qquad
A_q=U_q A U_q^{-1}=A+q(t)(K+t^4L).
```

The right side is polynomial and has the same strict identity endpoint
at one. The rational gauge tends to identity at infinity, so the original
worst infinity pole and its rank-one leading coefficient are retained.
Its complete infinity Smith profile need not stay fixed. The condition
`q(1)=0` leaves eleven independent coefficients and allows simultaneous
poles at zero and one, beyond the previously classified regular-at-zero
family. This construction works over every characteristic-two field.

Exploratory profiles for the 198 cases `q=c t^a(t+1)^b`, with
`a>=0`, `b>=1`, `a+b<=11`, `c` nonzero in `F4`, found no improvement at
the inherited physical heights. This finite catalogue does not optimize
the full pencil or other heights. Its profiles and floating pressure
estimates are retained in
[`joint-two-place-pencil-search.json`](../evidence/joint-two-place-pencil-search.json).
The optional search requires the existing NumPy/SciPy and native-oracle
dependencies:

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_joint_two_place_pencil.py
```

Adding `--random-trials N` also tests `N` arbitrary F4 polynomials in the
eleven-parameter pencil, with a reproducible seed.

## 4. Replay

```sh
make verify-joint-zero-prefix
```

The checker reconstructs the all-field contradiction and the pencil
identities. An independent dense F4 implementation checks 13,125 entries
of the quadratic equation assembly. It has no SMT dependency. The finite
profile catalogue above is exploratory and is not rerun by this gate.

Sources: [`check_joint_zero_prefix.py`](../tools/check_joint_zero_prefix.py)
and [`joint-zero-prefix.json`](../evidence/joint-zero-prefix.json).
This is ordinary proof with exact algebraic replay, not independent
mathematical review or new Lean formalization. Substantial OpenAI
assistant support.
