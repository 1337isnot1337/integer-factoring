# Constructing polynomial crossings by exact interpolation

## Result and scope

Prescribing local row divisibility at zero and infinity turns a crossing
construction problem into finite linear algebra. The method allows
non-Gaussian tensor entries and arbitrary Laurent polynomial degrees.
The degree range is forced by the prescribed physical gauges; it is
not a numerical search cutoff.

This note gives the construction method and a complete finite example.
For the existing four-label Schur Gaussian seed, its useful physical
gauges and two specified row lattices leave exactly a two-dimensional
space of endpoint-preserving even perturbations. Both directions are
local Gaussian changes at the `v` leg. Without the parity restriction,
the complete space has dimension six: every solution is one of 64
explicit local changes at that same leg. All 64 retain the original
two Smith profiles. A relaxed pair of row lattices also produces an
explicit non-Gaussian crossing, which checks that the interpolation
method can leave the Gaussian class.

These statements concern the specified row lattices. They do not exclude
other row bases, other gauges, coefficient recombinations, or a better
crossing. No improved factoring exponent is established here.

The starting seed is the
[exact Schur Gaussian crossing](schur-gaussian-digit-reduction.md).
The pressure conventions are those of the
[arbitrary-type transfer](type-optimized-crossing-pressure.md).
The interpolation method and exact checks were developed with substantial
OpenAI assistant support.

## 1. Prescribed local row divisibility

Let `T(z;x,y,u,v)` be a tensor over `F_2[z,z^-1]`, with each physical
label in `0,...,q-1`. Require the crossing endpoint

```math
T(1;x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u}.                   \tag{1}
```

Choose physical height functions `h_(0,ell)` and `h_(infinity,ell)` at
the four legs. For the physical tuple `lambda=(x,y,u,v)`, write their
sums as `h_0(lambda)` and `h_infinity(lambda)`. The gauged local entries
are

```math
z^{h_0(\lambda)}T(z;\lambda),\qquad
\xi^{h_\infty(\lambda)}T(\xi^{-1};\lambda),\quad \xi=z^{-1}.
                                                               \tag{2}
```

Flatten the first tensor in one crossing cut, and the second in the
other. At each end choose two matrices `U_e,V_e` that are invertible
over the local power-series ring. Choose nonnegative integer row targets
`tau_(e,i)`. Impose

```math
\operatorname{ord}_e (U_e M_e(T)V_e)_{i,j}\ge\tau_{e,i}
\qquad\text{for every }i,j.                                \tag{3}
```

Condition (3) is a useful rank certificate. It factors the transformed
matrix as `diag(z^tau_i) B(z)` with `B` regular at that end. Consequently
its sorted Smith orders dominate the sorted row targets. Equivalently,
every truncated block-Toeplitz rank is bounded by
`sum_i max(0,precision-tau_i)`. The usual coefficient transfer can use
the row-target profile directly.

The unit matrices preserve Smith orders. Their specified jets are enough
to check (3), since only finitely many coefficients precede the targets.
They can be taken to be polynomial matrices with invertible constant
terms; no numerical approximation to a formal inverse is needed.

## 2. The complete finite parameter space

Because every row target is nonnegative, (3) makes both gauged matrices
regular. Their unit inverses then make the original gauged entries
regular as well. Hence every nonzero monomial `z^k` in the entry
`T(z;lambda)` must satisfy

```math
-h_0(\lambda)\le k\le h_\infty(\lambda).                   \tag{4}
```

Thus the endpoint and row constraints have a finite parameter space,
even if arbitrary Laurent degrees were allowed at the outset.

Suppose a known endpoint tensor `T_*` lies in these entry windows. The
difference `T-T_*` vanishes at `z=1`. For each physical tuple, every such
difference has a unique expression as a linear combination of

```math
(1+z)z^k E_\lambda,
\qquad -h_0(\lambda)\le k<h_\infty(\lambda),               \tag{5}
```

where `E_lambda` is supported only on that tensor entry. To prove
completeness, divide the Laurent polynomial difference by `1+z`.
This division preserves its smallest exponent and lowers its largest
exponent by one. The result has exactly the monomial range in (5).

Insert the basis (5) into (3). Each forbidden coefficient is a linear
functional over `F_2`. Solving the resulting affine system therefore
either constructs a crossing or proves that none exists for those
specific row lattices and heights. Extra linear symmetries, such as
even fermion parity, can be imposed by removing the corresponding
entry variables.

### Why the two local expansions stay finite

At zero, a basis vector in (5) becomes
`z^(k+h_0(lambda))(1+z)`. At infinity it becomes
`xi^(h_infinity(lambda)-k-1)(1+xi)`. Both leading exponents are
nonnegative by (5). No truncation loses an unknown negative-power term.
For row `i`, only coefficients below `tau_(e,i)` enter the system.

## 3. A complete example for the four-label seed

Let `G(z)` be the four-label Gaussian tensor in the Schur digit note.
Use states `0=empty`, `1=plus`, `2=minus`, `3=both`, and these physical
height vectors:

| Leg | At zero | At infinity |
| --- | --- | --- |
| `x` | `0,0,0,0` | `0,1,1,2` |
| `y` | `0,0,0,0` | `0,1,1,2` |
| `u` | `0,1,1,2` | `0,0,0,0` |
| `v` | `0,-1,-1,-2` | `0,2,2,4` |

Use the row cut at zero and the column cut at infinity. Both exact Smith
profiles are

```math
\sigma=(0,1,1,2,3,3,4,4,4,4,5,5,6,7,7,8).                 \tag{6}
```

The checker constructs deterministic Smith row and column changes through
precision `20`. Their constant terms are invertible. They transform
the two seed matrices to the diagonal in (6) modulo `z^20` or `xi^20`.
Use these particular unit changes in (3), with row targets `sigma`.

There are `512` even-parity endpoint-preserving variables in (5).
The exact binary constraint matrix has rank `510`. Its kernel therefore
has dimension two.

### The entire fiber is explicit

Introduce the eight physical Grassmann variables, in the order
`x+,x-,y+,y-,u+,u-,v+,v-`. Over `F_2` their square-zero multiplication
has no sign ambiguity. Interpret `G` as its coefficient-generating
polynomial in these variables.

Every even Laurent polynomial crossing satisfying the prescribed row
lattices is exactly one of the four tensors

```math
\boxed{\quad
G(z)\left(1+(a z^2+b z^3)(1+z)\theta_{v,+}\theta_{v,-}\right),
\qquad a,b\in\mathbb F_2.
\quad}                                                     \tag{7}
```

The checker verifies that the two displayed perturbations lie in the
kernel and are linearly independent. Combined with rank `510`, this
proves completeness of (7); it is not an inference from testing the
four tensors alone.

These are local operations at the `v` leg. If `N_v` sends the empty
state to the doubly occupied state and kills the other states, then
`N_v^2=0`, and the corresponding local operator is `I+f(z)N_v`.
It is a unit at both gauged ends: at zero its coefficient is
`z^-2 f(z)`, and at infinity it is `xi^4 f(xi^-1)`.
Both are regular for the two polynomials in (7). Thus all four tensors
have exactly the original Smith profiles.

The factor in (7) adds a quadratic fermion edge within the `v` boundary,
so all four tensors remain Gaussian. Consequently an even
non-Gaussian perturbation cannot retain these particular row lattices.
This is a statement about the chosen lattices, which are stronger data
than the two unordered Smith profiles.

### Removing the parity restriction completely

The same finite-window argument includes every physical tensor entry.
There are now `1,024` Laurent endpoint-preserving variables. The
constraint matrix has rank `1,018`, so the entire unrestricted fiber
has dimension six.

Let `E_(a,b)` be the four-by-four matrix unit that sends physical state
`b` to state `a` at the `v` leg. All solutions have the following form:

```math
T(z)=L_v(z)G(z),                                         \tag{8}
```

where the local map is

```math
\begin{aligned}
L_v(z)=I
 &+(1+z)z\left(
 aE_{1,0}+bE_{2,0}+cE_{3,1}+dE_{3,2}\right)\\
 &+(1+z)(ez^2+fz^3)E_{3,0},\\
 &\hspace{22mm}a,b,c,d,e,f\in\mathbb F_2.
\end{aligned}                                            \tag{9}
```

Here `L_v` acts only on the `v` index of the tensor. The checker verifies
that the six displayed perturbations are linearly independent and lie
in the full kernel. The exact rank computation then proves that (9)
lists every solution, including all possible Laurent degrees.

Put `N_v=L_v-I`. Every nonzero entry of `N_v` increases occupation, so
`N_v^3=0` and

```math
L_v^{-1}=I+N_v+N_v^2.                                    \tag{10}
```

At zero, a one-occupation step in (9) acquires the factor `z^-1`;
at infinity it acquires `xi^2`. Thus `(1+z)z` is regular after both
gauges. A two-occupation step acquires `z^-2` and `xi^4`, respectively,
so both remaining coefficients are regular as well. The inverse in
(10) is regular at both ends. Also `L_v(1)=I`, which preserves the
crossing endpoint.

The `v` leg lies on the column side of both specified flattenings.
These regular local maps and their regular inverses therefore preserve
both Smith profiles. The checker verifies the local map, its inverse,
the resulting tensor, and both exact profiles for every one of the
64 parameter choices.

The even members are precisely `a=b=c=d=0`, recovering the four
tensors in (7). Thus the fixed-lattice classification extends beyond
even or Gaussian perturbations: every admissible tensor differs from
the seed by one of the local leg maps (9). This remains a statement
about these particular row lattices and physical gauges.

## 4. A non-Gaussian constructive control

The method itself does not require Wick identities. For example,
retain only target `3` at row number `3` of the zero-cut Smith basis,
with rows numbered from zero. Set every other row target to zero at
both ends. The affine system is solvable.

The checker constructs an explicit endpoint-preserving tensor in this
larger class and verifies a failure of a quartic Wick identity. It also
checks the requested row divisibility directly after applying the two
unit changes. Its coefficient changes are recorded in the evidence as
pairs `(mask,k)`, meaning addition of `(1+z)z^k` to that entry.

This control distinguishes failure of the stronger fixed-lattice search
from an implementation that only generates Gaussian tensors. It is not
an improved pressure certificate.

## 5. Reproduction

Run the standard-library checker:

```sh
python3 tools/check_polynomial_crossing_interpolation.py \
  --output evidence/polynomial-crossing-interpolation.json
```

It verifies the two unit basis changes, both the `512`-variable and
`1,024`-variable constraint ranks, and all spanning perturbations.
It checks all `64` unrestricted fiber members, including their local
map inverses and both Smith profiles. The original even-fiber checks
and an additional odd-parity member give `90` independent
block-Toeplitz rank checks. It checks all Gaussian Wick identities in
the even fiber and a non-Gaussian control outside it. The relaxed
interpolation example is constructed by the same exact linear solver.

See the [checker](../tools/check_polynomial_crossing_interpolation.py)
and [verification output](../evidence/polynomial-crossing-interpolation.json).
