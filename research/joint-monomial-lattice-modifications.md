# Discrete joint lattice modifications of the five-letter seed

A joint coordinate gauge `U(t)=diag(t^p_i)` may be singular at zero.
It can therefore change both weighted Smith profiles and evade the
[formal polynomial unit theorem](joint-polynomial-formal-rigidity.md).
For the five-letter reference, every integer gauge retaining polynomial
entries and the original worst infinity pole reduces to one of exactly
12 effective tensors. Both profiles of each tensor are weakly
prefix-dominated by the reference. None improves its cost at any scale
at the inherited physical heights.

The finite classification covers every integer charge in this class;
it does not bound charges by an experimental window. Section 4 also
classifies every energy-block projector at every positive lattice depth:
under the stated conditions only four tensor-preserving projectors remain.
Multiple interacting subspaces and changes allowing a worse pole stay open.
The strongest complete factoring score remains `1/2`.

## 1. Entrywise valuations give exact difference constraints

Use `A(t)=T(t)S`, `A(1)=I`, integer joint energies `H_i`, and
`C(z)=z^8 A(z^-1)D_H(z)` from
[joint energy conjugations](joint-energy-conjugations.md). For any
integer vector `p`, define

```math
A_p(t)_{rc}=t^{p_r-p_c}A(t)_{rc}.                       \tag{1}
```

The endpoint remains exact, since `U(1)=I`. For each nonzero polynomial
entry let `l_rc` and `d_rc` be its least and greatest degrees.
Polynomiality of (1) is exactly

```math
p_c-p_r\le l_{rc}.                                    \tag{2}
```

The weighted infinity entry becomes
`t`-gauge conjugation of `C`, with exponent shift `-p_r+p_c` at `z`.
Thus retaining regularity of `C_p`, equivalently no pole below `-8`
in `B_p`, is exactly

```math
p_r-p_c\le H_c-d_{rc}+8.                              \tag{3}
```

All bounds in (2)--(3) are nonnegative for the reference. Represent
`p_c-p_r<=w` by a directed edge `r -> c` of weight `w`. Opposite
zero-weight paths force equality of the corresponding charges.
Floyd--Warshall closure gives exactly these five equality groups:

| Group | Joint state indices `x+5y` |
| --- | --- |
| `G0` | `0` |
| `G1` | `1,2,3,5,10,15` |
| `G2` | `4,20` |
| `G3` | `6,7,8,11,12,13,16,17,18,24` |
| `G4` | `9,14,19,21,22,23` |

Let `q_i` denote the charge on `Gi`. After collapsing the groups,
the complete set of effective inequalities is

```math
0\le q_3-q_0\le1,\qquad |q_2-q_3|\le1,
\qquad 0\le q_4-q_1\le1.                              \tag{4}
```

The undirected support components are `{G0,G2,G3}` and `{G1,G4}`.
Adding one common charge on either component commutes with every entry
of `A` and leaves (1) unchanged. Normalize `q_0=q_1=0`. The remaining
solutions are precisely

```math
q_3=a\in\{0,1\},\quad q_2=b\in\{a-1,a,a+1\},
\quad q_4=c\in\{0,1\}.                               \tag{5}
```

There are `2*3*2=12` solutions and 12 distinct resulting tensors.
Equations (2)--(5) prove completeness for every integer charge.

## 2. Exact profiles compare the entire family

The [checker](../tools/check_joint_monomial_lattices.py) recomputes all
entry bounds, the shortest-path equality groups, every quotient edge,
the 12 tensors and all crossing endpoints. It computes both complete
weighted Smith profiles for every tensor with fraction-free polynomial
arithmetic, plus the two reference profiles: 26 profile computations.

For both places, write the orders increasingly. For every candidate
`p`, the exact prefix inequalities are

```math
\sum_{i=1}^k\nu_i(p)\le\sum_{i=1}^k\nu_i(0)
\quad(1\le k\le25).                                  \tag{6}
```

The hinge-sum proof in the
[constant-gauge comparison](joint-energy-conjugations.md#4-exact-profile-dominance-at-every-scale)
applies after a common shift making the orders nonnegative. Consequently
each candidate partition sum `sum s^(-nu_i)` is at least the reference
sum for every `s>=1`. The physical energies and restoration terms are
unchanged. Applying this to both places proves the claimed all-scale
cost comparison. Only the unchanged tensor has leading realignment
rank one; the other ranks range from two to four.

Every tensor and profile computation here is binary. Extending the
coefficient field in characteristic two preserves all ranks, minors,
valuations and the resulting profiles. This extends the same 12
tensors to those fields; it does not classify noncoordinate gauges
with new field parameters.

```sh
make verify-joint-monomial-lattices
python3 tools/check_joint_monomial_lattices.py --generate
```

The [evidence](../evidence/joint-monomial-lattices.json) retains the
quotient graph, every profile, leading rank and source hashes. It uses
only the Python standard library and has no floating-point step.

## 3. Noncoordinate modifications have an exact interface

The next discrete family replaces coordinate projections by an arbitrary
energy-block idempotent `P`, with `P^2=P` and `Q=I-P`. Set

```math
U(t)=tP+Q,\quad U(t)^{-1}=t^{-1}P+Q.
```

The complete transformation is

```math
A_U=PAP+QAQ+tPAQ+t^{-1}QAP.                            \tag{7}
```

It is polynomial if and only if `Q A(0) P=0`: the image of `P` must
be invariant under the cold operator. Since `P` preserves energy,

```math
C_U=PCP+QCQ+z^{-1}PCQ+zQCP.                           \tag{8}
```

The no-new-pole condition is exactly `P C_0 Q=0`: the kernel of `P`
must be invariant under the leading infinity operator. When it holds,
the new leading coefficient is

```math
(C_U)_0=PC_0P+QC_0Q+PC_1Q.                            \tag{9}
```

Thus complementary graded subspaces with these two invariant-subspace
conditions give an explicit noncoordinate search space. A general such
projection can change the zero-place profile too; preserving it must
not be assumed. The next section classifies this space when the leading
rank is also retained.

## 4. Every projector and every positive depth

**Theorem.** Over any characteristic-two field, let `P` be any
energy-block idempotent, `Q=I-P`, and `k>=1` an integer. If

```math
U(t)=t^kP+Q
```

gives a polynomial seed and introduces no infinity pole below `-8`,
then for `k>=2` the transformed seed is unchanged. For `k=1` the same
conclusion holds if the leading realignment has rank one.
In either case there are exactly four possible projectors: zero,
identity, and the two coordinate projections onto the neutral support
components in section 1. Their ranks are `0,25,13,12` respectively.
The dimensions 13 and 12 refer to joint-state subspaces.

### 4.1 One-step equations over an arbitrary field

At `k=1`, the necessary and sufficient equations are

```math
P^2=P,\quad (I-P)A_0P=0,\quad PC_0(I-P)=0,             \tag{10}
```

together with leading rank one for (9). Here `A_j=[t^j]A(t)`.
Using the last equation of (10) in characteristic two rewrites (9) as

```math
C'_0=C_0+[P,C_0]+PC_1+PC_1P.                          \tag{11}
```

The zero-energy state is a singleton block. On that state `P` is a
scalar idempotent, hence zero or one. Therefore the `(0,0)` entry of
the realignment `D=R(C'_0)` is exactly one. Leading rank one is
equivalent to the complete set of pivot-minor equations

```math
D_{rc}+D_{r0}D_{0c}=0\qquad(1\le r,c\le24).            \tag{12}
```

Expand all 69 energy-block entries of `P` in a binary basis with
arbitrary field parameters `alpha_i`. Equations (10) have degree two;
(12) has degree at most four. Treat repeated variables as formal
monomials: `alpha_i^2`, `alpha_i alpha_j` and higher powers are distinct.

The checker first reduces each entry in (11) by the binary row space
of (10). This preserves its value at every solution. Forming (12)
from those reduced entries is therefore an equivalent condition in
the presence of (10). It then performs the same all-field procedure
as the [quadratic/Frobenius proof](joint-energy-polynomial-dressings.md#62-exact-quadratic-elimination-for-degrees-1-through-17):
extract linear row-space consequences and pure-square consequences,
take the valid Frobenius roots of the latter, restrict to their common
linear kernel, and reconstruct the original equations in the new basis.
Quartic monomials are never replaced by binary-field values.

The exact chain is

```text
69 -> 13 -> 2.
```

The first stage has 1,191 distinct formal monomials and row rank 486;
it supplies 52 independent pure-linear and 20 pure-square consequences.
Their combined rank is 56. The next stage has 54 monomials and row
rank 35; its consequences have combined rank 11. Every solution over
every field in scope survives both restrictions.

### 4.2 Two jets handle every larger depth without a rank hypothesis

For general `k`, polynomiality is exactly

```math
QA_jP=0\quad(0\le j<k),                               \tag{13}
```

and no new infinity pole is exactly

```math
PC_jQ=0\quad(0\le j<k).                               \tag{14}
```

For `k>=2` it suffices to use the **necessary** subset with `j=0,1`
in both (13) and (14), together with `P^2=P`. This is only a quadratic
system, independent of `k`, and requires no leading-rank equation.
Exact linear/pure-square row-space elimination gives

```text
69 -> 33 -> 10 -> 3 -> 2.
```

This second terminal space is the same two-dimensional space as the
one-step quartic calculation. Thus one fixed two-jet proof covers
every `k>=2`; no cutoff or extrapolation from finitely many depths is
used.

### 4.3 The terminal space has exactly four idempotents

In both calculations the terminal basis is the pair of literal
coordinate projections `E,F` onto
`G0 union G2 union G3` and `G1 union G4`. Exact multiplication verifies

```math
E^2=E,\quad F^2=F,\quad EF=FE=0,\quad E+F=I,
\quad [E,C_j]=[F,C_j]=0\quad\text{for every }j.         \tag{15}
```

They preserve energy, so they commute with `A(t)` as well. Write any
surviving matrix as `P=alpha E+beta F`. The remaining idempotence
equations are `alpha^2=alpha` and `beta^2=beta`. Over a field these
force `alpha,beta in {0,1}`. Thus the only solutions are `0,E,F,I`,
and all four satisfy every original condition and leave the tensor
unchanged at every depth. This proves both necessity and sufficiency.

The [projector checker](../tools/check_joint_projector_lattices.py)
reconstructs both systems, both entire elimination chains, every
terminal identity (15), and 20,000 independent dense `F4` entry
controls for the assembled quadratic formulas. Its
[certificate](../evidence/joint-projector-lattices.json) retains all
linear and square-root consequences, terminal bases and source hashes.
The initial finite-field SAT experiments are not dependencies.

```sh
make verify-joint-projector-lattices
python3 tools/check_joint_projector_lattices.py --generate
```

### 4.4 One step at infinity needs no cold-place condition

In fact the one-step result has a stronger version. Drop the cold
condition `(I-P)A_0P=0` entirely. Idempotence, `PC_0(I-P)=0` and
leading rank one in (11) already force the same four neutral
projections. Repeating the formal-monomial elimination without the
cold equations gives the exact all-field chain

```text
69 -> 14 -> 2.
```

The checker replays this third system too. Thus the infinity condition
alone closes a larger rational family: for any scalar rational
function `f(t)=t+O(1)` at infinity, set `U=f(t)P+(I-P)`. Its conjugation
formula is

```math
C_U=PCP+QCQ+f(z^{-1})PCQ+f(z^{-1})^{-1}QCP.            \tag{16}
```

No new pole forces `PC_0Q=0`, and the leading coefficient is exactly
(9), independent of every lower term of `f`. Hence leading rank one
forces `P` to be one of the four neutral projections, and (16) is the
original tensor. This covers `f(t)=t-a` at **every** finite point `a`
in every field in scope, even before polynomiality is required.
Negative unit depths follow by interchanging `P,Q` and removing a
common scalar factor, which cancels from conjugation.

It does not cover scalar functions of degree zero at infinity. For
example `f(t)=(t-a)/(t-b)`, `a!=b`, gives a gauge tending to `I` at
infinity. Polynomiality requires the two exact finite-place conditions

```math
QA(a)P=0,\qquad PA(b)Q=0.                             \tag{17}
```

If `a,b` are both nonzero and different from one, the gauge is regular
and invertible at zero and one. It therefore preserves the entire
zero-place Smith profile and the strict crossing endpoint; its
leading infinity coefficient is also unchanged. These conditions give
a distinct rational-gauge construction target. They do not establish
the existence of a nontrivial projector or an improved profile.
The subsequent [single-evaluation algebra theorem](joint-evaluation-algebra.md)
closes this entire scalar-rational projector family over every
characteristic-two field, including higher multiplicities and
irreducible numerator/denominator factors, under regularity at zero,
one and infinity. Interacting projector directions remain outside it.

This theorem leaves interacting sequences of different lattice
subspaces, rational changes of degree zero at infinity, one-step changes
allowing a higher leading rank, and changes allowing a worse pole
outside its scope. This work used substantial OpenAI assistant support
and retains the Gaussian reference attribution. No new Lean work,
independent review, improved coefficient-rank bound or factoring
exponent is claimed.
