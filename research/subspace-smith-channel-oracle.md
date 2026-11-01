# Optimize regular channel bases with linear residue spaces

The complete [DVR channel oracle](dvr-adaptive-smith-channels.md)
enumerates primitive coefficient jets and then solves a vector-matroid
basis problem. The jet enumeration can be replaced by binary linear
algebra. This makes the local oracle polynomial in the number of
channels, the precision, and the explicitly declared lattice family.

The key is to allow a channel's normalized rows to lie in a declared
target lattice, instead of requiring that they generate that lattice
exactly. A charged valuation can likewise be smaller than the channel's
actual common valuation. Both choices preserve the exact physical
factorization and its all-length rank bound. With these safe choices,
the feasible residues for each target and charged valuation form a
linear subspace.

The implementation handles rank-sixteen channels at precision two
using systems with at most 32 binary variables. It reconstructs every
selected rational factorization entry exactly. The preceding companion
certificate `32053/2500` is retained; a stronger factoring exponent is
not claimed.

## 1. The admissible residue spaces

Let `O=F_2[[tau]]`. At a source lattice `M`, flatten the complete
physical action as in the
[state-adaptive construction](state-adaptive-smith-channels.md). Write

$$A_M=UR,\qquad
U_i=\tau^{\nu_i}\ell_i,\quad 1\le i\le r.\tag{1}$$

Here `U_i` is a column, every `ell_i` is regular, and the normalized
columns `ell_i mod tau` are linearly independent. Minimum-valuation
Schur elimination supplies such a factorization: its ordered pivot
rows give a triangular unit minor in the normalized columns. The
right factors `R_i` are regular scalar functions of the opposite
physical half.

Fix a finite list of full-tail row lattices

$$\tau^P O^D\subseteq N\subseteq O^D,\qquad P\ge1,$$

and a nonnegative charged order `a`. Define

$$
S_{a,N}=\left\{v\bmod\tau:
 v\in O^r,\quad
 \text{every auxiliary row of }\tau^{-a}Uv\text{ lies in }N
 \right\}\subseteq\mathbb F_2^r.\tag{2}
$$

The unnormalized condition is `Uv in tau^a N`, after reshaping the
column into its auxiliary rows. It is closed under addition and scalar
multiplication in `O`. Its image modulo `tau` is therefore a binary
linear subspace. A nonzero vector in this space has a primitive
regular lift `v`.

We need only `0<=a<=nu_max`. For a nonzero `v`, independence of the
normalized columns gives the exact valuation identity

$$
\mathrm{ord}_\tau(Uv)
 =\min_i(\nu_i+\mathrm{ord}_\tau v_i).\tag{3}
$$

For primitive `v`, the right side is at most `nu_max`. Thus larger
charged orders cannot have nonzero residues. Orders below `nu_min`
can be retained: they are valid undercharges, and adding them costs
only additional linear systems.

The exact-target DVR family is contained in (2). If its normalized
channel generates `N`, it certainly lies in `N`, at its actual order.
The present family also includes larger targets and lower charged
orders. There is no requirement that these extensions improve an
optimized rate on every crossing.

## 2. A finite linear system with at most `r(P+1)` variables

Put `delta_i=nu_i-a` and

$$w_i=\tau^{\delta_i}v_i.$$

By (3), regularity of `tau^-a Uv` forces every `w_i` to be regular.
Membership in `N` depends only on

$$\sum_i\ell_i w_i\pmod{\tau^P},\tag{4}$$

because the target contains the full tail `tau^P O^D`. The coefficient
restrictions and residue `b_i=v_i mod tau` are as follows.

| Value of `delta_i` | Allowed `w_i mod tau^P` | Residue `b_i` |
|---|---|---|
| `delta_i<0` | Any `P`-jet | `0` |
| `0<=delta_i<P` | Divisible by `tau^delta_i` | Coefficient of `tau^delta_i` in `w_i` |
| `delta_i>=P` | Zero | An independent invisible bit |

For a visible coefficient of `tau^j` in `w_i`, its lift into `v_i`
is the monomial `tau^(j-delta_i)`. An invisible residue bit is lifted
by the constant polynomial one. All these exponents are nonnegative.

Represent `N/tau^P O^D` by a binary basis in the `DP` coefficient
coordinates. For each auxiliary row in (4), reduction modulo this
basis gives homogeneous linear equations. Their joint kernel is
computed by Gaussian elimination while retaining preimages in the
coefficient variables. Project this kernel onto the residue coordinates
`b_i`, and perform another elimination with retained preimages. The
result is a basis of `S_(a,N)` and a polynomial lift for every basis
vector.

There are at most `rP` visible variables and `r` invisible bits.
No primitive jets or nonzero residue vectors are enumerated. If the
flattened left column has `J` auxiliary rows, at most `JDP` binary
constraint coordinates are needed before elimination.

### Exactness of the truncation

Any regular feasible vector `v` gives allowed jets and invisible bits
in the table, so its residue lies in the computed image. Conversely,
each kernel vector gives the indicated polynomial lift. Its normalized
left column has every row in `N` modulo `tau^P`, and the remaining
tail is also in `N`. Thus this lift satisfies the exact power-series
membership condition, not just an approximate one.

The polynomial representatives have degree less than

$$P+\nu_{\max}-\nu_{\min}.$$

This is the same finite precision bound underlying the exhaustive DVR
oracle. Linear elimination now avoids enumerating all coefficients
inside that bound.

## 3. Greedy optimization on labeled subspaces

Fix a positive potential `h_N` and `s>1`. Every nonzero residue in
`S_(a,N)` has the available cost

$$c_{a,N}=s^{\epsilon L-a}h_N,\tag{5}$$

where `epsilon=0` at zero and `epsilon=1` at infinity, and `L` is
the degree of the local crossing.

Conceptually, make a labeled vector-matroid ground set containing
every nonzero residue from every space, with cost (5). A feasible
channel basis is a selection of `r` linearly independent residues.
It is unnecessary to materialize this exponential ground set.

Sort the labeled subspaces by cost. Keep a basis `B` of the residues
already selected. For the next space `S`, scan its supplied basis
vectors and retain a vector whenever it increases the span of `B`.
This adds exactly

$$\dim(B+S)-\dim B$$

independent vectors, each with the space's cost and its retained
polynomial lift. At every cost threshold, the resulting span equals
the span of all ground-set elements of at most that cost. The usual
weighted matroid greedy argument therefore proves that the chosen
basis has minimum total cost among all admissible labeled residues.

This argument also explains why an arbitrary basis of each subspace
is sufficient. If the subspace extends the current span, at least one
of its basis vectors extends it; repeated extension realizes the
whole rank increment. Linear combinations inside the same subspace
cannot add another independent direction afterward.

The oracle is polynomial in `r`, `P`, the valuation range, and the
written source and target data. If the declared family has `s_0`
states, constructing all its actions requires at most
`s_0^2 (nu_max+1)` linear systems, with the state-dependent maximum
order used in the actual sum. This is polynomial in the *explicit*
state family. Enumerating every full-tail lattice for arbitrary
dimension and precision is a separate cost and can be exponential.

## 4. Recover an exact physical factorization

Let the selected polynomial lifts be the rows of `S`. Their residues
form a basis, so `S` is invertible over `O`, and `S^-1` is regular.
The algebraic identity is

$$UR=(US^T)((S^{-1})^T R).\tag{6}$$

For selected column `j`, divide `(US^T)_j` by its charged power
`tau^a_j`. Its reshaped auxiliary rows lie in the chosen target
`N_j`, so they have regular coordinates `H_j` in that target's basis.
Writing `R'_j` for the corresponding row of `(S^-1)^T R`, we obtain

$$B_M K=\sum_j\tau^{a_j}H_j B_{N_j}R'_j.\tag{7}$$

The left factors depend only on the first physical half. The right
factors are regular scalar functions of the second half. A source
state's choice is fixed by its potential and never depends on the
opposite physical input. Consequently (7) satisfies exactly the
state-adaptive all-length coefficient-rank theorem.

An undercharged left factor may have additional common valuation.
Ignoring this extra vanishing gives an upper bound in the rank
argument. Enlarging a target likewise preserves regular membership.
Neither extension approximates or changes the local polynomial
crossing itself.

The nonlinear potential map is therefore computed exactly by the
subspace greedy oracle. Its value is the minimum over regular
unimodular channel bases with the declared target/order choices.
The finite-policy and positive-potential reasoning from the
[matroid construction](matroid-adaptive-smith-channels.md) applies
without change.

## 5. Implementation and independent controls

The engine is
[subspace_smith_channels.py](../tools/subspace_smith_channels.py).
Its kernel elimination retains original coefficient preimages;
its residue elimination retains polynomial lifts. `choose` returns
a transfer matrix and a fixed policy. `verify` reconstructs the
regular inverse change of basis, every target membership, and every
entry of (6) as an exact rational-function identity.

Run:

```sh
python3 tools/check_subspace_smith_channels.py
```

The [checker](../tools/check_subspace_smith_channels.py) first uses
the independent exhaustive DVR oracle on the binary companion at
precisions one and two, at both endpoint cuts. It closes those exact
actions under target-lattice inclusion and order undercharging, using
exact Hermite membership. The resulting labeled residue sets agree
with the spans returned by the new linear systems. Thus the test
compares the complete action families, not only a selected policy.

For each of the resulting 28 source states, the checker enumerates
every independent residue basis and verifies the new oracle's exact
minimum for a supplied positive potential. Every returned polynomial
subspace-basis lift is independently checked for membership.

The preceding ten-state companion potentials still certify

$$\beta_\infty\le266/25,\qquad
\beta_0\le241/200,\qquad
\beta_\infty\beta_0\le32053/2500.$$

This preserves the earlier rate while eliminating the jet-enumeration
cost. It does not assert a better optimum for this particular crossing.

Finally, the same companion is evaluated on an alphabet of size four.
All ten states have sixteen physical channels. At precision two,
240 linear systems use at most 32 variables each. The checker
reconstructs 10,240 local identity entries and 5,120 selected target
memberships without enumerating a primitive jet or a nonzero residue
vector in this rank-sixteen control. This establishes an executable
larger-alphabet search interface. The
[evidence file](../evidence/subspace-smith-channels.json) records the
systems, action comparisons, policies, and verification counts.

The proposed subspace formulation and the independent implementation
were developed collaboratively during the current research pass.
Substantial OpenAI assistant support was used. No Lean formalization
or stronger complete factoring score is claimed.
