# Does the binary floor survive unequal supporting weights?

## Research question and changed hypothesis

The [universal binary Smith floor](universal-binary-smith-floor.md) excludes
balanced two-place pressure for every binary strict polynomial tensor and
all physical heights. Its displayed proof averages the two endpoint
profiles with equal weights. The [new uniform axis window](uniform-polynomial-axis-window.md)
excludes the extremes, but does not settle interior unequal supporting
weights. Can the asymmetric recursion yield a binary witness despite the
balanced obstruction?

The next diagnostic minimizes an explicit relaxation of *all* binary
endpoint-prefix inequalities and the shared second minor, at arbitrary
interior weights. A negative value is a necessary spectral target only;
no tensor is inferred. A clear negative (less than -1e-6) requires fresh
high-precision readback and an actual polynomial realization. If only
near-zero cold or axis limits are found, preserve the finite result and
stop the numerical diagnostic; do not infer a global binary theorem from
optimizer success or extend its last decimals.

## Complete relaxed profiles

Normalize both four-channel profiles to have first order zero. Let D>=0
be the physical-ground slack, and let 2X,2Y>=0 be the two physical
excitation energies, including the thermal scale. Set

```
S=X+Y,  m=min(X,Y),  C=4(D+S),  K=2(D+m).
```

Every endpoint profile p satisfies

```
p_1=0, p increasing,
p_2<=2D+2m, p_2+p_3<=3D+2S, sum p<=C,
r_2+s_2<=K.
```

These are the full first through fourth endpoint prefixes, plus the
shared two-minor. At fixed p_2=u, the smallest relaxed partition sets
sum p=C and balances the last two orders subject to the third prefix:

```
P(u)=(0,u,v,C-u-v),
v=min((C-u)/2,3D+2S-u).
```

The feasible upper bound is
`U=min(C/3,(3D+2S)/2,2D+2m)`. The partition decreases as u increases
on [0,U], so the coupled optimum uses r_2+s_2=K. This is possible because
K<=2U. Optimize u in `[max(0,K-U),min(U,K)]`, with the other profile
P(K-u). At norm-one weights a,b>0 the resulting lower relaxation is

```
D+a log(sum exp(-P(u)/a))+b log(sum exp(-P(K-u)/b))
 -log(1+exp(-2X))-log(1+exp(-2Y)).
```

No coefficient-rank equality, algebraic realization or uniform transfer
is supplied by these profiles alone.

## Bounded diagnostic registration

Use twelve fixed outer starts: `(D,X,Y)=(0,0,0),(1,1,1),(1,2,1),(4,1,1)`,
each at log tangent -2,0,2. Outer bounds are D,X,Y in [0,16] and log
tangent in [-4,4]. Each fit has at most 500 SLSQP iterations and 15 seconds;
the complete diagnostic has 180 seconds. For every outer query, minimize
the one-dimensional convex profile allocation on its full interval.
Retain all start/end coordinates, best queried profiles and statuses.
The only search objects are necessary real profiles, not tensor entries.
The complete factoring target is a realized negative critical margin
with full costs; score remains `1/2`.


## Completed bounded diagnostic

The [full attempt ledger](../evidence/asymmetric-binary-prefix-relaxation.json)
retains all twelve starts, best and last queried profiles, counts, source
hash and optimizer outcomes. All twelve reported local convergence after
**7,853 complete relaxation queries** in 1.78 seconds. No time cap was
reached. No margin below zero was found; the smallest was approximately
`9.5676e-13`, at D=0 with excitation energies about 30.17 and 30.69 and
weights about 0.99538 and 0.095996.

These near-zero cold values meet the registered stopping condition. In
fact profiles with D=0, X=Y tending to infinity and an interior allocation
make every term tend to zero, so such a numerical boundary is expected.
No high-precision negative-witness replay is triggered, and no actual
tensor search follows from this diagnostic. The outer convergence flags
do not prove a global nonnegative bound or an asymmetric binary theorem.
The one-dimensional convex allocation was also computed numerically.
Intermediate successful outer queries were counted, not retained; each
trial's best and last complete query remains explicit.

A continuation requires a new analytic inequality or a realizable
interior tensor mechanism. Increasing these outer cutoffs is not justified
by this finite result. The complete factoring score remains `1/2`.

This reuses the credited binary normalization, shared minors and the
asymmetric recursion criterion. Directed by Cinna Davis / 1337isnot1337
with substantial OpenAI assistant support; no new Lean or novelty claim.
