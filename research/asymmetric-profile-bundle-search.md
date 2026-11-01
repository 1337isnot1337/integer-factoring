# Fixed-type bundles for complete asymmetric crossing pressure

**Status:** a numerical optimizer comparison registered on 2026-11-01.
It introduces no new crossing, exact contraction or factoring exponent.
The complete score remains **1/2**. It reuses the fixed-type proximal
method in [the core bundle search](../tools/search_rectangular_core_bundle.py)
and [the coordinate adapter](pascal-coordinate-bundle-search.md).

The short mixed-graph SLSQP run attempted 51 free-height fits; only three
converged, 47 reached the iteration limit and one rejected a native query.
That motivates a separate comparison in raw height coordinates. The
registered two-hour coefficient walk is left unchanged and runs separately.

For fixed positive supporting coefficients `a,b`, use raw heights `h`
and intrinsic scale `t>0`. Let `r,s` be the complete weighted profiles and
`Gx,Gy` the two physical energy lists. Define

```
f(h,t) = a log sum exp(-r/a) + b log sum exp(-s/b).
```

Every determinantal prefix is the minimum of affine weighted-minor
functions, hence concave in `(h,t)`. For a nonincreasing probability
vector `p`, summation by parts expresses `sum p_i r_i` as a nonnegative
combination of these prefixes. The Gibbs variational formula therefore
expresses each log-profile pressure as a supremum of convex functions.
Thus `f` is convex in the raw variables at a fixed angle.

For physical types `piX,piY` of positive total entropy `H`, the fixed-type
rate is `(f+piX.Gx+piY.Gy)/H`. The common-type root is the **infimum** of
these rates. It need not be convex in the heights. At its Gibbs type,

```
grad f = H * grad(common-type rate) - physical_energy_vector.
```

Only this log-profile gradient supplies bundle planes. Root gradients
are never treated as global convex planes. Alternating fixed-type steps
and Gibbs updates is still a nonconvex numerical search, with no global
optimality claim. At a fixed raw point, changing only the supporting angle
uses the already complete profile pair and charges both physical energies.
All planes are discarded after an angle change.

[Independent controls](../tools/check_asymmetric_profile_bundle.py) check
99 log-profile derivatives, 168 directed held-out plane inequalities and
24 scaled/raw root identities on a nonbinary literal tensor at three
angles. Maximum derivative error was below `7e-10`. The
[control record](../evidence/asymmetric-profile-bundle-controls.json)
contains every point and plane. These are floating checks, not exact
convexity certificates or global pressure bounds.

## Registered comparison

The [manifest](../evidence/asymmetric-profile-bundle-manifest.json) freezes
three tensors and initial points: the reference, the completed selected
free-height best, and the best point of a recorded incomplete walk snapshot.
The ongoing walk is not represented as complete. Give each tensor 600
seconds, with 100 epochs and 100 steps per epoch as ceilings. The native
oracle must return every profile order; reject all unresolved queries.

Raw heights lie in `[-256,256]`, scale in `[exp(-9),4]` and normalized
heights in `[-64,64]`, enforced by linear inequalities. Log tangent lies
in `[-2,2]`. A master uses at most 200 SLSQP iterations. Gibbs types update
after six serious steps and at model stops. After an epoch, a bounded
one-dimensional angle fit includes both endpoints. Every third epoch
restarts near the best point using its recorded OS-random seed.

Every successful query retains its raw parameters, log-profile value,
subgradient, physical vector and full native details. Every new plane is
compared with every earlier point in both directions within the same
fixed-angle epoch; slack below `-2e-7` aborts the run. Rejected coordinates,
master failures and step statuses are recorded. Numerical stops and master
solutions do not certify an optimum. Any support below one requires a
fresh exact profile and rigorous pressure check before a recursion claim.

```
make verify-asymmetric-profile-bundle
python3 tools/search_asymmetric_profile_bundle.py --manifest evidence/asymmetric-profile-bundle-manifest.json --output /path/bundle.json
```

Credit: original rectangular/core and coordinate bundle methods, generic
asymmetric pressure, complete rectangular Smith backend and mixed-graph
constructor. Directed by Cinna Davis / 1337isnot1337 with substantial OpenAI
assistant support. No new Lean formalization or independent mathematical
acceptance is claimed.
