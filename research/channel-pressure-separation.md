# Certify the global channel-policy rate by linear separation

At a fixed rational scale `s>1`, the
[subspace channel oracle](subspace-smith-channel-oracle.md) can certify
lower bounds for **every** stationary policy in the declared state
model. It is unnecessary to enumerate policies or rely on the
convergence of a floating policy search.

For a proposed lower rate `beta`, the certificate is a nonzero
nonnegative vector `h` satisfying

$$F_s(h)\ge\beta h,\tag{1}$$

where each coordinate of `F_s` is the exact minimum cost of an
admissible channel basis at that source state. Feasibility of (1),
after normalizing `sum h=1`, is a rational linear-programming problem
with a polynomial-time separation oracle. The largest feasible rate
is exactly the minimum spectral radius over all policies, including
reducible policies.

This is an application of standard nonlinear Perron and independent-row
optimization ideas. For primary background, see
[Gaubert--Gunawardena](https://arxiv.org/abs/math/0105091v2) and
[Cvetkovic--Protasov](https://doi.org/10.1007/s10107-020-01585-z).
The proof below supplies the reducible finite-family argument needed
here. No new general spectral-optimization theorem is claimed.

For the ten-state binary companion at `s=5`, exact certificates now give

$$
12.8203599666
\le \min\rho(T_\infty)\,\min\rho(T_0)
\le 12.8212.\tag{2}
$$

Both minima range over the complete declared subspace-policy families.
This is a comparison inside a sufficient rank-certificate model. It
does not lower-bound the actual coefficient ranks or factoring time.

## 1. Independent row choices

Let `A_i` be the finite set of admissible nonnegative transition rows
at source state `i`. A complete policy chooses one row independently
from every `A_i`, giving a matrix in

$$\mathcal T=A_1\times\cdots\times A_d.$$

The regular coefficient lifts themselves can have many descriptions,
but the transition rows depend only on finitely many target/order
labels and their multiplicities. Thus `A_i` is finite. Its explicit
listing is unnecessary.

Define, for all `h>=0`,

$$F_i(h)=\min_{a\in A_i}a\cdot h.\tag{3}$$

The map `F` is continuous, order-preserving, positively homogeneous,
and coordinatewise concave. The labeled-subspace greedy algorithm
computes `F(h)` and a minimizing row at each state, even when some
coordinates of `h` vanish. Zero costs only produce ties and do not
invalidate matroid optimization.

The upper-potential identity from the earlier construction is

$$
\lambda_*:=\min_{T\in\mathcal T}\rho(T)
=\inf_{h>0}\max_i\frac{F_i(h)}{h_i}.\tag{4}
$$

The new useful form is

$$
\boxed{\quad
\lambda_*=
\max\{\beta:\exists h\ge0,\ \textstyle\sum_i h_i=1,
                       \ F(h)\ge\beta h\}.
\quad}\tag{5}
$$

Allowing zero coordinates is necessary in (5).

## 2. Proof of the lower-potential characterization

First suppose `h>=0` is nonzero and `F(h)>=beta h`, with `beta>=0`.
Every policy satisfies

$$Th\ge F(h)\ge\beta h.$$

Iterating gives `T^k h>=beta^k h`. The spectral-radius formula then
gives `rho(T)>=beta` for every policy. Hence every feasible `beta`
in (5) is at most `lambda_*`.

For the reverse direction, let `J` be the all-ones matrix and put

$$F_\epsilon(h)=F(h)+\epsilon Jh,\qquad\epsilon>0.$$

On the nonnegative unit simplex, `F_epsilon` is strictly positive.
Normalize it by the sum of its coordinates and apply Brouwer's
fixed-point theorem. This gives `h_epsilon>0`, of coordinate sum one,
and a positive scalar `lambda_epsilon` satisfying

$$F_\epsilon(h_\epsilon)=\lambda_\epsilon h_\epsilon.$$

Choose a policy `T_epsilon` whose rows attain the minimum in (3) at
this vector. Then

$$
(T_\epsilon+\epsilon J)h_\epsilon
=\lambda_\epsilon h_\epsilon.
$$

Since the eigenvector is positive, this eigenvalue is the spectral
radius of `T_epsilon+epsilon J`. Every other policy satisfies

$$
(T+\epsilon J)h_\epsilon
\ge\lambda_\epsilon h_\epsilon,
$$

and therefore has spectral radius at least `lambda_epsilon`. Thus

$$
\lambda_\epsilon
=\min_{T\in\mathcal T}\rho(T+\epsilon J).\tag{6}
$$

The policy family is finite. By continuity of the spectral radius,
the right side tends to `lambda_*` as `epsilon` decreases to zero.
The simplex is compact, so take a convergent subsequence of
`h_epsilon`, with limit `h>=0` and coordinate sum one. Continuity
gives

$$F(h)=\lambda_* h.$$

This vector is feasible in (5) at `beta=lambda_*`, completing the
proof. The limiting vector can have zero coordinates, even though
every perturbed eigenvector is positive.

## 3. The separating linear inequality

For fixed rational `beta`, condition (1) is equivalent to

$$
a\cdot h-\beta h_i\ge0
\quad\text{for every }i\text{ and every }a\in A_i,\qquad
h\ge0,\quad\sum_i h_i=1.\tag{7}
$$

These are linear constraints. Given a proposed rational `h`, check
the simplex constraints and then compute all minimum rows by the
subspace oracle. If `F_i(h)<beta h_i`, the returned minimizing row
`a` supplies the explicit violated inequality

$$a\cdot h-\beta h_i\ge0.$$

If no coordinate violates (1), every constraint in (7) holds. This
is strong separation for the implicit rational polytope.

At fixed rational `s`, the transition coefficients are sums of at
most `r` powers `s^(epsilon L-a)`. Their bit lengths are polynomial
in the written scale, degree, valuation range, and channel count.
The subspace oracle has polynomial bit cost in the explicitly
declared state family and these data. Standard rational-polyhedron
separation/optimization machinery therefore applies; see
[Groetschel--Lovasz--Schrijver, *Geometric Algorithms and Combinatorial Optimization*, Chapter 6](https://doi.org/10.1007/978-3-642-78240-4).
It includes lower-dimensional feasible sets, so zero coordinates do
not create a missing full-dimensionality assumption.

A rational upper bound for the rate is `max_i F_i(1)`: choose all
rows minimizing at the all-ones vector. Bisection between zero and
this bound, using (7), approximates the global rate to a requested
additive accuracy with polynomially many separation computations
through the standard ellipsoid reduction. This statement concerns
the fixed-scale rate value. It does not assert an implemented
polynomial iteration bound for the practical driver below, or a
global optimization over `s`.

## 4. Scope of a nonnegative lower witness

A global spectral witness need not certify growth from the chosen
initial state. For example,

$$T=\begin{pmatrix}1&0\\0&2\end{pmatrix},\qquad
h=(0,1)^T$$

satisfies `Th=2h`, and `rho(T)=2`. Words starting in the first state
and ending there grow at rate one. The lower witness is zero at
that initial state.

Accordingly, every certificate here explicitly concerns the spectral
radius on the complete declared state set. It is not promoted to a
boundary-specific word lower bound. The state set, physical cut,
precision, local crossing, and scale also remain part of its scope.
A different model or an actual-rank argument can give a smaller bound.

## 5. Exact certificate gate and practical search

[channel_pressure_separation.py](../tools/channel_pressure_separation.py)
provides an exact minimum-row oracle accepting nonnegative rational
potentials. It uses only subspace rank increments, so no residue
vectors or channel bases are enumerated. A verified lower certificate
records `h`, `F(h)`, all margins, and a minimizing transfer matrix.

The optional search driver keeps a finite set of discovered row
constraints and calls a floating linear-programming solver. Each
solution proposes a rational vector. The full exact oracle checks
that vector against (1); if it fails, minimizing rows add cuts.
Only a successful exact check returns `certified`. A floating
infeasibility message, repeated cut, or iteration-budget failure
returns `no_exact_certificate` and is inconclusive.

This is a practical cutting-plane driver, not an implementation of
the ellipsoid complexity proof. Floating arithmetic proposes the
witness and is not used to justify its lower bound. The default
replay checker does not import the linear-programming package.

For the companion at `s=5`, the exact lower and upper rates are:

| Endpoint | Lower bound for every policy | Upper bound from a verified policy |
|---|---:|---:|
| Infinity | `53197/5000 = 10.6394` | `266/25 = 10.64` |
| Zero | `1204989/1000000 = 1.204989` | `241/200 = 1.205` |
| Product | `64101799833/5000000000 = 12.8203599666` | `32053/2500 = 12.8212` |

The two endpoint policies are independent, so their minimum rates
multiply. The lower interval endpoint is a genuine global bound over
the declared policy families at this scale. It is not inferred from
a floating search's best observed value.

Replay the exact stored witnesses with:

```sh
python3 tools/check_channel_pressure_separation.py
```

To repeat the optional witness search:

```sh
python3 tools/check_channel_pressure_separation.py --search
```

The [checker](../tools/check_channel_pressure_separation.py) compares
60 minimum-row evaluations, including potentials with zero coordinates,
against exhaustive independent-basis enumeration. It reconstructs
1,280 entries of the two upper-policy identities. The reducible
two-state example checks a valid zero-coordinate witness and rejects
both an invalid positive witness and an overlarge rate. The
[evidence record](../evidence/channel-pressure-separation.json) separates
the floating search history from the exact rational certificates.

Substantial OpenAI assistant support was used. No Lean formalization
or stronger complete factoring score is claimed.
