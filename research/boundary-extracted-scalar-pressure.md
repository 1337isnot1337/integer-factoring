# Extracting a protected scalar line through bounded boundary words

## Result

The [protected scalar projection](protected-subtensor-pressure.md)
does not have to preserve the actual auxiliary terminal line.
It may preserve any reachable auxiliary line that the initial
boundary can observe. A finite combination of boundary words
extracts that line. The bounded extra lengths change only the
constant in an exponential word-pressure bound.

This extends the exact screening engine from the actual terminal
vector to every constant binary auxiliary line. It rejects a line
that is unreachable or unobservable, so a hidden scalar component
cannot supply a spurious lower bound.

The conclusion is about **uniform exponential full-word rates**,
or their limsup growth rates. The boundary words can have different
lengths, so a same-length pointwise lower bound is not asserted.
The existing restrictions remain: the theorem bounds regular-channel
pressure and does not lower-bound directly type-restricted coefficient
ranks or the complexity of integer factoring.

## 1. Reachability and a protected line

Work over `F=F_2(tau)` with local ring `O=F_2[[tau]]`. Let the
finite physical alphabet of auxiliary matrices be `K_lambda`,
and let the actual boundaries be `f,e`. Write

```math
\Gamma_t(\lambda_1,\ldots,\lambda_t)
=f^T K_{\lambda_1}\cdots K_{\lambda_t}e.
```

The matrices are regular after the usual uniform scalar shift.
The reachable column space and observable row space are

```math
\mathcal R=\mathrm{span}_F\{K_w e:w\text{ is a finite word}\},
\qquad
\mathcal L=\mathrm{span}_F\{f^T K_w:w\text{ is a finite word}\}.
```

Both spaces are found by finite linear-algebra closure. In dimension
`D`, the successive spans stabilize after at most `D` steps.

Let `v` be a nonzero vector satisfying

```math
v\in\mathcal R,
\qquad
\ell v\ne0\text{ for some }\ell\in\mathcal L.    \tag{1}
```

Choose regular physical projections `P,Q` at an endpoint. Suppose
the projected local action preserves the line `Fv` and has scalar
matrix `S` on it:

```math
\left(\sum_{a,b}P_{i,a}K(a,b)Q_{b,j}\right)v
=S_{i,j}v.                                       \tag{2}
```

The vector can be normalized to be regular and primitive. The
implemented screen currently enumerates only constant binary `v`,
for which this normalization is automatic.

## 2. Bounded boundary extraction

Reachability supplies a finite expression

```math
v=\sum_{j=1}^m c_j K_{w_j}e,
\qquad c_j\in F.                                  \tag{3}
```

The suffix words `w_j` need not have the same length. By visibility,
there is a fixed prefix word `u` with

```math
f^T K_u v\ne0.                                    \tag{4}
```

The prefix can be a single word because a nonzero linear combination
of the scalars `f^T K_w v` has at least one nonzero summand.

Multiply (3) by a fixed power of `tau` so that every coefficient
is regular. Apply `t` projected interior sites and then the prefix
row `f^T K_u`. By (2), the left side becomes a fixed nonzero scalar
times `S^tensor(t)`. The right side is a finite sum of regular
projections of the original tensors

```math
\Gamma_{t+|u|+|w_j|}.                              \tag{5}
```

Fixing a prefix or suffix physical word is a row and column selection
for either physical cut. It does not require a choice depending on
the opposite interior word. The projections on the interior tensor
as `P^tensor(t)` and `Q^tensor(t)`, so physical separability is
preserved.

Let `Phi_s(M)` be the weighted sum of finite Smith orders at scale
`s>1`. Finite Toeplitz ranks show both

```math
\Phi_s(A+B)\le\Phi_s(A)+\Phi_s(B)
```

and monotonicity under regular left and right multiplication.
Multiplication by a nonzero fixed scalar changes `Phi_s` only by
its fixed valuation factor. Therefore (3)--(5) imply

```math
\Phi_s(S)^t
\le\sum_{j=1}^m C_j\,Q_{t+L_j}(s),              \tag{6}
```

for constants `C_j>0`, bounded lengths `L_j=|u|+|w_j|`, and any
regular-channel word counts `Q_n(s)` for the actual crossing.

If `Q_n(s)<=C beta^n` uniformly in `n`, then (6) gives

```math
\Phi_s(S)^t\le C'\beta^t
```

with `C'` independent of `t`. Consequently

```math
\boxed{\quad\beta\ge\Phi_s(S).\quad}             \tag{7}
```

The same conclusion holds for the limsup exponential growth rate.
There is no need to pad the suffix words to a common length. A
finite collection of bounded shifts is enough for (7).

This is the reason an arbitrary nonzero line of a minimal auxiliary
representation may be tested. Every such line is reachable and
visible. For a nonminimal representation, both conditions in (1)
must be checked separately.

## 3. Exact implementation and an additional control

The function `boundary_scalar_projections` in
[`tools/scalar_terminal_projection.py`](../tools/scalar_terminal_projection.py)
computes the exact reachable and observable spaces over `F`. It
then examines all constant nonzero binary vectors within a declared
enumeration budget. For each eligible vector, the saturated-kernel
algorithm computes the maximal regular left projection satisfying
(2), with `Q=I`. The largest of its scalar Smith pressures is a
valid lower bound on a uniform full-word rate.

This search is complete for constant binary lines and regular left
projections. It is not complete for arbitrary `F`-lines or joint
left and right physical projections. Exceeding the line budget is
reported as inconclusive.

The new control starts with the elementary scalar crossing

```math
T(z)=C+z(E+C)
```

and sets

```math
K(z)=T(z)I_2+(1+z)C\,N,
\qquad
N=\begin{pmatrix}1&1\\1&1\end{pmatrix},
\qquad N^2=0.                                     \tag{8}
```

The actual boundaries are `f=e=(1,0)^T`. The matrix crossing has
the strict endpoint `K(1)=E I_2`. In the row cut, the actual terminal
line admits no nonzero left scalar projection: its residual contains
the full-rank physical identity tensor `C`.

Nevertheless, `v=(1,1)^T` is a scalar invariant line for every
physical symbol, because `Nv=0`. It is reachable and visible. Choose
the physical suffix `(x,y,u,v)=(0,1,0,1)`, where `C=1,E=0`, and put
`F=1+z`. At that symbol,

```math
K_{\rm suffix}=F\begin{pmatrix}0&1\\1&0\end{pmatrix},
\qquad v=e+F^{-1}K_{\rm suffix}e.
```

Thus the exact different-length identity is

```math
F\,\Gamma_t+\Gamma_{t+1,\,\mathrm{suffix}}
=F\,T^{\otimes t}.                               \tag{9}
```

It recovers the full four-dimensional scalar projection even though
the actual-terminal projection has dimension zero. At infinity in
the row cut, the recovered orders are `(0,1,1,2)` with declared
degree one, yielding lower rate `s+2+s^-1`.

The checker verifies every one of the `256` physical entries of
(9) at interior length two. It also checks a separate hidden scalar
component and rejects its line because the actual boundaries cannot
reach or observe it.

Run:

```sh
python3 tools/check_scalar_terminal_projections.py
```

Evidence is saved in
[`evidence/scalar-terminal-projections.json`](../evidence/scalar-terminal-projections.json).
This strengthens a scoped screening method; it does not improve the
currently proved factoring exponent.

## 4. A nilpotent vacuum extension cannot hide a single-letter obstruction

There is a useful stronger application to the four-letter reflection
seed from [the row-fiber checkpoint](triangular-auxiliary-row-fiber.md).
Let `T` be that scalar seed, and retain its stated endpoint gauges.
With the same nilpotent matrix `N` as in (8), set

```math
K^\circ(z;x,y,u,v)
=T(z;x,y,u,v)\left(I_2+\frac{1+z}{z}[x=0]N\right).
                                                        \tag{10}
```

Use the actual boundaries `f=e=(1,0)^T`. Allow an arbitrary finite
Laurent matrix perturbation `Delta` supported on just one nonzero
physical letter `x=a`, where `a` is `1`, `2`, or `3`. It may have
arbitrary auxiliary coefficients and degree. If a strict crossing
endpoint is required, impose `Delta(1)=0`.

The line `w=(1,1)^T` is still scalar on every unchanged physical
letter, since `Nw=0`:

```math
K(z;x,y,u,v)w=T(z;x,y,u,v)w\qquad(x\ne a).        \tag{11}
```

Unlike a hidden scalar component, this line is always reachable and
visible through the actual boundaries. The unchanged vacuum symbol
has `T(0,0,0,0)=z^4`, and therefore

```math
K_{(0,0,0,0)}e=z^4e+(z^3+z^4)w,
\qquad f^T w=1.                                  \tag{12}
```

Thus a suffix of length one and the empty suffix extract `w`.
Selecting only the physical rows with `x!=a` is regular at both
cuts. Equations (7), (11), and (12) import exactly the protected
scalar profiles already certified in
[the single-letter pressure theorem](protected-subtensor-pressure.md).
After restoring the endpoint monomial charges, every uniform
full-word regular-channel certificate for this entire family obeys

```math
\rho\ge
\begin{cases}
8/11,&a=1\text{ or }3,\\
5/7,&a=2.
\end{cases}                                      \tag{13}
```

All three floors exceed `1/sqrt(2)`. The conclusion includes
noncommuting auxiliary perturbations, state-dependent basis changes,
arbitrary finite lattice precision, and physical blocking within
the full-word pressure interface. As before, direct restriction to
a type before rank computation is outside this theorem.

One motivating family applies a regular joint physical and auxiliary
row operation to (10): for one allowed height-increasing pair `a,b`,
add `(1+z)z^k Q K^circ_b` or `(1+z)z^k K^circ_b Q` to the row with
`x=a`. Here `Q` is any constant auxiliary matrix and
`h(b)-h(a)<=k<0`. Since `h=(0,3,1,4)`, the destination `a` is
nonzero. This construction lies in (13), even if `Q` fails to
preserve `w` and the resulting auxiliary symbols do not commute.
Several row operations with the same destination remain covered;
operations affecting different nonzero letters need a separate test.

The checker clears the common Laurent pole by using `L=zK`. It adds
a noncommuting perturbation on each of the three possible letters,
recovers the twelve-dimensional scalar projections at both endpoints,
and checks all `110,592` length-two instances of the exact identity

```math
\Gamma_{t+1,\,\mathrm{vacuum}}+z^5\Gamma_t
=(z^4+z^5)(zT)^{\otimes t}\qquad(x_i\ne a).        \tag{14}
```

The recovered zero orders are `5+sigma_a`; the infinity orders are
`sigma_a` after its extra regularizing shift. Restoring that shift
cancels the common `z` used in the implementation. This is an exact
closure of a whole auxiliary family, not a claim based on a finite
numerical search.
