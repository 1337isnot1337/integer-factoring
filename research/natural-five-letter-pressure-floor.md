# A global balanced floor at the five-letter reference energy shape

Every scalar Laurent-polynomial strict crossing over any field, with
both normalized paired energy lists proportional by one common positive
scale to

```text
0,2,6,8,10,
```

has balanced two-place Smith parameter **rho >= 5/7 > 1/sqrt(2)**.
The two lists may be permuted independently. The statement allows
arbitrary polynomial degree, arbitrary splitting of these energies among
the eight physical height functions, and arbitrary regular local frames.
It therefore closes this entire balanced fixed-energy shape, beyond
conjugations or any bounded deformation family of the reference kernel.

This is an application of the existing
[cell-set caps](cell-minor-caps.md), [shared-prefix relaxation](shared-minor-profile-relaxation.md)
and [shared-minor all-scale certificate](shared-minor-smith-floor.md).
Those mechanisms and their contributor credit are retained. The new
contribution is the complete global certificate for this energy shape
and its consequences for the five-letter construction campaign.

The theorem concerns a sufficient crossing criterion. Other energy
shapes, asymmetric supporting lines, auxiliary modules and general
coefficient recombinations remain outside it. The factoring score is
still **1/2**.

## 1. Reduce to two nonnegative variables

Use the normalization of the cell-cap note. Let `d>=0` be the normalized
minimum-order deficit, and let `r,s` be the two nonnegative normalized
complete profiles, each with first entry zero. Write the two paired
energy lists as `c*E`, where `c>0` and `E=(0,2,6,8,10)`. At arbitrary
positive scale `t`, put

```math
u=td\ge0,\qquad v=tc>0,\qquad R=10/7.
```

Let `X_k` be the sum of the first `k` entries of `E`, and `Cell_k` the
sum of the `k` smallest of its 25 cell energies `E_i+E_j`.
For the increasing average profile `gamma=t(r+s)/2`, the existing
all-field minor inequalities give

```math
\sum_{i\le k}\gamma_i\le k u+v\,Cell_k\quad(1\le k\le25),
\qquad
\sum_{i\le k}\gamma_i\le(ku+vX_k)/2\quad(1\le k\le5),
\qquad\gamma_1=0.
```

Take the greatest convex minorant of all these prefix caps and write
its slopes as `gamma*(u,v)`. The cell-cap majorization argument gives

```math
G_R\ge L(u,v):=u+2\log\sum_i e^{-\gamma_i^*(u,v)}
             -2R\log\sum_j e^{-v E_j/R}.
```

Here `G_R` is the complete paired-type pressure at balanced parameter
`rho=R/2=5/7`. Positivity of `L` everywhere proves the claimed floor.
No independent optimization of the two actual profiles is presumed.

## 2. Three unbounded regions

**Large `u`.** Both partitions contain the zero channel and the physical
partition is at most five. Thus

```math
L(u,v)\ge u-2R\log5>0\qquad(u\ge5).
```

**Small `v`.** Each cap grows by at most `20*k*v` from its value at
zero, so the convex prefix ceiling grows by at most the same linear
function. Therefore

```math
L(u,v)\ge L(u,0)-40v.
```

At `v=0`, the existing five-letter shared-minor bound gives

```math
e^{u/2}\sum_i e^{-\gamma_i^*(u,0)}
\ge U_5(e^{u/8}),\qquad
U_5(x)=x^4+4/x+20/x^5\ge1029/100.
```

For completeness the checker reconstructs its exact Bernstein
certificate on `1<=x<=3` and its elementary tail on `x>=3`.
Consequently `L>0` for `0<=v<=1/10000`. The outward lower margin
after subtracting `40/10000` exceeds `0.0599`.

**Large `v`.** The second shared cap gives `gamma*_2<=u+v`. Hence

```math
u+2\log\sum_i e^{-\gamma_i^*}
\ge u+2\log(1+e^{-u-v})\ge2\log(1+e^{-v}).
```

The last inequality follows because its middle expression is
nondecreasing in `u>=0`. All four positive physical energies are at
least two. Using `log(1+x)>=x/2` for `0<=x<=1` and `log(1+x)<=x`,

```math
L(u,v)\ge e^{-v}-8R e^{-7v/5}>0\qquad(v\ge7).
```

The final sign follows from `e^(14/5)>8R=80/7`. The first five terms
of the positive exponential Taylor series already prove this strict
rational comparison.

## 3. Exact compact coverage

It remains to cover

```text
0 <= u <= 5,     1/10000 <= v <= 7.
```

For a rational rectangle `[ul,uh] x [vl,vh]`, monotonicity of every
prefix cap implies the rigorous lower bound

```math
ul+2\log\sum_i e^{-\gamma_i^*(uh,vh)}
  -2R\log\sum_j e^{-vl E_j/R}.
```

The [checker](../tools/check_natural_five_letter_floor.py) evaluates this
with 80-bit outward dyadic arithmetic, using the existing exponential
and logarithm remainder bounds. It certifies **2,456 rectangles**, of
maximum subdivision depth 19. Every lower endpoint is strictly positive;
the smallest recorded compact corner bound exceeds `0.00000522`.
This is a certified lower bound, not a numerical global minimum.

The [packet](../evidence/natural-five-letter-floor.json.xz) retains
every leaf path and lower bound, the tail certificates, exact source
hashes and controls. A separate replay reconstructs every rectangle and
rechecks every sign. It validates consistent splitting axes, prefix-free
binary paths and exact Kraft sum one, establishing full compact coverage
without relying on the generation queue.

```sh
make verify-natural-five-letter-floor
python3 tools/check_natural_five_letter_floor.py --output /tmp/natural-five-floor.json
```

## 4. Consequences for the retained kernels

Fifty independent fraction-free source orders reproduce both complete
reference profiles. Their first five combined prefix sums are exactly

```text
0,2,8,16,26,
```

which saturate all five shared-line bounds. In particular, every kernel
with this paired-energy shape satisfies

```math
\mu_1+\mu_2+\nu_1+\nu_2\le2.
```

The [graded zero construction](graded-zero-lattice.md) has normalized
zero prefix `8,10`. Pairing it with the old infinity prefix `-8,-7`
would make the left side three. Thus that particular optimistic pairing
is impossible for **any** strict crossing at these paired energies,
independently of which gauge family produced it. Its actual, poorer
opposite profile remains valid.

The packet also checks the known negative relaxation point with paired
energies `(0,0,2,8/3,3)`. It remains strictly negative at `R=10/7`.
This control demonstrates why changing the energy shape is a real
remaining direction and why the theorem must not be stated for every
five-letter crossing.

The full default verification retained at
[a2d894d](../evidence/verify-a2d894d.json) predates this certificate.
This certificate has its separate exact generation and replay status.
Substantial OpenAI assistant support was used for the application,
certificate implementation and text. No new Lean formalization or
stronger factoring algorithm is claimed.
