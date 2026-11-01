# Realizing the first admitted high-axis Gaussian profile

## Research question and stop conditions

The [shared-minor profile budgets](multiaxis-gaussian-profile-floor.md)
exclude the common-scale Gaussian criterion through eight axes, but admit
a relaxed nine-axis profile below `1/3`. Can an explicit homogeneous
Cayley kernel realize a sufficiently favorable *complete* profile there?
The mechanism is a common crossing on `2^9` legs, with nine rank cuts and
score `log_9(1/rho)`. It changes the number of regrouping axes, rather
than refining a blocked two-cut seed. The target is `rho<1/3`, plus the
complete effective multiaxis transfer.

The concrete kernel generalizes the credited three-axis construction.
For every nonempty axis subset S set
`c_S=product_(a in S) x_a^(D/|S|)`, with `D=lcm(1,...,d)`, and put
`c_empty=sum_(S nonempty)c_S`. In the characteristic-two augmentation
algebra use `f=sum c_S epsilon_S`, `epsilon_a^2=0`, and take its Cayley
matrix on cube vertices. Its diagonal is zero and its all-one endpoint
is the antipodal matching. All entries have total degree D. Its
principal-Pfaffian tensor is therefore an actual common Gaussian
crossing, with the full occupied-mode degree cost retained.

At an axis, multiplication by
`sum_(S contains a) c_S epsilon_(S minus a)` is its crossblock in a
constant augmentation basis. Its square is `x_a^(2D) I`; complete Smith
orders consequently occur in pairs summing to 2D. This additional
constraint still permits a numerical nine-axis relaxed rate about
0.329045, below `1/3`. That number is a diagnostic target, not a bound
for the displayed kernel.

The first diagnostic checks the literal crossblock profiles at a fixed
finite list of dimensions and coefficient-field specializations. Such
specialization can increase determinantal valuations, so favorable
specialized profiles are **not** a generic-profile proof. If even the
optimistic specialized profiles miss the charged target, stop this
kernel. If they pass, the next obligation is a symbolic factorization
or a complete determinantal-divisibility certificate over the other
variables, followed by a rational outward pressure certificate. No
realized nine-axis profile or stronger score is claimed in advance.


## Completed diagnostic and an exact nine-axis floor

The [specialization ledger](../evidence/high-axis-cayley-specializations.json)
retains every complete order at dimensions three through nine. No promised
pivot was unresolved. The nine-axis matrix has order 256 and degree 2520;
its specialized numerical minimum is about `0.421401383`, above `1/3`.
That numerical minimum alone is not a global certificate. The following
minor proof closes this particular kernel at every scale and frequency.

In the eight-variable augmentation algebra, multiplication by
`q=sum_(i<j) epsilon_i epsilon_j` has rank **120** over every
characteristic-two field. The associated alternating form has matrix
`I+J` of order eight, which is invertible because `J^2=0` and
`(I+J)^2=I`. A symplectic change of generators writes q as a sum of
four disjoint products. Multiplication by each two-generator product is
a four-dimensional square-zero complex with two-dimensional homology.
Tensoring the four explicit complexes gives total homology dimension
16. Thus its rank is `(256-16)/2=120`.

Choose a nonzero 120-by-120 minor of that multiplication matrix. The
[checker](../tools/check_high_axis_cayley.py) independently verifies a
literal binary minor and records its rows and columns. Nonzero determinant
terms all increase augmentation degree by two per entry. Consequently the
sum of selected row degrees minus selected column degrees is 240.

Now use the *same minor* in the axis crossblock, specializing all other
variables to one. An entry increasing augmentation degree by r has
monomial `x^(2520/(r+1))`. Every nonzero term of this determinant has
`sum r=240`. Strict convexity gives

```
sum 2520/(r+1) >= 120*2520/3 = 100800,
```

with equality exactly when every r is two. Its coefficient is precisely
the nonzero binary minor of multiplication by q. The determinant
therefore has **exact valuation 100800**. This proves, in the original
rational function field of the other variables,

```
Delta_120(a) <= 100800
```

at every axis. Specialization supplies an upper bound here, in the
correct direction: it cannot create a lower generic valuation. The
nonzero binary coefficient persists over every characteristic-two field.
The proof needs no independent replay of the 256 native orders.

For `0<t<1`, Jensen and monotonicity now give

```
log Z_a(t) >= 120*log(1+t^840)
           >= (256/3)*log(1+t^840).
```

Summing nine axes gives at least `768*log(1+t^840)`. Binary entropy
satisfies `H(p)+p*log u <= log(1+u)`. Set `u=t^840` and recall that
`mD=256*2520`. These identities imply

```
sum_a log Z_a(t) - mD*p*log t >= 768*H(p).
```

Dividing by `9*256*H(p)` proves **`R_9(p,t)>=1/3`** for every
`0<p<1` and every positive scale. It excludes a strict half-score
improvement by this kernel's unweighted common-scale degree-cap
calculation. It does not exclude other nine-axis kernels, different
physical codes, or tighter ranks after coefficient extraction.

The [exact receipt](../evidence/high-axis-cayley-controls.json) retains
the binary minor, its independent dense rank, all nine convexity checks,
512 endpoint coefficients, and 56 independent untruncated small orders.
The all-scale proof, rather than the finite numerical fits, is the
obstruction. Stop this balanced monomial kernel; increasing dimension
or refining its fitted rate is not the next research step.

```
make verify-high-axis-cayley
```

This generalizes the credited three-axis source in
[cube-homogeneous-gaussian.md](cube-homogeneous-gaussian.md), reusing its
augmentation algebra, Pfaffian rank construction and the multiaxis
transfer. Directed by Cinna Davis / 1337isnot1337 with substantial OpenAI
assistant support. The complete factoring score remains `1/2`.
