# Prescribed smooth descent can fail after every subexponential trial cutoff

## Result and changed hypothesis

The [prescribed prime-power obstruction](prime-power-descent-obstruction.md)
uses an elementary Bonferroni truncation under
`log Z log log Z=o(log D)`. That leaves trial division at the complete
positive-score budget outside its stated scope. A classical combinatorial
sieve removes this restriction. The geometric counterfamily and prescribed
exponent alignment remain Sol worker 1's construction; the change here is
the uniform sieve remainder estimate.

Fix `C>0`, `alpha>0` and

```math
F(D)=\sqrt D\exp\left(\frac{C\log D}{(\log\log D)^\alpha}\right).
```

Let `Z(D)>=3` be nondecreasing, and let `3<=Y(D)<=Z(D)`, with

```math
\log Z(D)=o(\log D).
```

Then there are arbitrarily large odd composite, non-perfect-power integers
N with no prime divisor at most `Z(N)`, and `N=3 mod 4`, for which the
canonical prime-three target at the **least** k satisfying
`3^k>ceil(F(N))` represents no `Y(N)`-smooth positive norm strictly below
`3^k`, even allowing imprimitive vectors.

This covers `Z(D)=exp(K log D/(log log D)^beta)` for every fixed `K,beta>0`:
trial division through that cutoff already spends a complete positive-score
runtime scale, up to polynomial bit costs. Raising the trial budget within
`D^o(1)` does not repair this prescribed one-shot descent statement.
The result does not exclude arbitrary preprocessing of the same cost,
adaptive/collective targets, full windows, or useful partial presentations.
It is not an infinite semiprime or fundamental-discriminant construction,
nor a lower bound on factoring. The complete score remains `1/2`.

## 1. The classical input, with its uniformity exposed

Use the fundamental lemma in
[Tao, Notes 4, Lemma 17 and Corollary 19](https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/).
For a fixed-dimensional admissible sieve with multiplicative density g,
there are lower weights in `{-1,0,1}`, supported on squarefree divisors
`d<=z^s`, whose main term is `V(z)(1+O(exp(-s)))`.
The constant is uniform when the Mertens product condition and dimension
are uniform. The resulting remainder is bounded by the sum of the absolute
divisibility errors through that level. The
[combinatorial-sieve construction and final inequality in Kedlaya, Theorem 12.2](https://kskedlaya.org/ant/chap-brun.html)
give the same input. These are classical sieve results, not new claims here.

Our densities are `g(2)=1/2`, `g(3)=1/3`, and `g(p)=2/p` for other
active primes. They satisfy `g(p)<=2/p` and `g(p)<=1/2`. Mertens' product
bound therefore gives sieve dimension two with an absolute product constant.
No coefficient of the linear forms, discriminant, prime q, or interval
length enters that constant. Choose one sufficiently large **fixed** s_0
so that the lower main term is at least `V/2`. No numerical s_0 is asserted
by the finite controls, and no sieve threshold is supplied as input advice.

## 2. An error supported on small products

Fix `m=3^k`, an integer `Z>=3`, a prime `q>Z`, and the progression

```math
a=a_0+q^2t,\qquad 4ma_0-1\equiv q\pmod{q^2}.
```

Let t range over X consecutive integers. At every prime p<=Z exclude
`p|a(4ma-1)`. There is one forbidden residue at p=2,3 and two otherwise.
For p>=5 the two residues are distinct, since the second linear form is
-1 at a=0. The progression step is a unit at every active prime. Thus,
for every squarefree product d of active primes, CRT gives

```math
A_d=X\frac{\nu(d)}d+r_d,\qquad
|r_d|\le\nu(d)\le2^{\omega(d)}\le\tau(d).
```

This exact bound is uniform even though m and q grow. Put
`W=(Z+1)^s_0` and `V=prod_(p<=Z)(1-nu(p)/p)`.
The fundamental lemma, with the densities set to zero above Z, gives

```math
S\ge\tfrac12 XV-\sum_{d<=W}|r_d|
 \ge\tfrac12 XV-W(1+\log W).                          \tag{1}
```

The last step uses the elementary divisor sum bound

```math
\sum_{d<=W}\tau(d)
 =\sum_{j<=W}\lfloor W/j\rfloor\le W(1+\log W).
```

Only squarefree active products enter the first sum; enlarging to all
integers only weakens the estimate. Unlike truncation by the number of
prime factors, the sieve weights restrict the **product** d. This is the
source of the improved error. Also `V` is bounded below by a fixed positive
multiple of `(1+log Z)^(-2)`. The predecessor's weaker elementary product
bound suffices for positivity as well.

If `X=m^(1-o(1))` and `log Z=o(log m)`, then `W=m^o(1)`,
`XV=m^(1-o(1))`, and (1) is positive. This is an existence estimate;
enumerating the long t interval is not claimed to be a fast factoring
sampler.

## 3. Apply the unchanged prescribed-target alignment

For large `m=3^k`, let D_0 solve `F(D_0)=m/2`. The predecessor proves
`D_0=m^(2-o(1))`, and for every integer `N in [D_0,2D_0]` eventually

```math
m/3<F(N)\le\lceil F(N)\rceil<m.
```

Set `Z_*=ceil(Z(ceil(2D_0)))`, choose `Z_*<q<=2Z_*` by Bertrand,
and impose the progression above. The integer a interval giving
`D_0<=4ma-1<=2D_0` has length `m^(1-o(1))`. Since
`log Z_*=o(log m)`, restricting it modulo q^2 still leaves a t interval
of length `X=m^(1-o(1))`. Equation (1) supplies an a for which both
a and N=4ma-1 have no prime divisor at most Z_*. Eventually `1<a<m`.

The imposed congruence gives `v_q(N)=1`, and `N>q`, so N is composite
and not a perfect power. It survives trial division through Z(N).
The reduced primitive form `(a,-1,m)` is the canonical norm-three ideal's
kth power, as proved in the predecessor. Its only positive norms below m
lie on y=0 and are divisible by rough a; hence none is Y(N)-smooth.
The displayed margin makes k exactly the prescribed least exponent.
No primality of a, unknown class-group index, or smoothness heuristic is used.

## Verification and attribution

The extension is an ordinary application of the sourced fundamental lemma
to the predecessor's exact CRT counts. It does not promote finite controls
to an asymptotic proof. `make verify-near-height-smooth-descent` replays the
existing exact progression/intersection controls, target arithmetic and
both contributed semiprime controls. The original elementary proof remains
available under its stronger cutoff assumption. The authorship and source
pins of those results are retained.

The [question and stopping test](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360248)
registered the previously uncovered full-trial-budget hypothesis before this
extension. Sol worker 1 supplied the original counterfamily and threshold
alignment; this continuation replaces only its sieve estimate, with OpenAI
Codex assistance. No independent mathematical review, new Lean, general
factoring lower bound or stronger complete factoring exponent is claimed.
