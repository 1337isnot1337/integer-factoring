# Contract synchronized Frobenius products jointly

## Research question and result

Can we avoid the representation growth caused by repeated one-variable
elimination? For an explicit product of synchronized Frobenius factors,
the answer is yes: sum all variables in one carry computation.

The earlier [relative Frobenius theorem](relative-frobenius-trace-contraction.md)
contracts one variable. Its output need not retain a small product form.
The construction here does not require that closure property. It keeps
correlations between variables in a joint coefficient table.

This gives a whole-network algorithm for the geometric-box family in
Section 5. It is not a representation of the complete factoring query.
In particular, an additive phase `zeta^(sxy)` is not automatically a
Frobenius factor. The complete factoring exponent remains `1/2`.

## 1. Explicit input and theorem

Let `E=F_(2^d)` be explicitly represented, and let `F_b` be its subfield,
where `b=2^w` and `w` divides `d`. Fix `t` dividing `w`, and put
`q=2^t`, `D=w/t`. There are `v>=1` summation variables.

A factor of row `rho`, with `0<=rho<t`, has the form

```math
A(U_1,\ldots,U_v)
 =\sum_{\mu}a_\mu
   \prod_{i=1}^v\prod_{j=0}^{D-1}
       (U_i^{2^\rho q^j})^{\mu_{i,j}}.
```

The finite list of terms and all nonnegative counts are explicit.
Let `delta(A)=max_mu sum_(i,j) mu_(i,j)`, with zero for a constant
or zero factor. Numerical degree bounds below are parameters of the
cost, not free binary-encoded exponents.

Consider an explicit product

```math
F(U)=c\prod_a N_a(U)\prod_b\iota(L_b(U)),\qquad
\iota(x)=x^{-1}\ (x\ne0),\quad\iota(0)=0.
```

Every numerator and denominator factor has one row. Rows can differ
between factors. For each variable, choose either `F_b` or `F_b^*`
as its summation domain. Define

```math
S=\sum_a\delta(N_a)+d\sum_b\delta(L_b).
```

**Joint contraction theorem.** The sum of `F` over the chosen product
domain has deterministic cost

```math
(S+2)^{O(vD)}\,\mathrm{poly}(L,d,w,v)
```

field operations and polynomial bit overhead per operation. Here `L`
is the written input length. Field construction is separate. The
algorithm does not enumerate the `b^v` variable assignments.

If all factors are constant, evaluate the constant product directly.
A full-field coordinate makes the sum zero in characteristic two.
A nonzero-field coordinate contributes the scalar one.

## 2. Remove inverse values without a pole assumption

For every `x in E`, including zero,

```math
\iota(x)=x^{2^{d+1}-3}
        =\prod_{h\in\{0\}\cup\{2,\ldots,d\}}x^{2^h}.
```

For `d=1`, the product has only the exponent-one factor. For larger
`d`, the exponent is positive and is congruent to `-1` modulo `2^d-1`.
Thus the identity holds at zero and at every nonzero element.

Replace each inverse by these `d` Frobenius powers. On `F_b`, reduce
every binary exponent modulo `w`. A factor of row `rho` becomes a
factor of row `(rho+h) mod t`. All its variables shift together.
Its ordinary degree in the formal Frobenius variables stays unchanged.
The resulting polynomial product `G` therefore has total such degree
at most `S`. In each ordinary variable, its degree is at most
`S*2^(w-1)`.

This replacement is an identity of values on the specified domain.
It does not assert an identity of rational functions outside that domain.
Zero denominators have the specified zero-preserving meaning throughout.

## 3. Cyclic carries compute the trace directly

Group the factors of `G` by their row. At row `rho`, multiply their
explicit polynomials in the `vD` formal symbols `Y_(i,j)`. Call the
result `P_rho`. Its coefficient table has at most `(S+1)^(vD)` positions.

An incoming carry vector `a` and outgoing carry vector `z` have weight

```math
[Y^h]P_\rho,\qquad h_{i,j}=2z_{i,j}-a_{i,j}.
```

Invalid count vectors have weight zero. Multiply these transitions for
`rho=0,...,t-1`. Unlike an exact coefficient query, a field trace has a
cyclic boundary. Retain entries satisfying

```math
z_{i,j}=a_{i,(j+1)\bmod D}.
```

For a full-field variable, also require `a_(i,0)>0`. For a nonzero-field
variable, allow `a_(i,0)=0`.

Here is why this cycle neither loses nor duplicates a monomial. Fix one
term choice per factor. For one variable, let `h_k` be its total count
at binary position `k`, before carries. Its ordinary exponent is
`K=sum_(k=0)^(w-1) h_k 2^k`. A cyclic carry sequence satisfies

```math
h_k+c_k=2c_{k+1},\qquad c_w=c_0.
```

Multiply by `2^k` and sum. All interior carries cancel, giving

```math
K=(2^w-1)c_0.
```

Thus a cycle exists only if `K` is divisible by `b-1`, and its initial
carry is uniquely `K/(b-1)`. Conversely, that integer initial carry
produces an integer sequence by binary addition. Integrality follows
successively from the remaining exponent and the displayed identity.
The final carry equals the initial carry. Each carry is nonnegative
and at most `S`: equivalently, it is the cyclically rotated weighted
sum of the counts, divided by `b-1`.

The cycle is zero exactly when `K=0`. Finite-field power sums retain
positive multiples of `b-1` on the full field. They retain all multiples,
including zero, on the nonzero field. The condition on `a_(i,0)` therefore
implements the correct domain for each variable. Apply this argument to
all variables together. Correlations remain in `P_rho`; no independence
between monomial choices is assumed.

Each term choice contributes exactly once, with its original coefficient,
if all required power sums are nonzero. Linearity proves the joint trace.
There are at most `X=(S+1)^(vD)` carry states. Coefficient tables and naive
multiplication of `t` transition matrices cost `(S+2)^O(vD)` times a
polynomial in the written input and `t`. This proves Section 1.

The cyclic boundary is specific to a trace. For an exact coefficient,
zero initial and terminal carries remain necessary. The implementation
retains that separate coefficient routine as a control. Direct cyclic
traces avoid both the extra overflow columns and the list of coefficient
queries used by the earlier one-variable algorithm.

## 4. Detect the largest common row modulus

Suppose the input factors instead list occurrences `U_i^(2^e)` in
explicit monomials. Form a bipartite incidence graph with variable
vertices `i` and factor vertices `f`. Each occurrence gives an edge
with integer label `e`. Keep parallel edges: different exponents of
one variable in one factor impose separate constraints.

A Frobenius change of variable `U_i=V_i^(2^sigma_i)` preserves both
summation domains. Synchronization modulo `t` requires

```math
\rho_f-\sigma_i=e\pmod t
```

on every edge. Assign integer vertex potentials along a spanning
forest. For every edge, record the residual between its label and
the potential difference. If these residuals are `Delta_e`, put

```math
t_{\max}=\gcd(w,\{\Delta_e\}_e).
```

A modulus `t` dividing `w` is feasible exactly when it divides every
residual. Necessity follows by summing constraints around cycles.
For sufficiency, reduce the forest potentials modulo `t`; all edge
constraints then hold. Hence `t_max` is the largest admissible modulus.
The construction also handles disconnected graphs and isolated vertices.
It preserves the counts in each monomial and requires no field enumeration.

A large `t_max` gives a small relative dimension `D`. A small value
can remove the algorithmic advantage. Synchronization alone is not a
factoring algorithm or a promise that the input has a small factor list.

## 5. A concrete whole-network application

Let `E=F_(b^2)`, with `b=2^w`, and define

```math
K_\gamma(u,v)=\sum_{a=0}^{b-1}(\gamma uv)^a,
\qquad u\in F_b^*,\quad v^{b+1}=1.
```

Take any bipartite graph with `p>=1` left variables, `r` right variables,
and an explicit tag `gamma` on each edge. Give each vertex an explicit multiplicative character weight,
`u_i^a_i` or `v_j^c_j`, with exponents reduced modulo `b-1` or `b+1`.
The task is to sum the weighted product of all edge kernels. This is a
joint network, not a sequence of independent star evaluations.

Choose `omega in E` outside `F_b`, and write `bar omega=omega^b`.
The norm-one group has the exact parametrization

```math
v=\frac{z+\omega}{z+\bar\omega}\quad(z\in F_b),
\qquad v=1\quad(z=\infty).
```

The finite map is injective because `omega != bar omega`. Its values
have norm one and never equal one, so adjoining infinity gives the
whole group. A binary basis supplies `omega` by testing its basis
vectors for membership in `F_b`; no search over the field is needed.

For every field element `x`, including `x=1`,

```math
\sum_{a=0}^{b-1}x^a=(1+x^b)\iota(1+x).
```

At `x=1`, both sides are zero. Substitute the finite parametrization
and use `u^b=u`. The edge kernel becomes

```math
\begin{aligned}
&\bigl(z+\omega+\gamma^b u(z+\bar\omega)\bigr)
 (z+\bar\omega)\\
&\quad\cdot\iota(z+\omega)
 \cdot\iota\bigl(z+\bar\omega+\gamma u(z+\omega)\bigr).
\end{aligned}
```

Both `z+omega` and `z+bar omega` are nonzero. The last denominator
can vanish; the displayed identity still holds with `iota(0)=0`.
Every factor is a polynomial of degree at most two in row zero,
with `t=w` and `D=1`. At infinity the edge instead becomes
`(1+gamma^b u) iota(1+gamma u)`.

The boundary weights also have short factor representations. Expand each
binary exponent into its nonzero bits. A left bit contributes one factor
`u_i^(2^h)`, of row `h mod w`. At a finite right vertex, a bit contributes

```math
(z^{2^h}+\omega^{2^h})\,
\iota(z^{2^h}+\bar\omega^{2^h}).
```

Reduce the variable Frobenius exponent modulo `w`, but retain the actual
coefficient powers in `E`. At infinity the character value is one.
This uses `O(w)` factors per vertex, each of degree one.

Enumerate the `2^r` choices of infinity vertices, not the field
assignments. Each patch is one product covered by Section 1, with
at most `p+r` variables and
`S=O((w+1)*number_of_edges+w(w+1)(p+r))`. The domains
are nonzero fields on the left and full fields for finite right vertices.
Thus the complete cost is

```math
2^r(S+2)^{O(p+r)}
\,\mathrm{poly}(L,w,p,r).
```

This is polynomial in `S,w` to an exponent `O(p+r)`, rather than
enumeration of roughly `2^(w(p+r))` assignments.
The field remains an explicit input. Its construction cost is not hidden.

### Why boundary weights matter

The unweighted network has a simple degeneracy when its number of edges
is less than `b+1`: its sum is exactly one. Expand each edge as a digit
sum. The left character sums require the total digit sum to be divisible
by `b-1`. The right sums require divisibility by `b+1`. These odd numbers
are coprime. The total digit sum is less than `b^2-1`, so it is zero.
All edge digits then vanish, and the constant contribution is one.

Thus unweighted small graphs are controls, not evidence of a useful hard
instance. Boundary characters change these congruences and give nonconstant
answers. For example, on the complete two-by-two graph with `b=4`, choose
left charges `b-2` and right charges `b`. Only the two perfect matchings
survive. The sum is `gamma_00 gamma_11 + gamma_01 gamma_10`.
The checker includes this coupled cycle, not only graphs with leaves.

## 6. Verification and the next obligation

The implementation and exact controls are:

- [Joint carries and synchronization](../tools/joint_frobenius_carries.py).
- [Coefficient and trace controls](../tools/check_joint_frobenius_carries.py).
- [Whole geometric-box networks](../tools/joint_geometric_box_network.py).
- [Network controls](../tools/check_joint_geometric_box_network.py).

Run `make verify-joint-frobenius` to reproduce the
[joint trace receipt](../evidence/joint-frobenius-carries.json) and the
[weighted network receipt](../evidence/joint-geometric-box-networks.json).
These are exact finite controls and an ordinary written argument, not
Lean formalization or independent mathematical acceptance.

The controls compare with ordinary polynomial multiplication, direct
finite-field sums, and exhaustive small synchronization assignments.
They include coupled monomials, mixed domains, zero denominators,
parallel incidence constraints, and the degree-one field.

The next research obligation is an exact representation of a relevant
arithmetic query by a controlled number of these products. Existing
[carry Fourier elimination](carry-fourier-elimination.md) still leaves
`zeta^(sxy)` and bounded operand domains. This note does not turn
those phases into fixed Frobenius monomials. A general arithmetic
circuit is not an explicit product. Expanding one can erase the gain.

Continue this route only if a concrete representation preserves those
correlations at a charged cost. The present result bypasses repeated-star
closure for the stated network class; it supplies no new factoring bound.
