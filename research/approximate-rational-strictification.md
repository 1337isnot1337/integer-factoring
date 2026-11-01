# Globally type-removable rational poles have polynomial pressure limits

## Result and research question

Extra finite scalar poles do not enlarge the class of strict negative
pressure witnesses when every paired-type word is already globally
Laurent-regular over a finite field. Rounded paired gauges and Frobenius
substitution produce a strict polynomial representative on the same
alphabets and field, with exactly charged margin error `4Lt/p^m` tending
to zero. Every strict witness therefore survives at some finite power.
The achieved factoring score remains `1/2`.

The [existing strictification](typed-rational-crossing-strictification.md)
removes endpoint poles exactly but assumes no other finite poles. At
other places only weak valuation inequalities are available, and their
paired potentials may have denominators prime to the characteristic.
The new research question was whether rounding those potentials could
avoid this denominator obstruction without erasing a strict pressure gain.

The exact diagnostic below uses a characteristic-three kernel with
necessary half-integer potentials at an extra pole. It verifies strict
endpoints, complete compensated profiles and the `4/D` margin penalty.
No favorable tensor or numerical pressure search is claimed. Matrix-valued
boundary contractions and selected-type-only constructions remain outside
the reduction.

## 1. Statement

Let k be a fixed finite field of characteristic p. Let K(w) be a finite
scalar rational tensor on fixed paired alphabets A,B,B,A. Assume that,
for **every tensor length and every paired full-letter-type restriction**,
each word coefficient is a Laurent polynomial in w and has the crossing
endpoint at w=1. Thus all of its other finite poles are removable on
every such word, not only on a selected list of types.

Let P_j be the distinct monic irreducible factors of the entries'
denominators other than w and w-1, and put `L=sum_j degree(P_j)`.
For all sufficiently large `D=p^m`, there is an explicit **strict
Laurent-polynomial crossing K_D(z) over the same field and alphabets**.
For every choice of original physical heights, compensating heights on
K_D have complete weighted orders and restored physical energies

```
r_D=D*r-4L,     s_D=D*s,     G_X,D=D*G_X,     G_Y,D=D*G_Y,
```

where r is the infinity cut and s the zero cut. Consequently its critical
margin at scale t/D and unchanged supporting weights is exactly

```
M_D=M+4L*t/D.                                           (1)
```

Every **strictly negative** margin witness in this rational class thus
has a strict polynomial representative on the same alphabets and fixed
finite field. This is approximation with a charged error, not equality
of each finite-D pressure. It introduces no asymptotic factoring gain.
All construction parameters depend on the finite candidate, not on the
integer being factored.

As a corollary of the [uniform axis window](uniform-polynomial-axis-window.md),
no negative critical margin in this rational class exists when a supporting
weight is at most `exp(-2048*(A*B)^2)`. In general this argument gives
nonnegative rational margin, rather than strict positivity at a limit.
When L=0 the earlier exact strictification gives strict positivity.

## 2. Weak potentials at the other poles

For an entry lambda=(x,y,u,v), use the usual charge vector
`a_lambda=(e_x-e_v,e_y-e_u)`. At a closed place P_j write nu_lambda,j
for its integer valuation. The all-type Laurent-regularity assumption
is exactly the weak circuit condition

```
sum n_lambda*nu_lambda,j >= 0
whenever n_lambda>=0 integers and sum n_lambda*a_lambda=0.
```

Each balanced multiset is a word with paired types, and a product of
nonzero scalar entries has the sum of their valuations. There is no
cancellation between different entries in this scalar word coefficient.
Farkas' lemma, with rational data, therefore gives rational paired
potentials f_j,g_j satisfying

```
nu_lambda,j+f_j(x)+g_j(y)-g_j(u)-f_j(v) >= 0.             (2)
```

At the endpoint w=1 the existing strict circuit argument gives potentials
f_1,g_1 with strict inequality on every nonzero offmatched entry. Matched
entries are regular there with value one, and their charge is zero.

Unlike that strict endpoint system, (2) can force denominators not equal
to powers of p. This is why simply invoking the previous exact proof
at every pole would be unjustified.

## 3. Frobenius rounding and its complete correction

Set D=p^m and round D*f_j,D*g_j down to integer functions F_j,G_j.
For every entry, the rounding error of its four signed leg potentials
has absolute value less than four. Thus (2) implies

```
D*nu_lambda,j+F_j(x)+G_j(y)-G_j(u)-F_j(v) >= -4.
```

Round the endpoint potentials too. For large enough D their finitely
many strict inequalities remain positive, without any scalar correction
at one.

Since k is finite, coefficient Frobenius is an automorphism. Let Q_j,D
be obtained from P_j by applying its inverse mth Frobenius to each
coefficient. Then `P_j(z^D)=Q_j,D(z)^D`. These polynomials have the same
degrees as P_j, still over k, and have nonzero values at zero and one.
No coefficient-field extension or extension-state index is introduced.

Multiply K(z^D) by the paired physical-leg gauges

```
(1-z)^(F_1(x)+G_1(y)-G_1(u)-F_1(v))
* product_j Q_j,D(z)^(F_j(x)+G_j(y)-G_j(u)-F_j(v)),
```

and by the common normalized correction

```
C_D(z)=product_j (Q_j,D(z)/Q_j,D(1))^4.
```

Every extra pole is now removed. The strict endpoint inequalities remove
the pole at one and kill every offmatched entry there. Matched gauges
have exponent zero, and C_D(1)=1, so the endpoint is literally crossing.
Only powers of z can remain in denominators. Hence K_D is a strict
Laurent-polynomial tensor. A common monomial may make it polynomial with
its usual exact two-place scalar cancellation.

The correction has degree **4L**, independently of D or original pole
orders. It is a unit at zero. This fixed correction is the entire added
cost; it must not be omitted at infinity.

## 4. Profile and physical-cost readback

At zero, all gauge polynomials just used are units. Give each leg the
new height `D*h_0`. Its complete zero profile is D*s.
At infinity, give each leg the height

```
D*h_infinity + endpoint_leg_potential
             +sum_j degree(P_j)*extra_pole_leg_potential.
```

The four leg potentials are respectively F,G,-G,-F. They cancel the
infinity orders of the paired gauge factors. Their matched X,V and Y,U
sums vanish identically, so the restored physical energies are exactly
D times the original lists. The remaining common correction shifts
every infinity order by -4L. This proves all displayed profile identities
and (1), with every degree cost retained.

If M<=-eta<0 is an established witness, choose a finite p-power D large
enough for the endpoint and `4Lt/D<eta`. Then M_D<0. The potentials,
rounding, inverse coefficient Frobenius and polynomial operations are
all effective finite constructions. This proves the negative-witness
reduction for the full two-place criterion, rather than only one rank,
one leading coefficient, or determinant means.

## 5. Exact diagnostic and limits

The [F3 control](../tools/check_approximate_strictification.py) uses extra
pole w=2 and four offmatched charges (+,+),(-,-),(+,-),(-,+), with
valuations -1,+1,0,0. The weak inequalities force f=g=1/2. At w=1 each
offmatched entry also contains w-1, giving the all-type crossing endpoint.
Balanced products have zero extra-pole valuation by the exact identity
`nu+(a_X+a_Y)/2=0`, proving all-length regularity there.

For D=1,3,9,27, rounding H=floor(D/2) leaves an actual pole of order one
on one local entry. The corrected strict representative is constructed
in full. The [receipt](../evidence/approximate-strictification-controls.json)
records all its polynomials, 64 literal endpoint checks and all 40
complete orders. Those orders are independently computed from **690
untruncated minor determinants**, with no native series window. They
satisfy the exact D scaling and -4 infinity shift, hence the retained
penalties 4,4/3,4/9,4/27 per unit original scale.

The result requires all paired-type words to be globally Laurent-regular
and the scalar product rule. A selected sequence of typical types,
matrix-valued auxiliary contractions, and characteristic-zero kernels
are outside the argument. It does not say uncharged arbitrary poles are
free or automatically removable. It is a reduction to the still-open
interior polynomial crossing problem, not an obstruction to all pressure
witnesses or a stronger factoring exponent.

```
make verify-approximate-strictification
```

This extends the credited typed-rational strictification, its charge
circuits and paired gauge accounting. The new step is rounded weak
potentials at extra closed places, with a bounded correction and a
vanishing normalized pressure cost. Directed by Cinna Davis / 1337isnot1337
with substantial OpenAI assistant support. No new Lean or independent
mathematical acceptance is claimed.
