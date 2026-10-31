# Exact representation search for one factored candidate norm

The [smooth-norm relation interface](smooth-norm-relations.md) needs a
primitive smooth represented norm. For one supplied factored candidate
`n`, this note gives an exact deterministic decision and witness algorithm
in

```math
Y\,\mathrm{poly}(b)+2^{\lceil\omega(n)/2\rceil}\mathrm{poly}(b)
=Y\,\mathrm{poly}(b)+2^{O(b/\log b)}
```

bit operations, where the primes of `n` are at most `Y`, `gcd(n,Delta)=1`,
and the arithmetic inputs have `O(b)` bits. Supplied splitting roots remove
the `Y` term. A positive result includes a primitive pair `(x,y)` with
`Q(x,y)=n`; a negative result proves that this one candidate has no such
primitive representation. It says nothing about other candidate norms.

The algorithm uses the classical prime-ideal orientation characterization
and meet-in-the-middle. No new representation-theory claim is made. The
repository contribution is an explicit integer algorithm, polynomial
composition without a search for coprime leading coefficients, complete
bit-cost accounting, and independent exact controls.

Selecting a suitable candidate `n` from a short deterministic list remains
open. No improved complete factoring score follows; the score remains
`1/2`.

## 1. Polynomial composition by two extended gcds

The older [small-form package](../tools/binary_quadratic_forms.py) finds a
coprime equivalent leading coefficient by a finite scan. Its documentation
correctly declines a polynomial bit-cost claim for that scan. The new
[arithmetic module](../tools/quadratic_form_arithmetic.py) instead multiplies
the two ideal lattices directly.

Write the input forms as `(a,b,c)` and `(A,B,C)`, of the same negative
discriminant `Delta`, and put

```math
s=\Delta\bmod2,\qquad t=(s^2-\Delta)/4,
\qquad u=(b+s)/2,\qquad v=(B+s)/2,
\qquad w=s-u-v.
```

In the order basis `1,omega`, where `omega=(s+sqrt(Delta))/2` and
`omega^2=s*omega-t`, the two proper ideals have bases
`(a,0),(-u,1)` and `(A,0),(-v,1)`. Their product lattice has the four
generators

```math
(aA,0),\quad(-av,a),\quad(-Au,A),\quad(uv-t,w).
```

Let `d=gcd(a,A,w)`. It divides every coordinate: the only less immediate
one follows from `v == s-u mod d` and
`u^2-su+t=ac`. It is exactly the integer content because it is already
the gcd of the second coordinates.

The six two-by-two minors, after removal of the common factor `aA`, are,
up to sign,

```math
a,\ A,\ w,\ u-v,\ C,\ c.
```

Their gcd is one. Any common prime divisor would divide `a,c` and
`w-(u-v)=-b`, contradicting primitivity of the first form. Thus the product
lattice has index `aA`, and division by its content gives index
`a'=aA/d^2`.

Two extended gcds give coefficients `r_1,r_2,r_3` with

```math
r_1a+r_2A+r_3w=d.
```

The corresponding first coordinate is

```math
X=-r_1av-r_2Au+r_3(uv-t).
```

After dividing by `d`, this is a lattice vector `(x_0,1)`, where `x_0=X/d`.
A rank-two integer lattice of index `a'` with surjective second coordinate
has basis `(a',0),(x_0,1)`. Therefore the primitive product form has

```math
a'=aA/d^2,\qquad
b'\equiv-2x_0-s\pmod{2a'},\qquad
c'=(b'^2-\Delta)/(4a').
```

Center `b'` and reduce this form. The standard correspondence between
primitive forms and proper invertible ideals proves that its class is
the product class and that the form remains primitive. This formula
works when the input leading coefficients are not coprime, including
shared prime powers and ramified reduced representatives.

All these integers have polynomially many bits in the input length.
Extended gcd, the displayed formulas and reduction have polynomial bit
cost. There is no loop through residues modulo `a` or `A`.

## 2. Reduction with an explicit equivalence matrix

The module returns both the reduced form and a matrix `U` in `SL_2(Z)`
such that `Q_reduced(z)=Q_input(Uz)`.

Each centering step uses the shear `(x,y) -> (x+ky,y)`, and each swap uses
`(x,y) -> (-y,x)`. Their matrices are multiplied into `U` as the algorithm
runs. The final boundary sign convention is handled by the same proper
swap when the first and third coefficients agree.

There are polynomially many arithmetic steps. More concretely, consider
a centered form `(a,b,c)` with `|b|<=a` and `c<a`, so a decreasing swap is
performed. If `c>a/2`, the next centering shift is only `-1`, `0` or `1`.
Its new trailing coefficient is either `a` or `a-|b|+c`, both at least
`c`. Consequently no further decreasing swap occurs. Every decreasing
swap except possibly the last therefore reduces the leading coefficient
by at least a factor two. There are `O(log(a+1))` such swaps.

The centered coefficients have polynomial bit length in the original
coefficients and discriminant. The product of polynomially many shear
matrices also has polynomial bit length. This gives a direct bit-cost
justification, including the returned equivalence matrix.

## 3. Prime-power orientations

Suppose

```math
n=\prod_{i=1}^k p_i^{e_i},\qquad \gcd(n,\Delta)=1,
```

with distinct primes `p_i`. If one `p_i` is inert, no primitive norm-`n`
ideal exists, so there is no primitive representation. This remains true
when the inert prime has even exponent: the possible representations
then carry integer content and are nonprimitive.

Otherwise each prime splits into two conjugate ideals. A primitive ideal
of norm `n` must put all `e_i` copies of `p_i` on one side of its conjugate
pair; using both would introduce rational prime content. Thus the complete
candidate family has exactly two orientation choices per distinct prime.

The canonical prime root can be found by testing integers from `0` through
`p_i` for `r^2 == Delta mod 4p_i`. This costs at most `Y*poly(b)` over the
prime list when `p_i<=Y`. If roots are supplied, no search is needed.

For odd primes, ordinary Hensel lifting gives the unique root branch
modulo `p_i^e_i`, and the required parity gives the root modulo
`2p_i^e_i`. For `p_i=2`, start with the orientation modulo four. A root
`r` modulo `2^(j+1)` satisfying `r^2 == Delta mod 2^(j+2)` is lifted by
testing the next bit, equivalently adding

```math
\left(\frac{\Delta-r^2}{2^{j+2}}\bmod2\right)2^{j+1}.
```

This retains the chosen root modulo four and gives the next precision.
All root lifting costs are polynomial in the bit length.

Let `F_i` be the reduced form of norm `p_i^e_i` from one lifted branch.
The complete class condition is

```math
[Q]=\prod_{i=1}^k[F_i]^{\sigma_i},\qquad \sigma_i\in\{-1,1\}.
```

## 4. Deterministic meet-in-the-middle search

Split the `k` prime powers into two groups of sizes `floor(k/2)` and
`ceil(k/2)`. Enumerate their signed products with the polynomial composer.
For each canonical class retain one orientation mask; replacing duplicate
classes by a representative loses no possible continuation.

For a left class `L`, the needed right class is `L^(-1)Q`. Sort the right
classes lexicographically by their exact integer coefficients, and use
binary search. The implementation also sorts and removes duplicates at
each half-enumeration step. It does not rely on an average-case hashing
claim for the deterministic bound.

Each half has at most `2^ceil(k/2)` states. Sorting, class multiplication,
inversion and comparisons contribute only polynomial factors in `b`,
since `k<=b`. Time and memory are therefore

```math
2^{\lceil k/2\rceil}\mathrm{poly}(b).
```

If the match fails, all primitive orientations have been covered, so
there is no primitive representation of this candidate `n`.

## 5. Recovering an integer vector

A successful pair of masks chooses one lifted root at every prime.
Chinese remaindering gives `B_n` modulo `2n` with

```math
B_n^2\equiv\Delta\pmod{4n}.
```

At an odd prime the CRT modulus is `p_i^e_i`; the parity condition is
included separately. When two divides `n`, its CRT component instead has
modulus `2^(e_2+1)`. This is necessary to retain the stronger form
congruence, rather than checking only modulo `n`.

Reduce

```math
F=\left(n,B_n,\frac{B_n^2-\Delta}{4n}\right)
```

and keep its equivalence matrix `U`. The successful class search ensures
that its reduced form is the reduced target. If the target was already
reduced, then

```math
(x,y)=U^{-1}(1,0)=(U_{22},-U_{21})
```

is the required primitive representation. For a nonreduced input target,
compose this vector with the target's own reduction matrix. The program
checks `gcd(x,y)=1` and the exact equation `Q(x,y)=n` before returning it.

This closes the search for coordinates once the factored norm is supplied.
It uses no scan through the ellipse `Q(x,y)<=n`.

## 6. The logarithmic cost and its exact scope

An elementary bound suffices to control `k=omega(n)`. For an integer of
at most `b` bits, there are at most `sqrt(b)` prime divisors at most
`sqrt(b)`. Each larger prime contributes at least `(log_2 b)/2` bits to
the product, so there are at most `2b/log_2 b` of them. Hence

```math
\omega(n)\le\sqrt b+2b/\log_2 b=O(b/\log b).
```

Small input lengths are covered by an absolute constant. Together with
the previous sections this proves the stated `2^O(b/log b)` candidate
test, with the root-search cost separately charged. For the intended
`Y=exp(O(sqrt(b log b)))`, the orientation search is the larger displayed
asymptotic term.

This is a score-1-scale component, not a score-1 factoring algorithm.
To use it in the relation framework, one must still construct a sufficiently
short list of candidate smooth norms, guaranteed to include a represented
one for at least one member of the biased target family. The Fourier
argument supplies an existing norm without identifying that list.

Enumerating every smooth norm is not the same as the orientation search
for one norm. Nor does a word walk in the existing factor-base subgroup
supply the missing global ambiguous directions. Both distinctions remain
part of the arithmetic research contract.

## 7. Exact evidence and replay

Run:

```sh
python3 tools/check_factored_norm_search.py
```

The [record](../evidence/factored-norm-search.json) checks:

- `58,131` form products for `|Delta|<1000` against the independent older
  composition algorithm, including `15,625` products with nontrivial
  ideal content;
- all product-lattice indices and normalized lattice bases in those
  cases, directly through two-by-two minors and divisibility;
- `22,195` reduction matrices, including nonreduced shears, by substituting
  the matrices into the original form;
- `26,196` candidate-norm decisions across `132` discriminants against
  independent coordinate enumeration, including `23,154` negative cases
  and `11,109` inert-prime rejections;
- `500` nonreduced target recoveries;
- a `629`-bit discriminant and a norm with `24` distinct primes: the
  complete orientation family has `16,777,216` members, while the two
  retained halves have `4,096` states each;
- a `2,074`-bit discriminant and a `312`-bit norm containing `150` prime
  occurrences but only eight distinct primes, exercising high prime-power
  lifting rather than treating every occurrence as an orientation choice.

The implementation is [factored_norm_search.py](../tools/factored_norm_search.py),
and its new arithmetic dependency is
[quadratic_form_arithmetic.py](../tools/quadratic_form_arithmetic.py).
The finite controls support the explicit algorithm and its scope; the
all-input cost is proved above. The work builds on the smooth-norm
checkpoints integrated at `601735e0a5c5d3cc7a5f19b0450af246a2a9ed05`, uses
substantial OpenAI assistant support, and adds no Lean code.
