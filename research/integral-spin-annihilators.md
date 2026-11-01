# The integral bad-prime spin limit disperses its Smith budget

Clearing denominators before reduction gives twelve-letter spin-eleven
crossings at the odd bad primes `3,5,7,11`. The characteristic-three seed
is exactly an old three-letter seed with four spectator labels. At the
other primes, both cut matrices admit polynomial inverse witnesses

```
M(u) A(u) = (u^p-u)^2 I.                                  (1)
```

Thus every local Smith order is at most two, at every algebraic place,
while their charged degrees are `4,6,10`. Every bare two-place product has
`rho>=1`; the bare critical asymmetric supporting-line test fails too.
The complete factoring score stays **1/2**.

This bounded diagnostic asked whether nonsemisimple modular extensions
could supply long delayed channels. In these seeds they instead disperse
the budget among several roots. This ends their bare two-place lane; no
larger-spin, adjacent-prime or thermal-scale sweep follows. It does not
exclude all integral spin constructions, physical types, joint endpoint
inversions, multipoint opposite flags, or factoring algorithms.

## 1. Combine over the integers before reducing

Use the explicit [rational Casimir construction](modular-spin-crossings.md)
`R_n(u)` at `n=11`. Combine the rational projectors first. The actual least
common denominator of the resulting coefficients is `11!=39916800`; put

```
F(u) = 11! R_11(u).
```

The source calculation checks integrality, integer content one,
`F(0)=11! swap` and leading coefficient `I`. Reducing individual projectors
at a bad prime would be invalid. Direct specialization of the earlier
six-vertex parameter to one would lose the confluent information.

Reduce the integer coefficients at each listed prime, compute their actual
monic polynomial gcd `g_p`, and divide every entry exactly. The result
`F_p` still has leading coefficient `I`. Define

```
T_p(u;x,y,u',v) = F_p(u;x,y,v,u').                         (2)
```

Its infinity fibre is the literal strict crossing. The substitution
`u=z/(1-z)` with homogeneous multiplier `(1-z)^D` puts that endpoint at
one without increasing degree. The two cuts are `(x,y)|(u',v)` and
`(x,u')|(y,v)`.

| Prime | Removed gcd degree | Actual degree D | Outcome |
| ---: | ---: | ---: | --- |
| 3 | 9 | 2 | Exact three-letter source times spectator crossing |
| 5 | 7 | 4 | Both inverse witnesses have degree at most 6 |
| 7 | 5 | 6 | Both inverse witnesses have degree at most 8 |
| 11 | 1 | 10 | Both inverse witnesses have degree at most 12 |

For `p=3`, write every label `x=3X+x_0`. The full checked identity is

```
F_3(x,y,u',v) = (2 R_2)(x_0,y_0,u'_0,v_0) [X=U][Y=V].     (3)
```

After (2), the coarse factor is a spectator strict crossing. Its letters
are not free. This is a literal coefficient identity with the same spectral
variable and no hidden basis transformation.

## 2. A polynomial inverse controls every local order

**Lemma.** If square polynomial matrices satisfy `M A=Q I` and `M` has an
invertible leading homogeneous fibre, each finite local Smith order of
`M` is at most `ord(Q)` there. At infinity all orders are zero.

Work over the local DVR after extending the constant field and write
`M=U diag(t^eta_i) V`. The displayed identity gives an integral right
inverse up to `Q`. Each diagonal channel divides `Q`, so
`eta_i<=ord(Q)`. The infinity claim follows from the invertible leading
coefficient. This argument uses no finite jet precision.

Here `u^p-u` has exactly the `p` distinct base-field elements as roots,
with derivative `-1`. Equation (1) gives

```
0 <= eta_i <= 2 at base-field finite places;
eta_i = 0 at every other algebraic place and at infinity.   (4)
```

These statements persist under constant-field extension and projective
coordinate change, using the homogeneous degree `D` trivialization.
Constant paired physical bases preserve the profiles.

Charge conservation splits each cut into 23 disjoint square blocks.
The packet supplies `A` for every block. Generation uses a reversed formal
inverse at infinity; verification only multiplies the complete polynomial
matrices and checks (1) coefficient by coefficient. The source gcd divisions
and strict leading endpoint are checked separately.

## 3. The charged degree excludes the bare target

At any two selected places let `r_i,c_j` be the unweighted local orders.
There are `n=q^2=144` in each profile. By (4), `r_i+c_j<=4<=D`.
For every common thermal scale `s>1`,

```
B(s) = s^D (sum_i s^(-r_i))(sum_j s^(-c_j))
     = sum_(i,j) s^(D-r_i-c_j) >= n^2 = q^4,
rho_bare = log(B)/(4 log q) >= 1.                         (5)
```

For the [critical asymmetric supporting line](asymmetric-support-line-crossings.md),
take positive `a,b` with `a^2+b^2=1` and `t>0`. Without physical heights
or types, its margin obeys

```
D t + a log(sum_i exp(-r_i t/a))
    + b log(sum_j exp(-c_j t/b)) - 2 log q
 >= (D-4)t + 2(a+b-1)log q > 0.                          (6)
```

The required strict negative margin is impossible. This does not account
for the different coefficient bounds of charged types or multiple opposite
flags, and makes no assertion that their costs have been paid.

## 4. Reproduction and stopping point

Run `make verify-integral-spin-annihilators`. The
[checker](../tools/check_integral_spin_annihilators.py) reconstructs the
combined rational source, checks all four integral reductions and (3), then
verifies 138 block inverse identities in the
[packet](../evidence/integral-spin-annihilators.json.xz). The
[receipt](../evidence/integral-spin-annihilator-controls.json) pins the source
and records 6,936 polynomial inverse entries, 20,736 spectator entries and
82,944 strict endpoint entries. To regenerate:

```
python3 tools/check_integral_spin_annihilators.py --generate
```

The initial native Smith diagnostic discovered the witness. The conclusion
uses neither a numerical pressure optimization nor that diagnostic's jet
window. No larger parameter sweep is pending.

The existing rational construction and its source attribution are retained.
The new content is the integral bad-prime limit, exact spectator factor and
all-place polynomial inverse certificate. These are ordinary arguments with
substantial OpenAI assistant support; no new Lean formalization or independent
mathematical review is claimed.
