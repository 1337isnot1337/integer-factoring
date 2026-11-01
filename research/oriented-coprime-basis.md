# Gcd-free oriented ideal refinement preserves useful torsion

## Result and scope

Given finitely many odd unramified quadratic-ideal keys `(R_i,B_i)`, one
can refine their norms into pairwise coprime integer atoms and express
every old ideal as an exact oriented atom word. Substitution into genuine
relation rows preserves the old presentation's evaluated torsion image.
It can expose additional useful involutions. A four-witness control below
extracts 1009 from 1022117 although the old residual norms are all distinct
and their unit-definition presentation is torsion-free.

Factoring supplied integers into a coprime basis is classical; see
D. J. Bernstein, [*Factoring into coprimes in essentially linear time*](https://cr.yp.to/coprimes.html),
Journal of Algorithms 54 (2005), 1--30. The present implementation uses
elementary repeated-gcd refinement and claims polynomial bit cost. This
note supplies its quadratic-root orientation layer, literal ideal replay,
presentation substitution theorem, and strict factor-gain control.

This is a transformation of supplied finite witnesses. It gives no bound
on their search cost and does not change the complete factoring score `1/2`.

## 1. The contract

Fix a negative order discriminant Delta, congruent to zero or one modulo
four. An input key satisfies

$$R_i\ge1,\qquad\gcd(R_i,2\Delta)=1,\qquad
B_i^2\equiv\Delta\pmod{4R_i}.$$

Normalize B_i modulo `2R_i` and write

$$I(R_i,B_i)=\left(R_i,\frac{-B_i+\sqrt\Delta}{2}\right).$$

The output consists of pairwise coprime integers `d_a>1`, valid roots C_a,
and signed integer words satisfying

$$R_i=\prod_a d_a^{|w_{ia}|},\qquad
[I(R_i,B_i)]=\prod_a[I(d_a,C_a)]^{w_{ia}}. \tag{1}$$

The atoms can remain composite or be prime powers. A sole key of norm 39
remains an atom of norm 39. A repeated-power control retains atom 169.
Neither norm factorization nor residual primality testing is performed.

## 2. Norm refinement

Initially attach the i-th coordinate vector to R_i. A piece `(a,u)` means
that a contributes exponent u_i to the i-th input norm. The invariant is
`R_i=prod_(a,u) a^u_i`. For pieces `(a,u),(b,v)` with `g=gcd(a,b)>1`, replace
them by

$$ (g,u+v),\qquad(a/g,u),\qquad(b/g,v), \tag{2}$$

discarding factors equal to one. This preserves every input norm. Repeat
until the piece norms are pairwise coprime.

The potential `sum log_2(a)` over current pieces decreases by `log_2(g)`
at every replacement. It starts below the total norm bit length S. Thus
there are fewer than S replacements and at most S nontrivial pieces at
any time. All exponents are at most S because `a^u_i` divides R_i.
The simple pair-rescanning implementation uses `O(S^3)` gcd decisions on
`O(S)`-bit integers and polynomially many divisions and vector operations.
No essentially-linear bound is claimed for this implementation.

## 3. Root-orientation refinement

For one coprime atom d, choose a used input root as C modulo `2d`, and
choose its canonical conjugate. For every other used root B, compute

$$d_+=\gcd\left(d,\frac{B-C}{2}\right),\qquad
d_-=\gcd\left(d,\frac{B+C}{2}\right). \tag{3}$$

The halves are integral because the roots have the same parity. Every
prime-power factor of d divides `(B-C)(B+C)/4`. No odd prime divides both
halves: it would divide B, C and Delta. Therefore

$$d_+d_-=d,\qquad\gcd(d_+,d_-)=1. \tag{4}$$

If both factors are nontrivial, replace `(d,u)` by `(d_+,u),(d_-,u)`.
This preserves pairwise coprimality with the other atoms and every input
norm. Repeat until every used root is uniformly C or uniformly -C on
each atom. The sign gives `w_ia=+u_i` or `-u_i`.

Each proper orientation split increases the atom count and preserves its
product. There are fewer than S such splits, so the orientation stage
also has polynomial bit cost. All root reductions and multiplicities
have polynomial written size.

At discriminant -23, the nonconjugate keys `(39,17)` and `(39,35)` leave
one norm atom 39 before this stage. Orientation refinement splits it into
3 and 13, with canonical roots 1 and 9. The signed words become
`(-1,-1)` and `(-1,+1)`. A single global conjugation cannot describe this
mixed root pattern.

## 4. Exact ideal identities

At every prime dividing an atom, the unramified root lifts uniquely to
higher prime powers. Thus the selected sign specifies the same local
ideal orientation at every exponent occurring in an input norm. Since
atom norms are pairwise coprime, this proves the literal integral identity

$$I(R_i,B_i)=
\prod_{w_{ia}>0}I(d_a,C_a)^{w_{ia}}
\prod_{w_{ia}<0}\overline{I(d_a,C_a)}^{-w_{ia}}. \tag{5}$$

The correctness proof uses the existence of prime factorizations; the
algorithm computes gcds, divisions and ideal products. The checker
verifies (5) by full two-dimensional HNF multiplication.

Since `I(d,C) conjugate(I(d,C))=(d)`, equation (5) gives (1). As a
fractional-ideal equation with signed powers, the corresponding formula is

$$I(R_i,B_i)=
\left(\prod_{w_{ia}<0}d_a^{-w_{ia}}\right)
\prod_a I(d_a,C_a)^{w_{ia}}. \tag{6}$$

The rational-integer factor in (6) is explicit. Equality of ideal classes
does not silently replace a literal ideal identity.

## 5. Torsion-image preservation

Let `phi:Z^s -> G` evaluate the old fixed-base and residual generators,
and let L be the lattice of supplied genuine relation rows. Substitution
defines an integer map `T:Z^s -> Z^t`. If psi evaluates the new generators,
then (1) gives `phi=psi T`. Retaining all transported rows defines lattice
`T(L)` and a commuting map of the presented quotients. A torsion class
maps to a torsion class with the same image in G. Hence

$$\phi(\mathrm{Tor}(\mathbb Z^s/L))
\subseteq\psi(\mathrm{Tor}(\mathbb Z^t/T(L))). \tag{7}$$

Transported rows can become dependent. The existing relation-halving
extractor accepts independent rows, so one can retain a maximal rationally
independent subset. Its integral row lattice K may differ from T(L), but
their rational spans agree. For any genuine relation lattice A, the
evaluated torsion image is the evaluation of

$$\mathrm{Sat}(A)=(A\otimes\mathbb Q)\cap\mathbb Z^s. \tag{8}$$

Equal rational row spans have equal saturations. Thus this independent-row
selection preserves the evaluated torsion image exactly. It does not
assert equality of the integral presentation lattices. The implementation
retains all transported rows; the checker records the independent subset
used for extraction separately.

In particular, useful old involutions remain in the new torsion image.
The exact relation-halving routine can expose them without a group order
or determinant. Additional torsion images can appear.

## 6. Four distinct residual norms expose a factor

Take `N=1022117`, `Delta=-4N`, and base prime forms for 2, 3 and 7.
Four supplied powers have these exact reduced forms:

| Power | Reduced form | Residual key |
| --- | --- | --- |
| `P_7^9` | `(418,98,2451)` | `(209,98)` |
| `P_7^10` | `(561,452,1913)` | `(187,78)` |
| `P_3^21` | `(874,814,1359)` | `(437,60)` |
| `P_7^45` | `(391,-78,2618)` | `(391,78)` |

Each leading norm is coprime to N and no selected target is an involution.
The residual norms are distinct, so the witnesses contain no equal or
conjugate residual-key pair. Their old four relation rows have a signed
identity block in the residual columns. The old quotient is free of rank
three, so its presented torsion image is trivial.

Gcd refinement yields atoms 11, 17, 19 and 23, with canonical roots 10,
10, 16 and 14. Their signed class words are

| Residual class | Atom word |
| --- | --- |
| `J_(209,98)` | `P_11 P_19^-1` |
| `J_(187,78)` | `P_11^-1 P_17` |
| `J_(437,60)` | `P_19^-1 P_23` |
| `J_(391,78)` | `P_17 P_23^-1` |

The original power identities therefore imply

$$P_2^{-2}P_3^{22}P_7^{64}
=J_{209}J_{187}J_{437}^{-1}J_{391}^{-1}=1. \tag{9}$$

All atom terms cancel. The even base relation is `(-2,22,64)`. Its half
`(-1,11,32)` evaluates to the ambiguous form `(1009,0,1013)`, extracting
the factor 1009. This is a strict gain for the supplied witness set.
The initial discovery used a longer finite power scan and removed
unnecessary witnesses; it is not a held-out timing comparison.

## 7. Replay and limits

- [Refinement and row substitution](../tools/oriented_coprime_basis.py).
- [Independent HNF and factor controls](../tools/check_oriented_coprime_basis.py).
- [Exact receipt](../evidence/oriented-coprime-basis.json).

```sh
python3 tools/check_oriented_coprime_basis.py \
  --output evidence/oriented-coprime-basis.json
```

The receipt checks 290 literal ideal identities, 29 mixed-orientation
batches, repeated prime powers, retained composite atoms and six invalid
inputs. It reconstructs the four power witnesses, verifies all old and new
relation identities, and extracts 1009 both from the transported system
and from the single explicit base relation.

The bit cost is polynomial in the supplied norms, roots, discriminant and
relation rows. The theorem supplies no bound on the number of residual
witnesses needed for a useful dependency or on the rejected powers before
them. Those costs remain part of a complete factoring algorithm.

## 8. More emitted rows on actual selector prefixes

The [prefix comparison](../tools/oriented_cofactor_prefixes.py) retains every
query from the existing round-robin prime-power scan, beginning at each
generator's norm threshold. It includes smooth norms and all repeated
keys. At fixed checkpoints it compares exact-key cancellation with the
oriented-atom refinement on precisely the same supplied rows.

Write the old relation matrix as `[A R]`, with the fixed-base columns
first. If W substitutes atom words for residual generators, the refined
matrix is `[A RW]`. Every old cancellation `cR=0` also satisfies `cRW=0`.
Thus the rational spaces of emitted base relations obey

$$\{cA:cR=0\}\subseteq\{cA:cRW=0\}. \tag{10}$$

To produce genuine integer relations, the implementation clears the
denominators of the entire row-combination certificate. Clearing only
the resulting row is insufficient. For example, eliminating the second
column of rows `(2,2),(0,1)` gives the rational combination
`row_2-row_1/2=(-1,0)`. The guaranteed integral relation is
`2 row_2-row_1=(-2,0)`. The checker includes this control, so an unchecked
half-relation cannot erase the order-two information being sought.

The emitted integer rows span the full rational cancellation space.
Equation (10), together with the saturation identity (8), proves that
the new base-only torsion image contains the old one. Every emitted row
retains its explicit integer combination of original rows.

Three matched finite prefixes give strict row and factor gains:

| Input N | Cutoff | Queries | Old base rank | Refined base rank | Refined factor |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1040399 | 5 | 16 | 1 | 2 | 1019 |
| 100160063 | 7 | 64 | 2 | 3 | 10007 |
| 1000036000099 | 31 | 320 | 0 | 3 | 1000003 |

The original `factor_with_cofactor_collisions` is run independently at
each exact query cap. It has the stated old rank and returns no factor
in all three cases. No earlier target involution or leading-norm gcd
yields a proper factor on these prefixes. An additional equal-outcome
control retains the old factor 1009 on input 1022117 after 16 queries.

On the 40-bit input, the original method retains 250 residual keys and
finds 70 dependent, trivial collision rows. Atom refinement yields 249
atoms and emits the three genuine base rows

$$
\begin{pmatrix}
78&-52&11\\
-37&23&52\\
-51&169&-53
\end{pmatrix}
$$

on primes `(5,7,19)`. Their exact extraction gives factor 1000003.
The 250 residual words contain 141 of support one, 107 of support two,
and two of support three; their maximum atom multiplicity is two.

The current prototype rebuilds the norm basis and uses rational row
elimination at each checkpoint. Its cumulative gcd-decision counts are
224, 3291 and 711903 for the three strict-gain controls. These are finite
query and row-supply improvements, with no wall-time speedup claim. Large
witness collections require more economical norm refinement and sparse
relation elimination.

The [independent prefix checker](../tools/check_oriented_cofactor_prefixes.py)
replays 416 literal powers with generic HNF composition, verifies every
emitted integer combination, checks its complete rational span, repeats
the exact extraction, and compares with the original selector. The
[retained receipt](../evidence/oriented-cofactor-prefixes.json) includes
all supplied witnesses and row-combination certificates.

```sh
python3 tools/check_oriented_cofactor_prefixes.py \
  --output evidence/oriented-cofactor-prefixes.json
```

This implementation and analysis used substantial OpenAI assistant support.
