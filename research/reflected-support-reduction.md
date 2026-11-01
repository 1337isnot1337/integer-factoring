# Far-reflected labels reduce to coordinate matching

## Result and scope

The high-energy reflected duals from
[fermionic-piece-polynomial.md](fermionic-piece-polynomial.md) have an
exact obstruction when the mode box is larger than every label part.
In the fixed-particle case they reproduce the elementary binary
coefficient construction. Allowing varying particle counts does not
remove the obstruction: a low-rank charge projection recovers the same
support condition.

Let `S` be any set of `d` distinct-odd partitions of one positive weight
`h`. Choose a power of two `M` satisfying

```math
M>\max\{k:k\text{ occurs in a label in }S\}.
\tag{1}
```

Use the unit dual `a_x=t^{x^*}`, where
`x^*={2M-k:k in x}`, and `b_v=t^v`. Let `A_s` be the complete
generalized-locality pieces. Put `L=max_(x in S) length(x)`.
For the row and column ranks `r_s,c_s`, respectively, we prove

```math
\max_s r_sc_s
\ \ge\ \frac{d^3}{(L+1)^2(2L+1)^2}.
\tag{2}
```

Thus if `d=2^lambda` and `log(L+1)=o(lambda)`, this particular reflected
family requires `rho>=3/4`. It cannot supply a rate below `1/sqrt(2)`.
The conclusion includes arbitrary label restrictions and varying particle
counts under (1). It does not cover smaller boxes where reflected modes
overlap the label modes. A finite negative control below shows that the
support conclusion can fail in that regime.

This is a restriction on a concrete construction, not on all crossing
families or factoring algorithms. It uses the already-proved
[coordinate-support theorem](coordinate-support-rank.md).

## 1. A matching with balanced particle counts has no mode changes

Use the four-block Pfaffian from
[fermionic-piece-polynomial.md, Section 3](fermionic-piece-polynomial.md).
Name its vertex groups `U,B,A,Y`; an `A` vertex associated with the
original label mode `x` has actual operator mode `2M-x`.

Condition (1) makes every `U--A` coefficient zero: the reflected
operator mode exceeds the output mode. Suppose an entry satisfies

```math
|x|+|u|=|v|+|y|,
\tag{3}
```

where the bars denote particle counts. Then the two vertex classes
`A union U` and `B union Y` have the same size. Since there is no edge
inside `A union U`, every perfect matching must pair these two classes
bijectively. In particular, it contains no `B--Y` edge.

All remaining nonzero edges have the following mode inequalities:

| Edge | Inequality in the original label modes |
|---|---|
| `U_u--B_b` | `u>=b` |
| `A_(2M-x)--Y_y` | `x>=y` |
| `A_(2M-x)--B_b` | `x>=b` |
| `U_u--Y_y` | `u=y` |

The first is the triangular creation coefficient. For the next two,
the binary-complement identities give

```math
D(2M-x,y)=[y\subseteq x],\qquad
D(b,2M-x)=[(b-1)\subseteq(x-1)].
\tag{4}
```

Every edge therefore has original left mode at least original right
mode, with left class `A union U`. Both classes have total original
weight `2h`. Summing their nonnegative differences forces **every
matched mode to be equal**. All surviving scalar coefficients are one.

At each physical mode, the only possible pairings are consequently

```math
x=v,\ y=u,
\qquad\text{or}\qquad
x=y,\ v=u.
\tag{5}
```

Encode a partition by its incidence vector in the odd modes. Equation
(5) is exactly the coordinate support with orientation mask zero.
It holds for every coefficient of the Pfaffian, even after cancellations.

## 2. Fixed particle count gives an exact elementary polynomial

Assume every label has `ell` particles. Then (3) holds for every entry.
For one binary incidence coordinate, write

```math
E=[x=v][y=u],\qquad R=[x=y][v=u],\qquad C=E+R.
```

If `P(w)=prod_k(R_k+w C_k)`, restricted to the label set, the entire
reflected family satisfies the entrywise Laurent-polynomial identity

```math
\sum_s A_s z^s=z^{-2M\ell}P(z^M).
\tag{6}
```

Here the product ranges over the finitely many odd-mode coordinates
used by the labels. No change of basis or approximation is involved.

To check the powers in (6), put `T=z^(2M)`. For an entry with support
(5), let `a` count modes occupied only by `A,B`, let `b` count modes
occupied only by `A,Y`, and let `c` count modes occupied by all four
groups. Then `a+b+c=ell`. Fixed particle counts also give `a` modes
occupied only by `U,Y` and `b` occupied only by `U,B`.

The exact-mode matchings from Section 1 contribute

```math
T^a(1+T)^b\bigl(T+(1+T)\bigr)^c=T^a(1+T)^b.
\tag{7}
```

For example, the `z` exponents of the `A--B` and `U--Y` edges sum to
`2M` times the number of `A--B` edges: their original-mode contributions
cancel because all four external weights are `h`. The analogous
`1+z` exponents also cancel, leaving `2M` times the number of
`A--Y` edges. This gives (7). The generalized locality offset is
`-2M ell`.

The elementary product is `w^(2a)(1+w)^(2b)`. Substitution `w=z^M`
and the characteristic-two identity `(1+w)^2=1+w^2` give (6).
In particular, only `ell+1` possible coefficient positions remain.
The larger physical modes supply no new tensors in this regime.

## 3. Varying particle counts: project the charge with a small rank cost

Define the entry mask

```math
F(x,y,u,v)=[|x|+|u|=|v|+|y|],\qquad B_s=A_s\mathbin{\circ}F,
\tag{8}
```

where `circ` is entrywise multiplication. Every crossing-diagonal entry
has `F=1`, so `sum_s B_s=E`. Section 1 proves that each `B_s` has
support (5).

The row flattening of `F` factors through the integer charge
`|x|-|y|=|v|-|u|`, with at most `2L+1` values. Its column flattening
factors through `|x|+|u|=|v|+|y|`, also with at most `2L+1` values.
The elementary rank inequality for an entrywise product therefore gives

```math
r(B_s)\le(2L+1)r(A_s),\qquad
c(B_s)\le(2L+1)c(A_s).
\tag{9}
```

This inequality follows directly by multiplying rank factorizations:
a sum of `r` outer products times a sum of `q` outer products is a sum
of at most `rq` outer products. It is valid over `F_2`.

There are at most `L+1` nonzero projected shifts. Indeed, for each
charge-balanced entry the exact-mode proof again gives a numerator
`T^a(1+T)^b`, with `a+b<=|x|`, and offset `T^(-|x|)`.
The cancellation of the original-mode exponents uses common weight,
not common particle count. Thus every projected entry is supported on

```math
s\in\{-2ML,-2M(L-1),\ldots,0\}.
\tag{10}
```

Apply the coordinate-support crossing theorem to the projected family:
some `s` has `r(B_s)c(B_s)>=d^3/(L+1)^2`. Equation (9) proves (2).

The argument also constrains later linear recombinations. If `N`
recombined tensors from this reflected span sum to the crossing, apply
the same mask to each. The resulting bound is

```math
\max_j r_jc_j\ge\frac{d^3}{N^2(2L+1)^2}.
\tag{11}
```

Consequently subexponentially many recombined pieces still cannot yield
`rho<3/4` when `log(L+1)=o(lambda)`.

## 4. Exact checks and the unexcluded regime

Run:

```sh
python3 tools/check_reflected_support_reduction.py \
  --output evidence/reflected-support-reduction.json
```

The checker compares (6) with the independent complete Wick evaluator,
checks the charge mask and both rank inequalities, and verifies the
projected crossing sum. It also checks two concrete exclusions:

- Without the charge mask, varying-particle families can have entries
  outside coordinate support even when (1) holds.
- With `M=8` and the four two-particle labels of weight 16, a central
  entry lies outside that support. Here the largest label part is 15,
  so (1) fails.

The second example identifies the remaining reflected-label route:
boxes with overlapping original and reflected modes permit the third
coordinate pairing and genuine changes of modes. The present theorem
does not bound that route.

The source operators are Anthropic's Appendix B construction at
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`. The generalized reflected
identity, this reduction, and the checks used substantial OpenAI
assistant support. No new Lean formalization is claimed.
