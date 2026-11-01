# Regrouping a d-axis cube as a two-axis crossing

## Result and scope

The [multi-axis transfer](multiaxis-cube-transfer.md) turns a `d`-axis
crossing family into a factoring recursion with exponent
`alpha = log_d(1/rho)`. This note observes that the same local family can
also be read at growth base two in every dimension `d >= 3`, with the same
pressure `rho` and with `lambda' = 2^(d-2) lambda`-bit superlegs.

> **Conditional regrouping lemma.**  Suppose a uniform `d`-axis family,
> `d >= 3`, has polynomially many pieces, exact factorizations of the
> same pieces at all `d` coordinate cuts, separate widths at most
> `2^(d-1) lambda`, and total width
> `r_i^1 + ... + r_i^d <= d 2^(d-1) rho lambda + O(log(lambda+2))`.
> Pair the `2^d` legs by the `d-2` hidden coordinates and use the other
> two coordinates as parent axes. Then the four superlegs form an exact
> `E_4` crossing with `lambda' = 2^(d-2) lambda`, the parent-pair total
> is at most `4 rho lambda' B + O(B log(lambda'+2))`, and each separate
> parent width is at most `2 lambda'`. Hence the relaxed two-axis
> transfer applies with `rho' = rho` and gives `alpha = log_2(1/rho)`.

The per-cube version in `ddb303e`/`00e00c9` remains rejected: choosing a
different cube cut after each summand is selected is not a parent-plane
factorization, and its uniform fixed-cut width claim conflicts with rank
subadditivity. The present lemma chooses the two parent axes **globally
for the complete piece assignment**, exactly as the established multi-axis
transfer chooses one global axis. It does not relabel individual
summands.

For `d=3` this is the global-pairing repair of worker 2's check-in. The
new content here is the dimension-independent arithmetic and the `d>=4`
target. No `d>=4` family is constructed, and no exponent above `1/2` is
established by this note alone. The [three-axis family](cube-homogeneous-gaussian.md)
and its [pressure certificate](cube-optimized-rational-type.md) are
inherited unchanged.

## 1. The relaxed two-axis interface

Use the exact two-axis geometry of [the multi-axis transfer](multiaxis-cube-transfer.md)
with `d=2`. The parent is a clique on `[g]^2`; the `B = binom(g,2)^2`
squares carry the four-label endpoint

```math
E_4(x,y,u,v)=[x=v][y=u].                         \tag{1}
```

The two coordinate cuts are `{X,Y}|{U,V}` and `{X,U}|{Y,V}`. The cited
theorem assumes, for every piece, exact factorizations at both cuts with
`r_i^1+r_i^2 <= 4 rho' lambda'` and separate caps `r_i^a <= 2 lambda'`.
Its proof uses the pair bound only after summing over the `B` cubes, and
uses the separate cap for the maximum child-edge weight. We therefore
state the weaker hypothesis actually needed.

> **Relaxed two-axis lemma.**  Suppose `E_4 = sum_i T_i` has polynomially
> many pieces and each *same* `T_i` factors exactly at both cuts with
> widths satisfying
>
> ```math
> 0 <= r_i^a <= 2 lambda',
> \qquad
> sum_(cubes) (r_(i(cube))^1 + r_(i(cube))^2)
>   <= 4 rho' lambda' B + O(B log(lambda'+2)).
>                                                               \tag{2}
> ```
>
> Then the two-axis grid step and the recursion of the cited theorem give
> the same contraction `U/V <= rho' + 1/g + O(log(lambda'+2)/lambda')`,
> hence the uniform exponent `alpha = log_2(1/rho')`.

**Proof.**  For a fixed assignment, choose the coordinate axis `a` with
the smaller total `sum_cubes r_(i(cube))^a`. The two totals sum to the
left side of (2), so the chosen total is at most `2 rho' lambda' B +
O(B log(lambda'+2))`. Every cube contributes its width `r^a` to exactly
one child edge, so the oriented volume satisfies

```math
U = 2 sum_(u<v) w_(uv)
  <= 4 rho' lambda' B + O(B log(lambda'+2)) + (other-edge term).
                                                               \tag{3}
```

The separate cap gives `w_(uv) <= V/g^2` exactly as in the cited proof.
The equitable grouping lemma uses only `U` and that maximum edge weight.
The rest of the recursion is unchanged. This proves the lemma. `square`

The audit's fixed-cut rank contradiction applies to a decomposition in
which **every** piece at one fixed cut has width at most
`2 rho' lambda' + O(log lambda')`. The relaxed hypothesis does not assert
that. Individual pieces may have width up to `2 lambda'`; only their
total over the assignment is bounded by (2). Rank subadditivity therefore
forces at least one large piece, as it must, without contradicting (2).

## 2. Global bundling in dimension d

Let the `d`-axis family act on the `2^d` legs indexed by
`v in {0,1}^d`, with endpoint

```math
E_(2^d)((x_v)_v) = prod_({v,bar v}) [x_v = x_(bar v)],
\qquad bar v = v xor (1,...,1).                   \tag{4}
```

Fix two distinct parent coordinates `a,b` and let `H` be the other `d-2`
coordinates. A superleg is indexed by the pair `(v_a,v_b)` and bundles
the `2^(d-2)` legs with those parent values. Each superleg carries

```math
lambda' = 2^(d-2) lambda                                \tag{5}
```

bits. The antipodal map flips all coordinates, so it sends superleg
`(0,0)` to `(1,1)` and `(0,1)` to `(1,0)`, and reverses the hidden
coordinate order. Reading the two antipodal superlegs in the opposite
internal order, as in the `d=3` pairing checker, gives

```math
E_(2^d) = [X=V][Y=U].                                   \tag{6}
```

The parent cuts `a` and `b` partition whole superlegs. Every hidden cut
splits all four superlegs into their internal halves. The hidden cuts are
not used as parent cuts. This is exactly the incidence that the
[per-cube audit](two-axis-cube-regrouping.md) requires; the difference is
that here the parent pair is chosen once for the whole assignment.

Each original factorization at a parent cut descends to a factorization
of the same bundled piece at the corresponding superleg cut, with the
same channel width. The pieces, their polynomial count, and their
uniform `2^(O(lambda)) = 2^(O(lambda'))` evaluation are unchanged.

### Width arithmetic

Write `S_i = r_i^1 + ... + r_i^d`. The `d`-axis interface gives

```math
S_(i(cube)) <= d 2^(d-1) rho lambda + O(log(lambda+2)).  \tag{7}
```

Choose the two parent axes with the smallest total widths. The remaining
`d-2` axes carry at least `(d-2)/d` of the total, so the parent pair
carries at most `2/d`:

```math
sum_(cubes) (r^a + r^b)
  <= (2/d) sum_(cubes) S_(i(cube))
  <= 2^d rho lambda B + O(B log(lambda+2))
  =  4 rho lambda' B + O(B log(lambda'+2)).               \tag{8}
```

The separate cap is `r_i^a <= 2^(d-1) lambda = 2 lambda'`. Thus (2)
holds with `rho' = rho`, and the relaxed lemma gives the exponent
`alpha = log_2(1/rho)`.

The parent axes are selected globally after the complete piece assignment
is known. Different axes may be selected for different assignments, as
the cited multi-axis transfer already permits. No summand is reindexed
independently of the others.

## 3. The three-axis instance and the four-axis target

For `d=3`, `lambda'=2 lambda`. The existing homogeneous Gaussian cube has
total three-cut width at most `12 rho lambda + O(log lambda)`, so (8) is
exactly `4 rho lambda' B + O(B log(lambda'+2))`. The [certificate](cube-optimized-rational-type.md)
gives

```math
rho < 312359221/500000000 = 0.624718442,
```

hence the conditional exponent

```math
alpha > log_2(500000000/312359221) = 0.6786...
```

This is the same conclusion as worker 2's global-pairing repair. The
existing one-parameter rational type has been numerically optimized to
about `0.624718442` (the exact certificate parameters already sit at that
minimum), so the registered `0.72` target cannot come from re-optimizing
this `d=3` family. Reaching `alpha > 0.72` needs `rho < 2^(-0.72) =
0.6070...`.

The [multi-axis profile floors](multiaxis-gaussian-profile-floor.md) give
`rho >= 21/40 = 0.525` for the common-scale four-axis profile criterion,
and the [four-axis shared-minor argument](four-axis-gaussian-profile-floor.md)
gives the weaker floor `1/2`. Neither excludes a realized `d=4` family
with pressure below `0.607`. Constructing such a family, with exact
four-cut factorizations and a certified profile, is the concrete next
task; the regrouping above then supplies the growth-base-two transfer at
the same pressure.

## 4. Checks

[`tools/check_two_axis_all_dimensions.py`](../tools/check_two_axis_all_dimensions.py)
checks

1. the endpoint identity (6) for every one of the `2^8` assignments at
   `d=3` and every one of the `2^16` assignments at `d=4`, for all parent
   pairs (768 and 393,216 cases);
2. the cut incidence: parent cuts split no superleg, hidden cuts split
   all four;
3. the exact width arithmetic (8) for `d=3,...,9`;
4. the inherited three-axis certificate, its negative control, and the
   exact one-particle profiles `(2,4,8,10)`;
5. a bounded `400 by 400` numerical grid for the displayed one-parameter
   pressure formula, reported as a diagnostic rather than a proof of
   global optimality.

The checker does not re-prove the Gaussian tensor factorization or the
two-axis transfer, and it constructs no `d>=4` family.

## 5. What is not claimed

- The current `main` still records `alpha=1/2`. This note is a conditional
  interface result and a `d>=4` research target, not an integrated score
  improvement.
- The rejected per-cube one-cut regrouping is not revived. Its hidden-cut
  choice and per-summand relabeling remain invalid.
- The `d=3` conditional conclusion depends on the relaxed two-axis lemma
  and on worker 2's global-pairing correction being accepted. No `d>=4`
  family, profile certificate, or complete machine cost is supplied here.
- The numerical pressure grid is not a global optimum certificate.
