# Closing the small characteristic-two pairing seeds

## Result and scope

The [all-field extension](all-field-pairing-barrier.md) completes the
three-pairing obstruction in every other characteristic. The stronger
characteristic-two constants proved here remain unchanged.

This strengthens the [pairing-type barrier](pairing-word-barrier.md).
For every fixed seed dimension `d>=2`, the full coefficient span of

```math
(rR+hH+cC)^{\otimes t},\qquad H=E+R+C,
```

cannot give a leading contraction below `1/sqrt(2)` on the complete word
domain `[d]^t`. Arbitrary collective recombinations, any piece count, all
characteristic-two extension fields, and crossing-preserving coordinate
symmetries are allowed. In addition to the earlier bound for `d>=5`, the
three remaining dimensions satisfy

| Seed | Forced leading `rho` | Approximation |
| --- | --- | ---: |
| 2 | `3/4` | 0.75 |
| 3 | `5/(4 log2 3)` | 0.788662191964... |
| 4 | `(2+log2 14)/8` | 0.725919365257... |

More precisely, every crossing decomposition has one piece with **both**
ranks at least `2^(b_d t-O_d(log^2(t+2)))`, where the three numerators
`b_d` are respectively `3/2`, `5/2`, and `1+(log2 14)/2`.
These are construction obstructions, not lower bounds on factoring and
not improved factoring exponents. New local tensors, arbitrary binary
restrictions, and unrestricted nonsymmetric tensor words remain outside
the assertion. The ordinary proof below has substantial OpenAI assistant
support and awaits independent review; no new Lean formalization is used.

Combining the estimates gives a universal lower contraction
`rho_*=(2+log2 14)/8` across all seed dimensions. This also holds if the
seed dimension grows: the preceding `d>=5` bound has an explicit uniform
loss `log2(t+1)`, while the three small dimensions have the losses above.
The old bound exceeds `rho_*` at `d=5` by exact logarithm comparison and
at `d>=6` because it exceeds `3/4`. Within the existing growth-base-two
contraction analysis, this pairing-type route has saving exponent at most
`log2(1/rho_*)=0.46211879...`; it cannot even match the half endpoint by
that analysis. This does not constrain a different recursion interface.

## 1. A nonzero Boolean coefficient forces rank

Let `A_m=F[N_1,...,N_m]/(N_1^2,...,N_m^2)`. Its monomial basis gives the
regular matrix of multiplication by any `f in A_m`. If the coefficient
`f_S` of a square-free monomial on a set `S` is nonzero, then

```math
rank(m_f)>=2^{m-|S|}.                                      (1)
```

Restrict the input basis to monomials `N_T` with `T` disjoint from `S`,
and the output basis to `N_(S union U)` with `U` disjoint from `S`.
An entry vanishes unless `T subset U`; all diagonal entries are `f_S`.
Order subsets by size. This is a triangular square submatrix of order
`2^(m-|S|)` with nonzero diagonal. Higher-degree terms cannot cancel it.

The local constant blocks in the preceding note are polynomial matrices
in independent square-zero `N_i`. Equation (1) therefore applies to each
tensor-product block. The key improvement is to **sum all compatible
blocks**, rather than retain only one block type.

## 2. Even dimensions

Put `L=binom(d,2)`. The local row blocks are one

```math
K=cI+(r+h)N,
```

`L` copies of `J=cI+(h+c)N`, and `d-2` scalar blocks `c`.
The column blocks exchange `r,c`. Write `alpha=r+h`, `beta=h+c`.
For a type with `a` special blocks, `j` off-diagonal blocks, and
`t-a-j` scalars, a coefficient with `k` special nilpotents and `l`
off-diagonal nilpotents has entry polynomial

```math
c^{t-k-l} alpha^k beta^l.                                 (2)
```

It is independent of the enclosing type as long as `a>=k,j>=l`.
If its linear functional is nonzero on a piece, (1) and direct-sum
additivity give the rank lower bound

```math
2^{-k-l} sum_{a>=k,j>=l,a+j<=t}
  multinomial(t;a,j,t-a-j) 2^a (2L)^j (d-2)^{t-a-j}.        (3)
```

Multiplicities here are integer counts of different blocks, never field
scalars. Formula (3) includes every compatible type of the same piece.

Set `t=2s`, where `s` is a power of two, and take the common polynomial
`p=(rc)^s`. Since `r=c+alpha+beta`, Frobenius gives

```math
p=c^{2s}+c^s alpha^s+c^s beta^s.                           (4)
```

These are three entry polynomials from (2). The column cut has the same
identity with `r,c` exchanged. Thus coefficient pairing with `p` defines
one actual matrix-entry functional in both cuts, and its value on the
crossing is `p(1,1,1)=1`.

For a piece detected by it, at least one term in (4) is nonzero. The first
forces full rank `d^(4s)`. For the second, (3) becomes

```math
2^{-s} sum_{a=s}^{2s} binom(2s,a) 2^a (d^2-2)^{2s-a}
 >= binom(2s,s)(d^2-2)^s.                                 (5)
```

For the third, the binomial success weight is `2L>=2`, rather than 2,
with the same total weight `d^2`. Its tail is no smaller: view the sum as
the number of words over an alphabet of size `d^2` with at least `s`
letters in a specified subset; enlarging that subset only adds words.
The first term's full rank also dominates (5). Consequently that same
piece has both ranks at least

```math
B_d(2s)=binom(2s,s)(d^2-2)^s
 >= [2 sqrt(d^2-2)]^{2s}/(2s+1).                           (6)
```

For `d=2`, the leading rank rate is `3/2`; for `d=4`, it is
`1+(log2 14)/2`. Dividing by `2 log2 d` proves the stated contraction
bounds. The seed-four comparison needs no numerical assumption:
`14^5>2^19` implies `rho_4>29/40`, and `2(29/40)^2>1`.

## 3. Seed three

For odd `d`, the blocks are a scalar `D=r+h+c`, `d-1` scalars `c`,
and `L` copies of `J=cI+(h+c)N`. With `a` special scalars and `j`
off-diagonal blocks, the degree-`l` nilpotent entry polynomial is

```math
D^a c^{t-a-l}(h+c)^l.                                    (7)
```

Nonzero detection again forces rank in every compatible block, totaling

```math
binom(t,a) 2^{-l} sum_{j=l}^{t-a} binom(t-a,j)
   (d-1)^{t-a-j}(2L)^j.                                  (8)
```

For `t=2s` and `s` a power of two, `r=D+(h+c)` gives

```math
p=(rc)^s=D^s c^s+c^s(h+c)^s.                              (9)
```

The column identity is the same with `r,c` exchanged. The first term of
(9), with `a=s,l=0`, has bound `binom(2s,s)(d^2-1)^s`.
The second has `a=0,l=s`. Its binomial success weight is
`2L=d(d-1)>=d-1`. At least half of all weighted words have at least half
their letters in the success subset, by pairing counts at `j` and `2s-j`.
Its bound is therefore at least `(d^2-1)^(2s)/2^(s+1)`.

For `d=3`, both alternatives are at least

```math
2^{5s}/[2(2s+1)].                                        (10)
```

The leading rank rate is `5/2`, giving `rho_3=5/(4 log2 3)`.
In particular `3^3<2^5` proves this is greater than `3/4`.

## 4. All word lengths, without a subsequence loophole

For an even general length `t`, split the positions into its nonzero
binary blocks `m_i=2s_i`, where each `s_i` is a power of two.
The full global type span is contained in the tensor product of the type
spans on these classes. Take the tensor product of their common
functionals from (4) or (9). It takes value one on the crossing.

Expand this functional as a sum of products of the detected entry
functionals. A nonzero value on a piece selects at least one such product
in each cut. On every compatible product of local block types it detects
a nonzero Boolean coefficient of degree equal to the sum of the detected
degrees. Equation (1) gives the product of their Boolean rank bounds.
Summing over all compatible type choices factorizes, giving at least
`prod_i B_d(m_i)` in each cut. The row and column selections may differ;
both have this same lower bound for the **same detected piece**.

There are `O(log(t+2))` binary classes, and their logarithmic denominator
losses total `O(log^2(t+2))`. If `t` is odd, fix the four coordinates of
one position to a single common label. This preserves the crossing and
does not increase either rank. The remaining tensor belongs to the
length-`t-1` type span: the local scalar pencil on that fixed entry is
`r+h+c`. This proves the all-length statement.

Together with the preceding `d>=5` result, every fixed pairing seed is
excluded. The same argument also permits separate pairing type spans
within `b=o(t)` position classes: apply the binary split inside each class.
The loss is `O_d(b log^2(2+t/b))=o(t)` by concavity, including the at most
one discarded position per class. This extension concerns pairing types;
it does not cover arbitrary invariant local tensors.

## 5. Checks and limitations

```sh
python3 tools/check_pairing_small_seeds.py
python3 tools/check_pairing_small_seeds.py --output evidence/pairing-small-seeds.json
make verify
```

The checker replays the Boolean triangular minors, the exact compatible
block sums and binomial inequalities, common polynomial identities in
both cuts, and binary-class losses. Direct Kronecker matrices independently
check small type-space ranks against the block formulas. These finite
checks support their exact stated scope; the arbitrary-size conclusion
uses the written proof. There is still no construction attaining a
factoring exponent greater than `1/2`.

The same evaluator exhausts the full `F_2` type span through word length
four, for seeds 2 through 6. The dimension is `binom(t+2,2)`: the three
local tensors are independent (each matching has an entry where the other
two vanish), so their ordered tensor words are independent, and distinct
types sum disjoint sets of words. Sorting all coefficient masks by rank
product and checking exact span membership therefore proves the finite
cutoffs, rather than just checking individual witnesses. For example,
seed two has cutoffs `4,36,224,2116` and seed three has
`16,276,8100,216240` at lengths 1 through 4. The last seed-three normalized
width is about 0.69884, below the half threshold at that finite length;
this is consistent with an asymptotic obstruction with a logarithmic
loss. Add `--max-word-length 5` for the larger optional full enumeration.
