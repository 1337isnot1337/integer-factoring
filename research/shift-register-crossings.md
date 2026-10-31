# Locally joined shift-register crossings

## Result and scope

A maximal binary linear-feedback register, with a local join at zero,
gives a regular full-cycle crossing in the
[cyclic-algebra constructor](cyclic-algebra-crossings.md). The physical
letters may use any balanced partition of the register coordinates.
This supplies explicit full cycles outside the one-layer controlled
permutation family.

The completed search covers all maximal feedback rules and all balanced
coordinate partitions through ten register bits: 8,186 systems at four,
eight, sixteen, and thirty-two physical letters. Every system has an
exact two-jet certificate `rho >= 17/24 > 1/sqrt(2)` for the bare
two-place criterion. A separate selected twelve-bit search also excludes
138 systems, without excluding all twelve-bit partitions. These results
do not establish a universal shift-register obstruction. The strongest
factoring exponent remains `1/2`.

## 1. A complete power-of-two cycle from a linear register

Let `m=2h` and write a binary register state as an integer `s<2^m`.
Choose a binary tap vector with top bit one, and define the invertible
linear map

```math
L(s)=((2s)\bmod 2^m)+\operatorname{parity}(s\mathbin{\&}t).
                                                               \tag{1}
```

Suppose its nonzero states form one cycle of length `2^m-1`.
The enumeration checks this condition directly, rather than depending
on a table of primitive polynomials. The zero state is its other cycle.
Let `e=2^(m-1)`. The top-tap condition gives `L(e)=1`.
Interchange the images of `0` and `e`:

```math
U(0)=1,\quad U(e)=0,\quad U(s)=L(s)\quad(s\ne0,e).           \tag{2}
```

This splices the two cycles into one full `2^m`-cycle. Thus, as a
permutation matrix in characteristic two, `I+U` is regular nilpotent.
No conjecture about register periods is used.

Choose an `h`-element coordinate set `A` and its complement `B`. Encode
the `A` bits as the first physical letter and `B` bits as the second,
using increasing bit order inside each letter. Then `q=2^h` and
`n=q^2=2^m`. Apply equation (2) of the cyclic-algebra note with `g(X)=X`.
That note proves the complete crossing endpoint, the column orders
`(0,n,...,n)`, and the uniform construction interface.

Complementing the bit partition conjugates `U` by the swap of its two
physical factors. Its row reshuffling becomes a transpose, up to fixed
index permutations. Both row Smith orders and all row jet ranks are
unchanged. Requiring `0 in A` therefore covers every balanced partition
up to this proved symmetry.

## 2. Linear cut ranks from small binary blocks

For any invertible linear map `L` on `A direct-sum B`, write its binary
matrix in blocks `[[P,Q],[R,S]]`. The exact row-reshuffling rank of its
permutation matrix over characteristic two is

```math
\operatorname{rank}\Phi(L)=2^{\operatorname{rank}P+
                                      \operatorname{rank}S}.   \tag{3}
```

To see this, the row and column variables of the reshuffling obey
`x=Pv+Qy`, `u=Rv+Sy`. Define binary linear maps

```math
F(x,y)=(x+Qy,Sy),\qquad G(u,v)=(Pv,u+Rv).
```

Its matrix entries are `[F(x,y)=G(u,v)]`. Rows with different attainable
values of `F` have disjoint supports; rows with the same value agree.
The rank over any field is therefore the number of elements in
`im(F) intersection im(G)`. The combined map is surjective, since
`x,u` are free. Its intersection dimension is
`rank F+rank G-m=rank P+rank S`, proving (3).

The correction in (2), as an operator on states, is the rank-one matrix

```math
\Delta=(e_0+e_1)(e_0+e_{2^{m-1}})^{\mathsf T}.              \tag{4}
```

Each vector in (4) is a product vector for every coordinate partition:
it varies in only one bit. Reshuffling their outer product consequently
also has rank one. In particular, the joined reshuffling rank differs
from the value in (3) by at most one. Every case in the complete search
checks (3) and this perturbation bound against literal binary matrices.
Equation (3) also makes selected larger searches inexpensive to specify.

## 3. Complete two-jet accounting

The first two row coefficients at normalized infinity are
`A0=Phi(U)` and `A1=Phi(U+U^2)`. Let `r=rank A0` and let `k` count
the order-one pivots, recovered by the block-Toeplitz matrix in the
cyclic-algebra note. Its exact pressure inequality gives

```math
B(s)\ge nr+(n-1)k\quad(s\ge1).                            \tag{5}
```

For a quick exclusion, `B(s)>=nr` sometimes suffices. The checker computes
`k` only when needed, using the exact test `(nr)^6 < q^17`. Every
recorded complete system then satisfies

```math
\bigl(nr+(n-1)k\bigr)^6\ge q^{17},\qquad
\rho\ge17/24>1/\sqrt2.                                   \tag{6}
```

An omitted `k` in the saved record means the zero-order rank already
proves (6); it does not mean zero order-one pivots. This is an obstruction
for the displayed bare profile criterion. Arbitrary physical gauges,
word-type optimizations, other feedback functions, and coefficient
recombinations are outside its scope.

## 4. Completed finite schedules

| Register bits | Physical alphabet | Maximal feedback rules | Partitions modulo complement |
| --- | --- | --- | --- |
| 4 | 4 | 2 | 6 |
| 6 | 8 | 6 | 60 |
| 8 | 16 | 16 | 560 |
| 10 | 32 | 60 | 7,560 |

For ten bits the smallest reported lower bound in (5) is `21500`, at
tap mask `516` and `A={0,1,3,4,7}`; this case has `r=17`, `k=4`.
Some larger zero-order ranks use the weaker bound `nr`, as described
above. The complete finite claim rests on every exact comparison (6),
rather than on a numerical minimum.

The separate twelve-bit schedule enumerates all 144 maximal rules and
66,528 partitions to select the 138 cases where both exponents in (3),
for `L` and `L^2`, are at most four. Literal joined matrices give
`r=k=17` in the minimum case, hence floor `139247`, still above `64^(17/6)`.
Only those 138 cases receive a two-jet exclusion in this schedule.
The other 66,390 twelve-bit partitions are not claimed excluded here.

## 5. Reproduction and status

```sh
make verify-shift-register
python3 tools/check_shift_register_crossings.py --max-bits 10 \
  --focused-twelve --output /tmp/shift-register-crossings.json
```

Both commands use Python's standard library. The small controls construct
complete polynomial tensors at four and eight letters, check every
endpoint entry, compare literal first jets with the constructor, and
recover complete row profiles by independent block-Toeplitz ranks.
The larger catalogues enumerate the exact feedback and partition spaces
listed above. Every joined cycle is tested, both linear block-rank
identities are checked where claimed, and the final pressure comparisons
use only integers.

Saved output: [`shift-register-crossings.json`](../evidence/shift-register-crossings.json).
The constructor is a proposed seed source; none of these results improves
the factoring score. The regular-nilpotent transfer retains its original
attribution. This construction, rank proof, and search used substantial
OpenAI Codex assistance, without new Lean formalization.
