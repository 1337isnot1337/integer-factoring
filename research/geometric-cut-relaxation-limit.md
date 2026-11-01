# The geometric compression bounds still admit the dyadic scalar profile

## Result and scope

The abstract profile from
[weighted-cut-relaxation-limit.md](weighted-cut-relaxation-limit.md)
satisfies every subset constraint supplied by the integer and
half-integer geometric split certificates. Its cost remains
`2d(log2(d)+1)` along powers of two. Thus these strengthened scalar
constraints still do not give a universal exclusion of costs below
`d^sqrt(2)`.

This is a statement about a list of rank pairs. It does not construct
tensors with those ranks, control cancellation in partial sums, or
improve a factoring exponent.

For `d=2^h`, `h>=1`, retain both pairs

```math
(d2^j,d/2^j),\quad(d/2^j,d2^j),\qquad 0\le j\le h.
\tag{1}
```

The two copies of `(d,d)` are distinct. Their row and column totals
are both `2d^2+d-1`; their geometric-rank cost is `2d(h+1)`.

## A uniform upper bound on all geometric split certificates

For `m>=2`, the geometric certificates give a complement lower bound
`s q` whenever either of the following holds:

* `a,b` are positive integers with `a+b=m`, or
* `a,b` are half-integers at least `3/2`, with `a+b=m`;

and the positive integers `s,q` satisfy

```math
s+a q\le d,\qquad b s+q\le d.
\tag{2}
```

Integer rounding in the three-map block construction is exactly
captured by (2): for an integer right side, `a q<=d-s` is equivalent
to `ceil(a q)<=d-s`, and similarly for `b s`.

Multiply the first inequality in (2) by `s`, the second by `q`,
and add. For `m>=3`, this gives

```math
(m-2)sq\le d(s+q)-(s+q)^2\le d^2/4,
\qquad sq\le\frac{d^2}{4(m-2)}.
\tag{3}
```

Let `F_geo(d,m)` be the maximum certified product, allowing both
split families. A positive pair forces `a,b<=d-1`, since `s,q>=1`.
Consequently `F_geo(d,m)=0` for `m>=2d-1`.
The same upper bound (3) applies to the zero-padding envelope
`max_{k>=m} F_geo(d,k)`: its right side decreases with `k`.

For exact crossings, include the separate rank-one and rank-two
complement bound `d(d-1)/2`. The following integer function therefore
dominates all the scalar complement lower bounds under discussion:

```math
B_d(m)=
\begin{cases}
d^2,&m=0,\\
d(d-1)/2,&m=1,2,\\
\left\lfloor d^2/[4(m-2)]\right\rfloor,&3\le m\le2d-2,\\
0,&m\ge2d-1.
\end{cases}
\tag{4}
```

The function `B_d` is an upper envelope used to check this abstract
profile. It is **not** claimed to be a valid complement lower bound
for arbitrary tensors.

## Every subset passes the stronger envelope

Let an index subset `J` of (1) have summed column rank `m`.
We prove

```math
\sum_{i\notin J}r_i\ge B_d(m).
\tag{5}
```

For `m=0,1,2`, the proof in the original relaxation note applies:
the empty subset leaves at least `d^2`, while a nonempty eligible
subset has one member and leaves at least `d^2+d-1`.

Suppose `3<=m<d`. Choose the smallest power of two `c>m`.
The list contains `(d^2/c,c)`, which cannot belong to `J`.
For `m=3`, it has `c=4` and row rank `d^2/4`, meeting (4).
For `m>=4`, `c<=2m` and

```math
\sum_{i\notin J}r_i\ge d^2/c
\ge d^2/(2m)\ge d^2/[4(m-2)].
\tag{6}
```

For `d<=m<=2d-2` and `d>=4`, the backward pair `(d/2,2d)`
cannot belong to `J`. It gives

```math
\sum_{i\notin J}r_i\ge d/2
\ge d^2/[4(d-2)]\ge d^2/[4(m-2)].
\tag{7}
```

The middle inequality holds precisely for `d>=4`.
At `d=2`, there is no unhandled interval: the cases `m<=2`
were covered already, and `m>=3=2d-1` have zero envelope.
Finally, for `m>=2d-1`, (5) is immediate. These cases exhaust all
subsets. Swapping the paired orientations proves the same statement
for the opposite cut.

For every integer `s>=0`, choosing the minimizing subset term by
term in (5) yields

```math
\sum_i\min(r_i,s^2c_i)
\ge \min_{0\le m\le2d-1}\bigl(B_d(m)+s^2m\bigr).
\tag{8}
```

The right side dominates the corresponding budget obtained from
`F_geo` or its zero-padding envelope. Values beyond `2d-1` cannot
decrease the minimum, because `B_d` is zero there and `s^2m`
is nondecreasing. The other cut follows by symmetry.

## Exact controls

Run from the repository root:

```sh
python tools/check_geometric_cut_relaxation_limit.py \
  --output evidence/geometric-cut-relaxation-limit.json
```

The checker keeps the duplicate center and both endpoint pieces.
For each power-of-two dimension through `4096`, it runs a descending
0/1 knapsack over exact summed cut budgets through `2d-1`, in both
orientations. It verifies (5) for every reachable budget, rather than
sampling subsets. It also checks (8) for every integer `0<=s<=d`.
For dimensions through `64`, it separately enumerates the complete
integer and half-integer split family and checks the envelope (4),
including its vanishing endpoint and its zero-padding maximum.

The exact `d=4096` comparison is inherited unchanged:
`106496^985<4096^1393` and `1393^2-2*985^2=-1` imply
`106496<4096^sqrt(2)`. No tensor realizability follows from these
arithmetic checks.

The construction, proof development, and replay used substantial OpenAI assistant support.
