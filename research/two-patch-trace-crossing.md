# A complete two-patch trace crossing and its finite rank control

The [bilinear trace kernels](trace-bilinear-conservation-kernels.md)
give individual pieces with low ranks in both cuts. They do not themselves
sum to the crossing tensor. This note supplies a related complete polynomial
crossing with at most `2t` nonzero coefficients. Its raw coefficient rank sums fail
in all five tested widths, so it is a synthesis control rather than an
improved factoring bound.

## 1. An exact endpoint without an exceptional label correction

Work over a coefficient field of characteristic two. Let
`L=F_(2^t)`, let `alpha` have degree `t`, and define the binary linear
coordinate isomorphism

```math
\tau(z)=\bigl(\operatorname{Tr}(\alpha^jz)\bigr)_{0\le j<t}.
```

For kernel variables `S,w` in `L`, set

```math
e=S+w,\qquad B_0=Se,\qquad B_1=(S+1)e.
```

The polynomial kernel and its four-leg tensor are

```math
G_t(z;S,w)=\prod_{j=0}^{t-1}
  (1+z\tau_j(B_0))(1+z\tau_j(B_1)),
\qquad
T_t(z)(x,y,u,v)=[x+y=u+v]G_t(z;x+y,x+u).
\tag{1}
```

At `z=1`, the product is one exactly when every trace coordinate of both
`B_0` and `B_1` is zero. Invertibility of `tau` makes this equivalent to
`B_0=B_1=0`. Their sum is `e`, so this holds exactly when `S=w`.
Consequently

```math
\boxed{T_t(1)(x,y,u,v)=[x=v][y=u].}
\tag{2}
```

This proof includes `S=0` and `S=1`; no additional singular-block piece is
needed. In fact the degree is exactly `2t-1`. For a coefficient of degree
`2t` to be nonzero, both trace vectors would have to be all ones. Then
`B_0=B_1`, which forces `e=0` and then `B_0=B_1=0`, a contradiction.
For the reverse inequality, let `R=tau^(-1)(1,...,1)` and let `beta`
be the inverse image of any coordinate unit vector. Take
`e=beta`, `S=R/beta`, and `w=S+beta`. Then `B_0=R` and
`B_1=R+beta`, so the two trace weights sum to `2t-1`. This also
covers `t=1`. Thus there are at most `2t` nonzero coefficient tensors,
and they sum to the crossing tensor. The constant coefficient has kernel one and has rank
`d=2^t` in each cut, since each of the `d` conservation blocks is all ones.

For fixed `S`, the two trace vectors in (1) are linear in `w`, with affine
constants. For fixed `w`, they are likewise affine linear in `S`: the map
`S -> S^2` is binary linear. This gives a bounded-degree trace interface,
but does not alone bound the ranks of its intermediate coefficients.

## 2. Exact coefficient ranks and scope of the screen

For each pair `S,w`, write

```math
m(S,w)=|\tau(B_0)|+|\tau(B_1)|.
```

The entry polynomial is `(1+z)^m`. Its coefficient of `z^k` is the parity
of `binom(m,k)`, or equivalently the predicate that every one-bit of `k`
is a one-bit of `m`. Thus the complete kernel coefficients are explicit.

For a fixed conservation block, both cuts are binary XOR-circulant
matrices. The engine computes their exact ranks and sums over the block
labels. It caches identical kernels, uses exact unit and linear-term
shortcuts in the square-zero group algebra, and otherwise performs packed
binary elimination. The independent checker constructs the full four-leg
coefficient tensors at `t=2,3,4`, multiplies every entry's polynomial
factors directly, and compares both cut ranks.

Let `S_rank=sum_k sqrt(r_k c_k)`, and define the displayed finite ratios
by `rho_sum=log2(S_rank)/(2t)` and
`rho_max=log2(max_k r_k c_k)/(4t)`.

| `t` | Field polynomial (integer) | Integer lower bound on `S_rank` | `rho_sum` | `rho_max` |
|---:|---:|---:|---:|---:|
| 2 | 7 | 31 | 1.242181 | 0.875000 |
| 3 | 11 | 164 | 1.228846 | 0.920366 |
| 4 | 19 | 863 | 1.219735 | 0.944283 |
| 5 | 37 | 4070 | 1.199222 | 0.959799 |
| 6 | 67 | 19180 | 1.185634 | 0.969278 |

The integer lower bound is `sum_k floor(sqrt(r_k c_k))`, so its comparison
with `d^2` needs no floating-point estimate. It exceeds `d^2` at every
tested width. The complete rank profiles are in the evidence file.

These are raw coefficient tests for five specified fields and coordinate
bases. They do not prove failure of all recombinations, other trace bases,
or the all-size family. Conversely, the endpoint identity alone does not
establish an asymptotic contraction or a better factoring score.

## Reproduction

```sh
python3 tools/check_two_patch_trace.py --output /tmp/two-patch-trace-verification.json
python3 tools/two_patch_trace_probe.py --lengths 2 3 4 5 6 --output /tmp/two-patch-trace-profiles.json
```

The checker regenerates all five profiles, compares 42 full coefficient
cuts, and verifies 4,672 endpoint entries within the conservation support.
It uses only the Python standard library and the existing
`check_trace_bilinear_kernel.py` field implementation. The evidence records
source hashes and exact integer ranks.
