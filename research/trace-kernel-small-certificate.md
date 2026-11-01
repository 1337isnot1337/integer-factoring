# A smaller binary trace-kernel counterexample with exact support

The [trace-kernel construction](trace-bilinear-conservation-kernels.md)
has an explicit binary instance on `d=2^19=524288` labels with

```math
r=c=67108736,\qquad \operatorname{rank}K=94183.
```

Consequently

```math
rc=4503582447517696
<4650663930232832=d(\operatorname{rank}K)^2.
\tag{1}
```

This reduces the dimension of the companion note's symbolic counterexample
from `2^33` to `2^19`. The large tensor is not materialized. Its ordinary
kernel rank is certified by a complete exact coefficient calculation and
the companion note's Vandermonde rank identity. Within that specific
trace-product family, nineteen is the smallest possible extension degree
for a counterexample to this inequality. This is a statement about a
proposed rank obstruction, not a new factoring score or a crossing synthesis.

## 1. The fully specified instance

Define

```math
L=\mathbf F_2[X]/(X^{19}+X^5+X^2+X+1),\qquad \alpha=X\bmod p(X),
```

and put

```math
f(z)=\prod_{j=0}^{6}\operatorname{Tr}_{L/\mathbf F_2}(\alpha^j z),
\qquad K(S,w)=f(Sw).
\tag{2}
```

The defining polynomial is irreducible. Exact division by all 511 monic
constant-one polynomials of degrees one through nine certifies this:
any proper factorization of a degree-nineteen polynomial would have a
factor in that finite list.

The seven trace forms are independent. The companion convolution proof
therefore gives both cut ranks exactly

```math
(2^{19}-1)2^7=67108736.
\tag{3}
```

## 2. Complete multiplicative coefficient support

Let `n=2^19-1`. Expand the successive products as functions on `L^*`,
reducing exponents modulo `n`. With `C_0(0)=1` and every other initial
coefficient zero, the exact recurrence is

```math
C_{j+1}(e)=\sum_{i=0}^{18}
 C_j(e-2^i\bmod n)\,\alpha^{j2^i},
\qquad 0\le j<7.
\tag{4}
```

This needs only two arrays of `n` field elements and field arithmetic.
The implementation evaluates every slot, with no sampling of support.
For this instance it obtains

```math
C_7(e)\ne0\quad\Longleftrightarrow\quad
1\le\operatorname{popcount}(e)\le7.
\tag{5}
```

The complete support counts after each multiplication are:

| Trace factors | Nonzero coefficients |
| ---: | ---: |
| 1 | 19 |
| 2 | 190 |
| 3 | 1,159 |
| 4 | 5,035 |
| 5 | 16,663 |
| 6 | 43,795 |
| 7 | 94,183 |

The evaluation matrix on `L^*` is invertible, and the kernel is
`V diag(C_7) V^T`. Thus (5) proves the exact rank

```math
\operatorname{rank}K=\sum_{j=1}^7\binom{19}{j}=94183.
\tag{6}
```

The computation is finite evidence for this fully specified instance.
It does not assert that every defining polynomial or every trace-product
length has full support in all the lower weight layers.

## 3. Minimal extension degree in the stated family

For every `1<=k<t`, cyclic binary carries in the companion proof give
the universal support bound

```math
\operatorname{rank}K\le U(t,k):=\sum_{j=1}^{k}\binom tj.
```

The two cut ranks are always `(2^t-1)2^k`. Therefore a necessary condition
for this family to violate `rc>=d rank(K)^2` is

```math
((2^t-1)2^k)^2<2^t U(t,k)^2.
\tag{7}
```

Exact integer comparisons exclude (7) for all 153 pairs with
`2<=t<=18` and `1<=k<t`. At `t=19`, the pair `k=7` satisfies it,
and (5) realizes the upper support bound. This proves minimality of
nineteen within the specified family. It says nothing about the smallest
counterexample among other conservation kernels.

## 4. Replay and independent checks

```sh
python3 tools/check_trace_kernel_support.py \
  --output /tmp/trace-kernel-support-t19.json
```

The driver compiles the C++17 coefficient recurrence and reconstructs all
524,287 coefficient slots. It checks the full weight support, all 94,183
nonzero Frobenius relations, the irreducibility certificate, and (1).
An independent Python implementation agrees on fifteen small rank profiles
and every coefficient in six small fields. Five full polynomial evaluations
are compared with the direct trace formula, including two positive values
obtained by solving the trace equations. The record includes source hashes
and a hash of the complete coefficient array.

These checks certify the small finite calculation alongside the separate
all-size rank proofs. A competitive crossing decomposition remains open.
