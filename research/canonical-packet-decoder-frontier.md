# Exact canonical decoder thresholds at small packet lengths

## Result and scope

When the binary row tensor is fixed and known, the following precisions
are necessary and sufficient, uniformly over all such row tensors, to
determine the original contraction from the complete odd-center local
polynomials of a canonical binary column tensor.

| Packet length | Original parity | Full integer count |
|---:|---:|---:|
| `2` | `3` bits | `4` bits |
| `4` | `12` bits | `13` bits |

These are exact finite separation theorems. They strengthen
[the visibility results](canonical-packet-upper-frontier.md): two distinct
canonical columns with the same recorded local data have the same requested
output. An exhaustive decoder exists by listing all binary columns. This
does not provide an efficient all-size decoder or a better factoring bound.

There is also an ordinary all-size lower construction. For `L=2^a>=4`,
two canonical binary columns and the same binary row tensor have original
counts zero and two, while their local data agree modulo

```math
2^{L(a+1)}.
```

Thus uniform full-count decoding in the fixed Cauchy coordinates needs
at least `L(a+1)+1` bits. This sharpens the parity lower witness by one
bit for the full count. It remains a local-data statement; input-dependent
coordinate changes or direct use of the tensors are separate possibilities.

The Cauchy packet is the one in
[OpenAI's quantitative contraction proof](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`. The
[alternating witness note](canonical-alternating-packet.md) supplies the
all-size local-content formula used below.

## 1. Differences make the reconstruction statement exact

Fix a binary row tensor `A`. If `B_1,B_2` are binary column tensors, their
difference `D=B_1-B_2` has entries in `{-1,0,1}`. Conversely every such
ternary tensor is a difference of two binary tensors. The packet
contraction is linear in the column tensor, so equal local data modulo
`2^q` are equivalent to vanishing local data for `(A,D)` modulo `2^q`.

At `L=4`, exact exhaustive checking establishes both implications:

```math
J(A,D)\equiv0\pmod{2^{12}}
 \quad\Longrightarrow\quad \sum_XA(X)D(X)\equiv0\pmod2,
\tag{1}
```

```math
J(A,D)\equiv0\pmod{2^{13}}
 \quad\Longrightarrow\quad \sum_XA(X)D(X)=0.
\tag{2}
```

Here `J` contains all 116 possible coefficients at the four odd centers.
The total finite domain has `2^16*3^16=6^16=2821109907456` pairs.
For (1), exactly `(6^16-2^16)/2=1410554920960` pairs have odd original
value. The zero row tensor need not be searched because its original
value is zero.

The checker uses three fixed modular linear combinations of the original
coefficient matrices as necessary filters. For each binary row tensor,
split the ternary column into two blocks of eight entries. Each half has
`3^8=6561` possibilities. Hashing their three projected sums enumerates
every pair that can pass all local coefficients. Every matching pair with
the relevant original value is then checked against all original matrices.
Duplicate hash keys retain every half-column. The projections change only
the runtime, so the exclusion is deterministic and exhaustive.

Positive controls at precision 11 in odd mode and precision 12 in nonzero
mode are returned and replayed in Python. At `L=2`, all `6^4=1296` pairs
are checked directly. The maxima for odd and nonzero original differences
are contents 2 and 3, respectively.

## 2. An all-size full-count collision

Let `X_*` be the alternating state with ones in all even positions.
Modify the alternating row tensor by ignoring its last, odd-positioned
bit. Its multilinear polynomial is

```math
\widetilde A(z)=
 \prod_{x\text{ even}}z_x
 \prod_{\substack{x\text{ odd}\\x<L-1}}(1-z_x).
\tag{3}
```

This is a binary row tensor of degree `L-1`, supported on the two states
`X_*` and `X_* xor e_(L-1)`. Let

```math
Q_L(T)=\widetilde A(W(T)X_*),\qquad
K_L=L(a+1)-1.
```

The omitted last-row factor has content zero at positive odd centers.
At negative odd centers it has content exactly one. Indeed, its constant
term is `1-w_t`, where

```math
w_t=(-1)^t\frac{\binom{2L-1}{L}}{\binom{2L-1}{L+t}}.
```

For negative odd `t`, both binomial coefficients are 3 modulo 4, so
`w_t=-1 modulo 4`. The higher local coefficients of every Cauchy entry
are even at that center. Thus the omitted factor has content one.
The binomial congruence follows from

```math
(1+z)^{2L}\equiv1+2z^L+z^{2L}\pmod4
```

by division by `1+z`: the coefficients below `L` alternate `1,-1`,
and the coefficient at `L` is 3 modulo 4.

The occupied row factors alone contribute at least `K_L`, as proved in
the alternating note. Dropping this one unoccupied factor therefore leaves

```math
\mathrm{cont}_2 Q_L(t+2Lz)\ge K_L
```

at every odd center. Equality remains at `t=L-1`.

Put `chi(X)=(-1)^(|X|)`. For any polynomial in the Boolean coordinates of
degree less than `L`, its signed sum against `chi` vanishes. Applying
this elementary finite-difference identity to (3) gives

```math
\sum_X\chi(X)\widetilde A(W(T)X)=0.
\tag{4}
```

Since `L/2` is even, `chi(X_*)=1`. Therefore

```math
D(X)=\chi(X)-2\mathbf1_{X=X_*}
```

has only the values `-1,1`. Split it as `D=B_+-B_-`, where both columns
are binary indicators. Equations (3)--(4) imply

```math
F_{A,B_+}(T)-F_{A,B_-}(T)=-2Q_L(T).
\tag{5}
```

The original counts are respectively zero and two: both row-support
states have `D=-1`. The complete local difference has minimum content
`K_L+1`, with equality at the positive endpoint. Thus the two columns
give exactly the same local data modulo `2^(K_L+1)` but different counts.

At `L=4`, this yields content 12 and proves the necessity of 13 bits in
the table. The ordinary alternating point witness against the zero column
proves the necessity of 12 bits for parity. The small exhaustive upper
certificates therefore meet both lower constructions.

## 3. What remains for an evaluator improvement

These finite results show that the entry-height restriction changes the
required precision. As a control, the checker supplies an explicit
canonical 22-bit column for the four-wire alternating point row. Its
original value is one and its minimum odd-center content is 22. The same
row's binary columns obey (1). The 22-bit example does not establish an
extra 22-bit loss beyond creation precision.

The current exhaustive decoder may list `2^(2^L)` binary columns. An
all-size separation bound alone would not remove that cost. The useful
next obligation is a uniform decoder that exploits the canonical entry
constraints with the allowed tensor-entry access cost, or a reconstruction
identity that proves the same bounds without enumerating the entries.

## 4. Reproduction

```sh
python3 tools/check_canonical_packet_decoder.py \
  --output /tmp/canonical-packet-decoder.json
```

The checker uses the Python standard library and an optional C++17
compiler. It records both exhaustive exclusions, both positive controls,
the explicit full-count collision, and source hashes in
[canonical-packet-decoder.json](../evidence/canonical-packet-decoder.json).
No floating-point optimization or numerical rank is used.
