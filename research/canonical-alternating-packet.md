# An all-size precision witness with binary packet tensors

## Result and scope

There is an explicit family of actual Boolean Cauchy packets with binary
tensor entries whose original contraction is one, but whose complete local
polynomials at every odd center vanish modulo

```math
2^{K_L},\qquad K_L=L(1+\log_2 L)-1.
\tag{1}
```

This holds for every power of two `L>=4`. The entries are already their
canonical representatives at creation precision one. Thus no later choice
of rational lifts is involved. Equality of local content holds at `t=L-1`,
so (1) is the exact minimum over the odd centers.

The result upgrades the four-wire binary example in
[the packet-realization note](cauchy-packet-realization.md) to an all-size
`Theta(L log L)` family. It concerns the odd centers, or retained valuation
level `u=1`. It does not supply the `Theta(L^2)` loss of the abstract tapered
witness with `M=L,u=1`, and it does not improve the factoring exponent.

The deformation and the canonical-creation convention are from
[OpenAI's quantitative contraction proof](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`. The proof below is an
ordinary exact argument; its numerical replay is supplementary.

## 1. The binary packet

Write `L=2^a`, `a>=2`, and `n=L/2`. Put

```math
C(T)=\prod_{q=-L}^{L-1}(T-q),\qquad
W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},\qquad 0\le x,y<L.
\tag{2}
```

Let `X` have a one in each even position and a zero in each odd position.
Both tensors are the indicator of this one state:

```math
A(Y)=\mathbf1_{Y=X},\qquad B(Y)=\mathbf1_{Y=X}.
\tag{3}
```

The deformed packet contraction is therefore exactly

```math
P_L(T)=
 \prod_{\substack{0\le x<L\\x\text{ even}}}(W(T)X)_x
 \prod_{\substack{0\le x<L\\x\text{ odd}}}
       \bigl(1-(W(T)X)_x\bigr).
\tag{4}
```

Since `W(0)=I`, `P_L(0)=1`. The description (3) is uniform and its entries
are zero or one. In particular, these are canonical one-bit input tensors.

For a nonzero polynomial over `Q`, let its two-adic content be the minimum
valuation of its coefficients. Contents add under multiplication. This
follows by dividing out their contents and reducing the resulting primitive
polynomials modulo two; their product is nonzero over `F_2[z]`.

## 2. A logarithmic-derivative identity

The following elementary local identity provides the extra logarithm.
For `n=2^b`, `b>=1`, define

```math
E_n(V)=\prod_{\substack{1\le r<n\\r\text{ odd}}}(V^2-r^2),
\qquad H_n(V)=\frac{E_n'(V)}{E_n(V)}.
\tag{5}
```

There is a power series `U(Z)` with coefficients in `Z_2` and odd constant
coefficient such that, for two-adically even `V`,

```math
H_n(V)=-nV\,U(V^2).
\tag{6}
```

To prove this, put `N=n/2` and

```math
S_k=\sum_{\substack{1\le r<n\\r\text{ odd}}}r^{-2k},\qquad k\ge1.
```

Inversion permutes the odd residue classes modulo `n`, and

```math
\sum_{\substack{1\le r<2^b\\r\text{ odd}}}r^{2k}
 \equiv 2^{b-1}\pmod {2^b}.
\tag{7}
```

For completeness, (7) follows by induction on `b`: pair an odd residue
`r<2^(b-1)` with `r+2^(b-1)`. Their even powers agree modulo `2^b`.
The base `b=1` is immediate. Hence `S_k/N` is integral and odd. Pairing
the roots `r,-r` in (5), then expanding each unit denominator, gives

```math
H_n(V)
 =2V\sum_{r\text{ odd},\ 0<r<n}\frac1{V^2-r^2}
 =-nV\sum_{j\ge0}\frac{S_{j+1}}N V^{2j}.
\tag{8}
```

This proves (6). The series converges coefficientwise after substituting
`V=V_0+2Lz` with `V_0` even. Its `U` factor is an integral power series
whose constant coefficient is odd, and therefore is a unit of `Z_2[[z]]`.

It follows that, on this disk,

```math
\mathrm{cont}_2 H_n(V_0+2Lz)
 =b+\min\{v_2(V_0),a+1\},
\tag{9}
```

where `v_2(0)=+infinity`. This use of content for a rational function is
unambiguous: its denominator has odd constant coefficient, and multiplying
by that denominator preserves content. The same argument applies after
multiplying by `E_n`, which is a unit on the disk.

## 3. Exact contents of the occupied row factors

Fix an odd `t` in `{-L,...,L-1}`. Every nonzero difference `t-q` has
valuation at most `a`, whereas the dilation `2L` has valuation `a+1`.
Consequently,

```math
\mathrm{cont}_2\frac{C(t+2Lz)}{C'(0)}=a+1.
\tag{10}
```

Indeed, the factor with `q=t` contributes `a+1`. The remaining factors
contribute `v_2(C'(t))`. That valuation equals `v_2(C'(0))`: with
`r=t+L` and `m=2L-1`, it is
`v_2(r!)+v_2((m-r)!)=m-(a+1)`, because the binary expansions of `r` and
`m-r` are complementary.

For an even row `x`, the selected differences `y-x` form the progression
`-x,2-x,...,L-2-x`. Therefore

```math
(W(T)X)_x
 =\frac{C(T)}{C'(0)}\,
   H_n\bigl(T+x-(n-1)\bigr).
\tag{11}
```

Combining (9)--(11), with `b=a-1`, gives the exact formula

```math
\mathrm{cont}_2 (W(t+2Lz)X)_x
 =2a+\min\{v_2(t+x-n+1),a+1\},\qquad x\text{ even}.
\tag{12}
```

All the other row factors in (4) have nonnegative content. One direct
check is that each individual `W_xy(t+2Lz)` is integral: when `t=y-x`,
its content is zero, and otherwise its content is
`a+1-v_2(t-y+x)>=1`. Taking sums and subtracting from one preserves
integrality.

As `x` runs through the even positions, `t+x-n+1` runs through `n`
consecutive even integers. Each has valuation at least one. Among any
`n=2^(a-1)` consecutive integers, exactly `n/2^j` are divisible by `2^j`
for `1<=j<=a-1`. This gives

```math
\sum_{x\text{ even}}
 \min\{v_2(t+x-n+1),a+1\}
 \ge n+\sum_{j=1}^{a-1}\frac n{2^j}
 =L-1.
\tag{13}
```

The occupied row factors contribute at least `2an+(L-1)=K_L`.
The remaining factors cannot reduce content, so

```math
\mathrm{cont}_2 P_L(t+2Lz)\ge K_L
\quad\hbox{for every odd }t.
\tag{14}
```

At `t=L-1`, each odd-row factor has constant coefficient one and content
zero. The even quantities in (13) are `n,n+2,...,3n-2`, whose valuations
sum to `L-1`; no clipping changes this sum. Thus equality holds in (14).

## 4. What precision this rules out

Suppose a reconstruction procedure receives only the local coefficients at
the odd centers, each reduced modulo `2^K_L`. The packet (3) produces the
same supplied local data as the zero contraction, but their values at zero
differ modulo two. Those data alone therefore do not determine the target
bit, even within the class of canonical binary packet tensors.

This is a statement about the local-data reconstruction interface. It is
not a lower bound on algorithms that inspect or use the tensors by other
means. The family also leaves the other valuation levels uncontrolled and
does not establish the stronger wire-scale tapered precision obstruction.

## 5. Exact replay

```sh
python3 tools/check_canonical_alternating_packet.py \
  --output /tmp/canonical-alternating-packet.json
```

The checker forms the shifted root product with unbounded integers and
divides it exactly by each linear factor. It checks (12) from the actual
row-polynomial coefficients, sums their contents, and checks every odd
center for `L=4,8,16,32,64`. For `L=4,8`, an independent rational packet
implementation also expands the complete contraction. The replay checks
the original value one, the sharp endpoint, the exceptional excluded
length two, and an even-center control.

The stored result is
[canonical-alternating-packet.json](../evidence/canonical-alternating-packet.json).
These finite checks do not replace the all-size proof above.
