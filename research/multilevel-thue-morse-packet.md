# Canonical binary packet collisions at several valuation levels

## Result and scope

Let `L=2^a`, with `a>=2`, and choose `1<=u<=a`. There are explicit
canonical binary tensors for one Boolean Cauchy packet such that their
original contraction is one and every complete local polynomial at

```math
-L\le t<L,\qquad v_2(t)<u
```

is divisible coefficientwise by

```math
2^{K_{a,u}},\qquad K_{a,u}=\frac{L}{4}(a-u+1).
\tag{1}
```

The construction uses a repeated Thue--Morse residue pattern. Its entries
are zero or one at creation precision one. The statement controls all
retained valuation levels below `u` simultaneously. For example,
`u<=a/2` gives a content lower bound at least `L a/8` on all those
levels. At `u=a`, it covers every nonzero center with `|t|<L`, with
the smaller bound `L/4`.

This is a reconstruction obstruction for the displayed local-data
interface. It does not improve the factoring exponent or exclude a
decoder that uses the tensors by another means. The known
[alternating construction](canonical-alternating-packet.md) remains
sharper at `u=1`: its exact minimum is `L(a+1)-1`. The new point is
simultaneous control of the other valuation levels with canonical
one-bit entries, rather than a favorable choice of high-precision lifts.

The deformation is the packet from
[OpenAI's quantitative contraction proof](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`. The proof below
uses its exact polynomial kernel. It does not assume an evaluator
runtime theorem beyond that definition.

## 1. The packet and a general residue-pattern lemma

Put

```math
C(T)=\prod_{q=-L}^{L-1}(T-q),\qquad
W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},\qquad 0\le x,y<L.
\tag{2}
```

Each entry is a polynomial, and `W(0)=I`. Let `m=2^u`, and choose
any subset `R` of the residues modulo `m`. Define the binary state

```math
X_x=[x\bmod m\in R],\qquad A=B=\mathbf1_{\{X\}}.
```

The exact deformed contraction is

```math
P(T)=\prod_{X_x=1}(W(T)X)_x
     \prod_{X_x=0}\bigl(1-(W(T)X)_x\bigr),\qquad P(0)=1.
\tag{3}
```

For `s` modulo `m`, set

```math
D_R(s)=|\{r\in R:r+s\notin R\}|.
```

**Residue-pattern bound.** For every center `t` not divisible by `m`,

```math
\mathrm{cont}_2 P(t+2Lz)
\ge 2(a-u+1)\frac{L}{m}D_R(t\bmod m).
\tag{4}
```

Here the content of a nonzero rational polynomial is the minimum
two-adic valuation of its coefficients. Contents add under
multiplication, since the product of two nonzero polynomials over
`F_2` is nonzero after their powers of two are removed.

## 2. An arithmetic-progression sum

If `A` is a unit of `Z_2[[z]]`, `D` is an even integer, and `b>=0`,
then

```math
\sum_{k=0}^{2^b-1}\frac1{A+Dk}\in 2^b Z_2[[z]].
\tag{5}
```

The assertion at `b=0` is immediate. Pair the terms with indices
`k` and `k+2^(b-1)`. Their sum is

```math
\frac2{A+Dk}
-\frac{D2^{b-1}}{(A+Dk)(A+D(k+2^{b-1}))}.
```

The sum of the first terms is divisible by `2^b` by induction.
Every second term is divisible by `2^b` because `D` is even and
both denominators are units. This proves (5) coefficientwise.
It also holds with negative `D`.

On the disk `T=t+2Lz`, the numerator in (2) has exact content

```math
\mathrm{cont}_2\frac{C(t+2Lz)}{C'(0)}=a+1.
\tag{6}
```

The root `q=t` contributes `a+1`. The other factors have contents
`v_2(t-q)`. Their sum is `v_2(C'(t))=v_2(C'(0))`, since the
`2L` consecutive grid points are a complete residue system modulo
`2L`. Equivalently,
`v_2(C'(t))=v_2((t+L)!)+v_2((L-1-t)!)=2L-a-2`.

Fix an occupied row `x` with `(x+t) mod m` outside `R`.
For each occupied residue `r`, all denominators

```math
T+x-r-mk,\qquad 0\le k<L/m,
```

have the same valuation `j=v_2(t+x-r)<u`. Divide them by `2^j`.
Their initial value is a power-series unit and their step has
valuation `u-j>=1`. Apply (5), with `b=a-u`. Thus

```math
\sum_{k=0}^{L/m-1}\frac1{T+x-r-mk}
\in 2^{a-u-j} Z_2[[z]]
\subseteq 2^{a-2u+1} Z_2[[z]].
```

Sum over `r in R` and multiply by (6). The occupied row factor in
(3) therefore has content at least `2(a-u+1)`. This argument does
not cancel a vanishing denominator: the mismatch condition gives
`j<u` for every term.

Every other row factor has nonnegative content. Indeed, each
`W_xy(t+2Lz)` is integral: its content is zero if `t=y-x`, and
otherwise is `a+1-v_2(t-y+x)>=1`. There are exactly
`(L/m)D_R(t mod m)` mismatched occupied rows. Add their content
bounds to obtain (4).

## 3. The Thue--Morse pattern has no almost-periodic shift

Choose

```math
R_u=\{0\le r<2^u:\text{the binary digit sum of }r\text{ is even}\}.
```

Write `epsilon_u(r)=(-1)^(digit sum(r))`, extended periodically
modulo `m=2^u`, and use the normalized cyclic autocorrelation

```math
c_u(s)=\frac1m\sum_{r=0}^{m-1}\epsilon_u(r)\epsilon_u(r+s).
```

Splitting even and odd indices gives the exact recurrences

```math
c_u(2s)=c_{u-1}(s),\qquad
c_u(2s+1)=-\frac{c_{u-1}(s)+c_{u-1}(s+1)}2.
\tag{7}
```

The adjacent sums at level `u` are one half of a difference of
two correlations at level `u-1`; hence they lie in `[-1,1]`.
The base level is `c_1(0)=1,c_1(1)=-1`. Induction in (7)
therefore gives

```math
c_u(s)\le\frac12\qquad(s\not\equiv0\pmod {2^u}).
\tag{8}
```

For an even nonzero shift, use the induction hypothesis. For an
odd shift, use the lower bound `-1` for the adjacent sum in (7).
The pattern is balanced, so the numbers of `1 -> -1` and
`-1 -> 1` transitions under a cyclic shift agree. Consequently

```math
D_{R_u}(s)=\frac{m}{4}(1-c_u(s))\ge\frac m8.
\tag{9}
```

Combining (4) and (9) proves (1). The pattern is explicit and its
membership can be evaluated from the `u` low binary digits of a
position. The two tensor entries merely test equality to that
explicit state, so their description has no hidden advice.

## 4. The precise decoder consequence

The packet in (3) and the zero column tensor produce the same
complete local data modulo `2^K_(a,u)` on every disk with
`v_2(t)<u`. Their original values differ modulo two. Even knowing
the fixed row tensor does not remove this local-data collision.

Thus those data alone do not determine the requested bit at this
precision. Input-dependent changes of the packet coordinates,
additional disks, or direct use of the input tensors remain outside
this statement. It is not a factoring lower bound, and it does not
establish the stronger wire-scale tapered loss `Theta(L^2/u)`.

## 5. Exact finite controls

```sh
python3 tools/check_multilevel_thue_morse_packet.py \
  --max-length 32 --output evidence/multilevel-thue-morse-packet.json
```

The checker forms every shifted root product with integers, divides
it exactly by each original linear denominator, and reconstructs
every actual row polynomial of (3). It checks every retained disk
for every `u` through packet length 32. Polynomial contents are
computed from all row coefficients; no predicted row valuation is
substituted for the literal polynomial. Their sum is the exact
content of the full contraction by content additivity. Separate
cyclic correlations verify (7)--(9) through period 512.

The finite replay supplements the all-size proof above. The proof,
implementation, and finite certificate used OpenAI assistant support.
