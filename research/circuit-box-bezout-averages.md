# Count-independent precision for circuit-isometric Bezout pairs

## Result and scope

The complete compiler in [circuit-box-bezout.md](circuit-box-bezout.md)
allows arbitrary-degree integer circuits in the index maps. Its short
factorial ratios can nevertheless require many guard bits. For example,
at output precision six, a domain of length `2^80` requires 86 coefficient
bits in that implementation.

The normalized averages in
[circuit-box-averages.md](circuit-box-averages.md) remove this dependence.
For output precision `p >= 2`, every normalized divided difference can be
computed using at most

\[
 W_{\max}(p)=p+v_2((p-2)!)+1
            =2p-1-\mathrm{popcount}(p-2)<2p
\]

coefficient bits. At precision one, the answer is simply the odd unit one.
The bound does not depend on the number of roots, the valuation of that
number, or the formal degrees of the input circuits.

This gives the same complete normalized Bezout pair, with the same exact
optimal exponent and fewer than `4p` stored coefficients per polynomial.
At `m = 2^80`, the new implementation uses ten coefficient bits and returns
exactly the pair from the 86-bit implementation.

The result applies to the certified separated isometric domains already
defined in the earlier note. It does not evaluate a singular denominator
whose valuation depends on unknown factors. The complete factoring
exponent remains unchanged.

## The divided-difference problem

Let `z_0,z_1,...` be dyadic integers satisfying

\[
 z_i\equiv i\pmod2,
 \qquad
 v_2(z_i-z_j)=v_2(i-j)\quad(i\ne j).
\]

In the full compiler, the two parity classes have globally regular
rational-circuit descriptions. More explicitly,

\[
\begin{aligned}
 z_{2i}
 &=\frac{a-1}{2^h}
    +\frac{D_A}{2^h}\bigl(i+2g_A(i)\bigr),\\
 z_{2i+1}
 &=\frac{(N-b)/2^h-(D_B/2^h)(i+2g_B(i))}
         {b+D_B(i+2g_B(i))}.
\end{aligned}
\]

The circuits `g_A,g_B` use only integer constants, addition and
multiplication. The second denominator is odd at every integer input.
The established conditions are

\[
\begin{gathered}
 h\ge1,\quad a,b\equiv1\pmod{2^{h+1}},\quad
 v_2(D_A)=v_2(D_B)=h+1,\\
 N\equiv1+2^h\pmod{2^{h+1}}.
\end{gathered}
\]

These conditions prove the isometry property. No numerical list of nodes
is part of the input.

For `k >= 1`, let `f[z_0,...,z_k]` denote the divided difference of the
parity function. Define its normalized odd unit by

\[
 D_k=(-1)^{k+1}k!\,2^{1-k}f[z_0,\ldots,z_k].
\]

Write

\[
 E_L(T)=\prod_{i=0}^{L-1}(1-iT),\qquad
 Z_L(T)=\prod_{i=0}^{L-1}(1-z_iT),
\]

and

\[
 H_{k,j}=[T^j]\frac{E_{k+j}(T)}{Z_{k+1}(T)},\qquad
 B_{k,j}=\prod_{t=1}^{j}(k+t),\qquad B_{k,0}=1.
\]

The parity Mahler expansion gives

\[
 D_k=\sum_{j\ge0}(-2)^j\frac{H_{k,j}}{B_{k,j}}.
\tag{1}
\]

Every quotient `H_{k,j}/B_{k,j}` is dyadically integral. One direct proof
uses the identity

\[
 \frac{H_{k,j}}{B_{k,j}}
 =k!\,[z_0,\ldots,z_k]\binom{X}{k+j}.
\]

The binomial polynomial takes integral dyadic values on dyadic integers.
In its Lagrange divided-difference expression, the coefficient at node
`i` has valuation

\[
 v_2(k!)-\sum_{\ell\ne i}v_2(z_i-z_\ell)
 =v_2\binom{k}{i}\ge0.
\]

Thus terms with `j >= p` vanish modulo `2^p`. The previous implementation
computed the short numerators directly before dividing by `B_{k,j}`. Its
guard is

\[
 \max_{0\le j<p}\bigl(v_2(B_{k,j})-j\bigr)
 =\mathrm{popcount}(k)
   -\min_{0\le j<p}\mathrm{popcount}(k+j).
\]

At `k = 2^s-1`, this can be of order `s`. The new construction absorbs
that large factor into a known-length normalized average.

## One large factor, followed by a small denominator

Fix `1 <= j < p`, and set `M=k+1`. Select `s` in `0,...,j-1` so that

\[
 M'=M+s
\]

has maximal dyadic valuation among `M,...,M+j-1`. The maximizer is unique.
Two integers with the same maximal valuation would have a more divisible
integer between them.

Factor

\[
 B_{k,j}=M'R_j.
\]

For every `t != s` in the short interval,

\[
 v_2(M+t)=v_2(t-s).
\]

Consequently,

\[
 \ell_j:=v_2(R_j)
 =v_2(s!)+v_2((j-1-s)!)
 \le v_2((j-1)!).
\tag{2}
\]

The odd part of `R_j` can have many bits. Its dyadic valuation is bounded
by the short length `j`, independently of `k`.

Introduce the correction polynomial

\[
 C_j(T)=
 \prod_{r=M'}^{M+j-2}(1-rT)
 \prod_{i=M}^{M'-1}(1-z_iT).
\]

Empty products equal one. Its degree is at most `j-1`, and

\[
 \frac{E_{M+j-1}(T)}{Z_M(T)}
 =C_j(T)\frac{E_{M'}(T)}{Z_{M'}(T)}.
\tag{3}
\]

Only `s <= j-1` extra nodes occur in this correction.

## Integral normalized moments

For `l >= 1`, define

\[
 A_l(M')=\frac1{M'}\sum_{i=0}^{M'-1}(z_i^l-i^l).
\tag{4}
\]

These numbers are dyadically integral. To see this, write
`M'=2^r u` with `u` odd. In every consecutive block of `2^r` indices,
the isometry sends the indices to a permutation of all residues modulo
`2^r`. The two power sums in (4) therefore agree modulo `2^r`. Division
by the remaining odd factor is harmless.

The normalized-average engine computes (4) without allocating `r` guard
bits. Let

\[
 \Lambda_L(q)=\frac2L\sum_{i=0}^{L-1}q(i).
\]

The earlier theorem computes this integral quantity modulo `2^q` using
only `q` coefficient bits for any globally regular circuit. Let `a,b` be
the numbers of even and odd indices in the prefix of length `M'`.
Let `Lambda_e`, `Lambda_o` be the normalized averages of the corresponding
node powers, and let `Lambda_id` be the normalized average of `i^l` on
the full prefix. If `M'` is even, then

\[
 A_l(M')=\frac{\Lambda_e+\Lambda_o-2\Lambda_{id}}4.
\tag{5}
\]

If `M'` is odd, then

\[
 A_l(M')=
 \frac{a\Lambda_e+b\Lambda_o-M'\Lambda_{id}}{2M'}.
\tag{6}
\]

An empty odd prefix contributes zero. Formula (5) needs two extra average
bits. Formula (6) needs only one; the implementation uses two in both
cases. The divisibility follows from (4), rather than from a heuristic
cancellation test.

The rational node circuit must use an actual globally odd inverse gate
in this step. Replacing it prematurely by a pointwise approximation can
lose precision when the average divides by a highly even count. The
regular-circuit quotient evaluates the inverse inside the
normalized-average algebra and preserves the required functional.

## Weighted exponential formula

As formal power series over the dyadic rationals,

\[
 \frac{E_{M'}(T)}{Z_{M'}(T)}
 =\exp\left(M'\sum_{l\ge1}\frac{A_l(M')}{l}T^l\right).
\]

Because `C_j` has degree below `j`, equation (3) yields

\[
 \frac{H_{k,j}}{M'}
 =[T^j]C_j(T)
 \frac{\exp\left(M'\sum_{l=1}^j A_l(M')T^l/l\right)-1}{M'}.
\tag{7}
\]

Keeping the factor `2^j` from (1) is essential. Define the integral series

\[
 G(T)=\sum_{l=1}^j
  2^{l-1-v_2(l)}\frac{A_l(M')}{\mathrm{odd}(l)}T^l,
\]

where `odd(l)=l/2^v2(l)`. The exponent of two is nonnegative for every
positive `l`. After replacing `T` by `2T` in (7),

\[
 \frac{2^jH_{k,j}}{M'}
 =[T^j]C_j(2T)\,V(T),
\tag{8}
\]

where, through degree `j`,

\[
\begin{aligned}
 V(T)
 &=\frac{\exp(2M'G(T))-1}{M'}\\
 &=\sum_{r=1}^j
   2^{\mathrm{popcount}(r)}
   \frac{(M')^{r-1}}{\mathrm{odd}(r!)}G(T)^r
   \pmod{T^{j+1}}.
\end{aligned}
\tag{9}
\]

Here the identity `r-v2(r!)=popcount(r)` removes every nonunit factorial
division. Each displayed summand has an explicit factor of two. All
remaining divisions in (9) are by known odd integers.

Finally, divide the coefficient in (8) by `R_j`. This gives exactly the
weighted contribution `2^j H_{k,j}/B_{k,j}` needed in (1).

## Precision accounting

Set

\[
 W_j=p+\ell_j.
\]

Compute the numerator in (8) modulo `2^W_j`, then divide its known factor
`2^ell_j` and invert the odd part of `R_j`. The result is correct modulo
`2^p`.

The moments need only `W_j-1` bits. An error divisible by `2^(W_j-1)` in
any moment gives such an error in `G` and every power of `G`. The explicit
factor of two in every term of (9) raises it to `2^W_j`.

Equations (5) and (6) therefore require normalized averages to precision

\[
 (W_j-1)+2=p+\ell_j+1.
\]

Using (2) and `j <= p-1`, the maximum is at most

\[
 p+v_2((p-2)!)+1<2p.
\]

The degree of each normalized-average quotient is `O(p)`. All node
circuits can have unrestricted formal degree. The total work is
polynomial in the written circuits, `p`, and the bit lengths of the
counts and parameters. The algorithm evaluates only polynomially many
boundary nodes and never enumerates a long interval.

This precision statement concerns coefficients in the modular algebras.
Binary counts, short factorial-ratio products and exact power-sum
functionals still carry their ordinary bit cost. For example, a
power-sum functional can contain integers with `O(p log M')` bits. The
algorithm does not treat these exact integers as unit-cost objects.

## Full normalized Bezout pairs

The existing full compiler uses the normalized units `D_k` only through
a small evaluator interface. Its Newton-support argument, odd factorial
units and reversed root-product jets remain valid.

For

\[
 P(X)=\prod_{i=0}^{m-1}(X-A_i),\qquad
 Q(X)=\prod_{i=0}^{m-1}(N-B_iX),
\]

the same unique rational pair of degree below `m` is normalized by

\[
 r=h(2m-1)-\mathrm{popcount}(m-1).
\]

The new evaluator returns all coefficients of `U=2^rU_0` and `V=2^rV_0`
modulo `2^p`, with `UP+VQ=2^r` and the proved optimal exponent `r`.
Each output is stored using fewer than `4p` reversed coefficients in
`X-1`; all other coefficients vanish modulo `2^p`.

The older evaluator can be preferable for a particular input. Its guard
may already be zero. Choosing the smaller of its known precision bound
and the new bound gives an immediate adaptive variant. The reference
implementation deliberately exposes the new evaluator separately, so
its precision receipt can be inspected directly.

## Implementation and evidence

Files:

- [Compiler](../tools/circuit_box_bezout_averages.py).
- [Independent checker](../tools/check_circuit_box_bezout_averages.py).
- [Deterministic receipt](../evidence/circuit-box-bezout-averages.json).

Replay:

```sh
python3 tools/check_circuit_box_bezout_averages.py \
  --output evidence/circuit-box-bezout-averages.json
```

The receipt contains 672 literal rational normalized-moment comparisons,
216 exact weighted factorial-ratio comparisons, 108 exact normalized
divided differences, and 48 complete rational Bezout-pair comparisons
covering 288 coefficients. The exact checks use rational Lagrange
interpolation or finite rational power series, independently of the
normalized-average implementation.

Two large complete-pair controls use the input perturbation
`g(i)=i^(2^200+1)` on both domains at output precision six. They compare
every stored coefficient against the earlier compiler and check the
global identity independently at unit points in a quadratic extension.

| Domain length | Earlier coefficient precision | New coefficient precision |
| --- | ---: | ---: |
| `2^80+17` | 6 | 9 |
| `2^80` | 86 | 10 |

A separate closed-form control uses `z_i=5i+6`,
`k=2^4096-1`, and output precision eight. Scaling the arithmetic nodes by
five gives the exact normalized unit `D_k=5^(-k)`, which is five modulo
256. The old guard would require 4103 coefficient bits. The new evaluator
returns five using thirteen bits and does not enumerate any prefix.

These tests establish reproducible evidence for the new normalization
formula. The proof supplies its uniform scope. Neither the tests nor the
formula provide an unknown-factor pole count or a stronger complete
factoring exponent.
