# Polynomials realized by a Boolean Cauchy packet

## Result and scope

The Cauchy deformation used in OpenAI's tensor evaluator has a concrete
polynomial realization property. Let `L>=2` be a power of two. Every rational
polynomial of degree at most `L` that is two-adically integral on the complete
integer grid `{-L,...,L-1}` is the exact deformed contraction of one `L`-wire
packet. Its row tensor is affine, and its column tensor is supported on zero
and the `L` one-hot states. The construction uses polynomial-time rational
arithmetic and supplies both tensors explicitly.

This is an ordinary all-size construction theorem. It does not establish a
better factoring exponent. Its degree range is too small to realize the
full tapered-Hermite witness at the source's wire scale. Two additional
finite certificates below belong to the actual packet class, but neither
proves a new asymptotic precision lower bound.

The source deformation and its creation-precision convention are in
[OpenAI's quantitative contraction proof](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`.
The separate [tapered-Hermite witness in PR #31](https://github.com/1337isnot1337/integer-factoring/blob/b1ba452fa4a3d19fe6fea3822c53613f3b582104/research/tapered-hermite-precision-witness.md)
concerns a larger polynomial class and explicitly leaves actual contraction
realization open. The result here supplies a restricted positive realization
and keeps that remaining distinction explicit.

## 1. The packet and the claimed realization

Put

```math
S=\{-L,\ldots,L-1\},\qquad
C(T)=\prod_{q\in S}(T-q),\qquad
W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},\qquad 0\le x,y<L.
\tag{1}
```

The apparent denominator cancels one root of `C`, so every matrix entry is
a polynomial of degree `2L-1`. At zero, `W(0)=I`.
For Boolean states `X,Y`, write `l_0(z)=1-z` and `l_1(z)=z`.
The deformed contraction of row tensor `A` and column tensor `B` is

```math
F(T)=\sum_{X,Y\in\{0,1\}^L}A(Y)B(X)
             \prod_{x=0}^{L-1}l_{Y_x}((W(T)X)_x).
\tag{2}
```

The coefficient ring for the theorem is `Z_(2)`, the rationals with odd
denominators. Assume `deg f<=L` and `f(q) in Z_(2)` for every `q in S`.
We construct `A,B` over this ring with `F(T)=f(T)` as rational polynomials.

Let `m=2L-1`, and let `Delta^j f(-L)` be the ordinary forward difference.
Define

```math
a_x=(-1)^{L-1-x}\binom{L-1}{x},\qquad c=f(-L),
\tag{3}
```

and let `b_y` be the coefficients of

```math
b(z)=\frac1{\binom{2L-1}{L}}
 \sum_{j=1}^{L}(-1)^j\binom{2L-1}{j}
  \Delta^j f(-L)\,z^{j-1}(1-z)^{L-j}.
\tag{4}
```

This has degree at most `L-1`. Set

```math
A(Y)=c+\sum_x a_xY_x,\qquad
B(0)=1-\sum_y b_y,\qquad B(e_y)=b_y,
\tag{5}
```

and set every other column entry to zero. These are explicit tensor-entry
algorithms. All coefficients lie in `Z_(2)`: forward differences use only
the integral values `f(-L),...,f(0)`, and `binom(2L-1,L)` is odd. One way
to see the latter fact is that, over `F_2`,

```math
(1+z)^{2L-1}=1+z+\cdots+z^{2L-1}.
```

## 2. Proof of the identity

For `q in S`, put `w_q=C'(q)/C'(0)`. Every `w_q` is a two-adic unit:
the grid is a complete residue system modulo `2L`, so all products of
nonzero differences have the same two-adic valuation. Direct factorial
cancellation gives

```math
\frac1{w_q}
 =(-1)^q\frac{\binom{2L-1}{L+q}}{\binom{2L-1}{L}}.
\tag{6}
```

The polynomials `C(T)/(C'(q)(T-q))` form the Lagrange basis on `S`.
Consequently,

```math
f(T)=f(-L)+
 \sum_{q=-L+1}^{L-1}c_q\frac{C(T)}{C'(0)(T-q)},
\qquad c_q=\frac{f(q)-f(-L)}{w_q}.
\tag{7}
```

Let `g(k)=f(k-L)-f(-L)`, so `deg g<=L` and `g(0)=0`.
Newton expansion and the binomial identity

```math
\sum_{k=0}^{m}(-1)^k\binom mk\binom kj z^k
 =(-1)^j\binom mj z^j(1-z)^{m-j}
\tag{8}
```

show that

```math
\sum_{q=-L+1}^{L-1}c_qz^{q+L-1}
 =(1-z)^{L-1}b(z).
\tag{9}
```

The missing constant term before division by `z` is zero because `g(0)=0`.
The sign in (6) is `(-1)^k` for `k=L+q`, since `L` is even.
Also `Delta^j g(0)=Delta^j f(-L)` for `j>=1`, giving exactly (4).

Since

```math
\sum_x a_x z^{L-1-x}=(1-z)^{L-1},
```

equation (9) says `c_q=sum_(y-x=q) a_x b_y`.
Substitution in (7) gives

```math
f(T)=c+\sum_{x,y}a_xW_{xy}(T)b_y.
\tag{10}
```

The Boolean interpolation of the affine function in (5) is the same
affine function over every commutative ring. Thus summing over `Y` in (2)
gives `c+sum_x a_x(W(T)X)_x`. Finally,

```math
\sum_X B(X)=1,\qquad \sum_XB(X)X_y=b_y.
```

Summing over `X` proves (10), hence `F=f`.

The construction uses `O(L^2)` arithmetic operations after evaluating `f`
at the required integers. For coefficients or sampled values of height `H`,
the rational bit lengths grow by at most a polynomial in `L+H`; the explicit
binomial factors have `O(L)` bits. Standard rational arithmetic therefore
gives a uniform polynomial bit-time construction. It does not enumerate
the `2^L` entries of either tensor.

## 3. An exact family of local contents

The polynomial

```math
f_L(T)=\prod_{\substack{q\in S\\q\text{ odd}}}
                   \left(1-\frac Tq\right)
\tag{11}
```

has degree `L`, coefficients in `Z_(2)`, and value one at zero. Therefore
the theorem realizes it by an actual packet. For every odd center `t in S`,

```math
\min_k v_2\bigl([z^k]f_L(t+2Lz)\bigr)=2L-1.
\tag{12}
```

For an even center the same content is zero. To prove (12), separate the
factor whose root is `t`; its content is `log_2 L+1`. The other odd-root
differences are twice a complete set of nonzero differences in a grid of
length `L`. Their total valuation is

```math
(L-1)+v_2(r!)+v_2((L-1-r)!)=2L-\log_2L-2.
```

Here `0<=r<L` is the position of `t` among the odd roots. The identity
`v_2(r!)+v_2((L-1-r)!)=L-1-log_2 L` follows from the binary-complement
relation between `r` and `L-1-r`. Contents add under multiplication over
`Z_(2)[z]`, because reducing primitive polynomials modulo two preserves a
nonzero product. Adding the two displayed contributions proves (12).

This is a realizable precision phenomenon, but its scale is `Theta(L)`
on `L` real wires and only the odd centers are suppressed. It is not the
source's full `Theta(MD/u)` tapered pattern.

## 4. Two finite certificates in the actual packet class

The checker also verifies two examples independently of the affine theorem.

| Packet | Tensor entries | Retained centers | Verified content |
|---|---|---|---|
| `L=4` | Both tensors are the binary indicator of state `0101` | `-3,-1,1,3` | Respectively `12,14,13,11` |
| `L=8` | Row indicator of zero; explicit 256-entry column tensor | All nonzero `t` with `v_2(t)<3` | At least `33` at every center |

Both original contractions are odd. The first has value one and uses only
binary tensor entries, so its creation precision can be one. It establishes
only a four-wire, odd-center example. The second has value `494022085` at
zero and is supplied in [the finite certificate](../evidence/cauchy-packet-witness.json).
Every column entry is the canonical representative modulo `2^33`.
Its input precision is therefore 33; it does **not** establish an additional
33-bit loss beyond creation precision.

The certificate checker reconstructs (1), shifts the rational polynomials,
and reduces each coefficient using only odd denominators. It then performs
the literal polynomial products and sums in (2) with the row indicator fixed.
Every possible coefficient through degree `L(2L-1)` is checked, including
the trailing zero coefficients. No numerical matrix rank is used to validate
either witness.

## 5. Canonical lifting is a real restriction

A realization over `Z_(2)` cannot automatically be used as a witness for a
low-precision input to the source algorithm. The source fixes canonical
representatives when a tensor is created; later increases in precision do
not replace those representatives by more favorable rational lifts.

There is already a concrete difference at `L=2`. The construction for
`f(T)=1-T^2` gives, in binary-state order `0,1,2,3`,

```math
A=(-3,-4,-2,-3),\qquad B=(3,-3,1,0).
\tag{13}
```

Every coefficient around either odd center is divisible by eight. But
canonical representatives modulo two give

```math
\overline A=(1,0,0,1),\qquad\overline B=(1,1,1,0).
```

At `T=1`, the only nonzero entry of `W` is `W_01=-3`. Direct substitution
in (2) gives `F_canonical(1)=6`, which is not divisible by eight.
The original contraction remains odd. The difference comes from the chosen
lifts, not from an error in the exact realization identity.

A useful asymptotic obstruction for the actual evaluator would have to
control the packet geometry, every retained valuation level, the creation
precision, and the additional loss simultaneously. An improved evaluator
could instead exploit restrictions of that canonical tensor class. The
present construction and finite examples settle neither alternative.

## 6. Reproduction and evidence

```sh
python3 tools/check_cauchy_packet_polynomials.py \
  --output /tmp/cauchy-packet-polynomials.json
```

The [recorded replay](../evidence/cauchy-packet-polynomials.json) checks the
global rational identity for 28 polynomials through `L=16`, including
integer-valued binomial polynomials with even coefficient denominators.
Two cases are also evaluated by the full Boolean sum. The replay checks
all local contents in (12) for those four packet lengths, both finite
packet witnesses, and the canonical-lift counterexample. Invalid degree
and integrality inputs are rejected.

The proof and exact implementation used OpenAI assistant support. The
Cauchy deformation and the canonical-lift convention retain their source
attribution. No new Lean formalization or improved factoring score is claimed.
