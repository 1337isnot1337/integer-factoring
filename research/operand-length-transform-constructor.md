# Smaller multiplication instances from the two operand bit lengths

## Result and scope

Conditioning on the divisor's bit length halves the leading transform
volume in the [adaptive constructor](adaptive-transform-constructor.md).
If an `n`-bit integer is a product and its first factor has `ell` bits,
the second factor needs at most `n-ell+1` bits. The two transform inputs
can therefore share a total bit budget close to `n`.

**Theorem.** Let `N>=2` have binary length `n`. Let `P>=2` satisfy
`P^40<=n`. For every `1<=ell<=n` and divisor-prefix query
`0<=k<=n`, `0<=x<2^(n-k)`, one deterministic algorithm writes an
Appendix B instance `I_ell` such that

```math
\begin{aligned}
\mathrm{WFe}(P,\lambda,B_{\rm tab},I_\ell),\qquad
\mathrm{zval}(P,I_\ell)
 &=\#\{d:d\mid N,\ \mathrm{bitlength}(d)=\ell,
                 \ \lfloor d/2^k\rfloor=x\},\\
P^2\lambda&\le49n,\qquad B_{\rm tab}\le57n/P.
\end{aligned}                                                    \tag{1}
```

The complete construction costs `2^O(n/P)` bit operations, with
absolute constants. The same bounds hold for every `ell`. Uniformly
as `P` tends to infinity subject to the cutoff,

```math
P^2\lambda\le(18+o(1))n,
\qquad P B_{\rm tab}\le(20+o(1))n.                  \tag{2}
```

An unrestricted divisor-prefix count is the sum of these `n` clique
values. Polynomially many such calls recover the smallest prime factor.
The completed half-rate evaluator consequently retains its
`2^O(n/sqrt(log n))` bit bound. This improves the written constructor's
constants from `99,117` and its leading constants from `36,40`; it
does not change the established logarithmic exponent `1/2` or the score.

The proof is ordinary mathematics. It changes the input supports and
parameter choice, while retaining the exact composite-ring transform,
bounded carry chain, generic grouping, and Appendix B encoding.

## 1. Parameters and the asymmetric input supports

Let `h` be the least power of two at least `P^9`. Put

```math
\begin{gathered}
A=2h,\qquad D=2h^2,\qquad T=A^2=2D,\qquad G=3P^8,\\
L=T-A-1,\qquad w=\left\lceil\frac{n+1}{L}\right\rceil,
\qquad b=2^w,\qquad t=\log_2D,\\
m=D\left\lceil\frac{2w+t}{D}\right\rceil,
\quad M=2^m+1,\quad f=m+1,\quad c=w+t,\quad K=D(b-1).
\end{gathered}                                                    \tag{3}
```

Here `L` is only the denominator used to choose the digit width; it is
not a transform length. All parameters in (3) are independent of `ell`.
For the length-conditioned input, set

```math
d_p=\lceil\ell/w\rceil,
\qquad d_q=\lceil(n-\ell+1)/w\rceil.                 \tag{4}
```

The first sequence is supported on `0,...,d_p-1`, and the second on
`0,...,d_q-1`. Both are zero-padded to length `T`. The ceiling inequality
`ceil(a)+ceil(b)<=ceil(a+b)+1` gives

```math
d_p+d_q\le\lceil(n+1)/w\rceil+1
\le L+1=T-A,\qquad \min(d_p,d_q)\le D.              \tag{5}
```

Thus the convolution has no cyclic aliasing. It also has at most `D`
products contributing to any coefficient, including for highly
unbalanced input lengths.

The first factor has **exactly** `ell` bits. Its top digit has width
`r_p=ell-w(d_p-1)`, and belongs to
`[2^(r_p-1),2^r_p-1]`. Its remaining digits are unrestricted `w`-bit
digits. The second factor's top digit has width
`r_q=n-ell+1-w(d_q-1)` and belongs to `[0,2^r_q-1]`; its other digits
are unrestricted. Positions above the stated supports are fixed zeros,
not additional local settings. The combined logarithmic input-setting
count is exactly `(ell-1)+(n-ell+1)=n`.

For every divisor `p` of `N` with bit length `ell`, its complementary
factor satisfies

```math
0<q=N/p<2^n/2^{\ell-1}=2^{n-\ell+1}.                \tag{6}
```

Consequently no desired factor pair is lost. The exact leading-bit
constraint on the first factor makes the `ell` classes disjoint.

## 2. The complete three-stage transform for full input supports

The input gates may now contain more than `h` digits individually.
For each `r in {0,...,A-1}`, store all nonzero-support positions
`r+A i`, with `0<=i<A`, at that input gate. This is the only topology
change needed: the same two families of `A` input gates still each send
one residue to every middle gate.

For completeness, here are formulas that work for arbitrary length-`T`
sequences, so no theorem restricted to half-length input supports is
being reused. Put `omega=2^(m/D) mod M` and use the explicit inverse
of `T` from the adaptive constructor. For input sequences `p,q`, the
first family sends

```math
a_r(k_1)=\omega^{r k_1}
          \sum_{i=0}^{A-1}p_{r+A i}\omega^{A i k_1},
\qquad
b_r(k_1)=\omega^{r k_1}
          \sum_{i=0}^{A-1}q_{r+A i}\omega^{A i k_1}.
\tag{7}
```

The middle gate `k_1` computes

```math
\begin{aligned}
P_{k_1,k_2}&=\sum_r a_r(k_1)\omega^{A r k_2},\\
Q_{k_1,k_2}&=\sum_r b_r(k_1)\omega^{A r k_2},\\
e_{k_1}(j_1)&=\omega^{-k_1j_1}
 \sum_{k_2=0}^{A-1}P_{k_1,k_2}Q_{k_1,k_2}\omega^{-A k_2j_1}.
\end{aligned}                                                    \tag{8}
```

Finally the inverse gate `j_1` emits

```math
a_{j_1+A j_2}=T^{-1}
 \sum_{k_1=0}^{A-1}e_{k_1}(j_1)\omega^{-A k_1j_2}.
\tag{9}
```

Expanding (7) identifies `P_(k_1,k_2)` with the length-`T` Fourier
transform at frequency `k_1+A k_2`; the same holds for `Q`. Expanding
(8)--(9) gives its inverse transform, since the additional factor
`omega^(-A^2 k_2j_2)` equals one. The principal-root identities over
the composite ring were proved in the adaptive constructor from
`omega^D=-1`; they do not depend on either input support.

By (5), every integer coefficient satisfies

```math
0\le a_j\le\min(d_p,d_q)(b-1)^2
\le D(b-1)^2<2^{2w+t}\le2^m<M.                     \tag{10}
```

Hence the canonical residue in (9) is the exact integer coefficient.
The same carry bound `K=D(b-1)`, carry width `c`, zero initial and final
carries, digit checks, and clipped outgoing messages apply unchanged.
The rejection condition for a computed carry above `K` remains mandatory.

Fixing all input digits forces all middle settings, inverse settings,
coefficients, and carries uniquely. The carry equations accept exactly
when `pq=N`. Local prefix checks on the first factor use its fixed
zero extension through `n` bits. Equation (6) then proves the exact
length-conditioned divisor-count identity in (1), including `ell=1`
and `ell=n`.

## 3. Paired input blocks keep the same cluster bound

At residue `r`, let `h_p(r),h_q(r)` count the stored digits of the two
input sequences. They satisfy

```math
\begin{aligned}
h_p(r)+h_q(r)
&\le\lceil d_p/A\rceil+\lceil d_q/A\rceil\\
&\le\lceil(d_p+d_q)/A\rceil+1\le A.
\end{aligned}                                                    \tag{11}
```

The last step uses `d_p+d_q<=A^2-A`. Highest-digit restrictions can
only reduce the setting counts. Thus the two input gates at the same
residue together have at most `2^(Aw)` settings, even though one input
gate alone may have more than the previous `h` digits.

Use the same clustering map `r -> floor(rG/A)` for both input families,
the middle family, and the inverse family. Its pairing of equal
residues means (11) can be used inside each cluster. With
`u_A=floor(A/G)+2` and `u_T=floor(T/G)+2`, the previous bounds therefore
remain valid:

```math
\begin{aligned}
U&=A^2(8f+2c),& S&=A^2(w+4f+c),\\
Q&=7u_AAf+u_T(f+2c),&
W&=u_AA(w+3f)+u_T(f+c).
\end{aligned}                                                    \tag{12}
```

Here `U` bounds ordered wire volume, `S` bounds the sum of logarithmic
setting counts, `Q` bounds a cluster's weighted degree, and `W` bounds
its logarithmic setting count. In particular, the `Aw` term in `W`
is a bound on the two input families together, not a bound that is
charged twice.

## 4. Uniform constants, including the small group counts

Write `C=Dw`. The earlier elementary bounds still give
`A>=1024`, `D>=2^19`, `D^2<=n/2`, and `t/D<=19/2^19`. Moreover,

```math
\frac1{1-(A+1)/A^2}<\frac{1001}{1000},
\qquad
C\le\frac{1001}{2000}(n+1)+D<\frac{501}{1000}n.
\tag{13}
```

For the final inequality use `n>=2^40` and `D/n<=2^-20`. As in the
adaptive constructor,

```math
U\le36C+20Dt+16D^2,
\qquad S\le20C+10Dt+8D^2.                           \tag{14}
```

At `P=2`, `D^2/n<=1/4`, so (13)--(14) imply

```math
U<\frac{221}{10}n,\qquad S<\frac{121}{10}n.
\tag{15}
```

For `P>=3`, use `D^2/n<19/50`. At `P=3`, this follows from the exact
integer inequality `50*2^62<19*3^40`; at `P>=4`, it follows from
`D^2/n<64/P^4<=1/4`. Thus

```math
U<25n,\qquad S<14n\qquad(P\ge3).                    \tag{16}
```

The same generic grouping theorem supplies a coloring with pair width
at most `floor((P U+GQ)/P^3)` and group size at most `floor(2G/P)`.
Set

```math
\lambda=\max\left(1,\left\lfloor\frac{P U+GQ}{P^3}\right\rfloor\right),
\quad H=W\lfloor2G/P\rfloor,\quad B_{\rm tab}=H+1.
\tag{17}
```

The unchanged estimates from the adaptive constructor give

```math
P^2\lambda\le P^2+U\left(1+\frac1P+\frac{2627}{1000P^2}\right),
\qquad
PH\le2S\left(1+\frac{2627}{1000P}\right).
\tag{18}
```

At `P=2`, substitute (15), `P^2<=n`, and the exact inequalities

```math
1+\frac{221}{10}\left(1+\frac12+\frac{2627}{4000}\right)<49,
\qquad
\frac{242}{10}\left(1+\frac{2627}{2000}\right)<56.
\tag{19}
```

For `P>=3`, the factors in (18) decrease with `P`; (16) gives

```math
1+25\left(1+\frac13+\frac{2627}{9000}\right)<49,
\qquad 28\left(1+\frac{2627}{3000}\right)<56.
\tag{20}
```

Consequently `P^2 lambda<=49n` and `PH<=56n` in all cases. Since
`P<=n`, adding one to `H` proves `B_tab<=57n/P`.

As `P` grows, `(A+1)/A^2`, `D/n`, and `D^2/n` tend to zero uniformly
under the cutoff. Equations (3) and (13) give `C/n -> 1/2` and
`f/w -> 2`. Therefore

```math
U=(18+o(1))n,\qquad S=(10+o(1))n.
\tag{21}
```

Substitution in (18), and the same negligible rounding bounds as before,
proves (2).

## 5. Written instances and the full factoring cost

For each fixed `ell`, construct the tables by enumerating original
local settings, enforcing all internal equalities, matching the
external message blocks, and requiring padding bits to be zero.
Every table entry is a sum of at most `2^H` zero-or-one contributions.
The identical off-diagonal codec and empty-group convention give
the Appendix B predicate `WFe` and preserve the exact count.

All arithmetic parameters, modular operations, carries, prefix tests,
and highest-digit bounds have polynomial bit cost in `n`. Coloring
search and explicit table construction cost at most

```math
n^{O(1)}\left(P^G+P\,2^{(P-1)\lambda+H}+1\right)
\le n^{O(1)}\left(P^G+P\,2^{105n/P}+1\right)
=2^{O(n/P)}.                                             \tag{22}
```

This includes writing the complete binary integer tables and the unary
channel width. The inequality `G log2 P<=3P^9<=3n/P` pays for the
deterministic coloring search. No modular arithmetic is treated as
unit cost.

For any prefix, sum the exact `n` length-conditioned clique values:

```math
\#\{d:d\mid N,\lfloor d/2^k\rfloor=x\}
=\sum_{\ell=1}^n\mathrm{zval}(P,I_\ell).
\tag{23}
```

The summands correspond to disjoint bit-length classes. Integers in
(23) have `O(n)` bits, so summation has polynomial bit cost. There is
no need to encode a selector in the clique channels.

To obtain the required smallest prime factor, subtract the contribution
of divisor one from each queried prefix count. Start at the empty
prefix. At each of the `n` steps, test the child prefix with next bit
zero; choose it if its adjusted count is positive, and otherwise choose
the child with next bit one. The current prefix always contains a
divisor at least two, since `N` itself is one. The final integer is
therefore the smallest divisor at least two. It is prime: a composite
minimal divisor would have a smaller divisor at least two. For prime
`N`, this procedure returns `N`.

There are at most `n^2` length-conditioned evaluator calls. This fixed
polynomial overhead, and (22), preserve the established half-rate
factoring bound with the same source level schedule and its ordinary
small-input fallback. The improvement is in the concrete construction
bounds, not in the evaluator's contraction parameter.

## 6. Exact finite controls

Run

```sh
python3 tools/check_operand_length_transform.py \
  --output evidence/operand-length-transform.json
```

The checker implements the full-input first transform stage, compares
it with independent integer convolution, verifies exact carries,
tests length-class counts and their sum, and recovers the smallest
prime factor including prime targets. It tests partial highest digits
and inputs with more than `D` occupied positions in one operand.
Separate large-integer parameter controls check all rounding bounds
and both `P=2` and `P>=3` constant arguments. The actual large clique
tables are not materialized by those parameter controls.

The source transform and grouping ingredients are attributed in the
linked adaptive constructor. The operand-length split, paired-block
bound, constants, and new checker were developed with OpenAI assistant
support. No Lean changes are included.
