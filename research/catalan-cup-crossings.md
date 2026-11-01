# Full Catalan cup codes with coordinate duals

The crossing endpoint needs an injective primal frame and a left inverse;
it does not require inversion of the frame's self-pairing. Keeping that
distinction retains all noncrossing cup tensors at the fixed field parameters
where their Gram matrix degenerates. The resulting ordered braid family has
Catalan-many labels, strict endpoints, and an effective construction.

At q=1 it is defined over F2 even though the entire primal Gram matrix is
zero. At the previous order-five q in F16, it retains channels discarded
by the [height-three quotient](quantum-braid-path-quotient.md). These are
new paired encoders and crossing candidates for this repository. A favorable
complete two-cut rank sum and a stronger factoring score remain open.

## 1. A unitriangular physical minor

Let n=2r. For each noncrossing perfect matching P of n ordered positions,
put one cup on every pair. In the binary physical basis the cup has
coefficient1 on (1,0), coefficient q on (0,1), and zero otherwise.
Let F_P be the product of these cup coefficients. No quantum integer
denominator is used. There are L=C_r such tensors, where

```math
C_r=\frac1{r+1}{2r\choose r}.
```

Let w_P put a1 at each opener of P and a0 at each closer. This is a
Dyck word, with1 interpreted as an up step. Orienting any cup the other
way lowers the word's prefix height by two between its endpoints.
Consequently

```math
F_Q(w_P)\ne0\ \Longrightarrow\ h_P(j)\le h_Q(j)\quad(0\le j\le n),
\qquad F_P(w_P)=1.                                  \tag{1}
```

Order the Dyck paths by increasing total height, breaking ties by their
words. The square matrix T=(F_Q(w_P)) is upper triangular with diagonal
one: equal total height together with dominance forces equality. Thus
the full physical frame F is injective, including at a root of unity
or in characteristic two.

Let S select precisely the Dyck-word physical coordinates, in this order,
and define the decoder

```math
D=T^{-1}S.\qquad DF=I_L.                            \tag{2}
```

The inverse is obtained by triangular substitution. It is not a Gram
inverse or an invariant transpose dual. At q=1 every closed cup pairing
has a factor delta=q+q^-1=0, so F^T F=0 for n>0; (2) still holds.
This supplies a particularly direct nonorthogonal example over F2.

## 2. Strict crossing endpoints without a negligible quotient

Use the [ambient cup and braid conventions](quantum-block-braid-crossings.md):
e=q^-1 cc^T, B=qI+e, with q a nonzero element of F16, including q=1.
The elementary cup slide is

```math
B_0B_1(c\otimes I)=q(I\otimes c).
```

Sliding the n-by-n block braid past all n/2 cups on its first block gives

```math
B_{n,n}(F_P\otimes F_Q)
 =q^{n^2/2}(F_Q\otimes F_P).                         \tag{3}
```

This uses only cups, braid relations and planar composition. It is valid
on the entire cup frame; no self-pairing norm, height cutoff, f4 quotient,
or quantum-group classification is needed.

Let gamma=q^(-n^2/2) and define

```math
H(t)=\gamma(D\otimes D)
         \prod_{i\text{ in block word}}(qI+t e_i)(F\otimes F),
\qquad
G(t)=(\gamma+(1+\gamma)t)H(t).                       \tag{4}
```

Products act in the stated application order. Equation (2) gives
H(0)=gamma^-1 I, and (3) gives H(1)=E. Therefore G(0)=I and G(1)=E,
with E(P,Q,R,S)=[P=S][Q=R]. The degree is at most n^2+1. At q=1,
gamma=1 and B is the ordinary binary swap, so the family lies in F2[t]
and has degree at most n^2.

## 3. Exact link-state compilation and the decoder boundary

The noncrossing cups on 2n positions form another injective physical
frame. Each e_i acts on its matching basis by an exact local surgery:
if i and i+1 are paired, multiply that matching by delta; otherwise
pair i with i+1 and reconnect their former partners, with coefficient1.
This follows directly from the cup/cap identities. Different source
matchings can have the same target, so the implementation explicitly
accumulates their contributions.

Concatenated matchings P,Q inject the input frame into this global
basis. The output map is

```math
\mathcal B=(D\otimes D)F_{2n}.
```

It is computed by evaluating the global cups on concatenated Dyck words
and applying T^-1 tensor T^-1. Its columns at concatenated input matchings
are I. Its other columns generally do not vanish. Replacing this map
by selection of middle-height-zero paths would change the construction.

There is an explicit warning at the old root q=8: on eight physical
strands the fourteen cup tensors are independent, their Gram has rank13,
and f4 on the first four strands has a nonzero rank-one image in the
Gram radical. The four-strand product decoder sees that image. Thus
the earlier negligible-wall boundary argument cannot discard f4 here.
In particular a quotient-only [spectral transport](fused-spectral-transport.md)
formula must have its ambient residual restored before use with this decoder.

## 4. Binary alphabets and effective construction

Choose the least positive even n with C_(n/2)>=2^lambda and restrict F,D
to one common subset of exactly 2^lambda labels. Equation (2) and both
crossing endpoints survive that restriction. Consecutive Catalan numbers
have ratio less than four, so the full L is less than4*2^lambda.
Their standard central-binomial estimates give

```math
n=\lambda+O(\log(\lambda+2)).                       \tag{5}
```

The primal entries are explicit products. Triangular inversion, global
link-state compilation and dense coefficient/rank-factor construction
use 2^O(lambda) work: 2^n, L and C_n are all 2^O(lambda), and there are
O(n^2) coefficients. This is within the construction allowance discussed
in the [effective modular interface](bounded-spin-modular-interface.md).
No polynomial-time entry bound is asserted for the dense decoder.
The q=1 family needs no extension-field descent.

The number of labels is asymptotically larger than the previous fixed
height-three Fibonacci frame. This does not by itself improve the
crossing rank budget. The complete same-piece square-root rank sum
still needs a uniform estimate below the half-rate threshold.

## 5. Evidence and next question

The [compiler](../tools/catalan_cup_crossing.py) includes the primal minor,
coordinate decoder and full link-state contraction. The
[checker](../tools/check_catalan_cup_crossing.py) compares physical cup
actions and independent physical-word polynomial evolution with the
compiled coefficients, checks strict endpoints, and records both actual
cut ranks. The [receipt](../evidence/catalan-cup-crossings.json) pins sources.

```sh
make verify-catalan-cup-crossing
```

The positive next task is a complete graded decomposition exploiting
the retained radical and freedom in the decoder. Neither projecting away
that radical nor importing a quotient-only identity is justified for this
family. No stronger complete factoring score or new Lean proof is claimed.
