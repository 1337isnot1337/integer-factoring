# A floor from the first fifty q8 orders

## Result and reason to stop higher-jet refinements

The first fifty protected orders already exclude the balanced two-axis
target for the retained q8 degree and paired-type invoice. No returning
layer or new tensor-profile computation is needed.

**Protected-prefix theorem.** Suppose each cut has 4096 finite nonnegative
Smith orders, total order at most `56*4096`, and distinct selected orders
consisting of one zero, fourteen orders at most six and thirty-five
further orders at most twelve. With the same common-scale, unweighted
64-letter paired-type transfer and word degree `32m E_nu` as the
[fixed curve](fixed-curve-quartic-floor.md), its profile-envelope pressure
is at least `5/7>1/sqrt(2)` for every common type with positive entropy.
The selected orders need not be the entire low prefix or have those exact
values. The two cut profiles need not be identical.

The [spatial-mode construction](spatial-mode-colour-fusion.md) satisfies
these premises at q8. Therefore changing only its terms after order
twelve, retaining the endpoint and degree, cannot cross the target by
this profile criterion. The proposed higher-jet matching refinement was
stopped before constructing or searching a tensor. This theorem is not
an exclusion of other alphabets, changed early structure or degree costs,
weighted gauges, unequal scales, or stronger coefficient-rank estimates.
The complete proven factoring exponent remains **1/2**.

The [scope and coverage question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19364676)
were registered before the exact certificate. A preceding scalar numerical
diagnostic motivated the proof; no numerical minimum is used below.

## 1. A conservative dimension-and-degree tail

Let `Z_i(t)=sum_j t^(v_ij)`, with `0<t<1`, be either cut partition.
The selected orders contribute at least `1+14t^6+35t^12`. There are
4046 remaining orders. They are nonnegative and their sum is at most the
**entire** determinant-order cap, even if the selected orders are smaller
than their stated upper bounds. Thus their mean is at most

```math
\frac{56\cdot4096}{4046}=\frac{114688}{2023}<57.
```

Convexity of `v -> t^v` and monotonicity in v give, independently for
both cuts,

```math
Z_i(t)\ge 1+14t^6+35t^{12}+4046t^{57}=:L(t).          \tag{1}
```

In particular, this does not incorrectly subtract the maximal possible
selected-order sum when only upper bounds on those orders are known.

The cap and word degree retain their original physical costs. There are
four q8 positions on each cut face. Their average occupied count is
`4*7/8`, so a homogeneous degree eight per occupied old vertex gives
full determinant degree `8*(2*4*7/8)*4096=56*4096`.
Four parent legs, each pairing two old vertices, give word degree
`8*(4m E_nu)=32m E_nu`. These are the fixed curve's costs, not a reused
three-axis normalization.

## 2. One exact positive polynomial

Put `t=u^5`. It suffices to prove

```math
\bigl(1+14u^{30}+35u^{60}+4046u^{285}\bigr)^7
 -(1+7u^{56})^{20}>0,\qquad0<u<1.                    \tag{2}
```

The difference has a factor `u^30`. Its quotient P has degree 1965,
82 nonzero integer monomials and constant coefficient 98. All **4097**
coefficients of P in the degree-4096 Bernstein basis are strictly positive.
This proves (2) uniformly on the interval, without sampling.

Here is the exact certificate calculation. If `P=sum_j a_j u^j`, let

```math
M=\operatorname{lcm}_{a_j\ne0}\binom{4096}{j},\qquad
B_k=\sum_{j\le k}a_j\frac{M}{\binom{4096}{j}}\binom{k}{j}.
```

Every B_k is an integer and is positive. The checker independently
reconstructs all 4097 monomial coefficients, including the zeros above
degree 1965, through

```math
\binom{4096}{j}\Delta^j B_0=M a_j.
```

The [receipt](../evidence/q8-leading-prefix-floor.json) retains the sparse
polynomial, exact positive minimum, its index and the complete integer
coefficient-list hash. The minimum occurs at index 1689; the hash is
`0d2cefcc6e2ff096bcd110594b9aba6c95a8d0527495ef3b714ac8525883a135`.
No floating-point arithmetic is used. Taking positive seventh roots of
(2), then using (1), gives

```math
Z_i(t)>(1+7t^{56/5})^{20/7}.                         \tag{3}
```

## 3. Arbitrary paired physical types

Let nu be any common type on the 64 paired letters, with positive entropy
H(nu) and mean occupied count E_nu in [0,2]. Its two components can be
correlated. The ordinary entropy inequality for the paired-letter
partition function gives

```math
H(\nu)\le 2\log(1+7t^{56/5})-\frac{56}{5}E_\nu\log t.
```

For the charged two-cut envelope,

```math
R_\nu(t)=\frac{\log Z_0(t)+\log Z_1(t)-32E_\nu\log t}{4H(\nu)}.
```

Apply (3) to both partitions and multiply the entropy inequality by
20/7. The energy coefficient is `(20/7)*(56/5)=32`, exactly the word
degree above, so `R_nu(t)>5/7`. Taking the relevant infima preserves
the floor at least 5/7. Finally `(5/7)^2=25/49>1/2`, excluding the
growth-base-two crossing target. This is a limitation of the specified
profile envelope, not a lower bound for actual coefficient ranks or
factoring algorithms.

## Evidence, credit and the next changed assumption

```sh
make verify-q8-leading-prefix-floor
```

The [checker](../tools/check_q8_leading_prefix_floor.py) reproduces the
integer polynomial, all signs and all monomial reconstructions. It also
checks that the existing saved q8 profiles satisfy the stated premises;
it does not generate another tensor profile. The [full default receipt](../evidence/q8-leading-prefix-full-verification.json)
and [complete compressed log](../evidence/q8-leading-prefix-full-verification.stdout.txt.xz)
record exit-zero `make verify` in 1449.59 seconds at clean frozen
`e1256cccd50ac21837cf072b5568ade3f3ade358`, tree
`3fc59b478f9db2f606ae9f7b54fc728fe03c2f81`. A fresh marker probe
verified private /tmp, read-only host root, only the frozen checkout
writable and network isolation; tracked state stayed clean. The log
includes the new checker. Merge `ccd12e7` has the same tree. Later
archival prose does not claim a second replay. No new Lean or independent
mathematical review is claimed.

The leading-prefix premise comes from Sol worker 1's PR114, with its
Gaussian/colour/PR97 construction credit retained. The fixed curve and
doubled degree invoice come from Sol worker 2's PR111. Jensen's inequality,
entropy duality and Bernstein positivity are classical. This
protected-prefix composition and its exact scalar certificate were
developed by Sol worker 2 with substantial OpenAI Codex assistance.

Further work must change a premise that can affect the target: the
protected early channels, the physical degree/entropy invoice, or the
coefficient-rank method. More late-jet cancellation on this retained q8
core is stopped. The theorem does not decide the growing-alphabet
spatial-mode family.
