# Small reduced norms retain density and force local event separation

## Result and scope

Truncating the smooth leading norm slightly below the full height preserves
the proved covered-subgroup density scale and gives an exact local
correlation rule. For an unramified split target prime `p`, choose integers
`T>=1,R>=1` with

```math
R\le\sqrt D/4,\qquad 4p^T R^2<D,\qquad D=|\Delta|.
                                                               \tag{1}
```

Count a row-change event for the **truncated** predicate: a reduced leading
norm is `Y`-smooth, coprime to `Delta`, and at most `R`. Two such events
cannot occur within `T` steps. Consequently their short-range connected
correlations can be computed exactly; they are not independent.

If `log(sqrt(D)/R)=O(b/(log b)^alpha)`, `b=bits(D)`, the existing
[even-product proof](small-index-smooth-density.md) retains density
`exp(-O_alpha(b/(log b)^alpha))` on every subgroup of index at most
`D^kappa`, fixed `kappa<1/4`, at cutoff
`Y=ceil(exp((log b)^(1+alpha)))`. The truncated row-change events inherit
that scale. This density deduction uses the same principal split-prime
estimate as the original theorem. The ideal-norm separation itself is
elementary under the stated unramified-word contract.

This supplies local arithmetic content for the
[connected-correlation interface](connected-event-window-criterion.md),
not the required long-range estimate. Actual target indices and affordable
uniform hit windows remain open. No complete factoring-score change:
the strongest bound remains **1/2**.

The question and stopping test were
[registered in Discussion 51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19355801).
There is one fixed saved-data control, with `T=2`, not a new orbit or
cutoff search. The raw factoring selector still accepts every original
success; truncation defines an auxiliary flag for the proof.

## 1. A nonzero principal word has norm at least D/4

Use the weighted relation-lattice gap already recorded in
[orda-aska's graph/lattice note, Lemma 1](leading-power-window-graphs.md#1-the-graph-and-the-lattice),
source `a427442`, with the later unramified-base scope repairs. The argument
is reproduced here to make the hypotheses of the event deduction explicit.

Fix one orientation `P_q` for each unramified split rational prime `q`
in the factor base. Let `w` be a nonzero integral vector satisfying
`prod_q [P_q]^(w_q)=1` in the proper ideal class group. Then

```math
\boxed{\prod_q q^{|w_q|}\ge D/4.}                        \tag{2}
```

To prove it, form the integral ideal

```math
I_w=\prod_{w_q>0}P_q^{w_q}\prod_{w_q<0}\overline P_q^{-w_q}.
```

All its primes are coprime to the order discriminant. At each rational
prime it contains only one of the two split orientations, so it is
primitive: no rational integer greater than one divides it as an ideal.
It is proper, invertible and principal, say `(beta)`, with `beta` in the
order and norm `M=prod_q q^|w_q|>1`. If `beta` were rational, it would be
an integer and `(beta)` would not be primitive. Thus write

```math
\beta=(u+v\sqrt\Delta)/2,\quad u,v\in\mathbb Z,\quad v\ne0.
```

Then `M=(u^2+Dv^2)/4>=D/4`. This also applies to nonmaximal orders:
the ideals used are away from the conductor, and the displayed expression
holds for every element of the order. No factorization of the discriminant
or algorithm for computing a principal generator is assumed.

The unramified restriction is real. At `Delta=-20` the ramified class
of `(2,2,3)` has square one. Its nonzero exponent vector `w_2=2` has
norm four, with `4M=16<20`; the ideal square is the rational ideal `(2)`.
Known ramified-square relations must be removed or handled separately
before using (2). Inert primes also do not supply this oriented-word model.

## 2. Distinct decoded rows obey an exact weighted-distance bound

At two smooth phases `k,k+t`, choose the usual canonical oriented
words `v(k),v(k+t)`. They satisfy

```math
[g_p]^k=\prod_q[g_q]^{v_q(k)},\quad
\prod_q q^{|v_q(k)|}=A_k,\quad r(k)=k e_p-v(k).
```

The difference `w=r(k+t)-r(k)` is a principal word. Thus either the
integer rows are equal, or

```math
4\prod_q q^{|t\mathbf1_{q=p}-v_q(k+t)+v_q(k)|}\ge D.
                                                               \tag{3}
```

In particular, equality of the rows is forced by
`4p^t A_k A_(k+t)<D`, since the product in (3) is at most
`p^t A_k A_(k+t)`. This is a statement about a true ideal relation,
not an assertion that distinct rows are rationally independent.
For a boundary-crossing pair use the integer difference `t` and periodic
words; do not compare the row labels modulo the period.

## 3. The small-norm row-change events are T-separated

Let `f_R(k)` mark the truncated predicate in the opening paragraph.
Use the periodic event definition from the original row-change lemma:
`e_R(k)=1` precisely when `f_R(k)=1` and either `f_R(k-1)=0` or its
integer row changes. Equivalently, a same-row continuation has
`v(k)-v(k-1)=e_p`. This handles the periodic boundary despite row drift.

Suppose two truncated hits occur at `k,k+t`, `1<=t<=T`. Equation (1)
and (3) force their rows to agree, hence
`v(k+t)=v(k)+t e_p`. For every intermediate integer `0<=j<=t`, the
word `v(k)+j e_p` represents `[g_p]^(k+j)`. Its norm is at most `R`:
only the `p` coordinate changes, and the absolute value of that affine
coordinate is bounded by the larger endpoint absolute value.

The corresponding primitive ideal is already reduced. Indeed,
its leading norm `a<=R<=sqrt(D)/4`; after centering its middle coefficient
`c`, the trailing coefficient `(c^2+D)/(4a)` is greater than `a`.
CRT orientation and uniqueness of the reduced representative give exactly
the same intermediate canonical word. All intermediate phases are
therefore truncated hits with the same row. The second endpoint is not
an event. We have proved

```math
\boxed{e_R(k)e_R(k+t)=0\quad(1\le t\le T).}             \tag{4}
```

This argument uses unramified split words throughout. Proper-gcd successes
can terminate the actual search; no unknown decoded word is assigned to
them in this auxiliary flag.

Write `delta_R` for truncated hit density and `lambda_R` for its event
density on the target cycle. The existing fixed-row argument, now with
`K_p(R)=min{k:p^k>R}`, still gives

```math
\lambda_R\ge\delta_R/(2K_p(R)-1)>0.                      \tag{5}
```

For distinct time indices `j_1,...,j_r` in an interval of diameter at most
`T`, every joint product involving two or more flags is zero. In the
partition definition of the joint cumulant, only the singleton partition
survives. Thus the exact connected correlation is

```math
\kappa(e_R(x+j_1),\ldots,e_R(x+j_r))
   =(-1)^{r-1}(r-1)!\lambda_R^r.                         \tag{6}
```

For example covariance is `-lambda_R^2`, while the third cumulant is
`2lambda_R^3`. Repeated indices need the full partition formula; (6)
is explicitly for distinct indices.

## 4. Density survives the reduced norm budget

Set `E=log(sqrt(D)/R)`, and assume `E=o(b)` and `R<=sqrt(D)/4`.
Repeat the original even-color proof with the largest even integer
`L<=log R/log Y`. Its concentration argument, common character means,
and distinct-product injection are unchanged. The one-sided negative
Fourier error is now

```math
\log\gamma=-\tfrac12\log R
       +O\left(\log Y+\frac b{\log Y}\log b\right)
 =-\tfrac14\log D+\tfrac12E+o(b).
```

Every fixed subgroup-index range `d<=D^kappa`, `kappa<1/4`, is therefore
still covered. Its normalized count, using
`|G|=O(sqrt(D)log D)`, gives

```math
\frac{|S_{Y,R}\cap U|}{|U|}
 \ge\exp\left(-E-O\left(\log Y+\frac b{\log Y}\log b\right)\right).
                                                               \tag{7}
```

Here `S_(Y,R)` consists of the original recognizable smooth classes with
leading coefficient at most `R`, coprime to `Delta`. It is not a
smooth-pair quotient or an unreduced norm witness. The floor/even-length
loss is `O(log Y)` as in the source proof. The class-number and analytic
contracts, including their effective sufficiently-large thresholds, are
inherited from that note.

For fixed `alpha>0`, use `Y=ceil(exp((log b)^(1+alpha)))`. If
`T log p=O(b/(log b)^alpha)`, choose, for example,

```math
R=\min\left(\lfloor\sqrt D/4\rfloor,
       \left\lfloor\sqrt{\frac{D-1}{4p^T}}\right\rfloor\right).
```

Eventually `R>=1`, (1) holds, and
`E=O(1+T log p)=O(b/(log b)^alpha)`. Equations (5) and (7) give
`lambda_R>=exp(-O_alpha(b/(log b)^alpha))` for covered target subgroups.
For every `p<=Y`, one may uniformly take
`T=floor(b/(log b)^(1+2alpha))` for sufficiently large `b`, without
changing this density scale. This specifies a polynomial separation
length; it is not a subexponential empty-window bound.

In this regime `p^T` has `O(b)` bits, so the displayed integer radius is
computable by powering and integer square root at polynomial bit cost.
The raw selector need not compute the auxiliary radius or event process;
it still stops at any original accepted query. No event-count oracle or
uncomputed principal generator is used by the factoring algorithm.

## 5. Verification and remaining analytic obligation

```sh
make verify-small-norm-event-separation
```

The [checker](../tools/check_small_norm_event_separation.py) uses the saved
complete smooth words at the existing 198850-period generator-5/7 cycles.
It checks (3) for every nonzero row difference, using capped integer
products to avoid constructing huge powers. The single predeclared
`T=2` truncation checks periodic event separation, the row-density loss,
and the exact covariance/third-cumulant formulas. The actual ramified
square above guards against dropping the hypothesis. The
[receipt](../evidence/small-norm-event-separation.json) pins the upstream
row data semantically; the existing group/row/event checker is a prerequisite.
There is no new arithmetic orbit or parameter selection.

Separation alone cannot bound long empty gaps: equally dense separated
sets can cluster into only part of a cycle. For the
[finite-order criterion](connected-event-window-criterion.md), tuples
whose diameter exceeds `T`, and repeated-index terms, still require an
arithmetic signed-sum or tree estimate. Neither (4) nor (6) supplies it.
In addition, the actual target indices remain unknown. These are the next
substantial obligations; optimizing this polynomial separation length
alone cannot justify a complete factoring improvement.

Credit: orda-aska's graph/lattice Lemma 1, the primitive
short-product/injectivity argument in PR87, PR82's
row-change construction and density bridge, and the saved selector/Grok
data. The principal-ideal norm gap is classical; this note records its
specific event-correlation consequence and density tradeoff. Substantial
OpenAI Codex assistance; no new Lean or independent mathematical review.

The [clean frozen default replay](../evidence/small-norm-event-full-verification.json)
passes at exact source `b36cf5799aa4537b75e988ab2fc907ef4fe3620e` in
1216.39 seconds. Its [complete compressed log](../evidence/small-norm-event-full-verification.stdout.txt.xz)
is archived. The replay used a private `/tmp` namespace with read-only host
root and only its frozen checkout writable, avoiding concurrent workers'
fixed temporary outputs. Main's later unary inverse integration preserves
this checker byte for byte and receives separate focused/docs checks.
