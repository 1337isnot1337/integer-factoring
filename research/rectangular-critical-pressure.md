# Critical pressure, frozen cores, and transported rectangular starts

## Result and scope

The rectangular defect search can be judged directly at the target
\(R_*=\sqrt2\), without first optimizing its rate parameter. Define

~~~math
F_*(\theta)=\log B(\theta)-R_*\log(Z_X(\theta,R_*)Z_Y(\theta,R_*)).
                                                               \tag{1}
~~~

A negative value is exactly the strict pressure target
\(\rho<1/\sqrt2\). Comparing (1) with an explicit Gaussian control
at the same bulk gauges separates a potential defect benefit from
ordinary convergence of the interval family.

The two frozen numerical snapshots analyzed here have **positive**
critical pressure and positive excess over that Gaussian control.
Their three-state cores are almost concentrated on the low state.
An explicit penalty path to the all-low boundary has a Gaussian limit,
proved below. This motivates starts which populate the top core state,
instead of relying only on larger intervals near that boundary.

The accompanying tool computes the discriminator and seven transported
or core-activation starts for each snapshot. It makes no new profile
queries. Its numerical values are exploratory floating data, not exact
pressure certificates or a claim of an improved factoring exponent.

## 1. The physical partition and the Gaussian control

Use the [bounded-defect coordinates](rectangular-defect-factorization.md)
or their [full-core extension](rectangular-full-core-factorization.md).
Write \(\ell=\log s\). All supplied height parameters \(\theta_i\)
already include multiplication by \(\ell\).

The first four parameters are the two bulk charges at infinity and
the two at zero. Put

~~~math
c_\infty=\theta_0+\theta_1,\quad c_0=\theta_2+\theta_3,
\quad g_a=(4a-2M+2)\ell,
\quad g^X_a=g_a+c_\infty+c_0,
\quad g^Y_a=g_a+c_\infty-c_0.                            \tag{2}
~~~

The physical alphabets have four extreme states on \(X\), and three
on \(Y\). Let \(v=(v_l,v_t,v_b)\) be the sums of all extra core
weights on the matching \(X,V\) pair. Let \(u=(u_l,u_t)\) be the
sums on \(Y,U\). The four- and three-state core energies are

~~~math
\begin{aligned}
E_X&=(0,\ g^X_0+v_l,\ g^X_{M-1}+v_t,
            \ g^X_0+g^X_{M-1}+v_b),\\
E_Y&=(0,\ g^Y_0+u_l,\ g^Y_{M-1}+u_t).
\end{aligned}                                           \tag{3}
~~~

Every interior mode contributes its binary factor. Therefore

~~~math
\log(Z_XZ_Y)=
\sum_{a=1}^{M-2}\left[
 \log(1+e^{-g^X_a/R})+\log(1+e^{-g^Y_a/R})\right]
+\log\sum_j e^{-E_{X,j}/R}+\log\sum_j e^{-E_{Y,j}/R}.       \tag{4}
~~~

The comparison Gaussian uses the same \(M,\ell,c_\infty,c_0\),
the full binary alphabet at every mode, and no extra core weights.
Its complete profiles are explicit:

~~~math
\begin{aligned}
\log B_G={}&
\sum_{a=0}^{M-1}\sum_{\sigma=\pm1}
 \log(1+e^{-[(2a+1)\ell+\sigma c_0]})\\
&+\sum_{j=0}^{2M-1}
 \log(1+e^{-[(2j+1-2M)\ell+c_\infty]}),\\
\log(Z^G_XZ^G_Y)={}&
\sum_{a=0}^{M-1}\left[
 \log(1+e^{-g^X_a/R})+\log(1+e^{-g^Y_a/R})\right].
\end{aligned}                                           \tag{5}
~~~

These profiles follow from the exact diagonal Gaussian factors in the
bounded-defect note, including the opposite zero-place charge signs.
Define \(F_G\) by (1) with (5), and define the total defect excess

~~~math
\Delta_*=F_*-F_G.                                       \tag{6}
~~~

The core gauges and physical alphabets differ between the defect and
this control. Thus (6) is a comparison of complete pressure expressions,
not a claim to isolate one local interaction or a universal lower bound.
The actual score test is \(F_*<0\), not merely \(\Delta_*<0\).

## 2. Why the critical value is a useful discriminator

For fixed energies, let \(Z(R)=Z_X(R)Z_Y(R)\). The Gibbs identity gives

~~~math
\frac{d}{dR}\bigl[R\log Z(R)\bigr]
=H_X(R)+H_Y(R).                                         \tag{7}
~~~

This entropy is positive whenever at least one physical distribution
has more than one supported state. Hence \(F(R)=\log B-R\log Z(R)\)
is strictly decreasing. If its root is \(R_\theta=2\rho_\theta\),
then

~~~math
F(R_*)=\int_{R_*}^{R_\theta}(H_X(R)+H_Y(R))\,dR.          \tag{8}
~~~

This proves the equivalence of the sign test and the rate comparison.
It also gives the local estimate
\(\rho_\theta-R_*/2\simeq F(R_*)/[2(H_X+H_Y)(R_*)]\).
The tool checks the saved root independently before reporting this
linearized estimate.

For the near-critical interval gauges, the known entropy scale is
\(\Theta(M/\log M)\), and the Gaussian rate gap is
\(\Theta((\log M)^2/M^2)\); see the
[interval asymptotic](quantitative-interval-frontier.md).
Equation (8) converts this to Gaussian critical pressure of order
\(\log M/M\). A negative defect excess of constant size would dominate
that positive remainder and give a rate gain of order \(\log M/M\).
A smaller negative excess can also suffice if it dominates the Gaussian
critical pressure. The two finite observations below do not establish
any limiting sign or decay law.

## 3. What the saved starts show

The [frozen input records](../evidence/rectangular-critical-pressure-inputs.json)
preserve the exact floating arrays used for this calculation. Live search
files may later contain different best points.

| Modes | Saved \(\rho\) | \(F_*\) | Gaussian control \(F_G\) | Excess \(\Delta_*\) |
|---:|---:|---:|---:|---:|
| 8 | 0.7202236588 | 0.1183920638 | 0.0962846737 | +0.0221073901 |
| 16 | 0.7115826971 | 0.0693110041 | 0.0597644975 | +0.0095465065 |

Both excesses are positive. The smaller value at 16 modes is not evidence
for a negative limiting defect. At the saved roots, the \(Y\) low-state
probability is approximately \(0.9998\), and the top-state probability
is below \(10^{-7}\). The useful search question is whether a region
with substantial top-state weight can change the sign of (6), while
retaining the Gaussian bulk benefit.

### An explicit frozen-boundary limit

The mixed partial repair fixes the all-low Gaussian background exactly,
for every polynomial parameter and every bulk entry; this is proved in
[the mixed repair note](rectangular-f4-repair.md). The input low state is
unchanged by its code embedding. The selected all-low output row is
unchanged by the repair. The quartic contribution creates a top mode
and is killed by the all-low output projection. Therefore fixing the
low extreme state on all four physical legs leaves exactly a Gaussian
coordinate restriction on the bulk.

The full-core gauges admit a direct limit to this restriction. At both
valuation places, give every non-low extreme state a penalty \(L\)
relative to the low state, and let \(L\to\infty\). Constants on a
physical leg can be subtracted freely, so these penalties fit the
normalization in which the empty-state weight is zero.

After removing the common low-state monomial, a local matrix minor using
only retained rows and columns has valuation independent of \(L\).
Every other minor has a nonnegative integral multiple of \(L\) added
to its valuation. For sizes up to the rank of the retained matrix,
the minimum minor valuations eventually equal those of the retained
matrix. Above that rank they diverge. Thus the finite Smith orders
converge to the orders of the all-low restriction, and all other profile
terms vanish. The physical partition similarly keeps only its low core
letter. The pressure limit is the Gaussian restriction's pressure.

This is a fixed-\(M\), explicit penalty-path result. It is not a lower
bound for all Gaussian gauges, nor does it classify every path along
which the physical type becomes concentrated. It shows why the all-low
boundary itself does not retain the non-Gaussian interaction. A strict
advance through this defect family must use something beyond that
unchanged background.

## 4. Transport and core-activation formulas

Let a source point have \(M\) modes and logarithmic parameter \(\ell\).
To move to \(N\) modes while preserving both extreme base energies,
set

~~~math
\ell'=\ell\frac{M-1}{N-1}                               \tag{9}
~~~

and keep all other raw parameters unchanged. Equations (2)--(3) show
that every physical core energy is unchanged. This differs from keeping
the normalized heights \(\theta_i/\ell\) fixed.

A second start uses the zero-charge Gaussian critical minimizer
\(\ell'_G\), found by a one-dimensional golden-section search of (5).
Put \(d=(N-1)\ell'_G-(M-1)\ell\). The new low base energy decreases
by \(2d\), and the top increases by \(2d\). Compensate the summed
core weights by

~~~math
\Delta v=(2d,-2d,0),\qquad \Delta u=(2d,-2d).             \tag{10}
~~~

Distribute each change equally among the corresponding physical legs
and valuation places. This preserves all core energies while putting
the bulk at the new Gaussian critical scale.

Finally, choose a desired core distribution \(q\), with positive empty
state probability \(q_0\). Its energies relative to the empty state
at \(R_*\) are

~~~math
E_j^{\mathrm{target}}=-R_*\log(q_j/q_0).                  \tag{11}
~~~

Apply the difference from (3) through the same summed core weights.
The tool checks that its actual Gibbs probabilities recover the chosen
targets. Besides the two transports, it emits starts with a populated
\(Y\) vacuum, equal low/top \(Y\) populations, uniform \(Y\) core,
simultaneously populated top states, and uniform cores on both pairs.
These are attempts to enter different search regions, not predicted
improvements.

For the 16-to-32 transport, the starts activating both cores can require
normalized heights as large as about 45. A search bound of 48 admits
them. The other starts fit within 24. This is an optimizer-domain choice;
the mathematical family permits the stated finite weights.

## 5. Reproduction

```sh
python3 tools/analyze_rectangular_critical_pressure.py \
  --starts-directory /tmp/rectangular-critical-starts
```

The default inputs and output are the two evidence files linked above
and [rectangular-critical-pressure.json](../evidence/rectangular-critical-pressure.json).
The output includes all seed vectors, their physical Gibbs distributions,
input and tool hashes, and the maximum normalized height of each start.
The separate seed files can be passed directly to the existing search:

```sh
python3 tools/search_rectangular_core_pressure.py --M 32 \
  --height-bound 48 \
  --starts /tmp/rectangular-critical-starts/defect-M32-transported-starts.json \
  --output /tmp/rectangular-critical-search.json
```

Only the first command was run for this evidence. It verifies the saved
root equations, the energy-preserving transports, and the target Gibbs
probabilities. It does not call the native profile oracle. Any candidate
from the second command still needs complete exact profile replay and
a strict rational pressure comparison.
