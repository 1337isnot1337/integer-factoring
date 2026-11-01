# Global factor moments from a fixed-scale Bezout jet

## Result and computational boundary

For arbitrary odd `N`, the complete regularized Bezout pair on two
ordinary odd intervals admits a short, known normalization. No optimal
primitive exponent or explicit reciprocal matching permutation is needed
to state the target. After normalization, every nonfactor row vanishes
at the requested precision. The remaining coefficients encode a measure
supported exactly on the actual factor pairs, with unit weights.

This removes polynomial degree and unknown primitive normalization from
the global coefficient-query contract. A polynomial bit-time
postprocessor converts a **supplied globally correct jet** into ordinary
factor moments. Computing that supplied jet for arbitrary `N` remains
open. The code does not replace that task with a fast provider, and no
factoring-score improvement is claimed.

The [companion extraction note](displaced-bezout-extraction.md) gives
the conditioned Hankel decoder and the sharper optimal-scale criterion.
The construction here supplies their fixed normalization and moment
contract, including independently sized rectangles.

Files: [postprocessor](../tools/regularized_bezout_factor_moments.py),
[finite checker](../tools/check_regularized_bezout_factor_moments.py),
and [evidence](../evidence/regularized-bezout-factor-moments.json).

## 1. Independent power-of-two domains

Let `r=2^u`, `s=2^v`, where `u,v>=0`, and let

\[
 A_i=a_0+2i\quad(0\le i<r),\qquad
 B_j=b_0+2j\quad(0\le j<s),                     \tag{1}
\]

with positive odd first nodes. Fix a positive odd integer `N`. Choose
`H>=max(u,v)+1` such that every nonzero `N-A_i B_j` has absolute
value below `2^H`. The extreme products give such an `H` in polynomial
binary-input cost. Define

\[
 E_A=2r-2-u,\qquad E_B=2s-2-v.                  \tag{2}
\]

These are the constant derivative orders on the respective grids:

\[
 \begin{aligned}
 v_2\prod_{i'\ne i}(A_i-A_{i'})&=E_A,\\
 v_2\prod_{j'\ne j}(B_j-B_{j'})&=E_B.
 \end{aligned}                                                   \tag{3}
\]

For example, the first product is
`2^(r-1)(-1)^(r-1-i)i!(r-1-i)!`.
The binomial `binom(r-1,i)` is odd because `r-1` has all low
binary digits equal to one. This proves the first identity in (3);
the second is identical. Empty derivative products have order zero.

Work in the unramified quadratic extension with `omega^2+omega+1=0`.
For numerical precision `p>=1`, take

\[
 K\ge H+p,\quad z=N+2^K\omega,\quad
 P(X)=\prod_i(X-A_i),\quad Q_z(X)=\prod_j(z-B_jX).
                                                               \tag{4}
\]

The polynomials in (4) have no common root. Let their unique rational
Bezout pair satisfy `U_0 P+V_0 Q_z=1`, `deg U_0<s`, `deg V_0<r`.
The fixed cleared pair is

\[
 \widehat U=\omega2^L U_0,\quad
 \widehat V=\omega2^L V_0,\qquad
 L=K+E_A+E_B.                                  \tag{5}
\]

The exponent `L` is stored in binary. Its numerical size is not a
working-precision requirement in the postprocessing theorem.

## 2. Every Lagrange summand is integral

For each odd `a=A_i`, exactly one `B_j` satisfies
`B_j=N/a mod2^(v+1)`, because the column interval is a full set of
odd residues at that modulus. For all other columns, the difference
valuation equals the valuation of `B_j-B_l`. Hence

\[
 v_2(Q_z(a))=E_B+
 \begin{cases}
 K,&ab=N\text{ for a column }b,\\
 v_2(N-aB_j),&\text{otherwise},
 \end{cases}                                                     \tag{6}
\]

where the second case has order below `H`.
The unique closest residue can be found by one modular inversion and
an index computation; the identity does not enumerate the columns.

Lagrange interpolation gives

\[
 \widehat V(X)=\sum_{a\in A}
 \frac{\omega2^L}{Q_z(a)P'(a)}
                  \frac{P(X)}{X-a}.             \tag{7}
\]

Equations (3) and (6) show that every scalar in (7) is integral.
If the row contains no actual factor pair, its scalar has order at
least `K-H+1>p` and its entire term vanishes modulo `2^p`.

If `ab=N`, write the odd derivative units as

\[
 \alpha_a=2^{-E_A}P'(a),\qquad
 \beta_b=2^{-E_B}\prod_{b'\ne b}(b-b').          \tag{8}
\]

Then, to precision `p`, the scalar in (7) is the real rational unit

\[
 c_{a,b}=\frac1{a^{s-1}\alpha_a\beta_b}.          \tag{9}
\]

Indeed the resonant factor of `Q_z(a)` is `2^K omega`, and
the remaining product is
`a^(s-1) prod_(b'!=b)(b-b')` times a unit congruent to one.
Every relative perturbation has order at least `K-v>=p`.
The leading `omega` in (5) cancels the pole's unit.

We have therefore proved the complete polynomial identity

\[
 \boxed{\widehat V(X)\equiv
    \sum_{\substack{a\in A,b\in B\\ab=N}}
               c_{a,b}\frac{P(X)}{X-a}\pmod{2^p}.}               \tag{10}
\]

Every weight in (10) is a unit. There is no generic specialization
assumption: the factor-supported terms are retained precisely at the
actual integer `N`.

## 3. Both complete polynomials have a short reverse support

Use `Y=X-1` and reverse `P,Q_z,V,U` at their prescribed degrees
`r,s,r-1,s-1`. Every root shift `a-1` in (7) is even. A coefficient
`d` places below the highest degree in `P(X)/(X-a)` therefore has
valuation at least `d`. Since every scalar in (7) is integral,

\[
 v_2([T^d]\widehat V^*)\ge d.                   \tag{11}
\]

Both input reverse products have the same degree-wise divisibility:
their roots `A_i` and `z/B_j` are odd. Their constant coefficients
are units, so their inverse unit series preserve that divisibility.
The reversed Bezout identity gives

\[
 \widehat U^*=-\widehat V^*Q_z^*/P^*
                              \pmod{T^s}.       \tag{12}
\]

The reversed scalar in (5) has degree `r+s-1` and is outside this
range, including singleton domains. Formula (12) and (11) prove
the same bound for `widehat U^*`.

Thus the **complete** pair modulo `2^p` has at most `min(r,p)`
and `min(s,p)` possible highest powers of `X-1`. This statement
holds for arbitrary reciprocal matching and arbitrary actual factor
incidence. The pair need not be primitive; no search for its optimal
scale is required. At the precision in (4), both complete outputs
are real modulo `2^p`.

## 4. A supplied jet gives factor-supported ordinary moments

Dividing (10) by the known reverse product gives

\[
 \frac{\widehat V^*(T)}{P^*(T)}
 \equiv\sum_{ab=N}\frac{c_{a,b}}{1-(a-1)T}
 =\sum_{j\ge0}\nu_jT^j\pmod{2^p},\qquad
 \nu_j=\sum_{ab=N}c_{a,b}(a-1)^j.               \tag{13}
\]

The ordinary moments are then

\[
 \mu_j=\sum_{k=0}^j\binom jk\nu_k
             =\sum_{ab=N}c_{a,b}a^j.            \tag{14}
\]

To output `q` moments, only `min(r,p,q)` coefficients of the global
`widehat V^*` jet are needed. The remaining depths from `p` onward
are zero by (11). The truncated product

\[
 P^*(T)=\prod_{i=0}^{r-1}(1-(a_0-1+2i)T)
\]

is a unit progression with coefficientwise-even increment. The
existing interval-jet algorithm computes the required prefix in
polynomial bit time in `log r`, endpoint lengths, `p`, and `q`.
Unit-series division and (14) use polynomially many bounded integer
operations. Every power-of-two guard required by the interval-jet
engine is charged by that engine. There is no numerical `L+p`
precision request.

The implementation rejects a missing required coefficient. Its
degree-wise divisibility check is only a necessary consistency check;
it does not establish that an arbitrary supplied jet is the correct
global one. In particular, zero input coefficients are not accepted
as a proof that no factor exists. For example, `N=15` and both
domains `{1,3,5,7}` have two factor pairs, but both complete
cleared polynomials are zero modulo two. The companion decoder
charges enough precision to separate this cancellation from an
empty factor set; a single low-precision zero test is insufficient.

This is the precise new provider contract: compute the indicated
short global jet at fixed, explicitly known normalization. The
[supplied-moment decoder](displaced-bezout-extraction.md) handles the
remaining conditioned reconstruction. That contract is weaker in
output and normalization requirements than constructing an expanded
global polynomial or determining its exact primitive content.

The [one-coefficient continuation](leading-bezout-factoring.md) further
chooses sufficiently long ordinary interval padding. On every actual
factor row, the derivative units then agree modulo `2^(2n+1)`, where
`n` is the original input bit length. One leading coefficient becomes
a common unit times a positive sum of proper divisors. This removes
the factor-count bound and all but one requested coefficient from
the factoring decision interface. Its global scalar provider remains
the unresolved step.

## 5. Ordinary interval rectangles require only polynomially many blocks

For odd `N`, retain only odd nodes in each original positive interval.
Their odd index lengths need not be powers of two. Decompose each
length into consecutive blocks whose lengths are its nonzero binary
powers. There are at most `O(log(length+1))` blocks per axis.
Their Cartesian products form only a polynomial number of rectangles
of the form (1), with independent exponents `u,v`.

The rectangles are disjoint, so their factor-root polynomials can be
combined without duplicated rows. No exponentially long padding list
or splitting of the longer side into blocks of the shorter side's
size is necessary. Endpoint generation, `H`, both derivative scales,
and all provider parameters have polynomial binary-input cost.
This removes an avoidable equal-size restriction from the interface.

## 6. The remaining matching information is a reciprocal carry

For a column grid of length `s=2^v`, an odd row value `a` satisfies
the exact normalized product identity

\[
 F(a):=2^{-(2s-1)}\prod_{j=0}^{s-1}(N-a(b_0+2j))
   =a^s\,s!_{\rm odd}\binom{(N/a-b_0)/2}{s}.    \tag{15}
\]

The generalized binomial is integral over `Z_2`. Put
`q=(N/a-b_0)/2` and let `j_0` be its residue modulo `s`, chosen in
`0,...,s-1`. All nonmatching factors have fixed total order
`v2((s-1)!)`, so

\[
 v_2(F(a))=v_2(N-aB_{j_0})-v-1                 \tag{16}
\]

when the difference is nonzero, and `F(a)=0` exactly when that
matching difference is zero. Reducing (15) modulo two gives

\[
 F(a)\bmod2=\operatorname{bit}_v(q).             \tag{17}
\]

One can see (17) directly from (16): the normalized product is odd
exactly when `q-j_0` has order `v`. This is also the power-of-two
binomial digit identity.

Equations (15)--(17) eliminate a long row valuation sum and identify
the actual information carried up the matched residue tree. A row
query is cheap: modular inversion gives `j_0`, and one integer
difference gives its depth or exact zero. The unresolved step is
aggregating those row-dependent carries into the global coefficients
in (10). The shorter formulas do not turn that aggregate into a
free operation.

The earlier balanced and rectangular constructors compute complete
primitive jets when the domains are separated by known residue
classes. They do not yet supply the arbitrary-`N` provider in
Section 4. The reflected quadratic construction is another exact
testbed, whose matching is already a square test. Neither restriction
is imposed on the algebraic target proved here.

## 7. Exact controls

The dense reference forms the full regularized pair over
`Q(omega)` without consulting a factor list. The checker then uses
literal factor pairs to verify (10), every coefficient in (11),
the complete Bezout identity, and the decoded moments. Independent
padding choices give the same residues. Singleton, rectangular,
factor-containing, and no-incidence examples are included.

Separate exact controls check the normalized row product, its
valuation, and its reciprocal carry. A 301-bit domain count tests
only the polynomial postprocessor on a jet synthesized from a
prescribed support; it is explicitly not a provider for unknown
factor incidence.

Run:

```sh
python3 tools/check_regularized_bezout_factor_moments.py \
  --output evidence/regularized-bezout-factor-moments.json
```

The proofs and computations used substantial OpenAI assistant support.
No new Lean formalization or improved unconditional factoring bound
is claimed.
