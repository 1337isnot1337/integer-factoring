# Exact interval-weighted Kloosterman blocks over powers of two

## Positive result and factoring scope

An entire two-parameter family of additive frequencies can be summed in
polynomial bit time, with both ordinary interval weights retained. Let
`q=2^m`, let `N` be odd, and let `I,J` be signed integer intervals with
binary endpoints. Write

\[
 e_q(t)=\exp(2\pi i t/q),\qquad
 \widehat I(u)=\sum_{a\in I}e_q(-ua),\qquad
 K_q(u,v)=\sum_{x\bmod q\atop x\text{ odd}}e_q(ux+vx^{-1}).
 \tag{1}
\]

For an odd residue `c0 mod 2^h`, define the complete frequency block

\[
 G_{m,h}(c_0;I,J)=
 \sum_{u\bmod q\atop u\text{ odd}}
 \sum_{c\bmod q\atop c\equiv c_0\ (2^h)}
 \widehat I(u)\widehat J(cu)K_q(u,cuN).
 \tag{2}
\]

If

\[
 \max(1,\lfloor m/2\rfloor)\le h\le m,
 \tag{3}
\]

then `G` is an integer computed by at most eight quadratic roots and
sixteen affine interval-congruence counts. Every affine count takes
polynomial bit time by Euclidean floor sums. The block contains

\[
 2^{m-1}\,2^{m-h}
 \tag{4}
\]

frequency pairs. Neither those pairs, the interval elements, nor a
cyclotomic field basis are enumerated.

At `m>=2` and `h=floor(m/2)`, the block has a particularly useful exact
interpretation. Put

\[
 C_r(c_0,h)=\#\{(a,b)\in I\times J:
       ab\equiv N\pmod{2^r},\quad a\equiv c_0b\pmod{2^h}\}.
\]

Then

\[
 \boxed{G_{m,h}(c_0;I,J)
       =\frac{q^2}{2}\bigl(2C_m(c_0,h)-C_{m-1}(c_0,h)\bigr).}
 \tag{5}
\]

If every pair in the rectangle satisfies `|ab-N|<q/2`, both congruence
counts in (5) are exact integer factor counts. Consequently

\[
 \boxed{\frac{2G_{m,h}(c_0;I,J)}{q^2}
       =\#\{(a,b)\in I\times J:ab=N,\ a\equiv c_0b\pmod{2^h}\}.}
 \tag{6}
\]

This is a positive factor-count query for one ratio class. Restricting
the first interval and using binary search recovers a factor from any
class with a nonzero count, in polynomially many queries.

The complete outer sum still has `2^(floor(m/2)-1)` ratio classes at the
largest conductor. The resulting full evaluator has
`2^(m/2)*poly(m)` bit cost, with endpoint bit lengths included. This
does not improve the established factoring score `1/2`. The positive
advance is an exact collective frequency operation that pays for the
interval weights and identifies its factor content. A cheaper way to
locate or sum the nonzero ratio classes remains necessary.

## 1. Relation to earlier representations

The [inverse-section algorithm](dyadic-inverse-sections.md) already
counts affine fibres at half the modulus precision. The
[complete fractional-linear orbit average](dyadic-orbit-aggregation.md)
sums a different complete group family. The
[original Fourier gate note](fourier-multiplication-gates.md) gives
characteristic-two cancellation for complete reciprocal phases at odd
moduli, while explicitly retaining the cost of the interval weights.

Here the modulus is a power of two, and the calculation uses ordinary
complex additive characters only to prove an integer identity. The
implemented algorithm is entirely integral. It sums both a full unit
frequency orbit and a complete progression of frequency ratios. In
particular, evaluating a single unweighted Kloosterman sum is not being
counted as evaluation of (2).

Stationary phase for prime-power exponential sums is classical; see
Dąbrowski and Fisher, [*A stationary phase formula for exponential sums
over Z/p^m Z and applications to GL(3)-Kloosterman sums*, Acta Arithmetica
80 (1997), 1–48](https://doi.org/10.4064/aa-80-1-1-48). The special
one-variable calculation needed below is proved directly. No analytic
estimate or unproved cancellation bound is used.

## 2. At most eight stationary representatives

Put

\[
 j=\lfloor m/2\rfloor,\qquad k=m-j,
 \qquad
 \mathcal R(c)=\{r:0\le r<2^k,\ r\text{ odd},\ r^2\equiv cN\pmod{2^j}\}.
 \tag{7}
\]

For every odd `u`,

\[
 K_q(u,ucN)
 =2^j\sum_{r\in\mathcal R(c)}
       e_q\bigl(u(r+cNr^{-1})\bigr).
 \tag{8}
\]

To prove this, write each odd `x mod q` uniquely as
`x=r+2^k t`, with `0<=r<2^k` and `0<=t<2^j`. Since `2k>=m`,

\[
 (r+2^kt)^{-1}\equiv r^{-1}-2^ktr^{-2}\pmod q.
\]

The inner sum over `t` is therefore

\[
 \sum_{t\bmod 2^j}e_{2^j}\bigl(ut(1-cNr^{-2})\bigr).
\]

It equals `2^j` precisely when `r^2=cN mod 2^j`, and otherwise
vanishes. This proves (8), including the small conductors.

An odd square has at most four roots modulo `2^j`. For `j>=3`, if `r,s`
are roots, one of `r-s,r+s` has valuation one; the other must be
divisible by `2^(j-1)`. Thus every root is congruent to `r` or `-r`
modulo `2^(j-1)`, giving at most four possibilities. The smaller
moduli can be checked directly. Passing from modulus `2^j` to the
representatives in (7) adds at most a factor of two, so

\[
 |\mathcal R(c)|\le8.
 \tag{9}
\]

Successively testing the two lifts of each surviving square root
constructs the complete list with `O(m)` operations on `O(m)`-bit
residues. In the exceptional notation `j=0`, the unique residue modulo
one is lifted and the odd representatives are retained.

## 3. Summing the interval weights and the ratio progression

Let

\[
 R=2^{m-h},\qquad c=c_0+2^h t\quad(0\le t<R).
\]

Condition (3) ensures that `c mod 2^j` is fixed. Hence the same at most
eight representatives `r` occur in (8) for every ratio in the block.
For an integer `E`, define the primitive dyadic character sum

\[
 \mathcal A_q(E)=\sum_{u\bmod q\atop u\text{ odd}}e_q(uE)
       =q\,1_{q\mid E}-\frac q2\,1_{q/2\mid E}.
 \tag{10}
\]

This is the full character sum minus the even-frequency sum.

Fix `r` and a physical pair `(a,b)`. Write `w=Nr^{-1} mod q`. The
phase in (2), after (8), is

\[
 u\bigl(r-a+c_0(w-b)+2^ht(w-b)\bigr).
\]

For every odd `u`, summing over `t` gives zero unless `w=b mod R`;
when this congruence holds it gives `R`. Equivalently,

\[
 \sum_{t\bmod R}\mathcal A_q(E+2^htD)
     =R\,1_{R\mid D}\,\mathcal A_q(E).
 \tag{11}
\]

This identity performs the entire ratio-progression sum. Its cost has
no factor of `R`.

Let `r0=r mod R` and `w0=w mod R`, using least nonnegative residues.
The nonzero terms in (11) require `b=w0 mod R`. Equation (10) then
also requires `a=r0 mod R`, since `R` divides `q/2`. Condition the
ordinary intervals on these progressions:

\[
 U_r=\{x\in\mathbb Z:r_0+Rx\in I\},\qquad
 V_r=\{y\in\mathbb Z:w_0+Ry\in J\}.
 \tag{12}
\]

These are again ordinary integer intervals. Their endpoints are
obtained by integer ceiling division. Define the integer

\[
 t_r=\frac{r+c_0w-r_0-c_0w_0}{R}.
 \tag{13}
\]

For a positive modulus `d`, put

\[
 L_d(c,t;U,V)=\#\{(x,y)\in U\times V:x+cy\equiv t\pmod d\}.
 \tag{14}
\]

Expanding both interval weights in (2), then applying (8)–(13), gives
the complete evaluation formula

\[
 \boxed{G_{m,h}(c_0;I,J)=2^jR
  \sum_{r\in\mathcal R(c_0)}
   \left(qL_{2^h}(c_0,t_r;U_r,V_r)
       -\frac q2L_{2^{h-1}}(c_0,t_r;U_r,V_r)\right).}
 \tag{15}
\]

The formula uses integer roots, inverses modulo `q`, interval
conditioning, and at most sixteen affine counts. No numerical roots
of unity are constructed.

## 4. Exact affine counts and their bit cost

For an interval of length `n=ad+r`, `0<=r<d`, every residue modulo
`d` occurs `a` times and the excess residues form one cyclic interval
of length `r`. Decompose both intervals in (14) this way. Because
`c` is odd and `d` is a power of two, multiplication by `c` permutes
the residue classes. The contributions from full periods are

\[
 a_Ua_Vd+a_Ur_V+a_Vr_U.
 \tag{16}
\]

Each excess cyclic interval splits into at most two ordinary residue
intervals. For each of the resulting at most four pairs, count the
indices whose affine residue lies in an ordinary residue interval.
The prefix formula is

\[
 \#\{0\le i<n:(Ai+B\bmod d)<T\}
 =n+S_d(n,A,B)-S_d(n,A,B+d-T),
 \tag{17}
\]

where

\[
 S_d(n,A,B)=\sum_{i=0}^{n-1}\left\lfloor\frac{Ai+B}{d}\right\rfloor.
\]

Signed `A,B` are handled first by Euclidean division. With their
remainders between zero and `d-1`, transpose the lattice region below
the line to exchange the residual slope and modulus. Repeating gives
the ordinary Euclidean floor-sum algorithm. It takes `O(log d)`
divisions, with all endpoint and count sizes retained. This proves
the polynomial bit bound for (14) and hence for (15).

If all original endpoints have at most `E` bits, every integer used in
(15) has `O(E+m+log(1+|N|))` bits. The exact output has a polynomial
number of bits as well: the literal definition bounds its magnitude
by `q^2 R |I||J|`. School integer arithmetic already proves the
claimed bit complexity. No unit-cost convention for large integers
or exponential-size algebraic numbers is needed.

## 5. Why the minimum blocks are ratio-filtered factor counts

Assume `m>=2` and `h=j=floor(m/2)`. Then `R=2^k`, so for a physical
pair the row condition in (12) selects at most one representative,
namely `r=a mod R`. The column condition says

\[
 ab\equiv N\pmod R.
 \tag{18}
\]

Together with the stationary equation `r^2=c0*N mod 2^j`, this is
equivalent to the additional ratio condition

\[
 a\equiv c_0b\pmod{2^j}.
 \tag{19}
\]

All relevant `a,b,r,c0` are odd. For the phase
`E=r-a+c0(N/r-b)`, multiplication by `r` gives

\[
 rE\equiv(r-a)(r-c_0b)+c_0(N-ab)\pmod q.
 \tag{20}
\]

The first product vanishes modulo `q`: its two factors are divisible
by `2^k` and `2^j`, and `j+k=m`. Multiplication by an odd unit does
not change either divisibility test in (10), so

\[
 \mathcal A_q(E)=\mathcal A_q(N-ab).
\]

Moreover, every pair that contributes to (10) satisfies (18), since
`q/2` is divisible by `R` for `m>=2`. Therefore the auxiliary condition
(18) can be dropped after (10) is applied. Since `2^jR=q`, equation
(15) becomes exactly (5).

For (6), an integer difference of absolute value strictly below
`q/2` is divisible by `q/2` or by `q` only when it is zero. A bound
on all rectangle differences is obtained from the four corner
products, because `ab-N` is affine in each variable separately.
Thus the no-wrap condition can be certified from binary endpoints.

If the count in (6) is positive, bisect `I` and evaluate the same block
on each chosen prefix. The ratio precision and the original certified
`m` remain fixed. Positivity identifies a half containing a factor
pair. After logarithmically many bisections, the row `a` is known;
exact integer division recovers `b=N/a`. This does not identify a
successful ratio class before the queries begin.

## 6. Complete Fourier inversion and the remaining outer cost

Let `Q=2^M`, and let `A_odd,B_odd` be the numbers of odd integers in
the two intervals. Fourier inversion on the inverse graph gives

\[
 \#\{(a,b)\in I\times J:ab\equiv N\pmod Q\}
 =\frac1{Q^2}\sum_{u,v\bmod Q}
       \widehat I(u)\widehat J(v)K_Q(u,vN).
 \tag{21}
\]

If `u,v` have different two-adic valuations and their joint conductor
is at least four, the Kloosterman sum vanishes. Indeed, after their
common power of two is removed, one coefficient is odd and the other
even. Pairing `x` with `x+q/2` gives opposite phases. If both are
divisible by `2^(M-m)`, reduction to conductor `q=2^m` contributes
the factor `2^(M-m)`.

The zero frequency and the conductor-two terms combine to
`2Q*A_odd*B_odd`. At every conductor `m>=2`, the surviving primitive
frequencies are uniquely `(u,cu)` with both `u,c` odd. Partition the
ratios by `c mod 2^h_m`, where `h_m=max(1,floor(m/2))`. Equation (21)
therefore has the entirely integral form

\[
 \boxed{\#\{ab\equiv N\ (Q)\}
  =\frac1{Q^2}\left(
   2Q A_{\mathrm{odd}}B_{\mathrm{odd}}
   +\sum_{m=2}^M2^{M-m}
      \sum_{c_0\bmod2^{h_m}\atop c_0\text{ odd}}
           G_{m,h_m}(c_0;I,J)\right).}
 \tag{22}
\]

The number of outer blocks is

\[
 \sum_{m=2}^M2^{h_m-1}=O(2^{M/2}).
 \tag{23}
\]

Streaming (22) needs polynomial bit space. Its worst-case time is the
same exponential scale as the earlier affine-fibre evaluator. In
particular, the `2^(3m/2)` frequency pairs inside one block cannot be
used to claim a faster global count while omitting (23).

The [inverse-window audit](circuit-box-inverse-windows.md) explains why
the straightforward physical residue tree can have exponentially many
live nodes even on an empty rough-semiprime rectangle. Formula (15)
does combine frequency families rather than visiting those nodes, but
the number of minimum ratio classes remains exponential. In the
no-wrap setting (6), the new concrete target is a collective query or
search over the sparse nonzero ratio classes. A fast point query and
sparsity alone do not provide that search.

## 7. Exact controls

The [implementation](../tools/circuit_box_kloosterman_blocks.py) includes
the block evaluator and the complete reference inversion (22). The
[independent checker](../tools/check_circuit_box_kloosterman_blocks.py)
and [receipt](../evidence/circuit-box-kloosterman-blocks.json) record:

- 1,000 affine interval counts against literal pairs, with signed
  endpoints and coefficients;
- 99 complete square-root lists and 320 stationary-phase identities
  in exact cyclotomic coefficient vectors;
- 288 opposite-parity frequency cancellations;
- 240 whole frequency blocks against independent complete ratio and
  physical-node scans;
- 48 raw weighted sums evaluated literally in
  `Z[zeta]/(zeta^(q/2)+1)`, checking both the integer value and the
  vanishing of every nonconstant coefficient;
- 600 exact identities (5), 400 full rectangle counts, and 149
  no-wrap integer-factor rectangle counts.

Four large controls have conductor exponents `127,256,513,1024`.
Each uses a singleton row and column in the matching ratio class.
The three targets differ from their product by `0,q/2,q/4`; the
independent exact-congruence answers in (5) are respectively
`q^2/2,-q^2/2,0`. An adjacent wrong ratio class gives zero. At
`m=1024`, one evaluated block contains `2^1535` frequency pairs.
The implementation still constructs at most eight roots and makes
at most sixteen affine count queries.

```sh
python3 tools/check_circuit_box_kloosterman_blocks.py \
  --output evidence/circuit-box-kloosterman-blocks.json
```

These controls establish reproducibility of the displayed algorithm.
They do not change its charged outer complexity or the factoring score.
