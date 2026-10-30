# Exact compression of the crossing coefficients

Status: a proved integer-algebra construction and a quantitative consequence for
the recursion. This file is not a Lean verification of a modified factoring
machine. The companion checker verifies the finite identities and small
coefficient examples with exact integers.

The upstream sources are pinned to
`IsaacRudich/SubexponentialDeterministicFactoring` commit
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`, in particular
[`Cross.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Cross.lean),
[`Pieces.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Pieces.lean),
and [`Level.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Level.lean).

## The integer identity

Write

\[
 E(x,y,u,v)=[x=v][y=u],\qquad R+C=2E,
\]

for upstream's eight-square crossing identity on three-bit blocks. Rows are
indexed by `(x,y)` and columns by `(v,u)`, so `E` is the identity matrix. The
row matrix of `R` is `AA^T`, where the eight columns of `A` are the functions
`psi_k(x,y)`. Their supports are disjoint, each support has eight entries of
magnitude one, and hence `A^T A=8I`.

Define

\[
 S=8E-R=3R+4C.
\]

`S` has an integer factorization with 56 channels, not 64. This is stronger
than saying that its rational matrix rank is 56: no division is needed.

For each XOR class `k`, write `s_x=sgn(x,x xor k)`; `s_0=1`. For `a` in
`{1,...,7}`, set

\[
 L_{x,a}=[x=a]-[x=0]s_a,\qquad
 H_{a,v}=8[a=v]-s_a s_v.
\]

Then

\[
 \sum_{a=1}^{7} L_{x,a}H_{a,v}=8[x=v]-s_xs_v.
\]

If `x!=0`, only `a=x` contributes. If `x=0`, the sum is
`-8[v!=0]s_v+7s_v`, which is 7 at `v=0` and `-s_v` otherwise.
Attach the indicator `x xor y=k` to the left factor and `v xor u=k` to the
right factor and sum over all eight `k`. This proves the claimed
56-channel factorization of `S` over the integers, and therefore over every
commutative ring.

`formal/CrossingComplement.lean` checks the full identity in the upstream
three-bit definitions. Its computable `complementFac` uses six channel bits:
three for the XOR class and three for a padded seven-way index. One index
per class is zero on both sides, leaving 56 potentially active labels. Lean
checks its value is `8E-R=3R+4C`. This local factorization does not implement
the multi-block compressed coefficient decoder or the factoring evaluator.

## Compressing every coefficient

Let `t=ceil(lambda/3)` and let `A_i` be upstream's coefficient of `T^i` in
the blockwise product `prod_s (R_s+T C_s)`. Define the scaled pieces
`B_i=4^t A_i`. The scaling removes all denominators:

\[
 B_i=[T^i]\prod_{s=1}^{t}\big((4-3T)R_s+TS_s\big).
\]

For a subset `Q` of blocks with `j=|Q|<=i`, its coefficient is

\[
 w_{t,i,j}=\binom{t-j}{i-j}4^{t-i}(-3)^{i-j}.
\]

Thus

\[
 B_i=\sum_{Q:\,|Q|\le i}w_{t,i,|Q|}
       \prod_{s\notin Q}R_s\prod_{s\in Q}S_s.
\]

`R` has eight channels and `S` has 56, giving the exact upper bound

\[
 r_i=\sum_{j=0}^{i}\binom tj8^{t-j}56^j
       \le 16^t7^i.
\]

The last inequality follows from `7^j<=7^i` and
`sum_{j<=i} binom(t,j)<=2^t`. An explicit channel encoding uses `t` bits
for `Q`, three bits per block for its XOR class, and `i` base-seven
digits for the nonpivot labels of the selected blocks. Unused digits are
fixed to one specified value. Channels with `|Q|>i` have zero weight.

Column-wise the roles of `R` and `C` are exchanged. Reversing the
coefficient index gives the bound `r_{t-i}<=16^t7^(t-i)`. Consequently the
two binary factor widths can be chosen with

\[
 c_i^{\rm row}+c_i^{\rm col}
 \le 8t+\lceil i\log_2 7\rceil
          +\lceil(t-i)\log_2 7\rceil
 \le (8+\log_2 7)t+2.
\]

The factors can also use the trivial width `2 lambda` when smaller, exactly
as upstream does. All channel functions can be evaluated in polynomial time
in `lambda`: scan the blocks, evaluate the fixed-size factors, and compute
the one binomial coefficient by exact integer arithmetic. Intermediate
integer bit lengths are `O(t)` apart from the input and channel indices.

The new pieces satisfy

\[
 \sum_{i=0}^{t}B_i=8^t[x=v][y=u],\qquad |B_i|\le 8^t.
\]

The multiplier in the counting recursion must therefore use `8^t` per
rectangle, replacing upstream's `2^t`. Leaving that multiplier unchanged
would give an incorrect algorithm. This replacement only changes a
constant in the bit-length bounds; it creates no need for modular division.

## Improved contraction

Set

\[
 \rho=\frac{8+\log_2 7}{12}=0.9006129101714669\ldots.
\]

For a `G by G` grid, the unchanged row/column construction gives

\[
 \min(U_{\rm row},U_{\rm col})
 \le G^3\lambda+\frac{G^4}{4}
       \max_i(c_i^{\rm row}+c_i^{\rm col})
 \le (\rho+1/G)G^4\lambda+3G^4.
\]

The last additive constant is conservative: `t<=lambda/3+2/3` gives
`(2rho+1/2)G^4<3G^4`.

Combined with a grouping whose pair widths satisfy
`K^2 Y <= U+GP/D`, with `P=G^3 lambda`, this yields

\[
 N'\le(\rho+1/G+1/D)N+3G^4,
 \qquad N=G^4\lambda,\quad N'=K^2\lambda'.
\]

A completely rational envelope is available from

\[
 7^{26}=9387480337647754305649
     <9444732965739290427392=2^{73}.
\]

It gives `rho <= 281/312`. Using the growing-group schedule whose
`log K_d=Theta(4^d)` and whose `1/D_d` errors are summable, the corresponding
log-saving powers are

\[
 \alpha=\frac{\log(1/\rho)}{\log4}
       =0.07551046810477115\ldots,
 \qquad
 \alpha_{\rm rational}=\frac{\log(312/281)}{\log4}
       =0.075487949\ldots.
\]

The recurrence and compiler-cost argument must be supplied separately to
turn this contraction into a factoring theorem. A finite matrix identity
alone is not a running-time proof.

## Sharpness for this coefficient family

Over the rationals, the channel count `r_i` above is the exact row rank.
On the eight-dimensional image of `A`, the pair `(R,C)` has eigenvalues
`(8,-6)`; on its 56-dimensional complement it has eigenvalues `(0,2)`.
For a tensor sector containing `j` complement positions, the coefficient
polynomial is `(8-6T)^(t-j)(2T)^j`. Its coefficient of `T^i` is nonzero
exactly when `j<=i`. Summing the sector dimensions proves the rank formula.

For even `t` and `i=t/2`, the summand with `j=i` gives

\[
 r_{t/2}^2\ge
 \binom{t}{t/2}^{2}\,448^t
 =2^{(8+\log_2 7)t-O(\log t)}.
\]

Therefore a uniform bound obtained only by refactoring these same
coefficients cannot lower the leading contraction below `rho`. A further
improvement would need to change the decomposition or exploit structure
outside the individual coefficient matrices. This is a limitation of this
coefficient method, not a lower bound for integer factoring.
