# Two coprime-basis passes attach every supplied orientation

This is a batched refinement of the
[supplied oriented-basis construction](oriented-coprime-basis.md).
It is an ordinary algorithmic consequence of
[Bernstein's coprime-basis and sparse-factorization algorithms](https://cr.yp.to/lineartime/dcba-20040404.pdf).
The classical norm engine is not implemented by this note. The supplied
data theorem gives no bound on the time needed to generate useful
class-group relations.

Let the input keys be `(R_i,b_i)`, with odd `R_i>1`,
`gcd(R_i,Delta)=1`, and `b_i^2=Delta mod 4R_i`. Normalize each root
modulo `2R_i`. Count every orientation key separately, even if norms
repeat, and write `S=sum_i bits(R_i)`. Unit keys contribute identity
words and can be removed first.

The first coprime-basis pass gives pairwise coprime `A_j>1` and sparse
positive exponents

\[
 R_i=\prod_j A_j^{e_{ij}}.
 \tag{1}
\]

For each used `A_j`, select one incident reference root. For every
incidence `e_ij>0`, compute

\[
 \beta_{ij}=b_i\bmod 2A_j,\qquad
 \delta_{ij}=\frac{\beta_{ij}-\beta_{0j}}2\bmod A_j.
 \tag{2}
\]

Both roots have the parity of `Delta`, so the division is exact.
Product/remainder trees compute all residues for row `i` with total
modulus length `O(bits(R_i))`: indeed
`sum_(j:e_ij>0) log A_j<=log R_i`. Pairwise coprimality of the moduli
`2A_j` is not required for a remainder tree. The normalized differences
have total bit length `O(S)`.

Apply a second coprime-basis pass to all `A_j` and all nonzero
`delta_ij`. Record zero differences as flags; omit units if desired.
Let its output components be `C_l>1`, with sparse factorizations

\[
 A_j=\prod_l C_l^{u_{jl}},\qquad
 \delta_{ij}=\prod_l C_l^{v_{ij,l}}\quad(\delta_{ij}>0).
 \tag{3}
\]

Discard components that occur in no `A_j`. Each retained component
belongs to exactly one parent `A_j`, because the first-pass bases are
pairwise coprime.

For a retained child `l` of `j`, the sign of row `i` relative to the
reference is plus if `delta_ij=0` or `v_ij,l>0`, and minus otherwise.
To prove this, the two half-differences of the roots have product zero
modulo `A_j`, while no prime dividing `A_j` divides both half-differences.
Thus their gcd with `A_j` selects whole prime-power components.
A second-pass component cannot mix selected and unselected primes,
since it also factors `delta_ij`. In fact, a nonzero `v_ij,l` must
satisfy `v_ij,l>=u_jl`.

Give `C_l` the reference root `beta_0j mod 2C_l`, using another remainder
tree within each parent. Replacing this root
by its canonical representative up to sign multiplies all corresponding
signs by one recorded column sign. The absolute input exponent is

\[
 |w_{il}|=e_{ij}u_{jl}.
 \tag{4}
\]

Unique unramified root lifting proves the exact integral-ideal identity:
use the reference ideal to this positive power for a plus sign, and
its conjugate for a minus sign. Signed exponents express the resulting
class-group word; a negative entry is not silently treated as an
integral inverse ideal.

Construct the words by traversing each original incidence and its own
parent's child list. This avoids a dense key-by-component scan. For
each row, the number of emitted entries is at most `log_2 R_i`, since
the positive norm product is still (1). Hence all incidence and
membership bookkeeping costs `O(S polylog S)` bit operations, including
indices and exponents.

Both coprime-basis inputs have `O(S)` integer bits. The two classical
passes and the remainder trees therefore use essentially linear
integer-data time. Root-congruence validation and reduced-form
construction add `(S+q) poly(n+log S)`, where `q` is the number of
keys and `n` bounds discriminant and norm bit lengths. No prime
factorization, residual primality test, or adaptive orientation-split
tree is required by this construction. Extra splitting induced by
unrelated differences only refines the valid basis.
