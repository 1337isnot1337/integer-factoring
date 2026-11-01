# Blocking and precision give monotone complete-channel hierarchies

The [complete subspace oracle](subspace-smith-channel-oracle.md) permits
two systematic extensions without worsening its optimal global rate:

* Group several consecutive sites and refactor their full physical
  action together.
* Add auxiliary lattice states, including all states at a finer
  full-tail precision.

The proof uses regular right inverses of the physical factors. They
show that inherited path columns generate the exact column lattice
of the new action. A subset of those columns is already a regular
basis of that lattice, so the complete oracle can retain their
targets and charged valuations.

These guarantees use the target-enlargement and undercharge freedom
of the subspace model. They do not establish monotonicity for the
earlier fixed-pivot or exact-target constructions.

There is a strict finite improvement on the binary companion. At the
same scale `s=5`, its recorded per-original-site product bounds are:

| Construction | State precision | Exact upper bound per original site |
|---|---:|---:|
| One site | 2 | `12.8212` |
| Two-site block | 2 | `<12.542` |
| Three-site block | 1 | `<12.386` |

The one-site complete model has the independent exact lower bound
`12.8203599666` at this scale. Both block certificates are strictly
better than every one-site policy in that declared model. These
rates remain above the separate elementary scalar certificate;
the complete factoring exponent is unchanged.

## 1. A column-lattice basis can be selected from generators

Work over `O=F_2[[tau]]`. Suppose a regular matrix is factored as

$$A=VR,$$

and `R` has a regular right inverse `S`. Then

$$V=AS,\qquad \mathrm{col}_O(A)=\mathrm{span}_O(V).\tag{1}
$$

Let the common column lattice have rank `r`. Choose any `O`-basis
`U` of it and write `V=UC`, with `C` regular. The columns of `C`
generate `O^r`. Therefore their residues span `F_2^r`; select `r`
residue-independent columns. The corresponding square submatrix
of `C` has unit determinant, so those selected columns of `V`
already form an `O`-basis of the lattice.

No arbitrary linear combinations of the generator columns are
needed for this selection. This matters because each column may
already carry a useful target lattice and valuation charge.
Selecting a subset retains those labels and cannot increase a
positive sum of their costs.

For a Schur factorization used by the channel oracle, the right
factors have a triangular unit minor on their pivot columns. Hence
they have a regular right inverse. A regular unimodular channel
basis change preserves that property.

## 2. Conditional path factors retain a right inverse

At a source state, let the first site's right factor be `R_0`,
with `r_0` rows and `n=q^2` physical columns. For its channel `j`,
the next state has a right factor `R_j`, with `r_j` rows and the
same `n` physical columns. The two-site right factor has entries

$$
(R_{\rm path})_{(j,k),(b_1,b_2)}
 =(R_0)_{j,b_1}(R_j)_{k,b_2}.\tag{2}
$$

Up to the fixed ordering of its columns, this is the product of
the block diagonal matrix with blocks `R_j` and `R_0 tensor I_n`.
Both factors have regular right inverses. Their product therefore
has one as well. The ranks `r_j` can vary, and the factors need
not be square.

Iterating proves the same statement for every block length. There
are at most `n^b` channel paths in a block of `b` sites. The right
factor maps physical-word coordinate vectors onto channel-path
coordinate vectors. Its regular right inverse is the algebraic
fact needed in (1).

The next-state choice must depend only on the channel path, as it
does in the existing state-adaptive construction. A factorization
chosen by inspecting the opposite physical word would not have
the conditional product structure in (2).

## 3. The complete block oracle dominates inherited policies

Fix a one-site policy with transition matrix `T(s)`, and fix a
positive potential `h` on its closed state family. Multiply its
local identities for `b` sites. A path contributes a left column
whose auxiliary rows lie in its final target after removing the
accumulated charged valuation

$$a_{\rm path}=a_1+\cdots+a_b.$$

These path columns, together with (2), factor the exact flattened
block action. By (1), they generate its column lattice. Choose an
`O`-basis subset of them. The selected columns keep their final
targets and charges, so their total potential cost is at most
the sum over all paths:

$$
(F_{b,s}(h))_M\le(T(s)^b h)_M.\tag{3}
$$

The complete subspace oracle permits this basis. Indeed, any two
`O`-bases of the same column lattice differ by a regular unimodular
matrix. Each selected column's coordinates in the block Schur
basis are primitive. Its actual valuation is consequently at most
the largest block Smith order, while its charged valuation is no
larger than its actual valuation. Thus its accumulated charge lies
inside the oracle's enumerated range. Its target is a declared
state, and exact generation of that target is unnecessary.

At infinity, the block degree is `bL`, so a path's weight is
`s^(bL-sum a_i)`, exactly the product of the one-site weights.
At zero the same equality holds with weights `s^(-sum a_i)`.
This proves (3) at both endpoint cuts.

If `Th<=beta h`, then (3) gives `F_b(h)<=beta^b h`. The block
oracle selects a fixed policy attaining its minimum row costs,
so it certifies a block rate at most `beta^b`. Taking infima over
positive potentials and one-site policies gives

$$\lambda_b(s)\le\lambda_1(s)^b.\tag{4}$$

More generally, applying the same argument to blocks of length
`a` gives

$$\lambda_{ab}(s)^{1/(ab)}\le\lambda_a(s)^{1/a}.\tag{5}
$$

Thus the per-site hierarchy is monotone along divisibility, for
example along repeated doubling. No assertion about every pair
of consecutive block lengths is needed. The two endpoint rates
can be treated separately and their bounds multiplied.

This result concerns the optimal complete oracle. A particular
floating policy search can miss that optimum and need not display
monotonicity until it finds or reuses an inherited certificate.

## 4. Adding states makes the optimal global rate nonincreasing

Suppose a new source lattice `M` is contained in an old source
lattice `N`. Write

$$B_M=W B_N$$

with `W` regular. Multiply the old state's local identity on the
left by `W`. In the flattened matrix, this is the regular row
operator `W tensor I`. The inherited right factor still has its
regular right inverse. Its new left columns therefore generate
the exact column lattice of the new source action.

Choose a basis subset as in Section 1. Each selected channel
retains its old target and its old nonnegative charge. Multiplying
by `W` can increase its actual valuation or shrink the span of its
auxiliary rows; the complete oracle permits both changes by
undercharging and target enlargement.

Consequently the new source has a valid policy row whose targets
are all old states. Extend any old policy this way at every added
state. In old/new state order, the resulting matrix is

$$\widetilde T=\begin{pmatrix}T&0\\B&0\end{pmatrix},\tag{6}$$

so its spectral radius is exactly `rho(T)`. The new states are
transient into the old graph. The minimum global policy rate
therefore cannot increase.

If the old family contains `O^D`, every new sublattice of `O^D`
satisfies the required containment. In particular, the complete
normalized full-tail family at precision `P+1` contains the one
at precision `P`, since a lattice containing `tau^P O^D` also
contains `tau^(P+1) O^D`. Hence

$$\lambda_{P+1}(s)\le\lambda_P(s)\tag{7}$$

for this complete subspace model.

For an elementary potential version, assign each added state the
height of a containing old state. Its selected inherited row has
cost at most that old state's row cost. Any old positive upper
potential therefore extends to an upper potential at the same
rate. This is the finite precision check used below.

Neither (6) nor (7) applies automatically when channel bases are
fixed, when only one Schur factorization is allowed, or when target
membership is replaced by exact target equality. Those restrictions
can remove the inherited basis used by the proof.

## 5. Exact companion block certificates

Let `K(z;x,y,u,v)` be the binary companion used by the preceding
adaptive notes. For a block of length `b`, encode each of the four
physical labels as a binary word of length `b`, least significant
digit first, and put

$$K^{[b]}(z;X,Y,U,V)
=K(z;x_1,y_1,u_1,v_1)\cdots K(z;x_b,y_b,u_b,v_b).$$

This is a polynomial local matrix of auxiliary dimension two,
alphabet size `2^b`, and degree at most `b`. Its first column at
`z=1` is `(E_block,0)^T`, so the same first-coordinate boundaries
give the required ordered crossing at every block length.

At `s=5`, use the ten precision-two states for the two-site block.
The integer potentials are

$$
\begin{aligned}
h_\infty={}&(1000000,512311,1000000,1000000,512311,1000000,512311,1000000,1000000,1000000),\\
h_0={}&(1000000,906557,645834,617226,895729,637956,895729,637956,540672,605830).
\end{aligned}
$$

Exact greedy selection and reconstruction certify

$$
\beta_\infty\le13904/125,\qquad
\beta_0\le707/500.
$$

Their block product is at most

$$\frac{2457532}{15625}=157.282048<(6271/500)^2.$$

Thus the per-original-site rate product is less than `12.542`.
This is strictly smaller than the
[global one-site lower certificate](channel-pressure-separation.md)
at `s=5`, which is `12.8203599666`.

The three-site block already improves the recorded upper bound at
precision one, with four states. Use

$$
h_\infty=(1000000,512311,1000000,1000000),\qquad
h_0=(1000000,895993,628166,625925).
$$

The simple exact rates

$$\beta_\infty\le1140,\qquad\beta_0\le5/3$$

give block product at most `1900 < (6193/500)^3`. Its
per-original-site bound is consequently less than `12.386`.
These rational rates deliberately leave slack; the proof does not
use the floating estimates that suggested the potentials.

The local-identity and positive-potential theorem supplies all
block lengths. The finite coefficient-rank controls are additional
checks, not an extrapolation to the asymptotic statement.

## 6. Reproduction

[block_subspace_channels.py](../tools/block_subspace_channels.py)
constructs the exact binary-polynomial block matrices. It also
implements unit-minor right inverses and conditional right-factor
products, including rectangular varying-rank inputs.

Run the shorter exact replay with:

```sh
python3 tools/check_block_subspace_channels.py
```

Include the rank-sixty-four three-site block with:

```sh
python3 tools/check_block_subspace_channels.py --maximum-block 3
```

The [checker](../tools/check_block_subspace_channels.py) verifies
every selected local rational identity and target membership. It
checks exact block domination and precision domination on supplied
positive potentials, and it checks regular right inverses of
conditional path factors. The complete literal coefficient tensors
at original lengths two, three, and four are compared with the
resulting block bounds. The
[evidence record](../evidence/block-subspace-channels.json) stores
both block policies, potentials, transfer matrices, and counts.

The larger alphabet makes the full block input grow exponentially
with block length. The theorem does not turn unbounded blocking
into a polynomial-size search. It does provide a rigorous way to
spend additional finite computation on a richer certificate family,
while retaining all earlier complete-oracle upper bounds.

Substantial OpenAI assistant support was used. No Lean formalization
or stronger complete factoring score is claimed.
