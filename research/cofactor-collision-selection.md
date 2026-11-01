# Exact relation selection from composite residual ideals

## Result and scope

The [selector](../tools/cofactor_collision_relations.py) extends smooth
leading-power selection by cancelling equal residual ideals. Their norms
may be composite. The selector neither factors those norms nor tests them
for primality. It retains a fixed prime-form base and obtains genuine
relations in that base.

The completed 92-bit and 124-bit controls use fewer leading-norm queries
than a smooth-only run with the same round-robin traversal. These finite
results do not prove a uniform short-window estimate. The accepted
factoring exponent remains one half. The
[admissible-quota theorem](cofactor-collision-quota.md) states the missing
uniform estimate precisely.

## 1. Exact partial words and residual keys

Let `N` be odd. Use discriminant `Delta=-N` when `N=3 mod 4`, and
`Delta=-4N` otherwise. Let `g_p` be the consistently oriented prime forms
in the fixed base. In the even-discriminant case include the ramified
class of norm two. Set

\[
 H=\left\lfloor\sqrt{|\Delta|/3}\right\rfloor.
\]

For a reduced primitive target `Q=(A,b,c)`, first handle a proper
`gcd(A,N)`. Otherwise `A<=H` and `gcd(A,N)=1`. If `P` is the product
of the distinct base primes, compute

\[
 S=\gcd\bigl(A,P^{\mathrm{bitlength}(A)}\bmod A\bigr),
 \qquad R=A/S.                                      \tag{1}
\]

This extracts the complete base-supported part: every exponent in `A`
is smaller than or equal to `bitlength(A)`. Thus `gcd(S,R)=1`.
The modular power uses repeated squaring; its exponent is the bit length,
not the integer `A`. Divide `S` only by the already known base primes.
The roots determined by `b mod 2p` give a signed base word `v`.

For `R>1`, define

\[
 b_R=\min(b\bmod 2R,-b\bmod 2R),\qquad
 \kappa=(R,b_R),\qquad
 J_\kappa=\left(R,b_R,\frac{b_R^2-\Delta}{4R}\right). \tag{2}
\]

Let `epsilon=+1` when the original root equals `b_R`, and `-1`
otherwise. The residual norm is odd and unramified. The formula is an
integral proper ideal form. Coprime ideal multiplication, or the Chinese
remainder theorem on the roots, gives the exact class identity

\[
 [Q]=g^v[J_\kappa]^{\epsilon},\qquad
 \sum_p |v_p|\log p=\log S.                         \tag{3}
\]

The ramified-two contribution has exponent one when it occurs. For an
even-discriminant primitive form with `N=1 mod 4`, writing `b=2d`
gives `Ac=d^2+N`; if `A` is even, this has exactly one factor of two.
Its class is self-inverse, so its recorded positive orientation suffices.

Equation (2) distinguishes residual ideals of the same norm. Matching
only `R` would be unsound. At `N=959`, the reduced forms `(15,1,16)`
and `(15,-1,16)` have the same key and opposite orientations, while
`(15,11,18)` has the same composite norm and a different key. A separate
control at `N=4969` uses `(70,2,71)` and its inverse: the base-supported
part is two and the residual norm is the composite integer 35.

## 2. Collision rows are genuine relations

If a target is `g^w`, (3) gives

\[
 g^{u}=J_\kappa,
 \qquad u=\epsilon(w-v).                            \tag{4}
\]

Keep the first normalized word `u_(kappa,0)` for each key. Each later
query with that key supplies the genuine relation

\[
 g^{u-u_{\kappa,0}}=1.                              \tag{5}
\]

A smooth target has `R=1` and supplies `w-v` directly. The construction
works equally for equal and opposite residual orientations. The target
words in the current implementation are `w=k e_p`; no discrete logarithm
or factorization of `R` is needed.

Keeping one reference does not lose any pair difference: `u_i-u_j`
equals `(u_i-u_0)-(u_j-u_0)`. There is also an integral presentation
argument. Add one temporary generator for each key and all equations
`u_i-J_kappa=0`. The first equation eliminates that generator with
coefficient one, leaving exactly the stored-reference differences.
The [quota note](cofactor-collision-quota.md) proves the corresponding
torsion and saturation statements.

## 3. Retaining independent original rows

The implementation uses exact rational elimination to select independent
original integer rows. It does not replace a relation by a rational
linear combination and then treat that combination as a genuine relation.
At most `m` rows are retained, where `m` is the number of base generators.

Let `L` contain every discovered relation and `L'` be the lattice of
retained rows. They have equal rational spans, hence equal saturation

\[
 T=(\mathbb QL)\cap\mathbb Z^m
   =(\mathbb QL')\cap\mathbb Z^m.                   \tag{6}
\]

For the class homomorphism `phi:Z^m -> G`, both lattices lie in its
kernel. The torsion subgroups of `Z^m/L` and `Z^m/L'` have the same
image `phi(T)` in `G`. Thus the complete two-primary relation extractor
has the same useful-involution image after this retention step. This
does not assert equality of integer lattices or their indices.

Every extraction verifies the retained relations by exact ideal
composition before applying two-primary halving. The selector triggers
extraction after a new mod-two dependency, every eighth retained row,
and at full rank; it also extracts once at the end if needed. There are
at most `m+1` extraction calls. The mod-two test schedules extraction;
it does not replace rational-rank selection or the exact extractor.

## 4. A deterministic finite search with all queries charged

For each base prime `p`, start at the least integer `K_p` with
`p^K_p>H`. Visit `K_p, K_p+1, ...` in round-robin order over the sorted
base. Each power advances by exact composition with the known prime
form. Each visited leading norm is charged, including rejected norms.

At residual cap `B`, decode a target only when `R<=B`. The special
choice `B=1` is the matched smooth-only traversal. A run ends at a
proper factor, full fixed-base relation rank, the query cap, or the
end of the finite power windows. Failure to find a factor is a result
of that bounded search and is not a primality certificate.

For `Q` visited targets, the modular extraction tests and form updates
cost `Q poly(n,m,log Q)`. Only admissible targets incur known-prime
division and residual-key storage. An ordered dictionary gives
deterministic logarithmic lookup overhead; the implementation uses the
language's dictionary. There are at most `B bitlength(B)` possible
keys by the root-count proof in the quota note, and at most `Q` are
encountered. Retaining witnesses therefore uses
`O(min(Q,B bitlength(B)) poly(n,m,log Q))` bits. Exact rational
elimination and relation extraction have polynomial bit cost in the
base dimension and row bit lengths, but can be significant in practice.

No uncharged smoothness, residual primality, residual factorization,
relation, or short-window oracle appears in these finite runs. A global
smoothness-density estimate alone would still not control the prescribed
power orbits. The quota theorem requires enough admissible queries along
each target orbit and explicitly leaves the time to collect them open.

## 5. Completed finite comparisons

Each row uses the same target-power order and extraction schedule. The
cofactor variants can stop at different retained ranks because their
relation spans develop differently.

| Input bits | Base cutoff | Residual cap | Norm queries | Retained keys | Independent collision rows |
|---:|---:|---:|---:|---:|---:|
| 92 | 1021 | 1 | 227,593 | 0 | 0 |
| 92 | 1021 | 1,042,441 | 92,062 | 1,053 | 52 |
| 92 | 1021 | 1,064,332,261 | 88,234 | 9,094 | 53 |
| 124 | 4093 | 1 | 3,553,722 | 0 | 0 |
| 124 | 4093 | 16,752,649 | 1,429,783 | 4,959 | 179 |
| 124 | 4093 | 68,568,592,357 | 1,400,316 | 51,108 | 180 |

The 92-bit integer is `4951760154835678088235319297`; the returned
proper factor is `2147483647`. The 124-bit integer is
`21267647932558653899591465697288388633`; the returned proper factor
is `2305843009213693951`. All six runs return verified divisors.
These are query-count comparisons. Additional decoding, storage, and
linear algebra mean that the query ratios are not wall-clock speedups.

## 6. Verification and reproduction

[check_cofactor_collision_relations.py](../tools/check_cofactor_collision_relations.py)
recomputes 2,766 small partial-ideal identities, exercises the composite
and ramified-two fixtures, and replays 13 small complete searches.
For every retained large-run witness, it independently reconstructs the
target power, splits its leading norm, verifies the canonical key and
orientation, checks the collision row, checks rational independence,
and verifies the class identity and returned divisor. The
[retained evidence](../evidence/cofactor-collision-relations.json)
contains the exact inputs, counters, relation rows, and witnesses.

The default checker verifies those retained certificates. Pass
`--replay-search` to repeat all rejected and accepted large queries and
compare the resulting receipts. Use `--regenerate --import-large PATH`
for each completed raw receipt to create the evidence file. The
[quota checker](../tools/check_cofactor_collision_quota.py) separately
verifies the finite key-count, strict-dominance, and saturation controls.

