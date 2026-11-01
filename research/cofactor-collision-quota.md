# A finite admissible quota for complete cofactor relations

## Result and remaining hypothesis

The [cofactor collision selector](../tools/cofactor_collision_relations.py)
cancels equal oriented residual ideals without factoring their norms.
This note gives an explicit sufficient number of admissible queries
**for each target generator**. Reaching that quota produces a full
rational-rank relation family. The cost of collecting that many
admissible leading norms remains unproved.

Use the fixed prime-form base from
[smooth leading-power selection](smooth-leading-power-selection.md).
Let `H=floor(sqrt(|Delta|/3))`, and suppose proper gcds with `N` are
handled when encountered. A reduced target form has leading norm
`A<=H`. Write

\[
 A=S R,
 \qquad S\text{ supported on the known base},
 \qquad \gcd(S,R)=1.
 \tag{1}
\]

A query is **admissible at residual cap B** if `R<=B`. The case `R=1`
is a fully smooth hit. For `R>1`, the residual norm is odd and coprime
to `2 Delta`, and its canonical key is

\[
 (R,\min(B_A\bmod2R,-B_A\bmod2R)),
 \tag{2}
\]

where `B_A` is the middle coefficient of the target form. Its
orientation is `epsilon=+1` or `-1` according to whether the original
root equals the canonical root or its negative.

For a target base prime `p`, define the exact integer thresholds

\[
 K_p=\min\{k:p^k>H\},\qquad
 L_p=\min\{l:p^l>H^2\}.
 \tag{3}
\]

For integer `B>=1`, put

\[
 K_B=B\,\mathrm{bitlength}(B).
 \tag{4}
\]

**Quota theorem.** Among more than `K_B L_p` admissible, distinct
exponents `k>=K_p` for a fixed target generator `g_p`, either a smooth
hit occurs, or two equal-key queries yield a relation strictly
diagonally dominant in the `p` column under column weights `log(q)`.
One successful quota for each base prime therefore produces a full
rational-rank relation family.

The exponents need not be consecutive. In a consecutive-power scan,
inadmissible queries do not count toward the quota and must still be
charged in the running time.

## 1. There are at most K_B residual keys

For an odd integer `R>1` coprime to `Delta`, a root of `Delta` modulo
an odd prime power has at most two choices. The Chinese remainder
theorem therefore bounds the root count modulo `R` by `2^omega(R)`.
The parity condition uniquely identifies the corresponding root
modulo `2R`. Passing to a root up to sign only reduces the count.

Since `2^omega(R)<=tau(R)`, the total number of keys, including a
possible smooth symbol, is at most

\[
 \sum_{R=1}^B\tau(R)
 =\sum_{d=1}^B\left\lfloor\frac Bd\right\rfloor
 \le B\sum_{d=1}^B\frac1d
 \le B\,\mathrm{bitlength}(B)=K_B.
 \tag{5}
\]

For the last inequality, group the positive integers into blocks
`[2^j,2^(j+1)-1]`. Each block contributes at most one to the harmonic
sum, and there are `bitlength(B)` nonempty blocks. This bound does not
require enumerating or factoring any residual norm.

## 2. A long gap or opposite orientation gives dominance

An admissible query at exponent `k` supplies

\[
 g_p^k=g^{v_k}h_{\kappa}^{\epsilon_k},
 \qquad
 \sum_q |v_{k,q}|\log q=\log S_k\le\log H.
 \tag{6}
\]

Here `h_kappa` is the canonical residual ideal class associated with
the key. For a nontrivial residual, normalize the exponent word to

\[
 u_k=\epsilon_k(k e_p-v_k),
 \qquad g^{u_k}=h_{\kappa}.
 \tag{7}
\]

Two equal-key queries `k_1<k_2` therefore give the genuine relation

\[
 r=u_{k_2}-u_{k_1}=\lambda e_p-v,
 \quad
 \lambda=\epsilon_{k_2}k_2-\epsilon_{k_1}k_1,
 \quad
 \sum_q|v_q|\log q\le2\log H.
 \tag{8}
\]

If the orientations differ, then
`|lambda|=k_1+k_2>=2K_p`, so `|lambda|log(p)>2log(H)`.
If the orientations agree and `k_2-k_1>=L_p`, the same strict
inequality follows from (3). In either case,

\[
 |r_p|\log p
   -\sum_{q\ne p}|r_q|\log q
 \ge |\lambda|\log p-\sum_q|v_q|\log q
 >0.
 \tag{9}
\]

A smooth query directly gives `r=k e_p-v_k`. Its corresponding
strict inequality uses `k log(p)>log(H)` and needs no collision.

## 3. Pigeonhole and full rank

Suppose there is no smooth hit and no pair covered by Section 2. For
each residual key, all of its queries have one orientation, and their
exponents lie in an interval of span at most `L_p-1`. Since the
exponents are distinct, that key occurs at most `L_p` times. Equation
(5) limits the total admissible count to `K_B L_p`, proving the quota
theorem.

Choose one dominant relation for each base prime. Multiplying column
`q` of their square relation matrix by the positive weight `log(q)`
makes it strictly row diagonally dominant. It is therefore nonsingular.
The rows are independent over the rationals. The weights are used
only in the proof; the selector does not need approximate logarithms.

The implementation stores the first normalized word `u_(kappa,0)`
for each key and retains differences `u_k-u_(kappa,0)`. Every pair
difference is an integer combination of two such stored-reference
differences. Consequently the selector's candidate span contains all
the dominant relations promised by the theorem, even if it did not
select those pairs explicitly.

## 4. Stored references and rational-basis retention preserve completeness

There is also an exact presentation interpretation. Temporarily add
one generator `h_kappa` per residual key, with relation
`u_k-h_kappa=0` for each partial query. The first relation for a key
has coefficient `-1` in its new generator, so eliminate that generator
over the integers. All remaining relations with that key become
`u_k-u_(kappa,0)=0`. Smooth relations are unchanged. Thus storing one
reference per key loses none of the torsion information in the full
partial-relation presentation.

Now let `L` be the integer lattice generated by every collision and
smooth relation, and let `L'` be generated by the final retained
independent original rows. Provided their rational spans agree, their
saturations agree:

\[
 S=(\mathbb Q L)\cap\mathbb Z^m
   =(\mathbb Q L')\cap\mathbb Z^m.
 \tag{10}
\]

For the homomorphism `phi:Z^m -> G` defined by the base classes,
both relation lattices lie in its kernel. The torsion subgroups of the
presented quotients are `S/L` and `S/L'`. Both have image exactly
`phi(S)` in the actual class group. Hence discarding rationally
dependent original rows preserves the useful-involution completeness
of [the two-primary relation extractor](../tools/relation_two_saturation.py).
It need not preserve the relation lattice itself or its index.
The quota conclusion applies once this final rational-span condition
has been established. Independence alone is not sufficient if a
streaming filter has deferred some candidate rows.

## 5. Cost and score conversion

For an `n`-bit input, `log H=O(n)` and `L_p=O(n)`. The sufficient
number of admissible queries per generator is

\[
 O(Bn\log(B+1)).
 \tag{11}
\]

Let `W(n)` be a uniform number of consecutive queries per generator
that guarantees the quota, or an earlier smooth hit or proper gcd.
Under the existing polynomial-size, generating factor-base contract,
the selector and exact relation processing have cost

\[
 W(n)\,\mathrm{poly}(n,\log W(n),\log(B+1)).
 \tag{12}
\]

Every rejected norm is included in `W(n)`. A bound
`W(n)<=exp(O(n/(log n)^alpha))` would transfer the same `alpha` to
the factoring reduction, with the existing preprocessing and
generation assumptions. Equation (11) alone is not such a bound.
For the deterministic complexity bound, residual keys can use an
ordered dictionary with logarithmic comparison cost. Retained row
entries have `O(n+log W)` bits; exact rational elimination has
polynomial bit cost in that size and the polynomial base dimension.

In particular, a statement about the global density of small residual
norms in the class group does not supply the required count along
each target power orbit. The missing hypothesis is now explicit:
accumulate the stated per-generator admissible quota within a uniform
improved query budget. No new score is claimed by this note.

## 6. Exact finite controls

[check_cofactor_collision_quota.py](../tools/check_cofactor_collision_quota.py)
checks root-key counts directly for small discriminants, exercises the
strict integer power threshold in (3), and verifies weighted dominance
using integer products for 1,000 bounded smooth-word controls. It also
checks the composite residual norm 15 with two distinct keys at
`N=959`, including both orientations of one key. The results are in
[cofactor-collision-quota.json](../evidence/cofactor-collision-quota.json).

The checker also runs 300 independent finite-abelian-group controls of
the final rational-span lemma. It constructs saturated lattices of
varying rank from unimodular row operations, adds dependent genuine
relations with different powers-of-two indices, and shuffles them.
The retained-row extractor agrees with the useful involutions in the
fully enumerated saturation image, allowing zero or one rejected top.
