# Source audit for scalar factoring recovery in a fixed prime field

## Result and scope

The fixed-prime recovery dependency used by
[all-field Pascal crossings](all-field-pascal-crossings.md) is a theorem
in [PR #31](https://github.com/1337isnot1337/integer-factoring/pull/31),
at the pinned revision
[`b1ba452fa4a3d19fe6fea3822c53613f3b582104`](https://github.com/1337isnot1337/integer-factoring/tree/b1ba452fa4a3d19fe6fea3822c53613f3b582104).
Its source file is
[fixed-prime-factor-transfer.md](https://github.com/1337isnot1337/integer-factoring/blob/b1ba452fa4a3d19fe6fea3822c53613f3b582104/research/fixed-prime-factor-transfer.md).
The source is not a local file in this checkout. This note pins and
audits that dependency rather than silently treating it as integrated.

The audit found the stated scalar recovery interface and its linear
network-size accounting consistent with the pinned primary definitions.
The published checker was replayed from the exact source bytes. Its
complete output, including both checker hashes, matched the pinned
evidence. This is an ordinary proof audit and a finite replay; it is
not a new Lean proof or a factoring exponent improvement.

## 1. The exact theorem required

Fix a prime \(p\). Suppose one uniform deterministic finite-alphabet
algorithm evaluates every explicitly written \(P\)-vertex clique over
\(\mathbb F_p\), with \(\lambda\) bits on every edge and arbitrary
field-valued vertex tables, in time \(E_p(P,\lambda)\). This cost
includes reading the written input. The source proves that complete
factorization of \(2\le N<2^n\), and hence its smallest prime factor,
has cost

~~~math
\mathrm{poly}_p(n)\left(
  2^{C_p n/P}+
  E_p\!\left(P,\left\lfloor\frac{20000p^2n}{P^2}\right\rfloor\right)
\right),\qquad P\ge2,\quad P^{40}\le n.                 \tag{1}
~~~

The characteristic is fixed before the input grows. The theorem requires
only scalar evaluations. It does not assume an evaluator over polynomial
rings, reconstruction of an integer count, or a generic isolation oracle.

### Nonzero counts for every composite

For the correctness proof only, write \(N=\prod_q q^{e_q}\). At each
level \(r=p^j\le n-1\), count the nontrivial ordered \((p-1)\)-tuples
whose product \(A\) satisfies \(A^r\mid N\). Their number is

~~~math
\prod_{q\mid N}\binom{\lfloor e_q/r\rfloor+p-1}{p-1}-1.  \tag{2}
~~~

If \(e_{\max}\ge p\), choose the largest power of \(p\) at most
\(e_{\max}\). All quotients in (2) lie between zero and \(p-1\),
and at least one is positive. For \(1\le e<p\),
\(\binom{e+p-1}{p-1}\) has \(p\)-adic valuation exactly one.
Thus (2) is \(-1\) modulo \(p\).

The remaining predicate counts lexicographically least representatives
of nonconstant cyclic orbits of ordered \(p\)-factor tuples. It removes
the orbit represented by \((1,\ldots,1,N)\). A nonconstant orbit has
size exactly \(p\), since \(p\) is prime. Its count is

~~~math
\frac{\prod_{q\mid N}\binom{e_q+p-1}{p-1}
      -[N\text{ is a perfect }p\text{th power}]}p-1.       \tag{3}
~~~

This division is an integer counting identity; the evaluator counts
the representatives directly. If every \(e_q<p\), then (3) is
\(-1\) modulo \(p\) when there are at least two distinct prime factors.
For \(N=q^e\), it is \(e^{-1}-1\) modulo \(p\), which is nonzero
for \(2\le e<p\). Prime inputs have no witnesses in either predicate.

Every witness contains a proper divisor greater than one. A nonzero
prefix count has a nonzero binary child because the two child counts
sum to the parent. At most \(pn\) queries recover a complete tuple.
Recursive splitting has at most \(2n-3\) nodes. The algorithm tries
all power levels and the cyclic predicate; the unknown factorization
is used only in the proof of completeness.

## 2. Network construction and quantitative dependencies

The free-output multiplication construction is in the pinned
[parity predecessor, Sections 4--7](https://github.com/1337isnot1337/integer-factoring/blob/b1ba452fa4a3d19fe6fea3822c53613f3b582104/research/parity-factor-transfer.md).
Its primary source is Rudich and Rousseau's reference at
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
The following definitions and theorems were inspected directly:

| Pinned source | Interface checked |
|---|---|
| [MulNet.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean) | The setting types `mS`, deterministic messages `mV`, and the fixed-output digit test in `mW`. |
| [MulCount.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean) | `good_pinned`, `val_fC`, and the carry-chain multiplication identity. |
| [Cluster.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Cluster.lean) | `mC_row_le`, the per-gate setting bounds, and total incident width. |
| [EngineGroup.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean) | `good_exists` for an arbitrary network fitting the supplied width and total bounds, with \(G=3P^8\). |

Removing the fixed-output digit comparisons exposes the uniquely forced
carry output. The proof of `good_pinned` uses equality messages and the
initial zero carry; it does not use the removed output comparisons.
`val_fC` identifies the coefficients with the integer product convolution.
Unrolling the carry recurrence therefore gives every bit of the full
product. The final carry vanishes because two \(C\)-bit inputs have
product below \(2^{2C}\). This does not require extending a theorem
whose stated target \(N\) is restricted to \(C\) bits.

For intended input width \(w\), the predecessor chooses a power of two
\(S\) with capacity \(C=4S^4\), where \(w\le C\le16w\). Put
\(A=2S\), \(f=8S^2+1\). The internal incident-width sum and the
setting-exponent sum are each at most \(10A^2f\). The \(4C\) exposed
bit ports and their incident edges add at most \(20C\) total mass.
Consequently one multiplier has total mass at most

~~~math
20A^2f+20C\le200C\le3200w,                            \tag{4}
~~~

and every vertex has mass at most \(512w^{3/4}\). Here mass means
incident bit width plus a bound on the base-two logarithm of its number
of local settings. Port-degree bounds include the stated copy edges.

For the root predicate, take \(m=\lceil n/r\rceil\), multiply the
\(p-1\) root coordinates, then perform \(j\) stages of \(p\)th
powering and multiply by a cofactor. The intended multiplier widths
sum to at most \(2p^2n\), since their stage widths form a geometric
series and \(rm<2n\). The cyclic predicate uses \(p-1\) successive
multiplications, of total intended width \(p(p-1)n/2\). All full
integer products are retained, so no modular truncation creates extra
witnesses.

Binary copy chains and deterministic three-state comparisons have
linear total mass for fixed \(p\). Strict comparison with each nonidentity
rotation enforces the orbit representatives in (3). The resulting
networks satisfy

~~~math
\sum_v(d_v+b_v)\le6500p^2n,
\qquad \max_v(d_v+b_v)\le1024p\,n^{3/4}.               \tag{5}
~~~

Greedy mass balancing into \(G=3P^8\) clusters makes every cluster's
mass at most
\(Q=\lceil6500p^2n/G\rceil+1024p\lceil n^{3/4}\rceil\).
The actual `good_exists` theorem gives pair width
\(\lfloor(P\cdot6500p^2n+GQ)/P^3\rfloor\) and at most
\(\lfloor2G/P\rfloor\) clusters per group. Under \(P^{40}\le n\),
these imply \(P^2\lambda\le20000p^2n\) and a per-group setting
exponent \(H\) with \(PH\le20000p^2n\).

Searching all \(P^G\) groupings and writing every clique table costs
\(2^{O_p(n/P)}\). For each entry, enumerate the original settings in
its group and test their deterministic messages. There is no additional
enumeration of internal channels. All arithmetic, copy, and comparison
settings extend each witness uniquely. This establishes the count
preservation and the written-input accounting in (1).

## 3. Connection to effective crossings

The local [scale-one transfer, Sections 1--3](scale-one-modular-transfer.md)
constructs a sum-product evaluator over a fixed finite field from
effective scale-one crossing pieces. Its circuit count is

~~~math
2^{e_\ell P^2\lambda}
\bigl(P(\lambda+1)\bigr)^{C_0P^2},
\qquad e_\ell=O(\rho^\ell),\quad P=K_\ell.             \tag{6}
~~~

All field operations have constant bit cost for fixed \(p\), and its
structural computation has fixed polynomial overhead. Thus it gives a
scalar evaluator of exactly the kind required in (1). This use does
not reinterpret the source's integer machine as a finite-field machine;
the local note constructs the sum-product evaluator separately.

Choose the largest level with \(P^{40}\le n\). The source schedule
satisfies \(n<P^{400}\) and \(2^\ell=\Omega(\log n)\).
For \(\rho\ge1/2\), the table construction term \(n/P\) and the
logarithmic overhead in (6) are \(O(n\rho^\ell)\). Combining (1)
and (6) gives

~~~math
T(n)=\exp\!\left(O\!\left(
\frac{n}{(\log n)^{\log_2(1/\rho)}}\right)\right).       \tag{7}
~~~

This implication still needs an effective crossing family at the claimed
\(\rho\). The all-field Pascal note supplies a constructor and an
explicit pressure criterion. It supplies no witness with
\(\rho<1/\sqrt2\), so (7) yields no new achieved factoring exponent.

## 4. Exact replay and provenance

At the pinned PR #31 revision, run

```sh
python3 tools/check_fixed_prime_factor_transfer.py \
  --output evidence/fixed-prime-factor-transfer.json
```

The checker imports the adjacent parity checker. Both were downloaded
from that exact revision without edits. Their source hashes are:

| File | SHA-256 |
|---|---|
| `tools/check_fixed_prime_factor_transfer.py` | `e1c49201238caad9f7c2664ecb6d9117f9beee0f72854473bf4e7b0d16f771c6` |
| `tools/check_parity_factor_transfer.py` | `ec1b7e496cdd51bf52885366b9b0d41d33a57afa0a843218c990a5f0425e5b5c` |

The replay matched the pinned evidence, including 262,136 completeness
cases, 4,092 complete factorizations, 290,536 nonzero-prefix queries,
49,679 composed-network mass controls, 248 constructor-cutoff controls,
and 5,460 deterministic comparison traces over characteristics 2, 3,
5, and 7. The comparison is recorded in
[fixed-prime-transfer-source-audit.json](../evidence/fixed-prime-transfer-source-audit.json).

Tuple enumeration is the checker's finite reference oracle. It is not
an implementation of the fast clique evaluator and does not prove an
asymptotic runtime by experiment. The all-input argument is the explicit
witness construction, unique-extension proof, and network accounting
audited above.
