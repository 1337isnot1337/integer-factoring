# A complete finite operator catalogue above a power of two

The [above-power reduction](above-power-identity-operator-reduction.md)
gives exact bounded cores at \(n=2^k+1\). The first completed search
uses \(n=33\), ramification \(P=1\), and the three terminal
directions \((30,31,32)\). The zero core additionally contains the
unchanged physical identity direction zero and auxiliary directions
\((0,1,31,32)\). All remaining Gaussian spectators and both physical
paired-type partitions are included.

The [driver](../tools/search_above_power_terminal_operators.py)
enumerates the same specified F4 projector catalogue used in the
three-direction operator searches: 674 idempotent matrices generate
2,660 paired occupation-sector patterns, with 16 color/top/endpoint
choices for each pattern. The total is **42,560 kernels**. Their
nonconstant entries use the Laurent support \(\{-1,0\}\); the
complete degree cost is included at infinity. This is a complete
enumeration of that defined catalogue, not of all three-direction
operator kernels.

All candidates completed and produced **1,111 distinct exact
profile pairs**. Each pair was evaluated on a 72-point scale/charge
grid. The 30 best grid profiles received a local continuous
refinement over the recorded bounds. The best critical pressure was
**0.032265836922881874**, attained by the Gaussian profile. Its
refined scale was approximately 0.1743004, with charge within
\(5\cdot10^{-5}\) of zero. No threshold crossing or same-gauge
cost improvement was detected. The smallest reported relative change,
about \(-7.1\cdot10^{-15}\), is floating-point noise.

The complete
[record](../evidence/above-power-operator-projectors-n33.json)
contains every profile class, its multiplicity and representative
coefficient row, grid and refinement results, and source hashes.
Five retained profile pairs were replayed from literal polynomial
kernels. The best zero profile was also replayed by expanding the
complete square active tensor before applying the Weyl transform.
The separate construction checker supplies the original physical
two-cut and endpoint controls.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/search_above_power_terminal_operators.py --n 33 --P 1 --one-height 1 --two-height 1 --top-height 1 --refinements 30 --output /tmp/above-power-operator-projectors-n33-replay.json
```

The exact profiles establish the evaluated finite family. The
floating pressure optimization supplies search evidence; it does not
prove a global obstruction for other Laurent supports, longer
intervals, different endpoint operators, or arbitrary gauge choices.
