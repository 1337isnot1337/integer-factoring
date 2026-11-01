# Affine combinations of filtered sixteen-letter reflection frames

## Result and scope

This checkpoint constructs and checks 211 distinct degree-96 polynomial
crossing tensors over `F4`, obtained from five filtered rank-one
directions on the sixteen-letter reflection seed. The completed search
found no improvement over its original balanced value
`0.716653659728187`. Its best asymmetric support value was
`1.0024570646034154`, above the required threshold one.

The exact construction and complete tensor reconstruction are separate
from the bounded numerical optimization. Neither the optimizer outcomes
nor the finite search establish a lower bound for the larger affine
family.

## 1. A three-parameter affine family from one direction

Let `G=u v^T` be a square-zero matrix in Bell coordinates, so `v^T u=0`.
The complete coupled physical action with parameter `f` is

```math
T(f)=T_0+fT_1+f^2T_2+f^3T_3+f^4T_4,\qquad T_0=T.
```

The [coupled-frame construction](coupled-rational-physical-frames.md)
retains all four physical-leg interactions before checking polynomiality.
For `f=w^-k`, the three nonzero field scalars satisfy `c^4=c`. Hence

```math
T(cf)-T
=c\,U_1+c^2U_2+c^3U_3,
\quad
U_1=fT_1+f^4T_4,\quad
U_2=f^2T_2,\quad
U_3=f^3T_3.                                      \tag{1}
```

The three-by-three coefficient matrix `(c,c^2,c^3)`, with
`c=1,omega,omega^2`, is invertible over `F4`. If all three complete
frames are polynomial, then each `U_j` is polynomial. Because all
three frames have the strict crossing endpoint, every `U_j` vanishes
at `w=1`. Therefore

```math
T+\lambda_1U_1+\lambda_2U_2+\lambda_3U_3,
\qquad \lambda_j\in F4,                           \tag{2}
```

is a polynomial tensor with the same strict endpoint.

Each term of positive `f` degree lowers the polynomial degree below
96. Thus the leading coefficient of `T` survives every combination
(2), and its exact degree remains 96. The
[paired-trace identity](reflection-paired-trace-invariant.md)
also survives (2), giving scalar polynomial content one.

The search tests all three scalar frames by their complete
coefficients before using (1). This is essential at higher poles:
admissibility of one scalar multiple alone does not justify extracting
the three polynomial directions.

## 2. The finite direction and pole catalogue

The [search archive](../evidence/filtered-reflection16-affine-search.json)
stores the full vectors `u,v` for five projective directions selected
from the natural physical filtration. These are explicit matrices on
the full sixteen-letter alphabet. No Gaussian one-particle
parameterization is assumed.

For directions 0, 2, and 3, all three nonzero scalar frames are
polynomial exactly at the checked positive poles 1 through 8.
For directions 1 and 4, the common checked range is 1 through 16.
The declared affine screen requests poles
`1,2,4,8,16,32,64` and retains only those accepted by the complete
three-scalar test.

There are 22 accepted direction/pole pairs. Testing the 63 nonzero
coefficient triples in (2) gives 1,386 parameter choices. After exact
tensor-hash deduplication, 211 tensors remain. Some nonzero triples
produce the original tensor because an extracted direction vanishes.
All such coincidences are retained in the parameter records rather
than interpreted as distinct constructions.

## 3. Completed numerical search

The worker began at `2026-10-31 22:41:53.591990 UTC` and finished at
`23:01:07.135661 UTC`, after `1153.543345591` measured seconds.
It completed its declared finite screen and refinement budgets.

| Quantity | Completed value |
|---|---:|
| Affine parameter choices | 1,386 |
| Exact duplicate tensors | 1,175 |
| Distinct complete tensors | 211 |
| Distinct natural two-cut profile vectors | 38 |
| Balanced refinements | 12 |
| Asymmetric refinements | 12 |
| Resolved objective calls | 899 |
| Unresolved objective calls | 0 |

Every natural profile uses both complete weighted cuts at zero and
infinity, the actual polynomial degree, and one common pair of physical
types. The initial screen optimizes its shared scale and Gibbs types.
The subsequent refinements optimize the physical gauges and common
types. The twelve balanced loops reported optimizer convergence. Each
asymmetric loop reached its 55-second trial deadline and retained its
best resolved objective call.

The best balanced tensor is exactly the original seed: direction zero,
pole one, coefficient triple `(0,0,1)`. Its replayed value is
`0.7166536597281871`. The best asymmetric tensor uses direction zero,
pole one, triple `(1,1,0)`, with digest

```text
8d851135d6b1917b6243b109bfaf05ae030eba66f5c4ccfff8a0ec26d5c15bbd
```

Its support value replays as `1.0024570646034154`.
These are resolved evaluations at recorded parameters. They are not
certified continuous minima.

## 4. Complete reconstruction and finite-place controls

The [replay checker](../tools/check_filtered_reflection16_affine_search.py)
checks every recorded search-source hash, reconstructs all 211 tensors,
and binds every refinement to its complete polynomial input. This
checks 13,828,096 polynomial entries. It also independently evaluates
the two recorded winners with fresh complete local profiles.
The [replay evidence](../evidence/filtered-reflection16-affine-replay.json)
records these checks.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/check_filtered_reflection16_affine_search.py \
  --output /tmp/filtered-reflection16-affine-replay.json
```

A separate [finite-place screen](../tools/screen_filtered_reflection16_affine_places.py)
checks both cuts at the two nontrivial `F4` points. Its
[complete record](../evidence/filtered-reflection16-affine-finite-places.json)
is bound to the final search archive:

| Ranks at `omega` | Ranks at `omega^2` | Tensors |
|---|---|---:|
| `(240,240)` | `(256,256)` | 48 |
| `(256,256)` | `(240,240)` | 48 |
| `(256,256)` | `(256,256)` | 115 |

Independent controls check 160 polynomial evaluations, four dense
256-by-256 rank calculations, and 1,048,576 reconstructed polynomial
entries. These are point-rank statements. A useful mixed-place
certificate still requires every local Smith order and the corresponding
weighted pressure.

The archive makes all five generators, pole tests, complete tensor
hashes, profile classes, optimizer parameters, and actual timings
available. Combinations of directions from different generators or
different pole orders are a separate search family.
