# Single-direction screens for the nonconstant fourth-root primitives

## Result and scope

The [fourth-root construction](source-slab-fourth-root-directions.md)
provides ten polynomial endpoint-zero slab primitives. This checkpoint
tests the seven nonconstant primitives around the actual degree-18
five-letter reflection seed. Their degrees are

\[
3,3,3,4,5,5,5.
\]

Every tested nonzero one-direction perturbation has a recorded cost
above the unchanged reference at both fixed gauges. This is a finite
coefficient screen with numerical scale optimization. It is not a
minimum certificate for the full affine family or for other gauges.

## 1. Exact scoped families

The full
[source record](../evidence/source-slab-fourth-root-directions-h16.json)
contains three degree-one constant-support primitives followed by the
seven primitives tested here. The constant supports were already
covered in earlier screens. The
[exporter](../tools/fourth_root_source_subfamilies.py) retains primitive
indices 3 through 9 and all five logical \(X\)-placements.

The [unshifted export](../evidence/source-root-nonconstant-unshifted-directions.json)
contains 35 directions. The
[shifted export](../evidence/source-root-nonconstant-shifted-directions.json)
also includes every multiplication by \(z^s\) that keeps the complete
degree at most 18. Its dimension over both \(\mathbf F_2\) and
\(\mathbf F_4\) is

\[
5\bigl(3\cdot16+15+3\cdot14\bigr)=525.
\]

The independent replays recompute these complete-coefficient ranks.
Every direction has zero endpoint, so adding a scalar multiple to
the strict reference preserves every endpoint entry. The fourth-root
operation is used as a separate exact tensor construction; no claim
that it preserves the original source-operator family is required.

## 2. Completed searches

Both runs use paired energies \((0,10,6,2,8)\) and
\((0,12,6,0,8)\). The unchanged reference costs are approximately
`0.7263479477044683` and `0.8335422390094076`.

The binary reference and binary direction coefficients make field
Frobenius an exact symmetry. Coefficients \(1\) and \(\omega\)
represent all nonzero single-direction choices over \(\mathbf F_4\).

| Family | Directions | Coefficient-orbit queries | Profile classes, natural / two-zero | Best nonzero recorded cost, natural / two-zero |
| --- | ---: | ---: | --- | --- |
| Unshifted | 35 | 70 | 38 / 38 | `0.8999489497` / `0.9164822863` |
| All allowed shifts | 525 | 1,050 | 577 / 580 | `0.7948361364` / `0.8939564928` |

Profile-class counts include the unchanged reference. The unshifted
pilot is a subset of the shifted screen, so these are not 1,120
distinct coefficient choices. No two-direction or larger-support
perturbation search is claimed here.

The complete records are the
[unshifted screen](../evidence/source-root-nonconstant-unshifted-search.json)
and [shifted screen](../evidence/source-root-nonconstant-shifted-search.json).
The [independent replay](../evidence/source-root-nonconstant-search-replays.json)
reconstructs 96 retained literal witnesses and certifies 200 complete
Smith profiles using regular polynomial frames. It verifies 62,500
endpoint entries across the repeated references and witnesses.

## 3. Reproduction

```sh
python3 tools/fourth_root_source_subfamilies.py \
  evidence/source-slab-fourth-root-directions-h16.json \
  --unshifted-output /tmp/source-root-nonconstant-unshifted-directions.json \
  --shifted-output /tmp/source-root-nonconstant-shifted-directions.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_saturated_dual_directions.py \
  evidence/source-root-nonconstant-shifted-directions.json \
  --pair-queries 0 --random-queries 0 --seconds 1200 \
  --output /tmp/source-root-nonconstant-shifted-search.json
OPENBLAS_NUM_THREADS=1 python3 tools/check_saturated_dual_search.py \
  evidence/source-root-nonconstant-unshifted-search.json \
  evidence/source-root-nonconstant-shifted-search.json \
  --output /tmp/source-root-nonconstant-search-replays.json
```

The retained searches completed their requested schedules. A time-limited
rerun that stops earlier must keep its incomplete status. Numerical
roots do not certify a stronger factoring exponent.
