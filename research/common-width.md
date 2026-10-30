# Checked common widths and one-level contraction

`formal/CommonWidth.lean` checks the integer arithmetic in equations
(8)--(9) of [the integration construction](integration.md). The statements
hold for arbitrary natural-number inputs under the stated hypotheses;
they are not limited to the finite Python examples.

This is an arithmetic interface for the proposed evaluator. The file does
not construct the row/column networks, prove their ordered-total bounds,
find a grouping, or implement the factoring machine.
The source reference remains commit
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1` of
`IsaacRudich/SubexponentialDeterministicFactoring`, as pinned in
`upstream.lock.json`. The formulas here use the rational `281/312` envelope.

## A common width for every branch

For positive `K,D`, put

```text
Umax = (281*G^4*lambda)//312 + G^3*lambda + 3*G^4
L = (D*Umax + G^4*lambda)//(D*K^2).
```

`common_width_iff` proves that `L` is exactly the largest integer `c`
satisfying `D*K^2*c<=D*Umax+G^4*lambda`. If a network satisfies
`U<=Umax`, `P<=G^3*lambda`, and a grouped pair satisfies
`D*K^2*c<=D*U+G*P`, `grouped_width_le_common` proves `c<=L`.
The resulting width depends only on `K,D,G,lambda`, so every piece branch
may use the same padding. This conclusion is conditional on the supplied
network and grouping inequalities, including for the stronger grouping
theorem that is currently handwritten.

`common_width_volume` proves the integral inequality

```text
312*D*(K^2*L) <= 281*D*G^4*lambda + 312*D*G^3*lambda
                  + 312*G^4*lambda + 936*D*G^4.
```

For positive `G,D`, this is precisely equation (9). Both floors round
down, as required by the pair guarantee and the volume upper bound.

## Executable schedule and replay

`tools.parameters.width_step` and `width_schedule`
implement these integer expressions. `schedule` returns levels from
bottom to top; `width_schedule` consumes that order and returns steps in
evaluation order, from top to bottom. Each step records the input/output
widths, volumes, and ordered-total upper bound. It computes parameters,
without constructing or evaluating clique tables.

```python
from tools.parameters import schedule, width_schedule

steps = width_schedule(schedule(4), top_width=2**4096 + 117)
assert [step.level.d for step in steps] == [3, 2, 1, 0]
assert all(a.output_volume == b.input_volume for a, b in zip(steps, steps[1:]))
```

Run the Python boundary and composition checks with:

```sh
python3 -m unittest discover -s tests -p test_common_width.py -v
bash tools/check_formal.sh
```

The tests cover exact/nonexact divisions on small synthetic levels and
a four-level composition with a 4,097-bit initial width. These
are arithmetic certificates, not factoring or timing benchmarks.

The next integration step is to establish the total, degree, and grouping
hypotheses for the actual transformed networks, and connect the padded
tables to the value-preserving generalized `nextT` operation.
