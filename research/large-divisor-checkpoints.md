# Exact divisor maxima through 32,768 input bits

## Completed checkpoints

The certified arbitrary-bound oracle and the independent prime-order
enumerator agree at both of these larger power-of-two bounds:

| Input bits `n` | `ceil(log2 D(n))` | Prime support of least maximizer | Main feasible profiles | Main search nodes | Independent terminal profiles |
|---:|---:|---:|---:|---:|---:|
| 16,384 | 1,398 | 1,354 | 5,825 | 261,598 | 38,707 |
| 32,768 | 2,568 | 2,511 | 11,109 | 528,359 | 199,393 |

Here `D(n)=max_(1<=m<2^n) tau(m)`. The complete exact values,
least maximizers, exponent vectors, parameter mesh, final windows,
and benefit certificates are in
[divisor-maximum-16384.json](../evidence/divisor-maximum-16384.json)
and
[divisor-maximum-32768.json](../evidence/divisor-maximum-32768.json).
The independent replay logs are
[divisor-audit-16384.json](../evidence/divisor-audit-16384.json)
and
[divisor-audit-32768.json](../evidence/divisor-audit-32768.json).

These checkpoints were computed and independently replayed before
the CLI regression below. The consolidation tool imports those
completed records and their hashes; it does not repeat their large
searches. The mathematical algorithm and certificate are described
in [the certified divisor oracle](certified-divisor-oracle.md).
This is exact preprocessing evidence, not a faster factoring result.

The 65,536-bit main search is a separate, larger experiment. It is
not included in this completed checkpoint because its independent
prime-order replay had not finished when these files were frozen.

## Decimal conversion fix

The least maximizers in the two completed records contain 4,933
and 9,865 decimal digits. The exclusive bounds have the same
respective decimal lengths. Python's default decimal conversion
guard therefore blocks both large `--below` arguments and JSON
serialization from `--bits` on interpreters with that guard.

The small accompanying patch changes the CLI entry point to call
`sys.set_int_max_str_digits(0)` when that API is available. The
call occurs before argument parsing. It permits the deliberately
arbitrary-size exact integer inputs and outputs of this tool.
The search, pruning, certificate arithmetic, and imported API
behavior are unchanged.

The patch is
[certified-divisor-cli.patch](../patches/certified-divisor-cli.patch).
The change is already applied in this checkout. The normal commands support
the larger records:

```sh
python3 tools/certified_divisor_oracle.py --bits 16384 --output result.json
python3 tools/certified_divisor_oracle.py --below LARGE_EXCLUSIVE_BOUND --output result.json
```

## Focused regression and reproduction

Run the completed-record and CLI controls with

```sh
python3 tools/check_large_divisor_checkpoint.py --bits 16384 32768
```

The checker constructs guarded and unguarded copies of the CLI in
a temporary directory. It uses each independently verified result
as the solver return value, checks that the CLI passes the exact
parsed bound to that solver, and compares the entire serialized
JSON output to the checkpoint. This isolates the decimal-conversion
boundary without rerunning the number-theoretic search.

For both input sizes, the controls verify:

- the unguarded `--below` path rejects the long decimal argument;
- the guarded `--below` path parses it and serializes the full result;
- the unguarded `--bits` path reaches output and fails at decimal
  serialization;
- the guarded `--bits` path serializes the same complete exact result.

The child processes clear any inherited `PYTHONINTMAXSTRDIGITS`
override. The fixture temporarily permits large integers only
while loading its verified result, then restores the default
guard before invoking the actual CLI entry point. Thus the
unguarded negative controls exercise the intended failure modes.

The combined evidence is
[large-divisor-checkpoint.json](../evidence/large-divisor-checkpoint.json).
It records source and audit hashes, the normalized patched source
hash, both failure controls, and both successful output comparisons.
Timings in the imported audit logs are observations, not part of
the proof or a claimed asymptotic speedup.
