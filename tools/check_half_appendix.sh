#!/usr/bin/env bash
set -euo pipefail

# Build the pinned Appendix B development independently of the factoring
# reference: the two upstream projects currently use different Lean versions.
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$root/.lake/appendix-source"
revision="16a5ce41304fcdaea4e0be7c132e1c361dc26da5"
logs="$root/build/half-appendix"
mkdir -p "$root/.lake" "$logs"
if [[ ! -d "$source_dir" ]]; then
  git clone --filter=blob:none --no-checkout https://github.com/anthropics/formal-math.git "$source_dir"
  git -C "$source_dir" checkout --detach "$revision"
fi
if [[ "$(git -C "$source_dir" rev-parse HEAD)" != "$revision" ]] ||
   [[ -n "$(git -C "$source_dir" status --porcelain --untracked-files=no)" ]]; then
  echo "Appendix source is not the clean pinned revision" >&2
  exit 1
fi
cd "$source_dir/eth"
lean --version | tee "$logs/toolchain.log"
lake exe cache get > "$logs/cache.log" 2>&1
lake build GeneralPieces.Main > "$logs/source-build.log" 2>&1

# The stage is below the source project so Lean can set a separate root while
# still finding the source project's built imports through `lake env`.
stage="$(mktemp -d "$source_dir/eth/.lake/half-bridge.XXXXXXXX")"
cp "$root/formal/HalfPad.lean" "$stage/HalfPad.lean"
cp "$root/formal-appendix/ExplicitPad.lean" "$stage/ExplicitPad.lean"
lake env bash -c '
  set -euo pipefail
  export LEAN_PATH="$1:$LEAN_PATH"
  lean --root="$1" -o "$1/HalfPad.olean" "$1/HalfPad.lean"
  lean --root="$1" "$1/ExplicitPad.lean"
' _ "$stage" > "$logs/bridge.log" 2>&1
python3 - "$logs/bridge.log" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text()
reports = re.findall(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", text, re.S)
expected = "IntegerFactoring.HalfAppendix.exists_explicit_padded_evaluator"
if "error:" in text or "sorryAx" in text or not any(name == expected for name, _ in reports):
    raise SystemExit("Missing or failed Appendix bridge axiom report")
allowed = {"propext", "Classical.choice", "Quot.sound"}
for name, axioms in reports:
    actual = {a.strip() for a in axioms.split(",") if a.strip()}
    if actual - allowed:
        raise SystemExit(f"Unexpected axioms for {name}: {sorted(actual)}")
print(f"Checked {len(reports)} axiom reports; logs: {sys.argv[1]}")
PY
