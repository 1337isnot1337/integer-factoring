#!/usr/bin/env bash
set -euo pipefail

# Fetch only the public, pinned reference; all new proofs stay in this checkout.
factoring_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
factoring_reference="$factoring_root/.lake/upstream-reference"
factoring_revision="c91cc3f4947caca8977e81a5d14244f73b5ae2b1"
factoring_logs="$factoring_root/build/formal"

mkdir -p "$factoring_root/.lake" "$factoring_logs"
if [[ ! -d "$factoring_reference" ]]; then
  git clone --filter=blob:none --no-checkout \
    https://github.com/IsaacRudich/SubexponentialDeterministicFactoring.git \
    "$factoring_reference"
  git -C "$factoring_reference" checkout --detach "$factoring_revision"
fi
if [[ "$(git -C "$factoring_reference" rev-parse HEAD)" != "$factoring_revision" ]]; then
  echo "Reference checkout differs from upstream.lock.json; refusing to use it." >&2
  exit 1
fi
if [[ -n "$(git -C "$factoring_reference" status --porcelain --untracked-files=no)" ]]; then
  echo "Reference checkout has tracked changes; refusing to use it." >&2
  exit 1
fi

cd "$factoring_reference"
lean --version | tee "$factoring_logs/toolchain.log"
lake exe cache get
lake build Solution > "$factoring_logs/upstream.log" 2>&1
for factoring_source in "$factoring_root"/formal/*.lean; do
  factoring_name="$(basename "$factoring_source" .lean)"
  lake env lean "$factoring_source" 2>&1 | tee "$factoring_logs/$factoring_name.log"
done
python3 "$factoring_root/tools/check_axioms.py" "$factoring_logs"
