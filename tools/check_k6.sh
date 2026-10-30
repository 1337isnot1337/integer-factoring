#!/usr/bin/env bash
set -euo pipefail

# Rebuild the variant from a fresh archive. Never reset an existing checkout.
factoring_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
factoring_cache="$factoring_root/.lake/k6-reference"
factoring_packages="$factoring_root/.lake/k6-packages"
factoring_reference="$factoring_cache/upstream.git"
factoring_revision="c91cc3f4947caca8977e81a5d14244f73b5ae2b1"
factoring_url="https://github.com/IsaacRudich/SubexponentialDeterministicFactoring.git"
factoring_patch="$factoring_root/variants/k6/reference.patch"
factoring_audit="$factoring_root/evidence/k6_source_audit.json"
factoring_logs="$factoring_root/build/k6"
factoring_build=""

mkdir -p "$factoring_cache" "$factoring_logs"
exec {factoring_lock_fd}>"$factoring_cache/verify.lock"
if ! flock --nonblock "$factoring_lock_fd"; then
  echo "Another k6 verification is using this package cache." >&2
  exit 1
fi

factoring_cleanup() {
  factoring_status=$?
  if [[ -n "$factoring_build" ]]; then
    if [[ "$factoring_status" == 0 && "${K6_KEEP_BUILD:-0}" != 1 ]]; then
      # This directory was created by this invocation. rm does not follow the
      # .lake/packages symlink, so the shared dependency cache remains intact.
      rm -rf -- "$factoring_build"
    else
      echo "Verification source tree retained at $factoring_build" >&2
    fi
  fi
  exit "$factoring_status"
}
trap factoring_cleanup EXIT

python3 - "$factoring_root/upstream.lock.json" "$factoring_patch" "$factoring_audit" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

lock = json.loads(Path(sys.argv[1]).read_text())
expected = {
    "repository": "https://github.com/IsaacRudich/SubexponentialDeterministicFactoring",
    "commit": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1",
    "lean_toolchain": "leanprover/lean4:v4.35.0",
    "mathlib_commit": "56a005855ea17e0b0958f73ec801ae140b0a1679",
}
if any(lock.get(key) != value for key, value in expected.items()):
    raise SystemExit("upstream.lock.json differs from the verified k6 inputs")
patch = Path(sys.argv[2])
if not patch.is_file() or patch.is_symlink():
    raise SystemExit("Missing regular file variants/k6/reference.patch")
audit = json.loads(Path(sys.argv[3]).read_text())
if audit["reference_commit"] != expected["commit"]:
    raise SystemExit("The source audit has a different reference commit")
if audit["reference_patch"]["path"] != "variants/k6/reference.patch":
    raise SystemExit("The source audit names a different reference patch")
if hashlib.sha256(patch.read_bytes()).hexdigest() != audit["reference_patch"]["sha256"]:
    raise SystemExit("reference.patch differs from the audited patch hash")
PY

if [[ ! -e "$factoring_reference" ]]; then
  git clone --bare --filter=blob:none "$factoring_url" "$factoring_reference"
fi
if [[ "$(git --git-dir="$factoring_reference" rev-parse --is-bare-repository)" != true ||
      "$(git --git-dir="$factoring_reference" remote get-url origin)" != "$factoring_url" ]]; then
  echo "The public reference cache has an unexpected repository or origin." >&2
  exit 1
fi
if ! git --git-dir="$factoring_reference" cat-file -e "$factoring_revision^{commit}"; then
  git --git-dir="$factoring_reference" fetch --depth=1 origin "$factoring_revision"
fi
if [[ "$(git --git-dir="$factoring_reference" rev-parse "$factoring_revision^{commit}")" != "$factoring_revision" ]]; then
  echo "The public reference cache does not contain the pinned commit." >&2
  exit 1
fi

factoring_build="$(mktemp -d "$factoring_cache/build.XXXXXXXX")"
git --git-dir="$factoring_reference" archive "$factoring_revision" | tar -x -C "$factoring_build"
# Without its own Git root, apply can interpret this archive as a subdirectory
# of the parent checkout and silently skip every patch path.
git -C "$factoring_build" init --quiet
if [[ "$(git -C "$factoring_build" rev-parse --show-toplevel)" != "$factoring_build" ]]; then
  echo "The fresh archive is not an isolated Git working tree." >&2
  exit 1
fi
git -C "$factoring_build" apply --check "$factoring_patch"
git -C "$factoring_build" apply "$factoring_patch"

# Verify the actual patched bytes before copying proofs or starting any build.
# A successful exit from git apply alone is not evidence that a patch applied.
python3 - "$factoring_audit" "$factoring_build" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

audit = json.loads(Path(sys.argv[1]).read_text())
build = Path(sys.argv[2])
changed = audit["changed_existing_lean_files"]
added = audit["added_reference_lean_files"]
if len(changed) != audit["changed_existing_lean_file_count"]:
    raise SystemExit("Inconsistent changed-source count in the audit")
if len(added) != audit["added_reference_lean_file_count"]:
    raise SystemExit("Inconsistent added-source count in the audit")
if len(changed) + len(added) != audit["reference_patch"]["changed_file_count"]:
    raise SystemExit("The audited source list does not cover the entire patch")
rows = changed + added + audit["protected_files"]
seen = set()
for row in rows:
    relative = Path(row["path"])
    if relative.is_absolute() or ".." in relative.parts or str(relative) in seen:
        raise SystemExit(f"Invalid or repeated audited source path: {relative}")
    seen.add(str(relative))
    source = build / relative
    if not source.is_file() or source.is_symlink():
        raise SystemExit(f"Missing regular audited source: {relative}")
    actual = hashlib.sha256(source.read_bytes()).hexdigest()
    if actual != row["variant_sha256"]:
        raise SystemExit(f"Patched source hash mismatch: {relative}")
print(f"Verified {len(changed)} changed, {len(added)} added, and "
      f"{len(audit['protected_files'])} protected source hashes against the audit.")
PY

factoring_check_fixed_inputs() {
  local factoring_path
  for factoring_path in Challenge.lean FactoringAux/Definitions.lean lean-toolchain lake-manifest.json; do
    if ! git --git-dir="$factoring_reference" show "$factoring_revision:$factoring_path" |
        cmp -s - "$factoring_build/$factoring_path"; then
      echo "Protected reference input changed: $factoring_path" >&2
      return 1
    fi
  done
}
factoring_check_fixed_inputs

# Copy only paths tracked by this repository's index. This includes staged
# additions, allowing a clean verification before the first commit.
python3 - "$factoring_root" "$factoring_build" "$factoring_logs/sources.json" <<'PY'
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

root, build, report = map(Path, sys.argv[1:])
prefix = Path("variants/k6/FactoringK6")
tracked = subprocess.check_output(
    ["git", "-C", str(root), "ls-files", "-z", "--", str(prefix)]
).decode().split("\0")
sources = sorted(Path(name) for name in tracked if name.endswith(".lean"))
if prefix / "Certified.lean" not in sources:
    raise SystemExit("Stage or commit variants/k6/FactoringK6/Certified.lean before verification")
hashes = {}
for relative in sources:
    source = root / relative
    if not source.is_file() or source.is_symlink():
        raise SystemExit(f"Proof source is not a regular file: {relative}")
    target = build / "FactoringK6" / relative.relative_to(prefix)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    hashes[str(relative)] = hashlib.sha256(source.read_bytes()).hexdigest()
config = build / "lakefile.toml"
parsed = tomllib.loads(config.read_text())
if any(lib.get("name") == "FactoringK6" for lib in parsed.get("lean_lib", [])):
    raise SystemExit("The reference patch must leave FactoringK6 library registration to this script")
with config.open("a") as stream:
    stream.write('\n[[lean_lib]]\nname = "FactoringK6"\nglobs = ["FactoringK6.+"]\n')
patch = root / "variants/k6/reference.patch"
report.write_text(json.dumps({
    "reference_commit": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1",
    "patch_sha256": hashlib.sha256(patch.read_bytes()).hexdigest(),
    "source_sha256": hashes,
}, indent=2) + "\n")
PY

# The cache can be seeded with a symlink to an existing package directory.
# Sources are checked before reuse; no dependency checkout is reset.
mkdir -p "$factoring_packages" "$factoring_build/.lake"
ln -s "$factoring_packages" "$factoring_build/.lake/packages"
factoring_check_packages() {
  python3 - "$factoring_build/lake-manifest.json" "$factoring_packages" "$1" <<'PY'
import json
from pathlib import Path
import subprocess
import sys

manifest = json.loads(Path(sys.argv[1]).read_text())
packages = Path(sys.argv[2])
complete = sys.argv[3] == "complete"
for package in manifest["packages"]:
    directory = packages / package["name"]
    if not directory.exists() and not complete:
        continue
    actual = subprocess.check_output(
        ["git", "-C", str(directory), "rev-parse", "HEAD"], text=True
    ).strip()
    changes = subprocess.check_output(
        ["git", "-C", str(directory), "status", "--porcelain", "--untracked-files=no"],
        text=True,
    ).strip()
    if actual != package["rev"] or changes:
        raise SystemExit(f"Dependency {package['name']} differs from the pinned clean source; refusing to reuse it")
PY
}
factoring_check_packages partial

cd "$factoring_build"
lean --version | tee "$factoring_logs/toolchain.log"
case "$(cat "$factoring_logs/toolchain.log")" in
  *"version 4.35.0,"*) ;;
  *) echo "Lean 4.35.0 is required." >&2; exit 1 ;;
esac
lake exe cache get 2>&1 | tee "$factoring_logs/cache.log"
factoring_check_packages complete
lake build Solution FactoringK6 2>&1 | tee "$factoring_logs/build.log"

# Run the endpoint source again, even when Lake considers its olean current.
lake env lean FactoringK6/Certified.lean 2>&1 | tee "$factoring_logs/Certified.log"
factoring_check_fixed_inputs
factoring_check_packages complete

python3 - "$factoring_build/FactoringK6/Certified.lean" "$factoring_logs" <<'PY'
import json
from pathlib import Path
import re
import sys

source = Path(sys.argv[1]).read_text()
logs = Path(sys.argv[2])
text = (logs / "Certified.log").read_text()
build = (logs / "build.log").read_text()
allowed = {"propext", "Classical.choice", "Quot.sound"}
expected = set(re.findall(r"^\s*#print\s+axioms\s+([A-Za-z0-9_.']+)\s*$", source, re.MULTILINE))
if not expected:
    raise SystemExit("Certified.lean must print the axioms of its final theorem")
if any(marker in text + build for marker in ("sorryAx", "uses 'sorry'", "error:")):
    raise SystemExit("An incomplete or failed proof appears in the verification logs")
reports = re.findall(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", text, re.MULTILINE)
reports += [(name, "") for name in re.findall(r"'([^']+)' does not depend on any axioms", text)]
found = {}
for declaration, values in reports:
    axioms = {name.strip() for name in values.split(",") if name.strip()}
    if axioms - allowed:
        raise SystemExit(f"Unexpected axioms for {declaration}: {sorted(axioms)}")
    found[declaration] = sorted(axioms)
if expected - found.keys():
    raise SystemExit(f"Missing fresh axiom reports: {sorted(expected - found.keys())}")
summary = {
    "protected_reference_files_unchanged": ["Challenge.lean", "FactoringAux/Definitions.lean"],
    "checked_declarations": [
        {"declaration": name, "axioms": found[name]} for name in sorted(expected)
    ],
}
rendered = json.dumps(summary, indent=2) + "\n"
(logs / "axioms.json").write_text(rendered)
print(rendered, end="")
PY

echo "Verified k6 machine and rate. Logs: $factoring_logs"
