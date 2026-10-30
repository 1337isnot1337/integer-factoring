"""Validate the completed k6 build and write receipts bound to its inputs.

Usage: python3 tools/check_k6_receipt.py ROOT BUILD LOGS

The build driver must call this only after both Lean commands have succeeded.
This checker independently requires Lake's success line, fresh endpoint
reports, and unchanged source bytes. An unsuccessful check removes previous
success receipts so they cannot be mistaken for results of the current run.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
import tempfile


REFERENCE_COMMIT = "c91cc3f4947caca8977e81a5d14244f73b5ae2b1"
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
ENDPOINTS = {
    "Factoring.K6.factoring_root_rate",
    "Factoring.K6.factoring_real_rate",
    "Factoring.K6.certified_improvement",
}
PROTECTED = {
    "Challenge.lean", "FactoringAux/Definitions.lean",
    "lean-toolchain", "lake-manifest.json",
}
PROOF_PREFIX = "variants/k6/FactoringK6/"
PATCH_PATH = "variants/k6/reference.patch"
LOG_NAMES = (
    "patch-check.log", "patch-apply.log", "toolchain.log",
    "cache.log", "build.log", "Certified.log",
)
RECEIPT_NAMES = ("axioms.json", "verification.json")
ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
PRINT_AXIOMS = re.compile(r"^\s*#print\s+axioms\s+([A-Za-z0-9_.']+)\s*$", re.MULTILINE)
AXIOMS = re.compile(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", re.MULTILINE)
NO_AXIOMS = re.compile(r"'([^']+)' does not depend on any axioms", re.MULTILINE)
BAD_LOG = re.compile(r"\bsorry(?:Ax)?\b|\berror\s*:|Some required targets logged failures", re.IGNORECASE)
BUILD_SUCCESS = re.compile(r"(?:^|\n)Build completed successfully \(\d+ jobs?\)\.\s*\Z")


class VerificationError(ValueError):
    """The supplied files do not constitute a successful k6 receipt."""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def invalidate_receipts(logs: Path) -> None:
    for name in RECEIPT_NAMES:
        (logs / name).unlink(missing_ok=True)


def read_regular(base: Path, name: str) -> bytes:
    relative = PurePosixPath(name)
    if relative.is_absolute() or ".." in relative.parts or str(relative) != name:
        raise VerificationError(f"Noncanonical or unsafe input path: {name!r}")
    path = base
    for part in relative.parts:
        path = path / part
        if path.is_symlink():
            raise VerificationError(f"Input must not traverse a symlink: {path}")
    if not path.is_file():
        raise VerificationError(f"Missing regular input file: {path}")
    return path.read_bytes()


def hash_map(manifest: dict, key: str) -> dict[str, str]:
    values = manifest.get(key)
    if not isinstance(values, dict) or not values:
        raise VerificationError(f"Missing nonempty source map: {key}")
    if any(not isinstance(name, str) or not isinstance(value, str)
           or re.fullmatch(r"[0-9a-f]{64}", value) is None
           for name, value in values.items()):
        raise VerificationError(f"Invalid SHA-256 map: {key}")
    return values


def verify_hash(base: Path, name: str, expected: str) -> bytes:
    data = read_regular(base, name)
    if sha256(data) != expected:
        raise VerificationError(f"Source hash mismatch: {base / name}")
    return data


def endpoint_reports(source: str, text: str, *, allow_other_reports: bool = False) -> list[dict]:
    printed = PRINT_AXIOMS.findall(source)
    if set(printed) != ENDPOINTS or len(printed) != len(ENDPOINTS):
        raise VerificationError("Certified.lean must print exactly the three required endpoints")
    reports = AXIOMS.findall(text) + [(name, "") for name in NO_AXIOMS.findall(text)]
    if allow_other_reports:
        reports = [(name, values) for name, values in reports if name in ENDPOINTS]
    found = {}
    for name, values in reports:
        if name in found:
            raise VerificationError(f"Duplicate endpoint axiom report: {name}")
        axioms = {value.strip() for value in values.split(",") if value.strip()}
        if axioms - ALLOWED_AXIOMS:
            raise VerificationError(f"Unexpected axioms for {name}: {sorted(axioms - ALLOWED_AXIOMS)}")
        found[name] = sorted(axioms)
    if set(found) != ENDPOINTS:
        raise VerificationError(
            f"Endpoint axiom reports differ from the required set; "
            f"missing={sorted(ENDPOINTS - found.keys())}, extra={sorted(found.keys() - ENDPOINTS)}"
        )
    return [{"declaration": name, "axioms": found[name]} for name in sorted(ENDPOINTS)]


def atomic_write(path: Path, data: bytes) -> None:
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".k6-receipt-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(data)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def verify_receipt(root: Path, build: Path, logs: Path) -> dict:
    root, build, logs = root.resolve(), build.resolve(), logs.resolve()
    invalidate_receipts(logs)
    try:
        manifest_bytes = read_regular(logs, "sources.json")
        manifest = json.loads(manifest_bytes)
        if not isinstance(manifest, dict) or manifest.get("reference_commit") != REFERENCE_COMMIT:
            raise VerificationError("Source receipt has the wrong pinned reference commit")
        patched = hash_map(manifest, "patched_reference_sha256")
        protected = hash_map(manifest, "protected_reference_sha256")
        sources = hash_map(manifest, "source_sha256")
        if set(protected) != PROTECTED:
            raise VerificationError("Source receipt does not cover exactly the protected reference inputs")
        if any(not name.startswith(PROOF_PREFIX) or not name.endswith(".lean") for name in sources):
            raise VerificationError("Proof source receipt contains a path outside FactoringK6")
        if PROOF_PREFIX + "Certified.lean" not in sources:
            raise VerificationError("Source receipt omits Certified.lean")

        patch = verify_hash(root, PATCH_PATH, manifest.get("patch_sha256", ""))
        destinations = [line[6:] for line in patch.decode().splitlines() if line.startswith("+++ b/")]
        if len(destinations) != len(set(destinations)) or set(destinations) != set(patched):
            raise VerificationError("Patched source receipt does not match the patch destinations")
        if "+++ /dev/null" in patch.decode().splitlines():
            raise VerificationError("Deleted reference inputs are unsupported by this receipt format")

        checked_build = {}
        for values in (patched, protected):
            for name, expected in values.items():
                verify_hash(build, name, expected)
                checked_build[name] = expected
        source_data = {}
        for name, expected in sources.items():
            verify_hash(root, name, expected)
            copied = "FactoringK6/" + name.removeprefix(PROOF_PREFIX)
            source_data[copied] = verify_hash(build, copied, expected)
            checked_build[copied] = expected
        actual_proofs = {path.relative_to(build).as_posix()
                         for path in (build / "FactoringK6").rglob("*.lean")}
        if actual_proofs != set(source_data):
            raise VerificationError("Built FactoringK6 source set differs from the captured source set")

        log_bytes = {name: read_regular(logs, name) for name in LOG_NAMES}
        texts = {name: ANSI.sub("", data.decode()) for name, data in log_bytes.items()}
        for name, text in texts.items():
            if BAD_LOG.search(text):
                raise VerificationError(f"Incomplete or failed verification in {name}")
        if any("Skipped patch" in texts[name] for name in ("patch-check.log", "patch-apply.log")):
            raise VerificationError("Git skipped a reference patch")
        if not BUILD_SUCCESS.search(texts["build.log"]):
            raise VerificationError("Build log has no terminal successful Lake build result")
        if not re.search(r"\bLean \(version 4\.35\.0,", texts["toolchain.log"]):
            raise VerificationError("Toolchain log does not identify Lean 4.35.0")
        declarations = endpoint_reports(
            source_data["FactoringK6/Certified.lean"].decode(), texts["Certified.log"]
        )
        # A fresh source build prints these reports itself. Requiring agreement
        # with the separate Lean pass prevents mixing an old endpoint log into
        # the new build's receipt merely because theorem names stayed the same.
        built_declarations = endpoint_reports(
            source_data["FactoringK6/Certified.lean"].decode(), texts["build.log"],
            allow_other_reports=True,
        )
        if declarations != built_declarations:
            raise VerificationError("Build and endpoint axiom reports disagree")
        summary = {
            "protected_reference_files_unchanged": sorted(PROTECTED),
            "checked_declarations": declarations,
        }
        summary_bytes = (json.dumps(summary, indent=2) + "\n").encode()
        verification = {
            "schema_version": 1,
            "status": "success",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "reference_commit": REFERENCE_COMMIT,
            "build_id": build.name,
            "source_manifest_sha256": sha256(manifest_bytes),
            "patch_sha256": sha256(patch),
            "build_source_sha256": checked_build,
            "root_proof_sha256": sources,
            "log_sha256": {name: sha256(data) for name, data in log_bytes.items()},
            "axioms_receipt_sha256": sha256(summary_bytes),
            "checker_sha256": sha256(Path(__file__).read_bytes()),
            "checked_declarations": declarations,
        }
        verification_bytes = (json.dumps(verification, indent=2) + "\n").encode()
        atomic_write(logs / "axioms.json", summary_bytes)
        atomic_write(logs / "verification.json", verification_bytes)
        return verification
    except BaseException:
        invalidate_receipts(logs)
        raise


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: python3 tools/check_k6_receipt.py ROOT BUILD LOGS")
    try:
        result = verify_receipt(*(Path(value) for value in sys.argv[1:]))
    except (OSError, ValueError) as error:
        raise SystemExit(f"k6 receipt verification failed: {error}") from None
    print(json.dumps({"status": result["status"], "checked_declarations": result["checked_declarations"]}, indent=2))


if __name__ == "__main__":
    main()
