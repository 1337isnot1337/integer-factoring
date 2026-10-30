"""Check the printed axiom reports from the formal verification script."""

from __future__ import annotations

import json
from pathlib import Path
import re
import sys


ALLOWED = {"propext", "Classical.choice", "Quot.sound"}
REPORT = re.compile(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", re.MULTILINE)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: python3 tools/check_axioms.py LOG_DIRECTORY")
    directory = Path(sys.argv[1])
    proof_dir = Path(__file__).resolve().parents[1] / "formal"
    summaries = []
    for proof in sorted(proof_dir.glob("*.lean")):
        log = directory / (proof.stem + ".log")
        text = log.read_text()
        reports = REPORT.findall(text)
        if not reports:
            raise SystemExit(f"No printed axiom report for {proof.name}")
        if "sorryAx" in text or "error:" in text:
            raise SystemExit(f"Incomplete or failed proof in {proof.name}")
        for declaration, values in reports:
            axioms = {name.strip() for name in values.split(",") if name.strip()}
            if axioms - ALLOWED:
                raise SystemExit(f"Unexpected axioms for {declaration}: {sorted(axioms)}")
            summaries.append({"declaration": declaration, "axioms": sorted(axioms)})
    if not summaries:
        raise SystemExit("No formal sources found")
    print(json.dumps({"checked_declarations": summaries}, indent=2))


if __name__ == "__main__":
    main()
