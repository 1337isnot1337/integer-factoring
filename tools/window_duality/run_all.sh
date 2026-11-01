#!/usr/bin/env bash
# Reproduce every number in duality.md. Needs python3 with cypari2 + numpy and the
# repo tools/ directory one level up (override with IF_TOOLS).
# Total runtime well under 5 minutes (about 20 s on the research box).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
mkdir -p data
$PY duality_scan.py --count 360 --lo 4 --hi 9 --seed 1 --out data/scan.jsonl > data/scan.log
$PY duality_scan.py --count 48 --lo 9 --hi 11 --seed 2 --u 2.5 3 4 --no-adversary --out data/scan_large.jsonl > data/scan_large.log
$PY analyze_scan.py data/scan.jsonl > output_scan.txt
$PY analyze_scan.py data/scan_large.jsonl > output_scan_large.txt
sha256sum output_*.txt > SHA256SUMS   # jsonl line order depends on worker scheduling; the outputs do not
