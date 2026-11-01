#!/bin/sh
# Reproduces data/*.jsonl (deterministic; appends, so delete old files first).
set -e
cd "$(dirname "$0")"
B="60 70 80 90 100 110 120 130 140 150 160"
python3 window_scan.py --bits $B --per-size 8 --schedule fixed4093 --targets 16 --hits 4 --cap 4000000 --jobs 8 --out data/scan_fixed4093.jsonl
python3 window_scan.py --bits $B --per-size 8 --schedule bsq --targets 16 --hits 4 --cap 4000000 --jobs 8 --out data/scan_bsq.jsonl
python3 window_scan.py --bits 80 100 120 --per-size 16 --first-index 100 --schedule fixed4093 --targets 100000 --hits 2 --cap 4000000 --jobs 8 --out data/scan_alltargets.jsonl
${PY:-python3} adversary.py --p 3 5 7 --bits 64 96 128 160 --cutoff 4093 --out data/adversary_cutoff4093.json
${PY:-python3} adversary.py --p 3 5 7 --bits 64 96 128 160 --cutoff 4093 --semiprime --out data/adversary_semiprime_cutoff4093.json
python3 window_scan.py --bits 90 100 110 --per-size 64 --first-index 1000 --schedule fixed4093 --targets 100000 --hits 2 --cap 4000000 --jobs 8 --out data/scan_tail.jsonl
python3 small_p.py data/scan_alltargets.jsonl data/scan_tail.jsonl > data/small_p.md
${PY:-python3} analyze.py > data/summary.md; python3 tail_check.py data/scan_alltargets.jsonl data/scan_tail.jsonl > data/tail_check.txt
