#!/usr/bin/env python3
"""Reanalyze PR75's saved data; this does not rerun its PARI scan."""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    folder=ROOT/'tools/window_duality';results=[]
    for input_name,expected_name in (('scan.jsonl','output_scan.txt'),
                                     ('scan_large.jsonl','output_scan_large.txt')):
        result=subprocess.run([sys.executable,str(folder/'analyze_scan.py'),str(folder/'data'/input_name)],
                              cwd=ROOT,capture_output=True,check=True)
        expected=(folder/expected_name).read_bytes()
        assert result.stdout==expected,expected_name
        results.append(dict(data=input_name,output=expected_name,sha256=sha256(expected).hexdigest()))
    for line in (folder/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split();assert sha256((folder/name).read_bytes()).hexdigest()==digest
    print(json.dumps(dict(status='PASS',scope='Byte-identical analysis of saved records; no PARI group scan.',results=results)))


if __name__=='__main__':main()
