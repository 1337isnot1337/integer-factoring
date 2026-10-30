#!/usr/bin/env python3
"""Run an exhaustive native correction screen in resumable exact-index chunks.

Each completed chunk is saved atomically. A restarted invocation skips only
chunks whose command and native source hash agree. Capped ranks are recorded
as lower bounds, not as an exact histogram above the stopping threshold.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--family", choices=("bilinear", "quadratic", "lookup"), required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--chunk", type=int, default=1000000)
    parser.add_argument("--max-new-chunks", type=int, default=8)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    if args.chunk <= 0 or args.start < 0 or args.stop <= args.start or args.max_new_chunks < 1:
        parser.error("require a positive chunk limit and a nonempty nonnegative index interval")
    args.directory.mkdir(parents=True, exist_ok=True)
    native_source = Path(__file__).with_name("explore_group_correction_rank.cpp")
    source_hash = hashlib.sha256(native_source.read_bytes()).hexdigest()
    completed, created = [], 0
    for start in range(args.start, args.stop, args.chunk):
        count = min(args.chunk, args.stop - start)
        command = [str(args.binary.resolve()), "--k", str(args.degree), "--family", args.family,
                   "--normalized", "--start", str(start), "--count", str(count)]
        if args.exact:
            command.append("--exact")
        destination = args.directory / f"chunk-{start:012d}-{count:012d}.json"
        if destination.exists():
            record = json.loads(destination.read_text())
            if record["command"] != command or record["native_source_sha256"] != source_hash:
                raise ValueError(f"incompatible checkpoint: {destination}")
        else:
            if created == args.max_new_chunks:
                break
            events = []
            with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
                assert process.stdout is not None
                for line in process.stdout:
                    event = json.loads(line)
                    if event["event"] in ("start", "improvement", "complete"):
                        events.append(event)
                if process.wait() != 0:
                    raise RuntimeError("native screen failed; incomplete chunk was not saved")
            if not events or events[-1]["event"] != "complete" or events[-1]["tested"] != count:
                raise RuntimeError("native screen did not finish its requested interval")
            if sum(number for _, number in events[-1]["histogram"]) != count:
                raise RuntimeError("native histogram does not cover the requested interval")
            record = {"command": command, "native_source_sha256": source_hash,
                      "exact_ranks": args.exact,
                      "rank_semantics": "exact" if args.exact else "rank or certified lower bound at the current stopping threshold",
                      "events": events}
            temporary = destination.with_suffix(".tmp")
            temporary.write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
            temporary.replace(destination)
            created += 1
        completed.append({"start": start, "count": count,
                          "best_rank_including_zero_baseline": record["events"][-1]["best_rank"],
                          "digest_fnv64": record["events"][-1]["digest_fnv64"]})
        print(json.dumps({"event": "chunk_saved", **completed[-1]}), flush=True)
    print(json.dumps({"event": "batch_complete", "new_chunks": created,
                      "completed_chunks": len(completed),
                      "next_index": completed[-1]["start"] + completed[-1]["count"] if completed else args.start,
                      "requested_stop": args.stop}), flush=True)


if __name__ == "__main__":
    main()
