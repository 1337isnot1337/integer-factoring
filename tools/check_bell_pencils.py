#!/usr/bin/env python3
"""Read back every saved Bell linear pencil and reconstruct every exact rank."""

import json
from pathlib import Path

from search_bell_pencils import profile_record, run


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    checked = 0
    for suffix in ("d2", "d4", "d8", "d16"):
        path = root / "evidence" / f"bell-linear-{suffix}.json"
        record = json.loads(path.read_text())
        lengths = [p["length"] for p in record["best"]["profiles"]]
        rebuilt = profile_record(tuple(record["best"]["values"]),
                                 record["bits_per_local_label"], lengths)
        assert json.loads(json.dumps(rebuilt)) == record["best"]
        if record["exhaustive_alphabet"]:
            replay = run(record["bits_per_local_label"], record["objective_length"],
                         record["exhaustive_alphabet"], 0, record["seed"])
            for key in ("candidate_evaluations", "distinct_histogram_pairs", "best"):
                assert json.loads(json.dumps(replay[key])) == record[key]
        checked += len(lengths)
    print(f"Rebuilt {checked} complete two-cut profiles and both exhaustive calibrations.")


if __name__ == "__main__":
    main()
