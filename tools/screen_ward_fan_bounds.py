#!/usr/bin/env python3
"""Complete rectangular upper profiles for the duplicated fan Ward code."""

import argparse
import hashlib
import json
from pathlib import Path

from check_ward_fan_family import family_record
from ward_channel_free_bounds import channel_free_score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--forms', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    family = family_record(args.forms)
    sources = (Path(__file__).name, 'check_ward_fan_family.py',
               'ward_channel_free_bounds.py', 'ward_rectangular_bounds.py')
    hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in sources}
    records = []
    for selected in (False, True):
        field = ('orientation_selected_central_labels' if selected
                 else 'central_grade_dimension')
        labels = family[field]
        if labels <= 1:
            raise ValueError('The profile requires more than one logical label.')
        record = channel_free_score(family['pairs'], family['pairs'],
                                    args.forms, labels, 0)
        record.update(family='duplicated fan columns', orientation_selected=selected,
                      source_sha256=hashes,
                      scope='Complete unrestricted-rectangle upper profile; no assertion about optimal actual ranks.')
        records.append(record)
        print(json.dumps({key: value for key, value in record.items()
                          if key not in ('profile', 'source_sha256')}, sort_keys=True), flush=True)
    args.output.write_text(json.dumps(records, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
