#!/usr/bin/env python3
"""Replay a source-guarded prepared search with sparse graph projections."""

from compact_oriented_power_search import _bind
from check_prepared_oriented_fresh_run import main as _reference_main
from projected_oriented_fresh_run import source_hashes
from projected_oriented_power_search import verify_projected_search


main = _bind(_reference_main, {
    'source_hashes': source_hashes,
    'verify_prepared_search': verify_projected_search,
})
assert main.__code__ is _reference_main.__code__


if __name__ == '__main__':
    main()
