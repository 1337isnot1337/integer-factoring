#!/usr/bin/env python3
"""Replay all compact witnesses and the complete final integer presentation."""

from compact_oriented_power_search import _bind
from compact_witness_oriented_fresh_run import source_hashes
from compact_witness_oriented_verifier import verify_compact_witness_search
from check_prepared_oriented_fresh_run import main as _reference_main

main = _bind(_reference_main, {
    'source_hashes': source_hashes,
    'verify_prepared_search': verify_compact_witness_search,
})
assert main.__code__ is _reference_main.__code__


if __name__ == '__main__':
    main()
