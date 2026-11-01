#!/usr/bin/env python3
"""Replay every power, ideal and direct integer unit substitution."""

from compact_oriented_power_search import _bind
from direct_unit_oriented_fresh_run import source_hashes
from direct_unit_oriented_power_search import verify_direct_unit_search
from check_prepared_oriented_fresh_run import main as _reference_main

main = _bind(_reference_main, {
    'source_hashes': source_hashes,
    'verify_prepared_search': verify_direct_unit_search,
})
assert main.__code__ is _reference_main.__code__


if __name__ == '__main__':
    main()
