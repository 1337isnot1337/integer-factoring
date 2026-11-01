#!/usr/bin/env python3
"""Run a declared power policy with compact reconstructible witnesses."""

from hashlib import sha256
from pathlib import Path

from check_compact_witness_oriented_power_search import SOURCES as MATCHED_SOURCES, query_logic_controls
from compact_oriented_power_search import _bind
from compact_witness_oriented_power_search import factor_with_compact_witnesses
from prepared_oriented_fresh_run import main as _reference_main
from projected_oriented_fresh_run import main as _failure_reporting_main

ROOT = Path(__file__).resolve().parents[1]
SOURCES = MATCHED_SOURCES + ('prepared_oriented_fresh_run.py',
    'check_prepared_oriented_fresh_run.py', 'projected_oriented_fresh_run.py',
    'compact_witness_oriented_fresh_run.py', 'check_compact_witness_oriented_fresh_run.py')


def source_hashes():
    return {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCES}


_main = _bind(_reference_main, {
    'source_hashes': source_hashes,
    'binding_controls': query_logic_controls,
    'factor_with_prepared_oriented_refinement': factor_with_compact_witnesses,
})
assert _main.__code__ is _reference_main.__code__
main = _bind(_failure_reporting_main, {'_main': _main})


if __name__ == '__main__':
    main()
