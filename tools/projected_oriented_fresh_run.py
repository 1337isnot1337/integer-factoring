#!/usr/bin/env python3
"""Run a declared policy with the sparse intermediate projection adapter."""

from hashlib import sha256
import json
from pathlib import Path
import resource
from time import perf_counter, process_time
import traceback

from compact_oriented_power_search import _bind
from check_projected_oriented_power_search import binding_controls
from prepared_oriented_fresh_run import SOURCES as PREPARED_SOURCES, main as _reference_main
from projected_oriented_power_search import factor_with_projected_oriented_refinement

ROOT = Path(__file__).resolve().parents[1]
SOURCES = PREPARED_SOURCES + ('sparse_projected_signed_graph.py',
    'projected_oriented_power_search.py', 'check_projected_oriented_power_search.py',
    'projected_oriented_fresh_run.py', 'check_projected_oriented_fresh_run.py')


def source_hashes():
    return {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCES}


_main = _bind(_reference_main, {
    'source_hashes': source_hashes,
    'binding_controls': binding_controls,
    'factor_with_prepared_oriented_refinement': factor_with_projected_oriented_refinement,
})
assert _main.__code__ is _reference_main.__code__


def main():
    started, cpu_started = perf_counter(), process_time()
    try:
        _main()
    except MemoryError as error:
        # Clear unwound large local arrays so the bounded failure can be
        # recorded. Traceback locations remain available to the caller.
        traceback.clear_frames(error.__traceback__)
        print(json.dumps({'event': 'failed', 'exception': 'MemoryError',
            'attempt_elapsed_seconds': perf_counter()-started,
            'attempt_cpu_seconds': process_time()-cpu_started,
            'maximum_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            'address_space_limit_bytes': resource.getrlimit(resource.RLIMIT_AS)[0],
            'complete_search_finished': False,
            'complete_raw_witness_record_saved': False}), flush=True)
        raise


if __name__ == '__main__':
    main()
