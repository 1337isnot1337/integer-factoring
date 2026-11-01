"""The frozen power search with prepared atom batches and small integer pivots."""

from batched_coprime_bucket_basis import BufferedBucketOrientedBasis
from compact_oriented_power_search import (
    _bind, _reference, _unit_after_verified_graph, verify_compact_search,
)
from small_pivot_integer_row_basis import SmallPivotIntegerRowBasis


_prepared = _bind(_reference, {
    'IncrementalOrientedBasis': BufferedBucketOrientedBasis,
    'RationalRowBasis': SmallPivotIntegerRowBasis,
    'eliminate_unit_residual_columns': _unit_after_verified_graph,
})


def factor_with_prepared_oriented_refinement(*args, **kwargs):
    """Preserve all witnesses and checks while replacing three exact components."""
    result = _prepared(*args, **kwargs)
    result['component_replacements'] = {
        'atom_tree': 'prepared coprime bucket remainders',
        'batch_size': 512,
        'bucket_size': 512,
        'row_selection': 'small-magnitude primitive residual pivots; original rows retained',
        'unit_elimination': 'compact verified unit substitution trace',
        'reference_search_body_unchanged': True,
        'reference_module_globals_unchanged': True,
        'graph_verified_before_compact_unit_pass': True,
        'partial_batch_flush_charged_in_transport_phase': True,
    }
    return result


verify_prepared_search = verify_compact_search
