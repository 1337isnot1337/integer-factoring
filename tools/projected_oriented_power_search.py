"""The exact prepared selector with a sparse intermediate graph projection."""

from check_sparse_oriented_power_search import verify_search as _reference_verifier
from compact_oriented_power_search import _bind
from prepared_oriented_power_search import (
    _prepared, _reference, factor_with_prepared_oriented_refinement as _prepared_wrapper,
)
from signed_residual_graph import verify_signed_graph_elimination
from sparse_projected_signed_graph import (
    eliminate_sparse_projected_signed_graph, verify_sparse_projected_signed_graph,
)
from unit_presentation_trace import verify_unit_presentation_trace


_projected = _bind(_prepared,
    {'eliminate_signed_residual_graph': eliminate_sparse_projected_signed_graph})
_wrapped = _bind(_prepared_wrapper, {'_prepared': _projected})
_unit_verifier = _bind(verify_unit_presentation_trace,
    {'verify_signed_graph_elimination': verify_sparse_projected_signed_graph})


def factor_with_projected_oriented_refinement(*args, **kwargs):
    result = _wrapped(*args, **kwargs)
    result['component_replacements']['graph_projection'] = 'sparse view with complete integer maps'
    result['component_replacements']['graph_scratch_released_before_verification'] = True
    return result


def verify_projected_search(data):
    graph = data['result']['graph_reduction']

    def verify_presentation(rows, presentation):
        schema = presentation.get('schema')
        if schema == 'integer-unit-substitution-trace-v1':
            return _unit_verifier(rows, graph, presentation, verify_graph=True)
        if schema == 'sparse-projected-signed-graph-v1':
            return verify_sparse_projected_signed_graph(rows, presentation)
        return verify_signed_graph_elimination(rows, presentation)

    return _bind(_reference_verifier,
        {'verify_signed_graph_elimination': verify_presentation})(data)
