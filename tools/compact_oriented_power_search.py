"""The frozen power search with exact compact presentation certificates."""

from types import FunctionType

from check_sparse_oriented_power_search import verify_search as _verify_reference
from primitive_integer_row_basis import PrimitiveIntegerRowBasis
from quotient_update_oriented_basis import QuotientUpdateOrientedBasis
from signed_residual_graph import verify_signed_graph_elimination
from sparse_oriented_power_search import factor_with_sparse_oriented_refinement as _reference
from unit_presentation_trace import eliminate_unit_residual_columns_compact, verify_unit_presentation_trace


def _unit_after_verified_graph(rows, graph, *, verify=True):
    # The unchanged search body has just constructed and verified this graph.
    assert graph['verification']['integer_presentation_isomorphism'] is True
    return eliminate_unit_residual_columns_compact(rows, graph, verify=verify, verify_graph=False)


def _bind(function, replacements):
    namespace = dict(function.__globals__)
    namespace.update(replacements)
    result = FunctionType(function.__code__, namespace, function.__name__,
                          function.__defaults__, function.__closure__)
    result.__kwdefaults__ = dict(function.__kwdefaults__ or {})
    return result


_compact = _bind(_reference, {
    'IncrementalOrientedBasis': QuotientUpdateOrientedBasis,
    'RationalRowBasis': PrimitiveIntegerRowBasis,
    'eliminate_unit_residual_columns': _unit_after_verified_graph,
})


def factor_with_compact_oriented_refinement(*args, **kwargs):
    """Retain the witness policy while replacing three proved exact components."""
    result = _compact(*args, **kwargs)
    result['component_replacements'] = {
        'atom_tree': 'exact leaf quotient updates',
        'row_selection': 'primitive temporary residuals; original relations retained',
        'unit_elimination': 'compact verified unit substitution trace',
        'reference_search_body_unchanged': True,
        'reference_module_globals_unchanged': True,
        'graph_verified_before_compact_unit_pass': True,
    }
    return result


def verify_compact_search(data):
    """Reuse the complete independent replay with a compact-proof dispatcher."""
    graph = data['result']['graph_reduction']

    def verify_presentation(rows, presentation):
        if presentation.get('schema') == 'integer-unit-substitution-trace-v1':
            return verify_unit_presentation_trace(rows, graph, presentation, verify_graph=True)
        return verify_signed_graph_elimination(rows, presentation)

    verifier = _bind(_verify_reference, {'verify_signed_graph_elimination': verify_presentation})
    return verifier(data)
