"""Bind two exact component replacements to the frozen power-search body.

The reference module and its globals stay unchanged, so archived receipts
retain their source guards. An isolated function namespace supplies the
quotient-update basis and sparse unit elimination. The search policy,
witness order, stopping rule, and verification calls are identical.
"""

from types import FunctionType

from quotient_update_oriented_basis import QuotientUpdateOrientedBasis
from sparse_oriented_power_search import factor_with_sparse_oriented_refinement as _reference
from sparse_unit_presentation import eliminate_sparse_unit_residual_columns


def _bind_exact_components():
    namespace = dict(_reference.__globals__)
    namespace['IncrementalOrientedBasis'] = QuotientUpdateOrientedBasis
    namespace['eliminate_unit_residual_columns'] = eliminate_sparse_unit_residual_columns
    function = FunctionType(_reference.__code__, namespace,
                            'factor_with_optimized_oriented_refinement',
                            _reference.__defaults__, _reference.__closure__)
    function.__kwdefaults__ = dict(_reference.__kwdefaults__ or {})
    return function


_optimized = _bind_exact_components()


def factor_with_optimized_oriented_refinement(*args, **kwargs):
    """Run the unchanged reference policy with exact cheaper components."""
    result = _optimized(*args, **kwargs)
    result['component_replacements'] = {
        'atom_tree': 'exact leaf quotient updates',
        'unit_elimination': 'sparse incidence with stable indices',
        'reference_search_body_unchanged': True,
        'reference_module_globals_unchanged': True,
    }
    return result
