"""Compact power witnesses with direct exact unit-presentation reduction."""

from compact_oriented_power_search import _bind
from compact_witness_oriented_power_search import factor_with_compact_witnesses as _reference
from compact_witness_oriented_verifier import verify_compact_witness_search as _reference_verifier
from direct_unit_presentation import eliminate_direct_unit_presentation, verify_direct_unit_presentation


def _identity_descriptor(rows, initial_size, *, generator_count, allow_extra_rows, verify):
    assert allow_extra_rows and verify
    assert 0 <= initial_size <= generator_count
    return {'schema': 'identity-original-sparse-presentation-v1',
            'initial_size': initial_size, 'generator_count': generator_count,
            'input_row_count': len(rows)}


def _reduce_direct(rows, descriptor, *, verify=True):
    assert descriptor['schema'] == 'identity-original-sparse-presentation-v1'
    assert descriptor['input_row_count'] == len(rows)
    return eliminate_direct_unit_presentation(rows, descriptor['initial_size'],
        descriptor['generator_count'], verify=verify)


_direct = _bind(_reference, {
    'eliminate_signed_residual_graph': _identity_descriptor,
    'eliminate_unit_residual_columns': _reduce_direct,
})


def factor_with_direct_unit_refinement(*args, **kwargs):
    result = _direct(*args, **kwargs)
    result['component_replacements']['graph_projection'] = 'identity descriptor; no graph maps constructed'
    result['component_replacements']['unit_elimination'] = 'direct unit trace with independent incidence replay'
    result['component_replacements']['reference_compact_query_body_unchanged'] = True
    return result


def verify_direct_unit_search(data):
    """Reuse every literal-power check with direct presentation verification.

The older verifier obtains its retained prefix length from a graph-row list.
A private constant-space range supplies only that length; it is neither saved
nor treated as a graph map. The original result is not modified.
    """
    record = data['result']
    descriptor = record['graph_reduction']
    if descriptor is None:
        proxy = None
    else:
        assert set(descriptor) == {'schema', 'initial_size', 'generator_count', 'input_row_count'}
        assert descriptor['schema'] == 'identity-original-sparse-presentation-v1'
        assert descriptor['initial_size'] == len(record['active_primes'])
        assert descriptor['generator_count'] == descriptor['initial_size'] + len(record['refinement']['atoms'])
        assert descriptor['input_row_count'] == record['unit_reduction']['input_row_count']
        proxy = dict(descriptor, old_row_combinations=range(descriptor['input_row_count']))

    def verify_identity(rows, presentation):
        assert presentation is proxy
        assert len(rows) == presentation['input_row_count']
        assert all(type(j) is int and 0 <= j < presentation['generator_count']
                   for row in rows for j in row)
        return {'identity_presentation_checked': True,
                'original_sparse_rows': len(rows), 'graph_maps_materialized': False}

    def verify_units(rows, presentation, trace, *, verify_graph=True):
        assert verify_graph
        verify_identity(rows, presentation)
        return verify_direct_unit_presentation(rows, presentation['initial_size'],
                                               presentation['generator_count'], trace)

    altered = {**data, 'result': {**record, 'graph_reduction': proxy}}
    verifier = _bind(_reference_verifier, {
        'verify_sparse_projected_signed_graph': verify_identity,
        '_unit_verifier': verify_units,
    })
    return verifier(altered)
