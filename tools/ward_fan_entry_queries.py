#!/usr/bin/env python3
"""Division-free physical Ward queries for duplicated fan incidence columns.

Values use a characteristic-two representation with XOR addition. Supply
the ring multiplication operation explicitly. It may be field multiplication,
polynomial multiplication, or multiplication in a commutative quotient ring.
"""


def collapse_parallel_sites(sites, multiply):
    ordinary, inserted = 1, 0
    for local_ordinary, local_inserted in sites:
        inserted = (multiply(inserted, local_ordinary)
                    ^ multiply(ordinary, local_inserted))
        ordinary = multiply(ordinary, local_ordinary)
    return ordinary, inserted


def fan_entry_query(forms, sites, multiply):
    """Evaluate [z_1...z_k] product of local Ward factors for the fan.

    The input order is two copies of spoke1, two copies of path12,
    two copies of spoke2, and so on. Each site supplies (m_i,b_i).
    The final scalar equals the physical entry of W(eta) M_P.
    Only two accumulated ring elements and one local pair are retained.
    """
    if forms < 1:
        raise ValueError('A duplicated fan has at least one form.')
    stream = iter(sites)

    def pair():
        try:
            first, second = next(stream), next(stream)
        except StopIteration as error:
            raise ValueError('The site stream ended inside the duplicated fan.') from error
        return collapse_parallel_sites((first, second), multiply)

    connected, pending = 0, 1
    ordinary, inserted = pair()
    connected, pending = (multiply(ordinary, connected) ^ multiply(inserted, pending),
                          multiply(ordinary, pending))
    for _ in range(1, forms):
        ordinary, inserted = pair()  # Path edge: the earlier pending component must survive.
        connected, pending = (multiply(inserted, connected),
                              multiply(ordinary, connected) ^ multiply(inserted, pending))
        ordinary, inserted = pair()  # Spoke edge: selecting it may close the pending component.
        connected, pending = (multiply(ordinary, connected) ^ multiply(inserted, pending),
                              multiply(ordinary, pending))
    try:
        next(stream)
    except StopIteration:
        return connected
    raise ValueError('The site stream contains more than 4*k-2 entries.')
