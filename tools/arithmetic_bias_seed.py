"""A coordinate-free, inverse-polynomial-gap seed given abelian generators.

This constructs exponent words only. Expander amplification is a published
algorithm used by the theorem, not implemented by this small seed routine.
"""
def powering_seed(generator_count, order_bit_bound):
    """Words for a lazy symmetric Cayley distribution when |G| <= 2**B."""
    if generator_count < 1 or order_bit_bound < 1:
        raise ValueError("positive generator count and order bit bound required")
    dimension, bound = generator_count, order_bit_bound
    zero = (0,) * dimension
    words = [zero] * (2*dimension*(bound+1))
    for i in range(dimension):
        for j in range(bound+1):
            word = tuple((1 << j) if k == i else 0 for k in range(dimension))
            words.extend((word, tuple(-value for value in word)))
    return words
