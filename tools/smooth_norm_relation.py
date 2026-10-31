"""Turn supplied smooth represented norms into short relation words.

An explicit pair of small forms and an exact balanced divisor are unnecessary
for the relation matrix.  Locating the supplied norm remains an open search
problem; see research/smooth-norm-relations.md.
"""

from math import gcd

from binary_quadratic_forms import reduce_form
from bounded_smooth_norm import bezout, smooth_factors


def norm_to_prime_word(form, x, y, smooth_bound, norm_bound=None):
    """Return the short prime-form word of a supplied primitive smooth norm.

    Coordinates are divided by their gcd first.  None denotes nonsmoothness
    or failure of the optional norm bound.  Invalid forms or noncoprime norms
    raise ValueError.  Each word entry uses the canonical root in [0,p].
    """
    a, b, c = form
    discriminant = b*b - 4*a*c
    common = gcd(x, y)
    if (discriminant >= 0 or a <= 0 or gcd(gcd(a, b), c) != 1
            or common == 0 or smooth_bound < 2):
        raise ValueError("a primitive positive form and a nonzero vector are required")
    x, y = x // common, y // common
    norm = a*x*x + b*x*y + c*y*y
    if gcd(norm, discriminant) != 1:
        raise ValueError("the primitive norm must be coprime to the discriminant")
    if norm_bound is not None and norm > norm_bound:
        return None
    factors = smooth_factors(norm, smooth_bound)
    if factors is None:
        return None
    v, negative_u = bezout(x, y)
    u = -negative_u
    assert x*v-y*u == 1
    middle = 2*a*x*u + b*(x*v+y*u) + 2*c*y*v
    assert (middle*middle-discriminant) % (4*norm) == 0
    word = []
    for prime in sorted(set(factors)):
        root = middle % (2*prime)
        canonical = min(root, 2*prime-root)
        sign = 1 if root == canonical else -1
        exponent = sign*factors.count(prime)
        assert (canonical*canonical-discriminant) % (4*prime) == 0
        prime_form = reduce_form((prime, canonical,
                                 (canonical*canonical-discriminant)//(4*prime)))
        word.append({"prime": prime, "root": canonical,
                     "exponent": exponent, "form": prime_form})
    length = sum(abs(entry["exponent"]) for entry in word)
    assert 1 << length <= norm
    return {"primitive_vector": (x, y), "norm": norm,
            "equivalent_middle": middle, "word": word, "word_length": length}


def dominant_relation_row(generator_index, generator_primes, sample_word,
                          norm_word, sample_bound, norm_word_bound):
    """Build the relation row after the group identity has been established.

    The caller supplies a witness that g_i^K * sample equals the norm class.
    This function checks the word bounds and exact strict diagonal dominance.
    It does not claim to verify that group identity from integer vectors alone.
    """
    size = len(generator_primes)
    if (len(sample_word) != size or not 0 <= generator_index < size
            or len(set(generator_primes)) != size):
        raise ValueError("invalid factor-base indexing or sample word")
    index = {prime: i for i, prime in enumerate(generator_primes)}
    short = [0]*size
    for entry in norm_word["word"]:
        if entry["prime"] not in index:
            raise ValueError("a norm prime is absent from the factor base")
        short[index[entry["prime"]]] += entry["exponent"]
    if (sum(abs(value) for value in sample_word) > sample_bound
            or sum(abs(value) for value in short) > norm_word_bound):
        raise ValueError("the supplied words exceed the certified bounds")
    diagonal = sample_bound + norm_word_bound + 1
    row = [left-right for left, right in zip(sample_word, short)]
    row[generator_index] += diagonal
    assert row[generator_index] > sum(abs(value) for i, value in enumerate(row)
                                     if i != generator_index)
    return row
