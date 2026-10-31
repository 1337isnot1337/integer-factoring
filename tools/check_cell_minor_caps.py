#!/usr/bin/env python3
"""Cell-set minor caps for the gauged two-place Smith relaxation.

Three independent parts, all exact:

1. Literal control.  Random scalar Laurent strict crossings over F_p with
   random integer heights of both signs are normalized exactly as in
   research/shared-minor-profile-relaxation.md section 1.  Their complete
   Smith profiles are computed by exact truncated power-series elimination,
   and every cell-set, shared-line and rectangular prefix cap is checked.
   Random cell-set minors are also expanded by exact interpolation and their
   endpoint value, degree and both orders are compared with the proof.

2. Relaxation verdicts.  The cell caps are added to the shared-line and
   rectangular caps, the exact rational greatest convex minorant is built,
   and outward dyadic bounds decide the sign of G_relaxed at the two targets
   of the predecessor note and at the sharper targets listed below.

3. A compact ternary region (q=3) on which the strengthened relaxation is
   certified positive by a monotone corner bound and a dyadic box cover.

Relaxation points are not tensors.  Nothing here is a crossing certificate
or a factoring improvement.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys

from charged_gaussian_floor import Dyadic


F = Fraction
PRIME = 10007


def require(condition, message):
    if not condition:
        raise SystemExit('FAILED: ' + message)


# ---------------------------------------------------------------------------
# Part 2 helpers: caps, convex minorant and outward pressure bounds.

def all_caps(d, ex, ey, cells=True, shared=True, rectangles=True):
    """Return {k: (cap, family)} for prefix sums S_k of gamma=(r+s)/2."""
    q = len(ex)
    ex = sorted(F(v) for v in ex)
    ey = sorted(F(v) for v in ey)
    require(len(ey) == q and ex[0] == 0 and ey[0] == 0 and d >= 0,
            'normalized data: equal alphabets, paired minima zero, d>=0')
    d = F(d)
    caps = {0: (F(0), 'empty'), 1: (F(0), 'r1=s1=0')}

    def put(k, value, family):
        # A rectangle cap equals the cell cap of its product set, so on ties
        # the cell-set label is kept: rectangles never bind strictly.
        if k not in caps or value < caps[k][0] or (
                value == caps[k][0] and family == 'cell set'):
            caps[k] = (value, family)
    px = [F(0)] + list(itertools.accumulate(ex))
    py = [F(0)] + list(itertools.accumulate(ey))
    if rectangles:
        for a in range(1, q + 1):
            for b in range(1, q + 1):
                put(a * b, a * b * d + b * px[a] + a * py[b], 'rectangle %dx%d' % (a, b))
    if shared:
        for k in range(2, q + 1):
            put(k, (k * d + min(px[k], py[k])) / 2, 'shared line')
    if cells:
        cell = sorted(x + y for x in ex for y in ey)
        total = F(0)
        for k in range(1, q * q + 1):
            total += cell[k - 1]
            put(k, k * d + total, 'cell set')
    return caps


def envelope(caps):
    """Greatest convex minorant of the cap points, as in the predecessor tool.

    Same algorithm as check_shared_minor_profile_targets.envelope, copied
    to avoid that module's optional scipy import chain.  Returns the sorted
    slopes gamma* and the hull vertices; both properties are re-checked.
    """
    hull = []
    for p in sorted(caps.items()):
        while len(hull) > 1 and ((hull[-1][1] - hull[-2][1]) * (p[0] - hull[-1][0])
                                 >= (p[1] - hull[-1][1]) * (hull[-1][0] - hull[-2][0])):
            hull.pop()
        hull.append(p)
    orders = []
    for a, b in zip(hull, hull[1:]):
        orders.extend([(b[1] - a[1]) / (b[0] - a[0])] * (b[0] - a[0]))
    require(orders == sorted(orders), 'hull slopes sorted')
    require(all(sum(orders[:k]) <= bound for k, bound in caps.items()), 'hull below caps')
    return tuple(orders), hull


def relaxed_orders(d, ex, ey, **families):
    caps = all_caps(d, ex, ey, **families)
    orders, hull = envelope({k: value for k, (value, _) in caps.items()})
    return orders, hull, caps


def point(arithmetic, value):
    value = F(value)
    return arithmetic.rational(value.numerator, value.denominator)


def log_at_least_one(arithmetic, value):
    """Outward bound for log(v), v>=1, by binary argument reduction."""
    require(value[0] >= arithmetic.scale, 'logarithm argument below one')
    shift = max(0, (value[0] // arithmetic.scale).bit_length() - 1)
    reduced = (value[0] >> shift, arithmetic.ceil_div(value[1], 1 << shift))
    tail = arithmetic.log_one_plus(arithmetic.subtract(reduced, arithmetic.one))
    log_two = arithmetic.log_one_plus(arithmetic.one)
    return (shift * log_two[0] + tail[0], shift * log_two[1] + tail[1])


def partition(arithmetic, values, root_two=False):
    total = arithmetic.zero
    for value in values:
        exponent = point(arithmetic, value)
        if root_two:
            exponent = arithmetic.divide_integer(
                arithmetic.multiply(exponent, arithmetic.root_two), 2)
        total = arithmetic.add(total, arithmetic.exp_negative(exponent))
    return total


def pressure(arithmetic, d, orders, ex, ey):
    """Outward interval for G_relaxed in (7) of the predecessor note."""
    raw = arithmetic.add(point(arithmetic, d), arithmetic.multiply(
        point(arithmetic, 2), log_at_least_one(arithmetic, partition(arithmetic, orders))))
    physical = arithmetic.multiply(arithmetic.root_two, arithmetic.add(
        log_at_least_one(arithmetic, partition(arithmetic, ex, True)),
        log_at_least_one(arithmetic, partition(arithmetic, ey, True))))
    return arithmetic.subtract(raw, physical)


def as_decimal(arithmetic, value):
    return float(F(value[0], arithmetic.scale)), float(F(value[1], arithmetic.scale))


def verdict(arithmetic, d, ex, ey, families=None):
    families = families or {}
    orders, hull, caps = relaxed_orders(d, ex, ey, **families)
    interval = pressure(arithmetic, d, orders, ex, ey)
    record = {
        'alphabet': len(ex), 'deficit': str(F(d)),
        'paired_energies_X': [str(F(v)) for v in ex],
        'paired_energies_Y': [str(F(v)) for v in ey],
        'families': families or 'cell+shared line+rectangle',
        'average_normalized_orders': [str(v) for v in orders],
        'convex_hull': [[k, str(v)] for k, v in hull],
        'binding_cap_families': sorted({caps[k][1].split(' ')[0] for k, _ in hull if k > 1}),
        'pressure_dyadic_interval': [str(v) for v in interval],
        'pressure_decimal_interval': as_decimal(arithmetic, interval),
    }
    return interval, record


# ---------------------------------------------------------------------------
# Part 1 helpers: exact literal control over F_p.

def series_inverse(series, precision):
    inverse = [0] * precision
    inverse[0] = pow(series[0], PRIME - 2, PRIME)
    for n in range(1, precision):
        total = 0
        for k in range(1, min(n, len(series) - 1) + 1):
            total += series[k] * inverse[n - k]
        inverse[n] = (-total * inverse[0]) % PRIME
    return inverse


try:  # numpy is optional; the pure-Python path gives identical results.
    import numpy
except ImportError:  # pragma: no cover
    numpy = None


def series_multiply(a, b, precision):
    if numpy is not None:
        product = numpy.convolve(numpy.array(a[:precision], dtype=numpy.int64),
                                 numpy.array(b[:precision], dtype=numpy.int64))[:precision] % PRIME
        return [int(v) for v in product] + [0] * (precision - len(product))
    out = [0] * precision
    for i, x in enumerate(a[:precision]):
        if x:
            for j, y in enumerate(b[:precision - i]):
                out[i + j] += x * y
    return [v % PRIME for v in out]


def valuation(series):
    for index, value in enumerate(series):
        if value:
            return index
    return None


def smith_orders(matrix, precision):
    """Sorted Smith orders at w=0 of a nonsingular polynomial matrix.

    Exact elimination over F_p[[w]] modulo w^precision.  Valid when
    precision exceeds ord_0(det), which the caller guarantees through an
    independent degree bound; an order >= precision raises an error.
    """
    rows = [[(list(entry) + [0] * precision)[:precision] for entry in row] for row in matrix]
    orders = []
    while rows:
        best = None
        for i, row in enumerate(rows):
            for j, entry in enumerate(row):
                v = valuation(entry)
                if v is not None and (best is None or v < best[0]):
                    best = (v, i, j)
        require(best is not None, 'precision too small for Smith orders')
        v, pi, pj = best
        pivot_row = rows.pop(pi)
        pivot = pivot_row[pj]
        unit_inverse = series_inverse(pivot[v:], precision - v)
        remaining = []
        for row in rows:
            entry = row[pj]
            if any(entry):
                factor = series_multiply(entry[v:], unit_inverse, precision - v)
                for col in range(len(row)):
                    if col != pj and any(pivot_row[col]):
                        shifted = series_multiply(factor, pivot_row[col], precision)
                        row[col] = [(a - b) % PRIME for a, b in zip(row[col], shifted)]
            remaining.append(row[:pj] + row[pj + 1:])
        rows = remaining
        orders.append(v)
    return sorted(orders)


def poly_eval(poly, point_value):
    total = 0
    for coefficient in reversed(poly):
        total = (total * point_value + coefficient) % PRIME
    return total


def determinant_mod(matrix):
    n = len(matrix)
    a = [row[:] for row in matrix]
    det = 1
    for c in range(n):
        pivot = next((r for r in range(c, n) if a[r][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            det = -det
        det = det * a[c][c] % PRIME
        inv = pow(a[c][c], PRIME - 2, PRIME)
        for r in range(c + 1, n):
            if a[r][c]:
                f = a[r][c] * inv % PRIME
                a[r] = [(x - f * y) % PRIME for x, y in zip(a[r], a[c])]
    return det % PRIME


def interpolate(xs, ys):
    """Coefficients (low first) of the unique polynomial of degree < len(xs)."""
    n = len(xs)
    coefficients = [0] * n
    for i in range(n):
        basis = [1]
        denominator = 1
        for j in range(n):
            if j != i:
                basis = [((basis[k - 1] if k else 0) - xs[j] * (basis[k] if k < len(basis) else 0)) % PRIME
                         for k in range(len(basis) + 1)]
                denominator = denominator * (xs[i] - xs[j]) % PRIME
        scale = ys[i] * pow(denominator, PRIME - 2, PRIME) % PRIME
        for k in range(n):
            coefficients[k] = (coefficients[k] + scale * basis[k]) % PRIME
    return coefficients


def polynomial_minor(entries, rows, columns):
    """Exact det of a polynomial submatrix over F_p via interpolation."""
    bound = sum(max((len(entries[r][c]) - 1 for c in columns if any(entries[r][c])), default=0)
                for r in rows)
    xs = list(range(2, bound + 3))
    ys = [determinant_mod([[poly_eval(entries[r][c], x) for c in columns] for r in rows]) for x in xs]
    coefficients = interpolate(xs, ys)
    while coefficients and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def random_crossing(rng, q, degree, density):
    """Laurent entries {exponent: coefficient} with T(1)=[x=v][y=u]."""
    tensor = {}
    for a in itertools.product(range(q), repeat=4):
        x, y, u, v = a
        target = int(x == v and y == u)
        entry = {}
        if target or rng.random() < density:
            low = rng.randint(-degree, degree)
            high = rng.randint(low, low + degree)
            for e in range(low, high + 1):
                if rng.random() < 0.6:
                    entry[e] = rng.randrange(PRIME)
            e = rng.randint(low, high)
            entry[e] = (entry.get(e, 0) + target - sum(entry.values())) % PRIME
            entry = {k: c for k, c in entry.items() if c}
        if entry:
            tensor[a] = entry
    return tensor


def monomial_crossing(rng, q, degree):
    """Only matched entries, T_{x y y x}=w^e: a tight family for the caps."""
    return {(x, y, y, x): {rng.randint(-degree, degree): 1}
            for x in range(q) for y in range(q)}


def tight_crossing(q, heights):
    """Matched monomials whose infinity orders all vanish.

    Then d=0 and the zero-place profile is the sorted list of cell
    energies, so every cell-set cap on that profile holds with equality.
    """
    hinf = heights[1]
    return {(x, y, y, x): {hinf[0][x] + hinf[1][y] + hinf[2][y] + hinf[3][x]: 1}
            for x in range(q) for y in range(q)}


def literal_case(rng, q, tensor, heights, cell_samples):
    legs = range(4)
    h0, hinf = heights
    keys = sorted(tensor)

    def gauged(place):
        out = {}
        for a in keys:
            shift = sum((h0 if place == 0 else hinf)[l][a[l]] for l in legs)
            out[a] = {(e if place == 0 else -e) + shift: c for e, c in tensor[a].items()}
        return out
    zero, infinity = gauged(0), gauged(1)
    mu0 = min(min(entry) for entry in zero.values())
    muinf = min(min(entry) for entry in infinity.values())

    def dense(laurent, offset):
        if not laurent:
            return [0]
        top = max(laurent) - offset
        poly = [0] * (top + 1)
        for e, c in laurent.items():
            poly[e - offset] = c
        return poly
    # Normalized matrices R (zero place) and its infinity version.
    R = {a: dense(zero[a], mu0) for a in keys}
    Rinf = {a: dense(infinity[a], muinf) for a in keys}
    H = [[h0[l][s] + hinf[l][s] for s in range(q)] for l in legs]
    EX = [H[0][s] + H[3][s] for s in range(q)]
    EY = [H[1][s] + H[2][s] for s in range(q)]
    bx = min(range(q), key=EX.__getitem__)
    by = min(range(q), key=EY.__getitem__)
    cl = [[H[l][s] - H[l][(bx if l in (0, 3) else by)] for s in range(q)] for l in legs]
    cX = [EX[s] - EX[bx] for s in range(q)]
    cY = [EY[s] - EY[by] for s in range(q)]
    d = EX[bx] + EY[by] - mu0 - muinf
    require(d >= 0, 'normalization gives d>=0')
    for a in keys:
        bound = d + sum(cl[l][a[l]] for l in legs)
        require(len(R[a]) - 1 <= bound, 'degree bound (2) on R')
        # infinity entry equals z^(d+sum c) R(1/z): reversal of R padded to bound
        padded = R[a] + [0] * (bound + 1 - len(R[a]))
        reversed_entry = padded[::-1]
        while len(reversed_entry) > 1 and reversed_entry[-1] == 0:
            reversed_entry.pop()
        ri = Rinf[a] + [0] * (len(reversed_entry) - len(Rinf[a]))
        require(reversed_entry == ri[:len(reversed_entry)] and not any(ri[len(reversed_entry):]),
                'infinity matrix is the reciprocal normalization')
        require(poly_eval(R[a], 1) == int(a[0] == a[3] and a[1] == a[2]), 'endpoint')
    rows_cut = [(x, y) for x in range(q) for y in range(q)]

    def row_matrix(source):
        return [[source.get((x, y, u, v), [0]) for (u, v) in rows_cut] for (x, y) in rows_cut]

    def col_matrix(source):
        return [[source.get((x, y, u, v), [0]) for (y, v) in rows_cut] for (x, u) in rows_cut]
    out = {}
    for name, (Mzero, Minf) in {'row cut at zero': (row_matrix(R), col_matrix(Rinf)),
                                'column cut at zero': (col_matrix(R), row_matrix(Rinf))}.items():
        precision = [1 + sum(max(len(e) for e in row) - 1 for row in M) for M in (Mzero, Minf)]
        r = smith_orders(Mzero, precision[0])
        s = smith_orders(Minf, precision[1])
        require(r[0] == 0 and s[0] == 0, 'normalized first orders vanish')
        caps = all_caps(F(d), sorted(cX), sorted(cY))
        cell = sorted(x + y for x in cX for y in cY)
        sx, sy = sorted(cX), sorted(cY)
        tight = {'cell_r': 0, 'cell_s': 0, 'line': 0}
        for k in range(1, q * q + 1):
            cap = k * d + sum(cell[:k])
            require(sum(r[:k]) <= cap and sum(s[:k]) <= cap, 'cell-set prefix cap')
            tight['cell_r'] += sum(r[:k]) == cap
            tight['cell_s'] += sum(s[:k]) == cap
            require(F(sum(r[:k]) + sum(s[:k]), 2) <= caps[k][0], 'combined prefix cap')
            if k <= q:
                line = k * d + min(sum(sx[:k]), sum(sy[:k]))
                require(sum(r[:k]) + sum(s[:k]) <= line, 'shared-line cap')
                tight['line'] += sum(r[:k]) + sum(s[:k]) == line
        orders, hull = envelope({k: v for k, (v, _) in caps.items()})
        gamma = [F(a + b, 2) for a, b in zip(r, s)]
        require(all(sum(gamma[:k]) <= sum(orders[:k]) for k in range(q * q + 1)),
                'actual averaged profile lies below the convex minorant')
        out[name] = {'r': r, 's': s, 'tight_prefixes': tight}
    # Direct minor checks: random cell sets in both cuts.
    minors = 0
    for _ in range(cell_samples):
        k = rng.randint(1, min(4, q * q))
        cells = rng.sample(rows_cut, k)
        cost = k * d + sum(cX[x] + cY[y] for x, y in cells)
        # row cut: rows (x,y) in S, columns (u,v)=(y,x)
        rows = [(x, y) for x, y in cells]
        entries = {(x, y): {(u, v): R.get((x, y, u, v), [0]) for (u, v) in rows_cut} for (x, y) in rows}
        m = polynomial_minor(entries, rows, [(y, x) for x, y in cells])
        require(m and poly_eval(m, 1) in (1, PRIME - 1), 'row cell minor is +-1 at the endpoint')
        require(len(m) - 1 <= cost, 'row cell minor degree <= kd+Cell(S)')
        # column cut: rows (x,u) in S, columns (y,v)=(u,x); infinity version
        entries = {(x, u): {(y, v): Rinf.get((x, y, u, v), [0]) for (y, v) in rows_cut} for (x, u) in cells}
        mi = polynomial_minor(entries, list(cells), [(u, x) for x, u in cells])
        entries0 = {(x, u): {(y, v): R.get((x, y, u, v), [0]) for (y, v) in rows_cut} for (x, u) in cells}
        m0 = polynomial_minor(entries0, list(cells), [(u, x) for x, u in cells])
        require(len(m0) - 1 <= cost, 'column cell minor degree <= kd+Cell(S)')
        recip = (m0 + [0] * (cost + 1 - len(m0)))[::-1]
        while recip and recip[-1] == 0:
            recip.pop()
        require(recip == mi, 'infinity cell minor equals z^(kd+Cell) m(1/z)')
        require(valuation(mi) <= cost, 'infinity cell-minor order <= kd+Cell(S)')
        minors += 2
    return {'deficit': d, 'paired_energies_X': cX, 'paired_energies_Y': cY,
            'profiles': out, 'cell_minors_checked': minors}


def literal_control(seed, cases):
    rng = random.Random(seed)
    summary = {'tensors': 0, 'complete_profiles': 0, 'cell_minors': 0,
               'tight_cell_prefixes': 0, 'tight_line_prefixes': 0,
               'nonzero_r_tail': 0, 'nonzero_s_tail': 0}
    digest = hashlib.sha256()
    for q, degree, density, family in cases:
        heights = tuple([[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2))
        if family == 'monomial':
            tensor = monomial_crossing(rng, q, degree)
        elif family == 'tight':
            tensor = tight_crossing(q, heights)
        else:
            tensor = random_crossing(rng, q, degree, density)
        record = literal_case(rng, q, tensor, heights, cell_samples=2 if q == 3 else 1)
        if family == 'tight':
            zero = record['profiles']['row cut at zero']
            cell = sorted(x + y for x in record['paired_energies_X'] for y in record['paired_energies_Y'])
            require(record['deficit'] == 0 and zero['r'] == cell and not any(zero['s']),
                    'tight family attains every cell-set cap on the zero profile')
            summary['tight_family_tensors'] = summary.get('tight_family_tensors', 0) + 1
        summary['tensors'] += 1
        summary['complete_profiles'] += 4
        summary['cell_minors'] += record['cell_minors_checked']
        for prof in record['profiles'].values():
            summary['tight_cell_prefixes'] += prof['tight_prefixes']['cell_r'] + prof['tight_prefixes']['cell_s']
            summary['tight_line_prefixes'] += prof['tight_prefixes']['line']
            summary['nonzero_r_tail'] += any(prof['r'][1:])
            summary['nonzero_s_tail'] += any(prof['s'][1:])
        digest.update(json.dumps(record, sort_keys=True, default=str).encode())
    summary['record_sha256'] = digest.hexdigest()
    return summary


# ---------------------------------------------------------------------------
# Part 3: compact ternary cover by the monotone corner bound.

def corner_lower_bound(arithmetic, box, q):
    """Lower bound of G_relaxed on a box of (d, X energies, Y energies).

    The caps are nondecreasing in d and every energy, so the hull orders at
    the upper corner dominate; Z(gamma*) is then smallest there.  The
    physical subtraction is largest at the lower corner, and d is smallest
    there.  Energies are sorted inside all_caps, so unsorted corners are
    harmless; the zero baseline energy on each side is fixed.
    """
    (dl, dh), xs, ys = box
    upper_x = [F(0)] + [h for _, h in xs]
    upper_y = [F(0)] + [h for _, h in ys]
    lower_x = [F(0)] + [l for l, _ in xs]
    lower_y = [F(0)] + [l for l, _ in ys]
    orders, _, _ = relaxed_orders(dh, upper_x, upper_y)
    raw = arithmetic.add(point(arithmetic, dl), arithmetic.multiply(
        point(arithmetic, 2), log_at_least_one(arithmetic, partition(arithmetic, orders))))
    physical = arithmetic.multiply(arithmetic.root_two, arithmetic.add(
        log_at_least_one(arithmetic, partition(arithmetic, lower_x, True)),
        log_at_least_one(arithmetic, partition(arithmetic, lower_y, True))))
    return arithmetic.subtract(raw, physical)[0]


def split(box):
    (dl, dh), xs, ys = box
    widths = [dh - dl] + [h - l for l, h in xs] + [h - l for l, h in ys]
    i = max(range(len(widths)), key=widths.__getitem__)
    intervals = [(dl, dh)] + list(xs) + list(ys)
    lo, hi = intervals[i]
    mid = (lo + hi) / 2
    out = []
    for part in ((lo, mid), (mid, hi)):
        new = list(intervals)
        new[i] = part
        nx = len(xs)
        out.append((new[0], tuple(new[1:1 + nx]), tuple(new[1 + nx:])))
    return out


def ordered_box_possible(box):
    """Keep only boxes meeting the fundamental domain of the symmetries.

    G_relaxed depends on each paired energy list only through its sorted
    order, and it is invariant under exchanging the X and Y lists.  Every
    point therefore has an image with both lists sorted and with a1<=b1,
    where a1,b1 are the first nonbaseline energies.  A box is discarded
    only if it misses that fundamental domain.  The corner bound remains
    valid on unsorted points because all_caps sorts.
    """
    _, xs, ys = box
    sorted_x = all(xs[i][0] <= xs[i + 1][1] for i in range(len(xs) - 1))
    sorted_y = all(ys[i][0] <= ys[i + 1][1] for i in range(len(ys) - 1))
    return sorted_x and sorted_y and xs[0][0] <= ys[0][1]


def compact_cover(arithmetic, q, d_max, e_max, margin, leaf_limit):
    root = ((F(0), F(d_max)), tuple((F(0), F(e_max)) for _ in range(q - 1)),
            tuple((F(0), F(e_max)) for _ in range(q - 1)))
    stack = [root]
    leaves = []
    threshold = point(arithmetic, margin)[1]
    while stack:
        box = stack.pop()
        if not ordered_box_possible(box):
            continue
        lower = corner_lower_bound(arithmetic, box, q)
        if lower >= threshold:
            leaves.append((box, lower))
            require(len(leaves) <= leaf_limit, 'compact cover leaf limit exceeded')
            continue
        stack.extend(split(box))
    return leaves


def audit_cover(leaves, q, d_max, e_max):
    """Independent cover audit: dyadic grid volume of sorted region covered."""
    # Every leaf is a product of dyadic subintervals of the root.  Check that
    # leaf volumes (restricted to the ordered region test) account for the
    # whole root by recursive re-splitting of the root against the leaf set.
    index = {}
    for box, _ in leaves:
        index[box] = True
    root = ((F(0), F(d_max)), tuple((F(0), F(e_max)) for _ in range(q - 1)),
            tuple((F(0), F(e_max)) for _ in range(q - 1)))
    stack = [root]
    visited = 0
    reached = set()
    while stack:
        box = stack.pop()
        visited += 1
        if box in index:
            require(box not in reached, 'leaf reached twice')
            reached.add(box)
            continue
        if not ordered_box_possible(box):
            continue
        require(visited < 50 * len(leaves) + 10, 'cover audit diverged: missing leaf')
        stack.extend(split(box))
    require(len(reached) == len(index), 'every saved leaf lies in the split tree')
    return visited


def ternary_cover(e_max, margin):
    d_max = F(4)
    arithmetic = Dyadic(60)
    leaves = compact_cover(arithmetic, 3, d_max, e_max, margin, 400000)
    visited = audit_cover(leaves, 3, d_max, e_max)
    # Independent re-evaluation of every leaf at a different precision.
    check = Dyadic(72)
    threshold = point(check, margin)[1]
    worst = None
    for box, _ in leaves:
        lower = corner_lower_bound(check, box, 3)
        require(lower >= threshold, 'leaf bound fails on re-evaluation')
        worst = lower if worst is None else min(worst, lower)
    # Beyond d_max: Z(gamma*)>=1 and each physical partition is at most 3.
    tail = check.subtract(point(check, d_max), check.multiply(
        check.multiply(point(check, 2), check.root_two), log_at_least_one(check, point(check, 3))))
    require(tail[0] > threshold, 'trivial large-deficit tail')
    return {'alphabet': 3, 'deficit_box': ['0', str(d_max)], 'energy_box': ['0', str(e_max)],
            'margin': str(margin), 'leaves': len(leaves), 'audit_nodes': visited,
            'smallest_leaf_lower_bound_decimal': float(F(worst, check.scale)),
            'large_deficit_tail_lower_bound_decimal': float(F(tail[0], check.scale)),
            'leaf_sha256': hashlib.sha256(json.dumps(
                [[str(v) for v in [b[0][0], b[0][1]] + [w for iv in b[1] + b[2] for w in iv]]
                 for b, _ in leaves]).encode()).hexdigest(),
            'statement': ('For q=3, every point with all paired energies in [0,%s] and any '
                          'd>=0 has G_relaxed >= %s.' % (e_max, margin))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--quick', action='store_true', help='smaller literal control')
    parser.add_argument('--skip-cover', action='store_true', help='omit the ternary compact cover')
    parser.add_argument('--ternary-energy', default='1', help='energy box edge for the ternary cover')
    args = parser.parse_args()
    arithmetic = Dyadic(120)
    result = {'scope': 'Necessary prefix inequalities for scalar Laurent strict crossings. '
              'Relaxation points are not tensors; nothing here is a crossing certificate '
              'or a factoring improvement.', 'dyadic_precision': 120}

    # Part 1: literal control.
    cases = []
    for trial in range(6 if args.quick else 24):
        cases.append((2, 2 + trial % 3, 0.5, 'random'))
        cases.append((3, 1 + trial % 3, 0.3 + 0.1 * (trial % 3), 'random'))
        cases.append((3, 3, 0, 'monomial'))
        cases.append((3, 0, 0, 'tight'))
    for trial in range(2 if args.quick else 4):
        cases.append((4, 1, 0.15, 'random'))
        cases.append((4, 2, 0, 'monomial'))
    summary = literal_control(2026103101, cases)
    require(summary['nonzero_r_tail'] > 0 and summary['nonzero_s_tail'] > 0,
            'literal control must exercise nontrivial profiles')
    require(summary['tight_cell_prefixes'] > 0, 'literal control must reach a tight cell cap')
    result['literal_control'] = summary
    print('literal control:', json.dumps(summary), flush=True)

    # Part 2: predecessor targets under the strengthened caps.
    targets = []
    q4 = (F(1), [F(0), F(0), F(2), F(8, 3)])
    q5 = (F(1), [F(0), F(0), F(2), F(8, 3), F(3)])
    for (d, e), expect in ((q4, 'excluded'), (q5, 'negative')):
        old, _ = verdict(arithmetic, d, e, e, {'cells': False})
        new, record = verdict(arithmetic, d, e, e)
        record['predecessor_pressure_dyadic_interval'] = [str(v) for v in old]
        require(old[1] < 0, 'predecessor target negative without cell caps')
        if new[0] > 0:
            record['status'] = 'Excluded by the cell-set caps: G_relaxed > 0.'
        elif new[1] < 0:
            record['status'] = 'Still a negative relaxation point with cell caps.'
        else:
            require(False, 'undecided sign at a predecessor target')
        if expect == 'excluded':
            require(new[0] > 0, 'q=4 predecessor target should be excluded')
        else:
            require(new[1] < 0, 'q=5 predecessor target should remain negative')
        targets.append(record)
        print(record['status'], record['alphabet'], record['pressure_decimal_interval'], flush=True)
    result['predecessor_targets'] = targets

    # Part 2b: sharper negative targets of the strengthened relaxation.
    # The data were located by a labelled floating-point search (see the
    # note); only the exact verdict below is claimed.
    sharper = []
    for d, e in ((F(7, 6), [F(0), F(0), F(19, 9), F(23, 9), F(37, 9)]),
                 (F(8, 7), [F(0), F(0), F(15, 8), F(27, 10), F(4), F(4)])):
        interval, record = verdict(arithmetic, d, e, e)
        require(interval[1] < 0, 'sharper target must have certified negative pressure')
        old, _ = verdict(arithmetic, d, e, e, {'cells': False})
        record['pressure_without_cell_caps_dyadic_interval'] = [str(v) for v in old]
        denominators = [F(v).denominator for v in [d] + e + list(map(F, record['average_normalized_orders']))]
        scale = 1
        for den in denominators:
            scale = scale * den // __import__('math').gcd(scale, den)
        record['common_denominator'] = scale
        record['status'] = ('Negative relaxation point satisfying every cell-set and shared-line '
                            'cap (rectangular caps are implied); no realizing tensor is known.')
        sharper.append(record)
        print('sharper target', record['alphabet'], record['pressure_decimal_interval'], flush=True)
    result['sharper_targets'] = sharper

    # Part 3: compact ternary region.
    if not args.skip_cover:
        cover = ternary_cover(F(args.ternary_energy), F(1, 100))
        result['ternary_compact_cover'] = cover
        print('ternary cover:', json.dumps(cover), flush=True)
    result['status'] = ('Cell-set caps passed the literal F_p control; the q=4 predecessor target is '
                        'excluded; the q=5 predecessor target and the sharper q=5 and q=6 targets '
                        'remain certified negative relaxation points'
                        + ('' if args.skip_cover else '; the compact ternary region is certified positive') + '.')
    args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
