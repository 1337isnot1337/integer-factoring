"""Optional native weighted Smith search over F4; Python uses only the stdlib.

Field coefficients are ``lo + omega*hi``, where omega**2 = omega + 1.
An entry is a pair of nonnegative polynomial bit masks. Dense entries use
index ``x + q*y + q**2*u + q**3*v``. On first use, a C++17 compiler builds
the adjacent source into a temporary, content-addressed cache.

Polynomial arithmetic is exact, but pivot weights and pressures use floats.
This is a SEARCH tool: a proposed improvement still needs an exact rational
replay. The 512-coefficient series window returns an error when it cannot
resolve a full profile. It never substitutes an incomplete profile.

``Oracle(values, q).evaluate(params, probabilities, want_spectra=True)``
returns ``(log_pressure, gradient, spectra)`` as Python values/lists. A search
program may convert these lists to NumPy arrays; neither NumPy nor SciPy is
needed here. Params are infinity heights, zero heights, then log(s)>0.
Each group of heights lists legs x,y,u,v and labels 1,...,q-1; label 0 is
fixed at height zero. Heights have ALREADY been multiplied by log(s).
Cut 1 uses C at infinity and R at zero; cut 0 reverses these choices.

``denominator_degree=d`` restores division by a common denominator of
degree d that is a unit at zero. This arithmetic adjustment does not prove
that a rational tensor is admissible. Supply its typed polynomiality proof
separately, and set endpoint_check=False for a noncrossing numerator.

The native library holds one active seed. Constructing another oracle makes
an older one inactive; call its activate() method before reusing it.
"""

from __future__ import annotations

import ctypes
import hashlib
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import threading


PRECISION = 512
_LIBRARY = None
_GENERATION = 0
_LOCK = threading.RLock()
_DOUBLE = ctypes.POINTER(ctypes.c_double)
_U64 = ctypes.POINTER(ctypes.c_uint64)
_I64 = ctypes.POINTER(ctypes.c_int64)


def _load_library():
    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    source = Path(__file__).with_suffix('.cpp')
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional native oracle needs a C++17 compiler; set CXX if needed.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC']
    try:
        cpu_flags = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu_flags = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu_flags:
        flags += ['-mpclmul', '-msse2']
    digest = hashlib.sha256(source.read_bytes() + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-native-f4-' + str(os.getuid()))
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = cache / (digest + '.so')
    if not target.exists():
        descriptor, temporary = tempfile.mkstemp(prefix=digest + '-', suffix='.so', dir=cache)
        os.close(descriptor)
        try:
            result = subprocess.run(compiler + flags + [str(source), '-o', temporary], capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError('Native F4 compilation failed:\n' + result.stderr[-12000:])
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    lib = ctypes.CDLL(str(target))
    lib.gf4_set_seed.argtypes = [ctypes.c_int, _U64, ctypes.c_int, ctypes.c_int]
    lib.gf4_set_seed.restype = ctypes.c_int
    lib.gf4_set_paired_seed.argtypes = [ctypes.c_int, _U64, ctypes.c_int]
    lib.gf4_set_paired_seed.restype = ctypes.c_int
    lib.gf4_cost_gradient.argtypes = [_DOUBLE, _DOUBLE, ctypes.c_int, _DOUBLE, _DOUBLE, _DOUBLE]
    lib.gf4_cost_gradient.restype = ctypes.c_int
    lib.gf4_prefixes.argtypes = [_DOUBLE, ctypes.c_int, _DOUBLE, _I64, _U64, _U64]
    lib.gf4_prefixes.restype = ctypes.c_int
    lib.gf4_paired_entry.argtypes = [ctypes.c_int] * 4 + [_U64]
    lib.gf4_paired_entry.restype = ctypes.c_int
    _LIBRARY = lib
    return lib


def _pack(values, expected, maximum_degree=PRECISION - 1):
    values = list(values)
    if len(values) != expected:
        raise ValueError(f'Expected {expected} dense tensor entries, got {len(values)}.')
    pairs = []
    for pair in values:
        if len(pair) != 2 or any(not isinstance(v, int) or v < 0 for v in pair):
            raise ValueError('Each entry must contain two nonnegative integer polynomial bit masks.')
        pairs.append(tuple(pair))
    bits = max((v.bit_length() for pair in pairs for v in pair), default=0)
    if bits > maximum_degree + 1:
        raise ValueError(f'Polynomial degree exceeds {maximum_degree}.')
    words = max(1, (bits + 63) // 64)
    mask = (1 << 64) - 1
    packed = (ctypes.c_uint64 * (2 * expected * words))(
        *((v >> (64 * k)) & mask for pair in pairs for v in pair for k in range(words))
    )
    return packed, words


def _logsumexp(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(v - maximum) for v in values))


class Oracle:
    def __init__(self, values, q, endpoint_check=True, denominator_degree=0):
        if not isinstance(q, int) or not 2 <= q <= 32:
            raise ValueError('The physical alphabet q must be an integer from 2 through 32.')
        self._initialize_dimensions(q, denominator_degree)
        self._packed, self._words = _pack(values, q ** 4)
        self._endpoint_check = bool(endpoint_check)
        self._bulk_modes = None
        self.activate()

    def _initialize_dimensions(self, q, denominator_degree):
        if not isinstance(denominator_degree, int) or denominator_degree < 0:
            raise ValueError('denominator_degree must be a nonnegative integer.')
        self.q, self.n = q, q * q
        self.stride = 4 * (q - 1)
        self.parameters = 2 * self.stride + 1
        self.denominator_degree = denominator_degree

    def activate(self):
        """Reload this retained seed after another Oracle has used the library."""
        global _GENERATION
        with _LOCK:
            lib = _load_library()
            _GENERATION += 1
            self._generation = None
            if self._bulk_modes is None:
                status = lib.gf4_set_seed(self.q, self._packed, self._words, int(self._endpoint_check))
            else:
                status = lib.gf4_set_paired_seed(self._bulk_modes, self._packed, self._words)
            if status:
                raise ArithmeticError(f'Native seed rejected (status {status}).')
            self._generation = _GENERATION
        return self

    def _check_active(self):
        if self._generation != _GENERATION:
            raise RuntimeError('Another seed replaced this oracle; call activate() before reuse.')

    def _params(self, params, cut):
        x = list(map(float, params))
        if len(x) != self.parameters or cut not in (0, 1):
            raise ValueError('Wrong parameter length or cut (expected 0 or 1).')
        if not all(math.isfinite(v) for v in x) or x[-1] <= 0:
            raise ValueError('Parameters must be finite and log(s) must be positive.')
        if max(map(abs, x[:-1]), default=0) / x[-1] > 1_000_000:
            raise ValueError('Normalized heights exceed the safe native integer-conversion range.')
        return x, (ctypes.c_double * len(x))(*x)

    def evaluate(self, params, probabilities, cut=1, want_spectra=False):
        x, native_x = self._params(params, cut)
        pi = list(map(float, probabilities))
        if len(pi) != self.q or any(not math.isfinite(v) or v < 0 for v in pi) or abs(sum(pi) - 1) > 1e-9:
            raise ValueError('probabilities must be a nonnegative finite distribution on q letters.')
        native_pi = (ctypes.c_double * self.q)(*pi)
        cost = ctypes.c_double()
        gradient = (ctypes.c_double * self.parameters)()
        spectra = (ctypes.c_double * (2 * self.n))()
        with _LOCK:
            self._check_active()
            status = _load_library().gf4_cost_gradient(native_x, native_pi, cut, ctypes.byref(cost), gradient, spectra)
        if status:
            raise ArithmeticError(f'Full profile unresolved (status {status}); increase precision or use exact replay.')
        value = cost.value - self.denominator_degree * x[-1]
        grad, orders = list(gradient), list(spectra)
        grad[-1] -= self.denominator_degree
        orders[:self.n] = [v + self.denominator_degree * x[-1] for v in orders[:self.n]]
        if not all(math.isfinite(v) for v in [value] + grad + orders):
            raise ArithmeticError('The floating pressure calculation overflowed.')
        return (value, grad, orders) if want_spectra else (value, grad)

    def prefixes(self, params, cut=1):
        """Return orders and attaining prefix-minor row/column masks (q<=8)."""
        if self.q > 8:
            raise ValueError('Prefix masks support q<=8 only.')
        x, native_x = self._params(params, cut)
        orders = (ctypes.c_double * (2 * self.n))()
        bases = (ctypes.c_int64 * (2 * self.n))()
        rows = (ctypes.c_uint64 * (2 * self.n))()
        cols = (ctypes.c_uint64 * (2 * self.n))()
        with _LOCK:
            self._check_active()
            status = _load_library().gf4_prefixes(native_x, cut, orders, bases, rows, cols)
        if status:
            raise ArithmeticError(f'Prefix profile unresolved (status {status}).')
        orders, bases = list(orders), list(bases)
        for i in range(self.n):
            orders[i] += self.denominator_degree * x[-1]
            bases[i] += self.denominator_degree * (i + 1)
        return orders, bases, list(rows), list(cols)

    def raw_evaluate(self, params, cut=1, want_spectra=False):
        """Return log(U*V) and its gradient before physical type restoration.

        This supports pair types or other externally supplied mean-height
        constraints without pretending that zero probabilities form a type.
        Any common-denominator restoration is still included.
        """
        x, _ = self._params(params, cut)
        cost, gradient, spectra = self.evaluate(x, [1 / self.q] * self.q, cut, True)
        cost -= sum(x[:-1]) / self.q
        gradient[:-1] = [v - 1 / self.q for v in gradient[:-1]]
        return (cost, gradient, spectra) if want_spectra else (cost, gradient)

    def types(self, params, cut=1):
        """Numerically optimize the common letter distribution by Gibbs duality."""
        x, _ = self._params(params, cut)
        weights = [0.] + [sum(x[self.stride * end + (self.q - 1) * leg + a - 1]
                                   for end in range(2) for leg in range(4))
                          for a in range(1, self.q)]
        _, _, spectra = self.evaluate(x, [1 / self.q] * self.q, cut, True)
        raw = _logsumexp([-v for v in spectra[:self.n]]) + _logsumexp([-v for v in spectra[self.n:]])
        minimum = min(weights)
        constant = minimum + raw
        def equation(rate):
            return constant - rate * _logsumexp([-(v - minimum) / rate for v in weights])
        low, high = 1e-12, 8.
        if equation(low) <= 0:
            raise ArithmeticError('No strictly positive type root at the supported numerical scale.')
        while equation(high) > 0:
            high *= 2
            if high > 1e12:
                raise ArithmeticError('The numerical type root could not be bracketed.')
        for _ in range(100):
            middle = (low + high) / 2
            if equation(middle) > 0:
                low = middle
            else:
                high = middle
        rate = (low + high) / 2
        normalization = _logsumexp([-v / rate for v in weights])
        pi = [math.exp(-v / rate - normalization) for v in weights]
        return rate / 4, pi, spectra, weights


class PairedOracle(Oracle):
    """Conjugate two copies with one shared core bit and independent bulk bits.

    one_copy_values has alphabet 2**(bulk_modes+1), with the core bit highest.
    The paired alphabet is 2**(2*bulk_modes+1). Only bulk_modes 1 and 2 are
    implemented. This initializer requires the two conserved parity blocks;
    a coefficient outside their four tensor components is rejected.
    """
    def __init__(self, one_copy_values, bulk_modes, denominator_degree=4):
        if bulk_modes not in (1, 2):
            raise ValueError('bulk_modes must be 1 or 2.')
        self._initialize_dimensions(1 << (2 * bulk_modes + 1), denominator_degree)
        self._packed, self._words = _pack(one_copy_values, (1 << (bulk_modes + 1)) ** 4, 255)
        self._bulk_modes = bulk_modes
        self.activate()

    def entry(self, x, y, u, v):
        if any(not isinstance(a, int) or not 0 <= a < self.q for a in (x, y, u, v)):
            raise ValueError('Letter outside the physical alphabet.')
        values = (ctypes.c_uint64 * 16)()
        with _LOCK:
            self._check_active()
            status = _load_library().gf4_paired_entry(x, y, u, v, values)
        if status:
            raise ArithmeticError(f'Entry query failed (status {status}).')
        return [sum(int(values[8 * plane + k]) << (64 * k) for k in range(8)) for plane in range(2)]
