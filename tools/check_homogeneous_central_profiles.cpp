// Exact F_2[Z] Pfaffian profiles for bounded finite homogeneous experiments.
// This is a finite checker, not an asymptotic factoring implementation.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__PCLMUL__)
#include <wmmintrin.h>
#endif

struct Poly {
  static constexpr int capacity = 16;
  std::array<uint64_t, capacity> data{};
  int used = 0;
  void normalize() { while (used && !data[used - 1]) --used; }
  int degree() const {
    return used ? 64 * (used - 1) + 63 - __builtin_clzll(data[used - 1]) : -1;
  }
  bool empty() const { return !used; }
  void flip(int exponent) {
    if (exponent < 0 || exponent >= 64 * capacity)
      throw std::runtime_error("polynomial capacity exceeded");
    data[exponent / 64] ^= uint64_t(1) << (exponent % 64);
    used = std::max(used, exponent / 64 + 1);
    normalize();
  }
  static Poly monomial(int exponent) { Poly p; p.flip(exponent); return p; }
  void add(const Poly &other) {
    for (int i = 0; i < other.used; ++i) data[i] ^= other.data[i];
    used = std::max(used, other.used);
    normalize();
  }
  void add_shifted(const Poly &other, int shift) {
    if (other.empty()) return;
    if (other.degree() + shift >= 64 * capacity)
      throw std::runtime_error("shift capacity exceeded");
    const int words = shift / 64, bits = shift % 64;
    for (int i = 0; i < other.used; ++i) {
      data[i + words] ^= other.data[i] << bits;
      if (bits) {
        const uint64_t high = other.data[i] >> (64 - bits);
        if (i + words + 1 < capacity) data[i + words + 1] ^= high;
        else if (high) throw std::runtime_error("shift capacity exceeded");
      }
    }
    used = std::min(capacity, std::max(used, other.used + words + int(bits != 0)));
    normalize();
  }
  Poly shifted(int shift) const { Poly p; p.add_shifted(*this, shift); return p; }
  int parity() const {
    int value = 0;
    for (int i = 0; i < used; ++i) value ^= __builtin_parityll(data[i]);
    return value;
  }
};

Poly multiply(const Poly &left, const Poly &right) {
  Poly answer;
  if (left.empty() || right.empty()) return answer;
  if (left.degree() + right.degree() >= 64 * Poly::capacity)
    throw std::runtime_error("product capacity exceeded");
  for (int i = 0; i < left.used; ++i) {
#if defined(__PCLMUL__)
    const __m128i a = _mm_set_epi64x(0, left.data[i]);
#endif
    for (int j = 0; j < right.used; ++j) {
#if defined(__PCLMUL__)
      const __m128i b = _mm_set_epi64x(0, right.data[j]);
      const __m128i c = _mm_clmulepi64_si128(a, b, 0);
      alignas(16) uint64_t halves[2];
      _mm_store_si128(reinterpret_cast<__m128i *>(halves), c);
#else
      uint64_t halves[2] = {0, 0};
      uint64_t remaining = right.data[j];
      while (remaining) {
        const int shift = __builtin_ctzll(remaining);
        halves[0] ^= left.data[i] << shift;
        if (shift) halves[1] ^= left.data[i] >> (64 - shift);
        remaining &= remaining - 1;
      }
#endif
      answer.data[i + j] ^= halves[0];
      if (i + j + 1 < Poly::capacity) answer.data[i + j + 1] ^= halves[1];
      else if (halves[1]) throw std::runtime_error("product capacity exceeded");
    }
  }
  answer.used = std::min(Poly::capacity, left.used + right.used);
  answer.normalize();
  return answer;
}

Poly divide_exact(Poly numerator, const Poly &denominator) {
  if (denominator.empty()) throw std::runtime_error("zero denominator");
  if (denominator.degree() == 0) return numerator;
  Poly quotient;
  const int degree = denominator.degree();
  while (numerator.degree() >= degree) {
    const int shift = numerator.degree() - degree;
    quotient.flip(shift);
    numerator.add_shifted(denominator, shift);
  }
  if (!numerator.empty()) throw std::runtime_error("inexact Pfaffian division");
  return quotient;
}

Poly binomial_poly(int exponent) {
  Poly answer = Poly::monomial(0);
  for (int shift = 1; exponent; exponent >>= 1, shift <<= 1)
    if (exponent & 1) answer.add(answer.shifted(shift));
  return answer;
}

Poly pfaffian(std::vector<std::vector<Poly>> matrix) {
  const int n = matrix.size();
  if (n % 2) return {};
  if (!n) return Poly::monomial(0);
  Poly previous = Poly::monomial(0), pivot;
  for (int k = 0; k < n; k += 2) {
    int found = k + 1;
    while (found < n && matrix[k][found].empty()) ++found;
    if (found == n) return {};
    if (found != k + 1) {
      std::swap(matrix[found], matrix[k + 1]);
      for (auto &row : matrix) std::swap(row[found], row[k + 1]);
    }
    pivot = matrix[k][k + 1];
    for (int i = k + 2; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        Poly numerator = multiply(pivot, matrix[i][j]);
        numerator.add(multiply(matrix[k][i], matrix[k + 1][j]));
        numerator.add(multiply(matrix[k][j], matrix[k + 1][i]));
        matrix[i][j] = divide_exact(numerator, previous);
        matrix[j][i] = matrix[i][j];
      }
    }
    previous = pivot;
  }
  return pivot;
}

using Partition = std::vector<int>;

int vacuum(const Partition &a, const Partition &b) {
  if (a.size() != b.size()) return 0;
  std::vector<uint64_t> rows(a.size());
  for (int i = 0; i < int(a.size()); ++i)
    for (int j = 0; j < int(b.size()); ++j)
      if (!(((a[i] - 1) / 2) & ((b[j] - 1) / 2))) rows[i] ^= uint64_t(1) << j;
  for (int j = 0; j < int(rows.size()); ++j) {
    int pivot = j;
    while (pivot < int(rows.size()) && !((rows[pivot] >> j) & 1)) ++pivot;
    if (pivot == int(rows.size())) return 0;
    std::swap(rows[pivot], rows[j]);
    for (int i = j + 1; i < int(rows.size()); ++i)
      if ((rows[i] >> j) & 1) rows[i] ^= rows[j];
  }
  return 1;
}

Poly profile(const Partition &a, const Partition &u,
             const Partition &b, const Partition &y,
             const std::vector<Poly> &binomials) {
  std::vector<std::pair<int, int>> vertices;
  for (int t : u) vertices.emplace_back(0, t);
  for (int t : a) vertices.emplace_back(1, t);
  for (int t : b) vertices.emplace_back(2, t);
  for (int t : y) vertices.emplace_back(3, t);
  std::vector<std::vector<Poly>> matrix(vertices.size(),
                                       std::vector<Poly>(vertices.size()));
  for (int i = 0; i < int(vertices.size()); ++i) {
    const auto [first, n] = vertices[i];
    for (int k = i + 1; k < int(vertices.size()); ++k) {
      const auto [second, m] = vertices[k];
      Poly value;
      if (first == 0 && second == 1 && !((m - 1) & ~(n - 1)))
        value = binomials[m].shifted(n);
      else if (first == 1 && second == 3 && (n & m) == 1)
        value = binomials[n];
      else if (first == 1 && second == 2 && (n & m) == 1)
        value = Poly::monomial(n);
      else if (first == 0 && second == 2 && !((m - 1) & ~(n - 1)))
        value = binomials[m];
      else if (first == 2 && second == 3 && (n & m) == 1)
        value = binomials[n].shifted(m);
      else if (first == 0 && second == 3 && n == m)
        value = Poly::monomial(m);
      matrix[i][k] = matrix[k][i] = value;
    }
  }
  return pfaffian(std::move(matrix));
}

int binary_rank(std::vector<uint64_t> matrix, int rows, int columns) {
  const int words = (columns + 63) / 64;
  std::vector<int> pivots(columns, -1);
  int rank = 0;
  for (int row = 0; row < rows; ++row) {
    uint64_t *current = &matrix[row * words];
    while (true) {
      int word = words - 1;
      while (word >= 0 && !current[word]) --word;
      if (word < 0) break;
      int column = word * 64 + 63 - __builtin_clzll(current[word]);
      if (pivots[column] < 0) {
        pivots[column] = row;
        ++rank;
        break;
      }
      const uint64_t *other = &matrix[pivots[column] * words];
      for (int i = 0; i <= word; ++i) current[i] ^= other[i];
    }
  }
  return rank;
}

int main(int argc, char **argv) try {
  const bool emit_profiles = argc > 1 && std::string(argv[1]) == "--emit-profiles";
  int j, h, energy, d;
  if (!(std::cin >> j >> h >> energy >> d)) throw std::runtime_error("missing input header");
  if (j < 1 || j > h || energy < 1 || d < 1 || d > 64)
    throw std::runtime_error("invalid bounded experiment parameters");
  const int degree_bound = 2 * energy + h + j;
  if (2 * degree_bound + 2 >= 64 * Poly::capacity)
    throw std::runtime_error("requested degree exceeds fixed checker capacity");
  auto read_domain = [&](int expected) {
    std::vector<Partition> domain(d);
    for (auto &part : domain) {
      int count, sum = 0;
      if (!(std::cin >> count) || count < 1 || count > 32)
        throw std::runtime_error("invalid partition size");
      int previous = -1;
      for (int k = 0; k < count; ++k) {
        int value;
        std::cin >> value;
        if (value <= previous || value % 2 != 1) throw std::runtime_error("nonregular partition");
        part.push_back(value); sum += value; previous = value;
      }
      if (sum != expected) throw std::runtime_error("wrong partition weight");
    }
    return domain;
  };
  const auto operators = read_domain(j), partners = read_domain(h), states = read_domain(energy);
  int largest = 0;
  for (const auto *domain : {&operators, &partners, &states})
    for (const auto &parts : *domain) for (int value : parts) largest = std::max(largest, value);
  largest = std::max(largest, degree_bound);
  std::vector<Poly> binomials(largest + 1);
  for (int value = 0; value <= largest; ++value) binomials[value] = binomial_poly(value);
  const int size = d * d, words = (size + 63) / 64;
  std::vector<std::vector<uint64_t>> matrices(degree_bound + 1,
                                             std::vector<uint64_t>(size * words));
  std::vector<Poly> emitted;
  if (emit_profiles) emitted.reserve(size * size);
  int minimum_residual_valuation = degree_bound + 1;
  std::vector<int> jet_indices;
  for (int k = 1; k <= degree_bound; k *= 2) jet_indices.push_back(k);
  std::vector<std::vector<uint64_t>> jets(jet_indices.size(),
                                         std::vector<uint64_t>(size * words));
  const auto started = std::chrono::steady_clock::now();
  uint64_t checked = 0;
  for (int x = 0; x < d; ++x) {
    for (int y = 0; y < d; ++y) for (int v = 0; v < d; ++v) for (int u = 0; u < d; ++u) {
      const Poly value = profile(operators[x], states[u], partners[v], states[y], binomials);
      if (value.degree() > degree_bound) throw std::runtime_error("profile degree bound failed");
      const int crossing = int(y == u && vacuum(operators[x], partners[v]));
      if (value.parity() != crossing)
        throw std::runtime_error("raw summed crossing failed");
      ++checked;
      if (emit_profiles) emitted.push_back(value);
      Poly residual = value, translated;
      if (crossing) residual.flip(energy + j);
      for (int word = 0; word < residual.used; ++word) {
        uint64_t remaining = residual.data[word];
        while (remaining) {
          const int bit = __builtin_ctzll(remaining), degree = 64 * word + bit;
          translated.add(binomials[degree]);
          remaining &= remaining - 1;
        }
      }
      if (!translated.empty()) {
        int word = 0;
        while (!translated.data[word]) ++word;
        minimum_residual_valuation = std::min(minimum_residual_valuation,
                                              word * 64 + __builtin_ctzll(translated.data[word]));
      }
      for (int k = 0; k < int(jet_indices.size()); ++k) {
        const int exponent = jet_indices[k];
        if ((translated.data[exponent / 64] >> (exponent % 64)) & 1) {
          const int row = x * d + y, col = v * d + u;
          jets[k][row * words + col / 64] ^= uint64_t(1) << (col % 64);
        }
      }
      for (int word = 0; word < value.used; ++word) {
        uint64_t remaining = value.data[word];
        while (remaining) {
          const int bit = __builtin_ctzll(remaining), degree = 64 * word + bit;
          const int row = x * d + y, col = v * d + u;
          matrices[degree][row * words + col / 64] ^= uint64_t(1) << (col % 64);
          remaining &= remaining - 1;
        }
      }
    }
    std::cerr << "operator " << x + 1 << "/" << d << " complete\n";
  }
  std::vector<std::array<int, 3>> ranks;
  int maximum_product = 0;
  double rank_sum = 0;
  auto column_rank_of = [&](const std::vector<uint64_t> &matrix) {
    std::vector<uint64_t> columns(size * words);
    for (int row = 0; row < size; ++row) for (int word = 0; word < words; ++word) {
      uint64_t remaining = matrix[row * words + word];
      while (remaining) {
        int bit = __builtin_ctzll(remaining), col = word * 64 + bit;
        int target_row = (row / d) * d + col % d, target_col = (col / d) * d + row % d;
        columns[target_row * words + target_col / 64] ^= uint64_t(1) << (target_col % 64);
        remaining &= remaining - 1;
      }
    }
    return binary_rank(std::move(columns), size, size);
  };
  for (int degree = 0; degree <= degree_bound; ++degree) {
    auto &matrix = matrices[degree];
    if (std::none_of(matrix.begin(), matrix.end(), [](uint64_t value) { return value != 0; })) continue;
    const int row_rank = binary_rank(matrix, size, size);
    const int column_rank = column_rank_of(matrix);
    ranks.push_back({degree - energy - j, row_rank, column_rank});
    maximum_product = std::max(maximum_product, row_rank * column_rank);
    rank_sum += std::sqrt(double(row_rank) * column_rank);
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::cout << "{\"operator_weight\":" << j << ",\"partner_weight\":" << h
            << ",\"state_energy\":" << energy << ",\"labels\":" << d
            << ",\"complete_crossing_entries\":" << checked
            << ",\"maximum_rank_product\":" << maximum_product
            << ",\"rank_sum\":" << std::setprecision(17) << rank_sum
            << ",\"finite_seed_threshold\":" << std::pow(double(d), std::sqrt(2.0))
            << ",\"minimum_residual_taylor_valuation\":" << minimum_residual_valuation
            << ",\"seconds\":" << seconds << ",\"nonzero_profile\":[";
  for (size_t i = 0; i < ranks.size(); ++i) {
    if (i) std::cout << ',';
    const auto &r = ranks[i];
    std::cout << '[' << r[0] << ',' << r[1] << ',' << r[2] << ']';
  }
  std::cout << ']';
  std::cout << ",\"power_two_residual_jet_ranks\":[";
  for (int k = 0; k < int(jet_indices.size()); ++k) {
    if (k) std::cout << ',';
    std::cout << '[' << jet_indices[k] << ','
              << binary_rank(jets[k], size, size) << ',' << column_rank_of(jets[k]) << ']';
  }
  std::cout << ']';
  if (emit_profiles) {
    std::cout << ",\"profiles_words\":[";
    for (size_t i = 0; i < emitted.size(); ++i) {
      if (i) std::cout << ',';
      std::cout << '[';
      for (int k = 0; k < emitted[i].used; ++k) {
        if (k) std::cout << ',';
        std::cout << emitted[i].data[k];
      }
      std::cout << ']';
    }
    std::cout << ']';
  }
  std::cout << "}\n";
  return 0;
} catch (const std::exception &error) {
  std::cerr << error.what() << '\n';
  return 1;
}
