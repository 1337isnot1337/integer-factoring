#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

// Exact ranks for bivariate, site-weighted two-pair crossings.
// The algebra rank engine is shared mathematically with two_state_run_profile.cpp.
using Polynomial = std::vector<uint64_t>;
using Wide = unsigned __int128;

bool coefficient(const Polynomial& f, unsigned s) {
  return (f[s >> 6] >> (s & 63)) & 1;
}
void toggle(Polynomial& f, unsigned s) { f[s >> 6] ^= uint64_t(1) << (s & 63); }
void set(Polynomial& f, unsigned s) { f[s >> 6] |= uint64_t(1) << (s & 63); }
Polynomial zero(int m) { return Polynomial(((1u << m) + 63) / 64); }
std::vector<unsigned> terms(const Polynomial& f) {
  std::vector<unsigned> answer;
  for (unsigned w = 0; w < f.size(); ++w) {
    uint64_t word = f[w];
    while (word) {
      unsigned j = __builtin_ctzll(word); word &= word - 1;
      answer.push_back(64 * w + j);
    }
  }
  return answer;
}
unsigned compress(unsigned support, unsigned retained) {
  unsigned result = 0, bit = 1;
  while (retained) {
    unsigned low = retained & -retained; retained ^= low;
    if (support & low) result |= bit;
    bit <<= 1;
  }
  return result;
}
Polynomial multiply(const Polynomial& a, const Polynomial& b, int m) {
  Polynomial result = zero(m);
  unsigned all = (1u << m) - 1;
  for (unsigned support : terms(a)) {
    unsigned available = all ^ support, part = available;
    do {
      if (coefficient(b, part)) toggle(result, support | part);
      if (!part) break;
      part = (part - 1) & available;
    } while (true);
  }
  return result;
}

uint32_t binary_rank(std::vector<uint32_t> rows) {
  std::array<uint32_t, 32> basis{};
  uint32_t answer = 0;
  for (uint32_t row : rows) while (row) {
    int pivot = 31 - __builtin_clz(row);
    if (!basis[pivot]) { basis[pivot] = row; ++answer; break; }
    row ^= basis[pivot];
  }
  return answer;
}

struct Eliminator {
  size_t words;
  std::vector<Polynomial> basis;
  uint32_t answer = 0;
  explicit Eliminator(unsigned columns): words((columns + 63) / 64), basis(columns) {}
  void insert(Polynomial row) {
    int top = int(words) - 1;
    while (top >= 0) {
      if (!row[top]) { --top; continue; }
      unsigned pivot = unsigned(64 * top + 63 - __builtin_clzll(row[top]));
      if (basis[pivot].empty()) { basis[pivot] = std::move(row); ++answer; return; }
      for (int w = 0; w <= top; ++w) row[w] ^= basis[pivot][w];
    }
  }
};

struct AlgebraRanker {
  uint64_t calls = 0, hits = 0, quadratics = 0, homogeneous = 0, gaussian = 0;
  std::unordered_map<std::string, uint32_t> cache;

  uint32_t direct(const Polynomial& f, int m) {
    unsigned n = 1u << m, all = n - 1;
    Eliminator matrix(n);
    for (unsigned s = 0; s < n; ++s) {
      Polynomial row = zero(m);
      unsigned available = all ^ s, part = available;
      do {
        if (coefficient(f, part)) set(row, part | s);
        if (!part) break;
        part = (part - 1) & available;
      } while (true);
      matrix.insert(std::move(row));
    }
    return matrix.answer;
  }

  uint32_t homogeneous_rank(const Polynomial& f, int m, int degree,
                            const std::vector<unsigned>& monomials) {
    ++homogeneous;
    unsigned n = 1u << m;
    uint32_t answer = 0;
    // Perfect pairing by complementary monomials identifies the two
    // transpose blocks with source degrees w and m-degree-w.
    for (int weight = 0; 2 * weight <= m - degree; ++weight) {
      std::vector<int> index(n, -1);
      std::vector<unsigned> sources;
      unsigned columns = 0;
      for (unsigned s = 0; s < n; ++s) {
        int pop = __builtin_popcount(s);
        if (pop == weight) sources.push_back(s);
        if (pop == weight + degree) index[s] = columns++;
      }
      Eliminator matrix(columns);
      for (unsigned s : sources) {
        Polynomial row(matrix.words);
        for (unsigned a : monomials) if (!(s & a)) {
          int position = index[s | a]; assert(position >= 0);
          set(row, unsigned(position));
        }
        matrix.insert(std::move(row));
      }
      answer += matrix.answer * (2 * weight == m - degree ? 1 : 2);
    }
    return answer;
  }

  uint32_t rank(const Polynomial& f, int m) {
    ++calls;
    auto monomials = terms(f);
    if (monomials.empty()) return 0;
    if (coefficient(f, 0)) return 1u << m;
    for (int i = 0; i < m; ++i) if (coefficient(f, 1u << i)) return 1u << (m - 1);
    unsigned used = 0, common = (1u << m) - 1;
    int degree = __builtin_popcount(monomials[0]);
    bool is_homogeneous = true;
    for (unsigned a : monomials) {
      used |= a; common &= a;
      if (__builtin_popcount(a) != degree) is_homogeneous = false;
    }
    if (common || used != (1u << m) - 1) {
      unsigned retained = used ^ common;
      int smaller = __builtin_popcount(retained);
      Polynomial reduced = zero(smaller);
      for (unsigned a : monomials) set(reduced, compress(a, retained));
      return (1u << (m - __builtin_popcount(used))) * rank(reduced, smaller);
    }
    std::string key(1, char(m));
    key.append(reinterpret_cast<const char*>(f.data()), f.size() * sizeof(uint64_t));
    auto found = cache.find(key);
    if (found != cache.end()) { ++hits; return found->second; }
    uint32_t answer = 0;
    if (is_homogeneous && degree == 2) {
      // Alternating quadratic forms have a Darboux basis. Each xy pair
      // doubles homology, giving rank (2^m - 2^(m-r))/2 for r pairs.
      std::vector<uint32_t> rows(m);
      for (unsigned a : monomials) {
        int i = __builtin_ctz(a); int j = __builtin_ctz(a ^ (1u << i));
        rows[i] ^= 1u << j; rows[j] ^= 1u << i;
      }
      unsigned paired_rank = binary_rank(rows); assert(!(paired_rank & 1));
      answer = ((1u << m) - (1u << (m - paired_rank / 2))) / 2;
      ++quadratics;
    } else {
      int first = -1, second = -1;
      for (int i = 0; i < m && first < 0; ++i)
        for (int j = i + 1; j < m; ++j)
          if (coefficient(f, (1u << i) | (1u << j))) { first = i; second = j; break; }
      if (first >= 0) {
        unsigned retained = ((1u << m) - 1) ^ (1u << first) ^ (1u << second);
        std::array<Polynomial, 4> block{zero(m - 2), zero(m - 2), zero(m - 2), zero(m - 2)};
        for (unsigned support : monomials) {
          unsigned which = ((support >> first) & 1) + 2 * ((support >> second) & 1);
          set(block[which], compress(support, retained));
        }
        assert(coefficient(block[3], 0));
        assert(!coefficient(block[0], 0) && !coefficient(block[1], 0) && !coefficient(block[2], 0));
        Polynomial g = multiply(multiply(block[3], block[1], m - 2), block[2], m - 2);
        for (unsigned w = 0; w < g.size(); ++w) g[w] ^= block[0][w];
        assert(!coefficient(g, 0));
        ++quadratics;
        answer = (1u << (m - 2)) + 2 * rank(g, m - 2);
      } else if (is_homogeneous) {
        answer = homogeneous_rank(f, m, degree, monomials);
      } else {
        ++gaussian;
        answer = direct(f, m);
      }
    }
    if (cache.size() >= 32768) cache.clear();
    cache.emplace(std::move(key), answer);
    return answer;
  }
};

std::string decimal(Wide value) {
  if (!value) return "0";
  std::string out;
  while (value) { out.push_back(char('0' + value % 10)); value /= 10; }
  std::reverse(out.begin(), out.end()); return out;
}

int main(int argc, char** argv) {
  assert(argc >= 5 && std::string(argv[1]) == "profile");
  int t = std::stoi(argv[2]); assert(1 <= t && t <= 18 && argc == t + 4);
  std::vector<unsigned> weights(t);
  unsigned total_weight = 0;
  for (int i = 0; i < t; ++i) {
    weights[i] = std::stoul(argv[i + 4]);
    total_weight += weights[i];
  }
  assert(total_weight <= 1000000);
  unsigned stride = total_weight + 1;
  std::vector<uint64_t> ranks((t + 1) * stride);
  AlgebraRanker ranker;
  unsigned n = 1u << t;
  for (unsigned active = 0; active < n; ++active) {
    int m = __builtin_popcount(active);
    uint64_t multiplicity = uint64_t(1) << (t - m);
    std::vector<unsigned> local_weights;
    for (int i = 0; i < t; ++i) if (active & (1u << i)) local_weights.push_back(weights[i]);
    std::vector<unsigned> subset_weight(1u << m), subset_size(1u << m);
    std::map<unsigned, Polynomial> coefficients;
    for (unsigned s = 0; s < (1u << m); ++s) {
      if (s) {
        unsigned low = s & -s;
        subset_weight[s] = subset_weight[s ^ low] + local_weights[__builtin_ctz(low)];
        subset_size[s] = subset_size[s ^ low] + 1;
      }
      unsigned key = subset_size[s] * stride + subset_weight[s];
      auto found = coefficients.find(key);
      if (found == coefficients.end()) found = coefficients.emplace(key, zero(m)).first;
      set(found->second, s);
    }
    for (const auto& item : coefficients)
      ranks[item.first] += multiplicity * ranker.rank(item.second, m);
    if (active % 1024 == 1023 || active + 1 == n)
      std::cerr << "active subsets " << active + 1 << "/" << n
                << " homogeneous " << ranker.homogeneous << " cache " << ranker.hits << std::endl;
  }
  Wide maximum_product = 0; long double rank_sum = 0;
  uint64_t integer_lower_bound = 0;
  unsigned pieces = 0;
  for (unsigned k = 0; k <= unsigned(t); ++k) for (unsigned w = 0; w <= total_weight; ++w) {
    uint64_t column = ranks[k * stride + w];
    uint64_t row = ranks[(t - k) * stride + total_weight - w];
    assert(bool(row) == bool(column));
    if (!row) continue;
    ++pieces;
    Wide product = Wide(row) * column;
    maximum_product = std::max(maximum_product, product);
    rank_sum += std::sqrt((long double)product);
  }
  auto output = [&](std::ostream& out) {
    out << "{\"t\":" << t << ",\"weights\":[";
    for (int i = 0; i < t; ++i) out << (i ? "," : "") << weights[i];
    out << "],\"pieces\":[";
    bool first = true;
    for (unsigned k = 0; k <= unsigned(t); ++k) for (unsigned w = 0; w <= total_weight; ++w) {
      uint64_t column = ranks[k * stride + w];
      if (!column) continue;
      uint64_t row = ranks[(t - k) * stride + total_weight - w];
      if (!first) out << ","; first = false;
      out << "[" << k << "," << w << "," << row << "," << column << "]";
    }
    out << "],\"nonzero_pieces\":" << pieces
        << ",\"maximum_rank_product\":" << decimal(maximum_product)
        << ",\"rank_sum_estimate\":" << std::setprecision(18) << rank_sum
        << ",\"maximum_product_rho\":" << std::log2((long double)maximum_product) / (4 * t)
        << ",\"rank_sum_rho\":" << std::log2(rank_sum) / (2 * t)
        << ",\"algebra_rank_calls\":" << ranker.calls
        << ",\"quadratic_reductions\":" << ranker.quadratics
        << ",\"homogeneous_ranks\":" << ranker.homogeneous
        << ",\"generic_ranks\":" << ranker.gaussian << ",\"cache_hits\":" << ranker.hits << "}\n";
  };
  std::ofstream out(argv[3]); assert(out); output(out); output(std::cout);
}
