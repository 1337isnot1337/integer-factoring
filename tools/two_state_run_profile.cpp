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

// Exact ranks in F_2[e_1,...,e_m]/(e_i^2), followed by the exact
// decomposition (J_2 + 2 J_1)^{tensor t} for the two local pairings.
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
  assert(argc >= 3);
  if (std::string(argv[1]) == "selftest") {
    int m = std::stoi(argv[2]); assert(1 <= m && m <= 4);
    AlgebraRanker ranker;
    uint64_t count = uint64_t(1) << (1u << m);
    for (uint64_t f = 0; f < count; ++f)
      assert(ranker.rank(Polynomial{f}, m) == ranker.direct(Polynomial{f}, m));
    std::cout << "{\"status\":\"passed\",\"polynomials\":" << count << "}\n";
    return 0;
  }
  assert(argc == 10 && std::string(argv[1]) == "profile");
  int t = std::stoi(argv[2]); assert(1 <= t && t <= 18);
  std::array<int, 2> vertex{std::stoi(argv[4]), std::stoi(argv[5])}; // S,C
  std::array<int, 4> edge{}; // SS,SC,CS,CC
  for (int i = 0; i < 4; ++i) edge[i] = std::stoi(argv[6 + i]);
  for (int value : vertex) assert(value >= 0);
  for (int value : edge) assert(value >= 0);
  int maximum_degree = t * std::max(vertex[0], vertex[1])
                     + (t - 1) * *std::max_element(edge.begin(), edge.end());
  struct Group { unsigned representative = 0; uint64_t count = 0; };
  std::map<std::vector<int>, Group> groups;
  unsigned n = 1u << t;
  for (unsigned active = 0; active < n; ++active) {
    std::vector<int> key;
    for (int j = 0; j < t;) {
      if (!(active & (1u << j))) { ++j; continue; }
      int start = j;
      while (j < t && (active & (1u << j))) ++j;
      int boundaries = (start > 0 ? 1 : 0) + (j < t ? 2 : 0);
      key.push_back(32 * boundaries + j - start);
    }
    std::sort(key.begin(), key.end());
    Group& group = groups[key];
    if (!group.count) group.representative = active;
    ++group.count;
  }
  std::vector<uint64_t> row(maximum_degree + 1), column(maximum_degree + 1);
  AlgebraRanker ranker;
  unsigned completed = 0;
  for (const auto& item : groups) {
    unsigned active = item.second.representative;
    int m = __builtin_popcount(active);
    uint64_t multiplicity = item.second.count * (uint64_t(1) << (t - m));
    std::vector<unsigned> expanded(1u << m);
    std::vector<int> positions;
    for (int i = 0; i < t; ++i) if (active & (1u << i)) positions.push_back(i);
    for (unsigned s = 1; s < expanded.size(); ++s) {
      unsigned low = s & -s;
      expanded[s] = expanded[s ^ low] | (1u << positions[__builtin_ctz(low)]);
    }
    for (int cut = 0; cut < 2; ++cut) {
      std::vector<Polynomial> coefficients(maximum_degree + 1, zero(m));
      for (unsigned s = 0; s < expanded.size(); ++s) {
        unsigned word = cut ? expanded[s] : (n - 1) ^ expanded[s]; // C positions
        int previous = word & 1, cost = vertex[previous];
        for (int i = 1; i < t; ++i) {
          int symbol = (word >> i) & 1;
          cost += vertex[symbol] + edge[2 * previous + symbol]; previous = symbol;
        }
        set(coefficients[cost], s);
      }
      auto& target = cut ? column : row;
      for (int degree = 0; degree <= maximum_degree; ++degree)
        target[degree] += multiplicity * ranker.rank(coefficients[degree], m);
    }
    ++completed;
    if (completed % 32 == 0 || completed == groups.size())
      std::cerr << "groups " << completed << "/" << groups.size()
                << " quadratic " << ranker.quadratics << " homogeneous " << ranker.homogeneous
                << " generic " << ranker.gaussian << " cache " << ranker.hits << std::endl;
  }
  Wide maximum_product = 0; long double rank_sum = 0;
  for (int d = 0; d <= maximum_degree; ++d) {
    assert(bool(row[d]) == bool(column[d]));
    Wide product = Wide(row[d]) * column[d];
    maximum_product = std::max(maximum_product, product);
    rank_sum += std::sqrt((long double)product);
  }
  auto output = [&](std::ostream& out) {
    out << "{\"t\":" << t << ",\"method\":\"two-state-run-types\",\"vertex\":["
        << vertex[0] << "," << vertex[0] << "," << vertex[1] << "],\"edge\":[";
    int lifted[9] = {edge[0], edge[0], edge[1], edge[0], edge[0], edge[1], edge[2], edge[2], edge[3]};
    for (int i = 0; i < 9; ++i) out << (i ? "," : "") << lifted[i];
    out << "],\"row_ranks\":[";
    for (int d = 0; d <= maximum_degree; ++d) out << (d ? "," : "") << row[d];
    out << "],\"column_ranks\":[";
    for (int d = 0; d <= maximum_degree; ++d) out << (d ? "," : "") << column[d];
    out << "],\"maximum_rank_product\":" << decimal(maximum_product)
        << ",\"rank_sum_estimate\":" << std::setprecision(18) << rank_sum
        << ",\"maximum_product_rho\":" << std::log2((long double)maximum_product) / (4 * t)
        << ",\"rank_sum_rho\":" << std::log2(rank_sum) / (2 * t)
        << ",\"run_types\":" << groups.size() << ",\"algebra_rank_calls\":" << ranker.calls
        << ",\"quadratic_reductions\":" << ranker.quadratics
        << ",\"homogeneous_ranks\":" << ranker.homogeneous
        << ",\"generic_ranks\":" << ranker.gaussian << ",\"cache_hits\":" << ranker.hits << "}\n";
  };
  std::ofstream out(argv[3]); assert(out); output(out); output(std::cout);
}
