#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Full local even-parity support: K_(x,g,h)(z)=[g=h]I+(1+z)B_(x,g,h).
// No simultaneous-bit-flip symmetry or XOR-circulant assumption is used.
using Row = std::array<uint64_t, 4>;
struct VectorHash {
  size_t operator()(const std::vector<uint64_t>& values) const {
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    for (uint64_t word : values) h ^= word + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
  }
};
struct Ranker {
  int n, words;
  size_t limit;
  uint64_t calls = 0, hits = 0, gaussian = 0;
  std::unordered_map<std::vector<uint64_t>, uint16_t, VectorHash> cache;
  explicit Ranker(int t): n(1 << t), words((n + 63) / 64) {
    limit = std::min<size_t>(65536, (1u << 25) / (8 * n * words + 64));
    cache.reserve(limit);
  }
  unsigned rank(const Row* rows) {
    ++calls;
    uint64_t nonzero = 0;
    std::vector<uint64_t> key; key.reserve(n * words);
    for (int i = 0; i < n; ++i) for (int w = 0; w < words; ++w) {
      nonzero |= rows[i][w]; key.push_back(rows[i][w]);
    }
    if (!nonzero) return 0;
    auto found = cache.find(key);
    if (found != cache.end()) { ++hits; return found->second; }
    std::vector<Row> basis(n);
    std::vector<uint8_t> present(n);
    unsigned answer = 0;
    for (int i = 0; i < n; ++i) {
      Row row = rows[i]; int top = words - 1;
      while (top >= 0) {
        if (!row[top]) { --top; continue; }
        int pivot = 64 * top + 63 - __builtin_clzll(row[top]);
        if (!present[pivot]) { present[pivot] = 1; basis[pivot] = row; ++answer; break; }
        for (int w = 0; w <= top; ++w) row[w] ^= basis[pivot][w];
      }
    }
    ++gaussian;
    if (cache.size() >= limit) cache.clear();
    cache.emplace(std::move(key), uint16_t(answer));
    return answer;
  }
};

struct Profile {
  std::array<unsigned, 8> matrices{};
  int boundary = 0, off_shift = 0;
  std::vector<unsigned> row, column;
  uint64_t maximum_product = 0;
  long double rank_sum = 0;
  std::vector<unsigned> flip_witness, order_witness, off_crossing_witness;
  std::vector<unsigned> identity_letters, paired_identity_letters;
  int order_swap = -1;
};

uint64_t evaluate_word(const std::vector<unsigned>& word, const Profile& p) {
  uint64_t dp[4] = {1, 0, 0, 1};
  for (unsigned q : word) {
    unsigned b = p.matrices[q]; bool diagonal = ((q >> 1) & 1) == (q & 1);
    uint64_t next[4] = {};
    for (int i = 0; i < (p.boundary ? 2 : 1); ++i) for (int j = 0; j < 2; ++j) {
      if (diagonal) next[2 * i + j] ^= dp[2 * i + j];
      for (int k = 0; k < 2; ++k) if (b & (1u << (2 * k + j))) {
        int shift = k == j ? 0 : p.off_shift;
        next[2 * i + j] ^= (dp[2 * i + k] << shift) ^ (dp[2 * i + k] << (shift + 1));
      }
    }
    std::copy(next, next + 4, dp);
  }
  return p.boundary ? dp[1] ^ dp[2] ^ dp[3] : dp[0];
}

void short_witnesses(Profile& p) {
  assert(3 * (1 + p.off_shift) < 64);
  for (unsigned q = 0; q < 8; ++q) if (((q >> 1) & 1) == (q & 1) && !p.matrices[q]) {
    p.identity_letters.push_back(q);
    if (!p.matrices[q ^ 4]) p.paired_identity_letters.push_back(q);
  }
  for (int length = 1; length <= 3; ++length) for (unsigned code = 0; code < (1u << (3 * length)); ++code) {
    std::vector<unsigned> word(length);
    bool off = false;
    for (int i = 0; i < length; ++i) {
      word[i] = (code >> (3 * i)) & 7;
      off |= ((word[i] >> 1) & 1) != (word[i] & 1);
    }
    uint64_t value = evaluate_word(word, p);
    if (off && value && p.off_crossing_witness.empty()) p.off_crossing_witness = word;
    if (p.flip_witness.empty()) {
      auto flipped = word; for (auto& q : flipped) q ^= 4;
      if (evaluate_word(flipped, p) != value) p.flip_witness = word;
    }
    if (p.order_witness.empty()) for (int i = 0; i + 1 < length; ++i) {
      auto swapped = word; std::swap(swapped[i], swapped[i + 1]);
      if (evaluate_word(swapped, p) != value) { p.order_witness = word; p.order_swap = i; break; }
    }
    if (!p.flip_witness.empty() && !p.order_witness.empty() && !p.off_crossing_witness.empty()) return;
  }
}

Profile evaluate(int t, const std::array<unsigned, 8>& matrices, int boundary, int off_shift, Ranker& ranker) {
  int n = 1 << t;
  int maximum_degree = t * (1 + off_shift); assert(maximum_degree < 64);
  std::vector<Row> rows((maximum_degree + 1) * n * n), columns((maximum_degree + 1) * n * n);
  Profile result; result.matrices = matrices; result.boundary = boundary; result.off_shift = off_shift;
  short_witnesses(result);
  for (int g = 0; g < n; ++g) for (int h = 0; h < n; ++h) for (int x = 0; x < n; ++x) {
    uint64_t dp[4] = {1, 0, 0, 1};
    for (int position = 0; position < t; ++position) {
      unsigned gi = (g >> position) & 1, hi = (h >> position) & 1, xi = (x >> position) & 1;
      unsigned b = matrices[4 * xi + 2 * gi + hi];
      uint64_t next[4] = {};
      for (int i = 0; i < (boundary ? 2 : 1); ++i) for (int j = 0; j < 2; ++j) {
        if (gi == hi) next[2 * i + j] ^= dp[2 * i + j];
        for (int k = 0; k < 2; ++k) if (b & (1u << (2 * k + j))) {
          int shift = k == j ? 0 : off_shift;
          next[2 * i + j] ^= (dp[2 * i + k] << shift) ^ (dp[2 * i + k] << (shift + 1));
        }
      }
      std::copy(next, next + 4, dp);
    }
    uint64_t polynomial = boundary ? dp[1] ^ dp[2] ^ dp[3] : dp[0];
    assert(__builtin_parityll(polynomial) == (g == h));
    while (polynomial) {
      int degree = __builtin_ctzll(polynomial); polynomial &= polynomial - 1;
      int u = x ^ h, y = x ^ g;
      rows[(degree * n + g) * n + x][u >> 6] |= uint64_t(1) << (u & 63);
      columns[(degree * n + h) * n + x][y >> 6] |= uint64_t(1) << (y & 63);
    }
  }
  result.row.resize(maximum_degree + 1); result.column.resize(maximum_degree + 1);
  for (int d = 0; d <= maximum_degree; ++d) {
    for (int block = 0; block < n; ++block) {
      result.row[d] += ranker.rank(rows.data() + (d * n + block) * n);
      result.column[d] += ranker.rank(columns.data() + (d * n + block) * n);
    }
    assert(bool(result.row[d]) == bool(result.column[d]));
    uint64_t product = uint64_t(result.row[d]) * result.column[d];
    result.maximum_product = std::max(result.maximum_product, product);
    result.rank_sum += std::sqrt((long double)product);
  }
  return result;
}

void output(std::ostream& out, const Profile& p, int t) {
  out << "{\"t\":" << t << ",\"boundary\":" << p.boundary << ",\"off_shift\":" << p.off_shift << ",\"matrices\":[";
  for (int i = 0; i < 8; ++i) out << (i ? "," : "") << p.matrices[i];
  out << "],\"row_ranks\":[";
  for (size_t d = 0; d < p.row.size(); ++d) out << (d ? "," : "") << p.row[d];
  out << "],\"column_ranks\":[";
  for (size_t d = 0; d < p.column.size(); ++d) out << (d ? "," : "") << p.column[d];
  auto array = [&](const char* name, const std::vector<unsigned>& values) {
    out << ",\"" << name << "\":[";
    for (size_t i = 0; i < values.size(); ++i) out << (i ? "," : "") << values[i];
    out << "]";
  };
  out << "]";
  array("bit_flip_witness", p.flip_witness); array("order_witness", p.order_witness);
  array("off_crossing_witness", p.off_crossing_witness);
  array("identity_letters", p.identity_letters); array("paired_identity_letters", p.paired_identity_letters);
  out << ",\"order_swap\":" << p.order_swap << ",\"maximum_rank_product\":" << p.maximum_product
      << ",\"rank_sum_estimate\":" << std::setprecision(18) << p.rank_sum
      << ",\"maximum_product_rho\":" << std::log2((long double)p.maximum_product) / (4 * t)
      << ",\"rank_sum_rho\":" << std::log2(p.rank_sum) / (2 * t) << "}\n";
}

int main(int argc, char** argv) {
  assert(argc >= 5);
  std::string mode = argv[1]; int t = std::stoi(argv[2]), boundary = std::stoi(argv[4]);
  assert(1 <= t && t <= 8 && (boundary == 0 || boundary == 1));
  Ranker ranker(t);
  std::ofstream out(argv[3]); assert(out);
  out.exceptions(std::ios::badbit | std::ios::failbit);
  if (mode == "profile") {
    assert(argc == 13 || argc == 14);
    int off_shift = argc == 14 ? std::stoi(argv[13]) : 0; assert(off_shift >= 0);
    std::array<unsigned, 8> matrices{};
    for (int i = 0; i < 8; ++i) { matrices[i] = std::stoul(argv[5 + i]); assert(matrices[i] < 16); }
    auto p = evaluate(t, matrices, boundary, off_shift, ranker);
    output(out, p, t); out.flush(); output(std::cout, p, t); return 0;
  }
  assert(mode == "screen" && (argc == 6 || argc == 7));
  int off_shift = argc == 7 ? std::stoi(argv[6]) : 0; assert(off_shift >= 0);
  std::ifstream candidates(argv[5]); assert(candidates);
  unsigned diagonal0, diagonal1, upper, lower, count = 0;
  uint64_t best_product = UINT64_MAX; long double best_sum = 1e100L;
  while (candidates >> diagonal0 >> diagonal1 >> upper >> lower) {
    assert((diagonal0 | diagonal1 | upper | lower) < 256);
    std::array<unsigned, 8> matrices{};
    for (int q = 0; q < 8; ++q)
      matrices[q] = ((diagonal0 >> q) & 1) | (((upper >> q) & 1) << 1)
                  | (((lower >> q) & 1) << 2) | (((diagonal1 >> q) & 1) << 3);
    auto p = evaluate(t, matrices, boundary, off_shift, ranker); output(out, p, t);
    best_product = std::min(best_product, p.maximum_product);
    best_sum = std::min(best_sum, p.rank_sum); ++count;
    if (count % 4096 == 0)
      std::cerr << "checked " << count << " boundary " << boundary
                << " best product " << best_product << " best sum " << double(best_sum) << std::endl;
  }
  out.flush();
  std::cout << "{\"checked\":" << count << ",\"boundary\":" << boundary
            << ",\"best_maximum_product\":" << best_product
            << ",\"best_rank_sum_estimate\":" << std::setprecision(18) << best_sum
            << ",\"rank_calls\":" << ranker.calls << ",\"gaussian_ranks\":" << ranker.gaussian
            << ",\"cache_hits\":" << ranker.hits << "}\n";
}
