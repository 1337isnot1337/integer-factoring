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

// K_(g,h)(z) = [g=h] I + (1+z) B_(g,h), with four independent
// binary 2x2 matrices. The endpoint is the exact binary-word crossing.
struct Ranker {
  int t, n;
  std::array<uint64_t, 6> lower{};
  std::unordered_map<uint64_t, uint8_t> cache;
  uint64_t calls = 0, hits = 0, gaussian = 0;
  explicit Ranker(int bits): t(bits), n(1 << bits) {
    for (int j = 0; j < t; ++j) for (int h = 0; h < n; ++h)
      if (!(h & (1 << j))) lower[j] |= uint64_t(1) << h;
    cache.reserve(1 << 18);
  }
  uint64_t translate(uint64_t word, int shift) const {
    for (int j = 0; j < t; ++j) if (shift & (1 << j)) {
      int d = 1 << j;
      word = ((word & lower[j]) << d) | ((word & ~lower[j]) >> d);
    }
    return word;
  }
  int direct(uint64_t kernel) const {
    std::array<uint64_t, 64> basis{};
    int answer = 0;
    for (int x = 0; x < n; ++x) {
      uint64_t word = translate(kernel, x);
      while (word) {
        int pivot = 63 - __builtin_clzll(word);
        if (!basis[pivot]) { basis[pivot] = word; ++answer; break; }
        word ^= basis[pivot];
      }
    }
    return answer;
  }
  int rank(uint64_t kernel) {
    ++calls;
    if (!kernel) return 0;
    if (__builtin_parityll(kernel)) return n;
    for (int j = 0; j < t; ++j)
      if (__builtin_parityll(kernel & ~lower[j])) return n / 2;
    auto found = cache.find(kernel);
    if (found != cache.end()) { ++hits; return found->second; }
    int answer = direct(kernel); ++gaussian;
    if (cache.size() >= (1u << 21)) cache.clear();
    cache.emplace(kernel, uint8_t(answer));
    return answer;
  }
};

struct Profile {
  std::array<unsigned, 4> matrices{};
  int boundary = 0;
  std::vector<unsigned> row, column;
  uint64_t maximum_product = 0;
  long double rank_sum = 0;
};

Profile evaluate(int t, const std::array<unsigned, 4>& matrices,
                 int boundary, Ranker& ranker) {
  int n = 1 << t;
  Profile result; result.matrices = matrices; result.boundary = boundary;
  std::vector<std::vector<uint64_t>> rows(t + 1, std::vector<uint64_t>(n));
  std::vector<std::vector<uint64_t>> columns(t + 1, std::vector<uint64_t>(n));
  for (int g = 0; g < n; ++g) for (int h = 0; h < n; ++h) {
    uint64_t dp[4] = {1, 0, 0, 1};
    for (int position = 0; position < t; ++position) {
      unsigned gi = (g >> position) & 1, hi = (h >> position) & 1;
      unsigned b = matrices[2 * gi + hi];
      uint64_t next[4] = {};
      for (int i = 0; i < (boundary ? 2 : 1); ++i) for (int j = 0; j < 2; ++j) {
        if (gi == hi) next[2 * i + j] ^= dp[2 * i + j];
        for (int k = 0; k < 2; ++k) if (b & (1u << (2 * k + j)))
          next[2 * i + j] ^= dp[2 * i + k] ^ (dp[2 * i + k] << 1);
      }
      std::copy(next, next + 4, dp);
    }
    uint64_t polynomial = boundary ? dp[1] ^ dp[2] ^ dp[3] : dp[0];
    assert(__builtin_parityll(polynomial) == (g == h));
    while (polynomial) {
      int degree = __builtin_ctzll(polynomial); polynomial &= polynomial - 1;
      rows[degree][g] |= uint64_t(1) << h;
      columns[degree][h] |= uint64_t(1) << g;
    }
  }
  result.row.resize(t + 1); result.column.resize(t + 1);
  for (int d = 0; d <= t; ++d) {
    for (int g = 0; g < n; ++g) {
      result.row[d] += ranker.rank(rows[d][g]);
      result.column[d] += ranker.rank(columns[d][g]);
    }
    assert(bool(result.row[d]) == bool(result.column[d]));
    uint64_t p = uint64_t(result.row[d]) * result.column[d];
    result.maximum_product = std::max(result.maximum_product, p);
    result.rank_sum += std::sqrt((long double)p);
  }
  return result;
}

void output(std::ostream& out, const Profile& p, int t) {
  out << "{\"t\":" << t << ",\"boundary\":" << p.boundary << ",\"matrices\":[";
  for (int i = 0; i < 4; ++i) out << (i ? "," : "") << p.matrices[i];
  out << "],\"row_ranks\":[";
  for (int i = 0; i <= t; ++i) out << (i ? "," : "") << p.row[i];
  out << "],\"column_ranks\":[";
  for (int i = 0; i <= t; ++i) out << (i ? "," : "") << p.column[i];
  out << "],\"maximum_rank_product\":" << p.maximum_product
      << ",\"rank_sum_estimate\":" << std::setprecision(18) << p.rank_sum
      << ",\"maximum_product_rho\":" << std::log2((long double)p.maximum_product) / (4 * t)
      << ",\"rank_sum_rho\":" << std::log2(p.rank_sum) / (2 * t) << "}\n";
}

int main(int argc, char** argv) {
  assert(argc >= 5);
  std::string mode = argv[1]; int t = std::stoi(argv[2]);
  assert(1 <= t && t <= 6);
  int boundary = std::stoi(argv[4]); assert(boundary == 0 || boundary == 1);
  Ranker ranker(t);
  std::ofstream out(argv[3]); assert(out);
  if (mode == "profile") {
    assert(argc == 9);
    std::array<unsigned, 4> matrices{};
    for (int i = 0; i < 4; ++i) { matrices[i] = std::stoul(argv[5 + i]); assert(matrices[i] < 16); }
    auto p = evaluate(t, matrices, boundary, ranker);
    output(out, p, t); output(std::cout, p, t); return 0;
  }
  assert(mode == "screen" && argc == 5);
  uint64_t best_product = UINT64_MAX; long double best_sum = 1e100L;
  for (unsigned code = 0; code < 65536; ++code) {
    std::array<unsigned, 4> matrices{};
    for (int i = 0; i < 4; ++i) matrices[i] = (code >> (4 * i)) & 15;
    auto p = evaluate(t, matrices, boundary, ranker);
    output(out, p, t);
    best_product = std::min(best_product, p.maximum_product);
    best_sum = std::min(best_sum, p.rank_sum);
    if (code % 4096 == 4095)
      std::cerr << "checked " << code + 1 << " boundary " << boundary
                << " best product " << best_product << " best sum " << double(best_sum) << std::endl;
  }
  std::cout << "{\"checked\":65536,\"boundary\":" << boundary
            << ",\"best_maximum_product\":" << best_product
            << ",\"best_rank_sum_estimate\":" << std::setprecision(18) << best_sum
            << ",\"rank_calls\":" << ranker.calls << ",\"gaussian_ranks\":" << ranker.gaussian
            << ",\"cache_hits\":" << ranker.hits << "}\n";
}
