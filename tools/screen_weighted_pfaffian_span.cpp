// Exact F2/F4 ranks for different weighted direct-sum blocks in the two cuts.
// Reuse the audited field elimination, but not the original input layout.
#define main unweighted_pfaffian_screen_main
#include "screen_pfaffian_span_f4.cpp"
#undef main

struct WeightedBlock {
  int size, stride;
  uint64_t multiplicity;
  std::vector<Matrix> generators;
};

std::string decimal128(unsigned __int128 value) {
  if (!value) return "0";
  std::string result;
  while (value) { result.push_back(char('0' + value % 10)); value /= 10; }
  std::reverse(result.begin(), result.end());
  return result;
}

int weighted_block_rank(const WeightedBlock &block, uint64_t code,
                        int cap, bool binary_only) {
  size_t words = size_t(block.size) * block.stride;
  Matrix lo(words, 0), hi(words, 0);
  for (size_t i = 0; i < block.generators.size(); ++i) {
    int scalar = int((code >> (2 * i)) & 3);
    if (scalar & 1)
      for (size_t j = 0; j < words; ++j) lo[j] ^= block.generators[i][j];
    if (scalar & 2)
      for (size_t j = 0; j < words; ++j) hi[j] ^= block.generators[i][j];
  }
  if (binary_only)
    return binary_rank_capped(std::move(lo), block.size, block.stride, cap);
  return field_rank(std::move(lo), std::move(hi), block.size, block.stride, cap);
}

int main(int argc, char **argv) {
  try {
    bool binary_only = false, all = false;
    uint64_t first = 1, count = 0;
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--field2") binary_only = true;
      else if (arg == "--all") all = true;
      else if (arg == "--start" && i + 1 < argc) first = std::stoull(argv[++i]);
      else if (arg == "--count" && i + 1 < argc) count = std::stoull(argv[++i]);
      else throw std::runtime_error("unknown argument: " + arg);
    }
    uint64_t labels, gate;
    int dimension;
    if (!(std::cin >> labels >> dimension >> gate))
      throw std::runtime_error("missing labels, coefficient dimension, and integer cost gate");
    if (labels < 2 || labels > 1000000000 || dimension < 1 ||
        dimension > (binary_only ? 20 : 13) || gate < 2)
      throw std::runtime_error("invalid dimensions");
    uint64_t total_rows = labels * labels;
    std::vector<WeightedBlock> cuts[2];
    for (auto &cut : cuts) {
      int blocks;
      if (!(std::cin >> blocks) || blocks < 1)
        throw std::runtime_error("invalid block count");
      uint64_t covered = 0;
      for (int b = 0; b < blocks; ++b) {
        WeightedBlock block;
        if (!(std::cin >> block.size >> block.multiplicity) ||
            block.size < 1 || !block.multiplicity)
          throw std::runtime_error("invalid weighted block");
        if ((unsigned __int128)block.size * block.multiplicity > total_rows - covered)
          throw std::runtime_error("weighted blocks exceed the complete cut dimension");
        covered += uint64_t(block.size) * block.multiplicity;
        block.stride = (block.size + 63) / 64;
        block.generators.resize(dimension);
        for (Matrix &matrix : block.generators) {
          matrix.resize(size_t(block.size) * block.stride);
          for (uint64_t &word : matrix) {
            std::string value;
            if (!(std::cin >> value)) throw std::runtime_error("truncated matrix input");
            word = std::stoull(value, nullptr, 16);
          }
        }
        cut.push_back(std::move(block));
      }
      if (covered != total_rows)
        throw std::runtime_error("weighted blocks do not cover a complete cut");
    }
    uint64_t limit = uint64_t(1) << ((binary_only ? 1 : 2) * dimension);
    if (!count) count = limit - first;
    if (!first || first >= limit || count > limit - first)
      throw std::runtime_error("invalid coefficient interval");
    unsigned __int128 product_cap = (unsigned __int128)gate * gate - 1;
    uint64_t tested = 0, eligible = 0;
    auto began = std::chrono::steady_clock::now();
    std::cout << "{\"event\":\"start\",\"labels\":" << labels
              << ",\"dimension\":" << dimension << ",\"field_size\":" << (binary_only ? 2 : 4)
              << ",\"integer_cost_gate\":" << gate << ",\"all_ranks\":" << (all ? "true" : "false")
              << ",\"code_start\":" << first << ",\"code_count\":" << count << "}" << std::endl;
    size_t blocks = std::max(cuts[0].size(), cuts[1].size());
    for (uint64_t index = first; index < first + count; ++index) {
      uint64_t code = index;
      if (binary_only) {
        code = 0;
        for (int i = 0; i < dimension; ++i) code |= ((index >> i) & 1) << (2 * i);
      }
      int last = (63 - __builtin_clzll(code)) / 2;
      if (((code >> (2 * last)) & 3) != 1) continue;
      ++tested;
      uint64_t ranks[2] = {0, 0};
      bool rejected = false;
      for (size_t b = 0; b < blocks && !rejected; ++b)
        for (int cut = 0; cut < 2; ++cut) {
          if (b >= cuts[cut].size()) continue;
          const auto &block = cuts[cut][b];
          uint64_t maximum = total_rows;
          if (!all && ranks[1 - cut])
            maximum = uint64_t(std::min<unsigned __int128>(total_rows, product_cap / ranks[1 - cut]));
          if (ranks[cut] > maximum) { rejected = true; break; }
          int cap = int(std::min<uint64_t>(block.size, (maximum - ranks[cut]) / block.multiplicity));
          ranks[cut] += block.multiplicity * weighted_block_rank(block, code, cap, binary_only);
          if (!all && (unsigned __int128)ranks[0] * ranks[1] > product_cap) {
            rejected = true; break;
          }
        }
      if (!rejected) {
        auto product = (unsigned __int128)ranks[0] * ranks[1];
        bool affordable = product <= product_cap;
        eligible += affordable;
        std::cout << "{\"event\":\"" << (all ? "rank" : "eligible") << "\",\"index\":" << index
                  << ",\"code\":" << code << ",\"row_rank\":" << ranks[0]
                  << ",\"column_rank\":" << ranks[1] << ",\"product\":" << decimal128(product)
                  << ",\"eligible\":" << (affordable ? "true" : "false") << "}" << std::endl;
      }
    }
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - began).count();
    std::cout << "{\"event\":\"complete\",\"projective_tested\":" << tested
              << ",\"eligible\":" << eligible << ",\"milliseconds\":" << milliseconds << "}" << std::endl;
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << std::endl;
    return 2;
  }
}
