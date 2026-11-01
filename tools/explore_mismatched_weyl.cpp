#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using Polynomial = unsigned __int128;
using Label = std::vector<int>;
using Blocks = std::unordered_map<uint64_t, std::vector<uint32_t>>;

int first_bit(Polynomial value) {
  uint64_t low = static_cast<uint64_t>(value);
  return low ? __builtin_ctzll(low)
             : 64 + __builtin_ctzll(static_cast<uint64_t>(value >> 64));
}

bool choose_odd(int n, int k) {
  return k >= 0 && k <= n && (k & ~n) == 0;
}

Polynomial multiply(Polynomial a, Polynomial b) {
  Polynomial result = 0;
  while (b) {
    result ^= a << first_bit(b);
    b &= b - 1;
  }
  return result;
}

int rank_slice(uint64_t *matrix, int rows, int words) {
  std::vector<int> pivots(rows, -1);
  int result = 0;
  for (int i = 0; i < rows; ++i) {
    uint64_t *row = matrix + static_cast<size_t>(i)*words;
    int word = words - 1;
    while (word >= 0) {
      while (word >= 0 && row[word] == 0) --word;
      if (word < 0) break;
      int bit = 63 - __builtin_clzll(row[word]);
      int column = 64*word + bit;
      if (pivots[column] < 0) {
        pivots[column] = i;
        ++result;
        break;
      }
      const uint64_t *pivot = matrix + static_cast<size_t>(pivots[column])*words;
      for (int j = 0; j <= word; ++j) row[j] ^= pivot[j];
    }
  }
  return result;
}

int main(int argc, char **argv) {
  int energy = argc > 1 ? std::stoi(argv[1]) : 63;
  int radix = 3;
  bool stats_only = false;
  std::string dump_prefix;
  for (int i=2;i<argc;++i) {
    std::string option=argv[i];
    if(option=="--stats-only") stats_only=true;
    else if(option=="--radix" && i+1<argc) radix=std::stoi(argv[++i]);
    else if(option=="--dump-prefix" && i+1<argc) dump_prefix=argv[++i];
    else return 2;
  }
  if (energy < 1 || energy > 127 || radix<2) return 2;
  auto started = std::chrono::steady_clock::now();
  std::vector<int> weights;
  for (int w=1;w<=energy;w*=radix) weights.push_back(w);
  std::vector<Label> labels;
  Label building(weights.size());
  std::function<void(int,int)> generate = [&](int index, int remaining) {
    if (index == 0) {
      building[0] = remaining;
      labels.push_back(building);
      return;
    }
    for (int value=0;value*weights[index]<=remaining;++value) {
      building[index] = value;
      generate(index-1,remaining-value*weights[index]);
    }
  };
  generate(static_cast<int>(weights.size())-1,energy);
  const int d = labels.size();
  if (static_cast<uint64_t>(d)*d > UINT32_MAX) return 3;
  std::vector<int> shifts, widths;
  int total_bits = 0;
  for (int w: weights) {
    shifts.push_back(total_bits);
    int bits=1;
    while ((1<<bits)<=2*(energy/w)) ++bits;
    widths.push_back(bits);
    total_bits += bits;
  }
  assert(total_bits <= 64);
  Blocks row_blocks, col_blocks;
  for (int a=0;a<d;++a) for (int b=0;b<d;++b) {
    uint64_t sum=0,difference=0;
    for (size_t i=0;i<weights.size();++i) {
      sum |= static_cast<uint64_t>(labels[a][i]+labels[b][i]) << shifts[i];
      difference |= static_cast<uint64_t>(labels[a][i]-labels[b][i]+energy/weights[i]) << shifts[i];
    }
    uint32_t code = static_cast<uint32_t>(a*d+b);
    row_blocks[difference].push_back(code);
    col_blocks[sum].push_back(code);
  }
  auto block_stats = [](const Blocks &blocks) {
    uint64_t entries=0, maximum=0;
    for (const auto &item: blocks) {
      uint64_t size=item.second.size();
      entries += size*size;
      maximum=std::max(maximum,size);
    }
    return std::pair<uint64_t,uint64_t>(entries,maximum);
  };
  auto row_stats=block_stats(row_blocks), col_stats=block_stats(col_blocks);
  std::cerr << "energy=" << energy << " labels=" << d
            << " row_blocks=" << row_blocks.size() << " col_blocks=" << col_blocks.size()
            << " entries_per_cut=" << row_stats.first << "," << col_stats.first
            << " max_block=" << row_stats.second << "," << col_stats.second << "\n";
  assert(row_stats.first==col_stats.first);
  if (stats_only) return 0;
  std::vector<Polynomial> translation_masks;
  for (int step=1;step<=energy;step*=2) {
    Polynomial mask=0;
    for (int i=0;i<=energy;++i) if(i&step) mask |= Polynomial(1)<<i;
    translation_masks.push_back(mask);
  }
  uint64_t crossing_checks=0;
  std::array<Polynomial,128> coefficient_pivots{};
  int coefficient_span_dimension=0;
  auto entry = [&](int a, int y, int v) {
    Polynomial answer=1;
    for (size_t j=0;j<weights.size();++j) {
      Polynomial local=0;
      for (int k=0;k<=std::min(labels[a][j],labels[y][j]);++k)
        if (choose_odd(labels[y][j],k) && choose_odd(labels[v][j],labels[a][j]-k))
          local |= Polynomial(1)<<(weights[j]*k);
      answer=multiply(answer,local);
      if (!answer) break;
    }
    int step=1;
    for (Polynomial mask: translation_masks) {
      answer ^= (answer&mask)>>step;
      step*=2;
    }
    int parity=__builtin_parityll(static_cast<uint64_t>(answer)) ^
               __builtin_parityll(static_cast<uint64_t>(answer>>64));
    assert(parity==(a==v));
    Polynomial residual=answer;
    while(residual) {
      int pivot=first_bit(residual);
      if(coefficient_pivots[pivot])residual^=coefficient_pivots[pivot];
      else {
        coefficient_pivots[pivot]=residual;
        ++coefficient_span_dimension;
        break;
      }
    }
    ++crossing_checks;
    return answer;
  };
  auto compute = [&](const Blocks &blocks, bool row_cut) {
    std::vector<uint64_t> ranks(energy+1);
    uint64_t finished=0;
    std::ofstream dump;
    auto write32=[&](uint32_t value) {
      for(int b=0;b<4;++b)dump.put(static_cast<char>((value>>(8*b))&255));
    };
    auto write64=[&](uint64_t value) {
      for(int b=0;b<8;++b)dump.put(static_cast<char>((value>>(8*b))&255));
    };
    if(!dump_prefix.empty()) {
      dump.open(dump_prefix+(row_cut?".row.bin":".col.bin"),std::ios::binary);
      if(!dump)throw std::runtime_error("cannot open block dump");
      dump.write("WBIN001\0",8);
      for(uint32_t value:{static_cast<uint32_t>(energy),static_cast<uint32_t>(d),
                         static_cast<uint32_t>(blocks.size()),static_cast<uint32_t>(energy<64?1:2),
                         static_cast<uint32_t>(weights.size()),static_cast<uint32_t>(row_cut)})write32(value);
      for(int value:weights)write32(value);
      for(const Label &label:labels)for(int value:label)write32(value);
    }
    for (const auto &item: blocks) {
      int n=item.second.size(), words=(n+63)/64;
      std::vector<std::pair<int,int>> pairs;
      for (uint32_t code:item.second) pairs.emplace_back(code/d,code%d);
      if(dump.is_open()) {
        write32(n);
        for(uint32_t code:item.second)write32(code);
      }
      size_t plane=static_cast<size_t>(n)*words;
      std::vector<uint64_t> matrices(static_cast<size_t>(energy+1)*plane);
      for (int i=0;i<n;++i) for (int j=0;j<n;++j) {
        int a=pairs[i].first;
        int y=row_cut?pairs[i].second:pairs[j].first;
        int v=row_cut?pairs[j].first:pairs[j].second;
        Polynomial polynomial=entry(a,y,v);
        if(dump.is_open()) {
          write64(static_cast<uint64_t>(polynomial));
          if(energy>=64)write64(static_cast<uint64_t>(polynomial>>64));
        }
        while (polynomial) {
          int degree=first_bit(polynomial);
          matrices[static_cast<size_t>(degree)*plane+static_cast<size_t>(i)*words+j/64] |= uint64_t(1)<<(j%64);
          polynomial &= polynomial-1;
        }
      }
      for (int k=0;k<=energy;++k)
        ranks[k]+=rank_slice(matrices.data()+static_cast<size_t>(k)*plane,n,words);
      if ((++finished%10000)==0) std::cerr << (row_cut?"row":"col") << " blocks " << finished << "/" << blocks.size() << "\n";
    }
    return ranks;
  };
  auto rows=compute(row_blocks,true), cols=compute(col_blocks,false);
  long double sum=0;
  int pieces=0;
  for (int k=0;k<=energy;++k) {
    assert((rows[k]==0)==(cols[k]==0));
    if(rows[k])++pieces;
    sum+=std::sqrt(static_cast<long double>(rows[k])*cols[k]);
  }
  long double target=std::pow(static_cast<long double>(d),std::sqrt(2.0L));
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  std::cout << std::setprecision(18) << "{\n  \"field\":2,\n  \"energy\":" << energy
            << ",\n  \"weight_radix\":" << radix
            << ",\n  \"labels\":" << d << ",\n  \"nonzero_pieces\":" << pieces
            << ",\n  \"coefficient_span_dimension\":" << coefficient_span_dimension
            << ",\n  \"crossing_entry_checks\":" << crossing_checks
            << ",\n  \"row_block_count\":" << row_blocks.size()
            << ",\n  \"col_block_count\":" << col_blocks.size()
            << ",\n  \"entries_per_cut\":" << row_stats.first
            << ",\n  \"row_ranks\":[";
  for(int k=0;k<=energy;++k)std::cout << (k?",":"") << rows[k];
  std::cout << "],\n  \"col_ranks\":[";
  for(int k=0;k<=energy;++k)std::cout << (k?",":"") << cols[k];
  std::cout << "],\n  \"rank_sum_decimal\":" << sum
            << ",\n  \"half_threshold_decimal\":" << target
            << ",\n  \"rho_decimal\":" << std::log(sum)/(2*std::log(static_cast<long double>(d)))
            << ",\n  \"seconds\":" << seconds
            << ",\n  \"scope\":\"Exact finite coefficient ranks; no asymptotic rank theorem or improved score.\"\n}\n";
}
